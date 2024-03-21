# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     James Krieger (jmkrieger@cnb.csic.es)
# *
# * Centro Nacional de Biotecnologia, CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************


"""
This module will provide ProDy mode editing tools.
"""
import os
import numpy as np

from pwem.objects import AtomStruct, SetOfNormalModes, SetOfPrincipalComponents, String

from pyworkflow.utils import glob, logger
from pyworkflow.protocol.params import (PointerParam, EnumParam, BooleanParam,
                                        FloatParam, IntParam, LEVEL_ADVANCED)

import prody
from prody2.protocols.protocol_modes_base import ProDyModesBase
from prody2 import fixVerbositySecondary, restoreVerbositySecondary

NMA_SLICE = 0
NMA_REDUCE = 1
NMA_EXTEND = 2
NMA_INTERP = 3

class ProDyEdit(ProDyModesBase):
    """
    This protocol will edit a SetOfNormalModes object to have more or fewer nodes
    """
    _label = 'Edit modes'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label='ProDy edit')

        form.addParam('modes', PointerParam, label='Input set of modes',
                      pointerClass='SetOfNormalModes',
                      help='The input modes can be a SetOfNormalModes '
                           'from an atomic model (true PDB) or a pseudoatomic model '
                           '(an EM volume compared into pseudoatoms)'
                           'or a SetOfPrincipalComponents.')

        form.addParam('edit', EnumParam, choices=['Slice', 'Reduce', 'Extend', 'Interpolate'],
                    default=NMA_SLICE,
                    label='Type of edit',
                    help='Modes can have the number of nodes decreased using either eigenvector slicing '
                    'or the slower but often more meaningful Hessian reduction method (aka vibrational subsystem '
                    'analysis; Hinsen et al., Chem Phys 2000; Woodcock et al., J Chem Phys 2008) for ProDy vectors. \n'
                    'The number of nodes can be increased by extending (copying) eigenvector values '
                    'from nodes of the same residue or by through-space thin plate splines interpolation')

        form.addParam('newNodes', PointerParam,
                      label='new nodes',
                      pointerClass='AtomStruct',
                      help='Atoms or pseudoatoms to use as new nodes.')   

        form.addParam('norm', BooleanParam, default=True, 
                      condition='edit==%d' % NMA_SLICE,
                      label='Normalise sliced vectors',
                      help='Elect whether to normalise vectors.')
                      
        form.addSection(label='Animation')
        form.addParam('doAnimation', BooleanParam, default=False,
                      label='Make animations for ContinuousFlex viewer')
        animCheck = 'doAnimation == True'
        form.addParam('rmsd', FloatParam, default=5,
                      condition=animCheck,
                      label='RMSD Amplitude (A)',
                      help='Used only for animations of computed normal modes. '
                      'This is the maximal amplitude with which atoms or pseudoatoms are moved '
                      'along normal modes in the animations. \n')
        form.addParam('n_steps', IntParam, default=10,
                      condition=animCheck,
                      expertLevel=LEVEL_ADVANCED,
                      label='Number of frames',
                      help='Number of frames used in each direction of animations.')
        form.addParam('pos', BooleanParam, default=True,
                      condition=animCheck,
                      expertLevel=LEVEL_ADVANCED,
                      label="Include positive direction",
                      help='Elect whether to animate in the positive mode direction.')
        form.addParam('neg', BooleanParam, default=True,
                      condition=animCheck,
                      expertLevel=LEVEL_ADVANCED,
                      label="Include negative direction",
                      help='Elect whether to animate in the negative mode direction.')

    # --------------------------- STEPS functions ------------------------------
    # This is inherited from modes base protocol
    def _insertAllSteps(self):
        modes = prody.parseScipionModes(self.modes.get().getFileName())
        self.nzero = len(np.nonzero(modes.getEigvals() < prody.utilities.ZERO)[0])

        super(ProDyEdit, self)._insertAllSteps(len(self.modes.get()), self.nzero)

    def computeModesStep(self):
        fixVerbositySecondary(self)
        
        self.inputStructure = self.modes.get().getPdb()
        modes = prody.parseScipionModes(self.modes.get().getFileName(),
                                        pdb=self.inputStructure.getFileName())

        oldNodes = prody.parsePDB(self.inputStructure.getFileName(), altloc="all")
        newNodes = prody.parsePDB(self.newNodes.get().getFileName(), altloc="all")

        nodesList = [oldNodes, newNodes]
        numAtomsArr = np.array([nodes.numAtoms() for nodes in nodesList])
        smaller = nodesList[np.argmin(numAtomsArr)]
        bigger = nodesList[np.argmax(numAtomsArr)]

        amap = prody.alignChains(bigger, smaller, match_func=prody.sameChid)[0]
        
        if self.edit == NMA_SLICE:
            self.outModes, self.atoms = prody.sliceModel(modes, bigger, amap, norm=self.norm)

        elif self.edit == NMA_REDUCE:
            modesPath = os.path.dirname(os.path.dirname(
                self.modes.get()._getMapper().selectFirst().getModeFile()))

            fromPrody = len(glob(modesPath+"/*npz"))
            if fromPrody:
                modes = prody.loadModel(glob(modesPath+"/*npz")[0])
                self.outModes, self.atoms = prody.reduceModel(modes, bigger, amap)
                zeros = bool(np.any(modes.getEigvals() < prody.utilities.ZERO))
                self.outModes.calcModes(modes.numModes(), zeros=zeros)
            else:
                logger.warn('ContinuousFlex modes cannot be reduced at this time. Slicing instead')
                self.outModes, self.atoms = prody.sliceModel(modes, bigger, amap, norm=self.norm)

        elif self.edit == NMA_EXTEND:
            self.outModes, self.atoms = prody.extendModel(modes, amap, bigger, norm=True)

        else:
            self.outModes, self.atoms = prody.interpolateModel(modes, amap, bigger, norm=True)

        prody.writePDB(self._getPath('atoms.pdb'), self.atoms)
        prody.writeScipionModes(self._getPath(), self.outModes, write_star=True)

        if isinstance(self.outModes, prody.PCA):
            self.nmdFileName = self._getPath('modes.pca.nmd')
        elif isinstance(self.outModes, prody.GNM):
            self.nmdFileName = self._getPath('modes.gnm.nmd')
        else:
            self.nmdFileName = self._getPath('modes.nmd')
        prody.writeNMD(self.nmdFileName, self.outModes, self.atoms)

        if isinstance(self.outModes, prody.GNM):
            self.gnm = True

        restoreVerbositySecondary(self)

    def createOutputStep(self):
        fnSqlite = self._getPath('modes.sqlite')

        inputClass = type(self.modes.get())
        nmSet = inputClass(filename=fnSqlite)
        nmSet._nmdFileName = String(self.nmdFileName)
        nmSet.setPdb(self.newNodes.get())

        self._defineOutputs(outputModes=nmSet)
        self._defineSourceRelation(self.newNodes, nmSet)
