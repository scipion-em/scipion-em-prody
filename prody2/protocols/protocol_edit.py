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

from pwem import *
from pwem.objects import AtomStruct, SetOfNormalModes, String
from pwem.protocols import EMProtocol

from pyworkflow.utils import *
from pyworkflow.protocol.params import (PointerParam, EnumParam, BooleanParam,
                                        FloatParam, IntParam, LEVEL_ADVANCED)

import prody
from prody.utilities import ZERO

from prody2.protocols.protocol_modes_base import ProDyModesBase

NMA_SLICE = 0
NMA_REDUCE = 1
NMA_EXTEND = 2

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

        form.addParam('edit', EnumParam, choices=['Slice', 'Reduce', 'Extend'],
                      default=NMA_SLICE,
                      label='Type of edit',
                      help='Modes can have the number of nodes decreased using either eigenvector slicing '
                      'or the slower but often more meaningful Hessian reduction method (aka vibrational subsystem '
                      'analysis; Hinsen et al., Chem Phys 2000; Woodcock et al., J Chem Phys 2008) for ProDy vectors. \n'
                      'The number of nodes can be increased by extending (copying) eigenvector values '
                      'from nodes of the same residue')

        form.addParam('modes', PointerParam, label='Input SetOfNormalModes',
                      pointerClass='SetOfNormalModes',
                      help='The input SetOfNormalModes can be from an atomic model '
                           '(true PDB) or a pseudoatomic model '
                           '(an EM volume compared into pseudoatoms).')

        form.addParam('newNodes', PointerParam,
                      label='new nodes',
                      pointerClass='AtomStruct',
                      help='Atoms or pseudoatoms to use as new nodes.')   

        form.addParam('norm', BooleanParam, default=True, 
                      condition='edit==%d' % NMA_SLICE,
                      label='Normalise sliced vectors',
                      help='Elect whether to normalise vectors.')
                      
        form.addSection(label='Animation')        
        form.addParam('rmsd', FloatParam, default=5,
                      label='RMSD Amplitude (A)',
                      help='Used only for animations of computed normal modes. '
                      'This is the maximal amplitude with which atoms or pseudoatoms are moved '
                      'along normal modes in the animations. \n')
        form.addParam('n_steps', IntParam, default=10,
                      expertLevel=LEVEL_ADVANCED,
                      label='Number of frames',
                      help='Number of frames used in each direction of animations.')
        form.addParam('pos', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Include positive direction",
                      help='Elect whether to animate in the positive mode direction.')
        form.addParam('neg', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Include negative direction",
                      help='Elect whether to animate in the negative mode direction.')

    # --------------------------- STEPS functions ------------------------------
    # This is inherited from NMA base protocol
    def _insertAllSteps(self):
        super(ProDyEdit, self)._insertAllSteps(len(self.modes.get()))
        
    #     # Insert processing steps
    #     self._insertFunctionStep('editModesStep')
    #     self._insertFunctionStep('createOutputStep')

    def computeModesStep(self):
        modes_path = os.path.dirname(os.path.dirname(self.modes.get()[1].getModeFile()))
        
        from_prody = len(glob(modes_path+"/*npz"))
        if from_prody:
            modes = prody.loadModel(glob(modes_path+"/*npz")[0])
        else:
            modes = prody.parseScipionModes(modes_path)

        self.inputStructure = self.modes.get().getPdb()
        structureEM = self.inputStructure.getPseudoAtoms()

        old_nodes = prody.parsePDB(self.inputStructure.getFileName(), altloc="all")
        new_nodes = prody.parsePDB(self.newNodes.get().getFileName(), altloc="all")

        nodes_list = [old_nodes, new_nodes]
        n_atoms_arr = np.array([nodes.numAtoms() for nodes in nodes_list])
        smaller = nodes_list[np.argmin(n_atoms_arr)]
        bigger = nodes_list[np.argmax(n_atoms_arr)]

        amap = prody.alignChains(bigger, smaller, match_func=prody.sameChid, pwalign=False)[0]
        
        if self.edit == NMA_SLICE:
            self.outModes, self.outAtoms = prody.sliceModel(modes, bigger, amap, norm=self.norm)

        elif self.edit == NMA_REDUCE:
            if from_prody:
                self.outModes, self.outAtoms = prody.reduceModel(modes, bigger, amap)
                zeros = bool(np.any(modes.getEigvals() < ZERO))
                self.outModes.calcModes(zeros=zeros)
            else:
                prody.LOGGER.warn('ContinuousFlex modes cannot be reduced at this time. Slicing instead')
                self.outModes, self.outAtoms = prody.sliceModel(modes, bigger, amap, norm=self.norm)

        elif self.edit == NMA_EXTEND:
            self.outModes, self.atoms = prody.extendModel(modes, amap, bigger, norm=True)

        prody.writePDB(self._getPath('atoms.pdb'), self.outAtoms)
        prody.writeScipionModes(self._getPath(), self.outModes, write_star=True)
        prody.writeNMD(self._getPath('modes.nmd'), self.outModes, self.atoms)

    def createOutputStep(self):
        fnSqlite = self._getPath('modes.sqlite')
        nmSet = SetOfNormalModes(filename=fnSqlite)
        nmSet._nmdFileName = String(self._getPath('modes.nmd'))
        nmSet.setPdb(self.newNodes.get())

        self._defineOutputs(outputModes=nmSet)
        self._defineSourceRelation(self.newNodes, nmSet)
