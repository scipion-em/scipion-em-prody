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
from pyworkflow.protocol.params import PointerParam, EnumParam, BooleanParam

import prody
from prody.utilities import ZERO

NMA_SLICE = 0
NMA_REDUCE = 1
NMA_EXTEND = 2
NMA_INTERP = 3

class ProDyEdit(EMProtocol):
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
        form.addSection(label='ProDy compare')

        form.addParam('edit', EnumParam, choices=['Slice', 'Reduce', 'Extend', 'Interpolate'],
                      default=NMA_SLICE,
                      label='Type of edit',
                      help='Modes can have the number of nodes decreased using either eigenvector slicing '
                      'or Hessian reduction (aka vibrational subsystem analysis; Hinsen et al., Chem Phys 2000; '
                      'Woodcock et al., J Chem Phys 2008). \nThe number of nodes can be increased by extending '
                      '(copying) eigenvector values from nodes of the same residue or by through-space '
                      'interpolation')

        form.addParam('modes', PointerParam, label='Input SetOfNormalModes',
                      pointerClass='SetOfNormalModes',
                      help='The input SetOfNormalModes can be from an atomic model '
                           '(true PDB) or a pseudoatomic model '
                           '(an EM volume compared into pseudoatoms).')

        form.addParam('newNodes', PointerParam,
                      label='new nodes',
                      pointerClass='AtomStruct',
                      help='Atoms or pseudoatoms to use as new nodes.')   

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('editModesStep')
        self._insertFunctionStep('createOutputStep')

    def editModesStep(self):
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
            self.outModes, self.outAtoms = prody.sliceModel(modes, bigger, amap)

        elif self.edit == NMA_REDUCE:
            if from_prody:
                self.outModes, self.outAtoms = prody.reduceModel(modes, bigger, amap)
                zeros = bool(np.any(modes.getEigvals() < ZERO))
                self.outModes.calcModes(zeros=zeros)
            else:
                prody.LOGGER.warn('ContinuousFlex modes cannot be reduced at this time. Slicing instead')
                self.outModes, self.outAtoms = prody.sliceModel(modes, bigger, amap)

        elif self.edit == NMA_EXTEND:
            self.outModes, self.outAtoms = prody.extendModel(modes, amap, bigger, norm=True)

        else:
            self.outModes, self.outAtoms = prody.interpolateModel(modes, amap, bigger, norm=True)

        prody.writePDB(self._getPath('atoms.pdb'), self.outAtoms)
        prody.writeScipionModes(self._getPath(), self.outModes, write_star=True)
        prody.writeNMD(self._getPath('modes.nmd'), self.outModes, self.outAtoms)

    def createOutputStep(self):
        fnSqlite = self._getPath('modes.sqlite')
        nmSet = SetOfNormalModes(filename=fnSqlite)
        nmSet._nmdFileName = String(self._getPath('modes.nmd'))

        # outputPdb = AtomStruct()
        # outputPdb.setFileName(self._getPath('atoms.pdb'))
        # nmSet.setPdb(outputPdb.get())
        nmSet.setPdb(self.newNodes.get())

        self._defineOutputs(outputModes=nmSet)#, outputStructure=outputPdb)
        # self._defineSourceRelation(outputPdb, nmSet)
        self._defineSourceRelation(self.newNodes, nmSet)
