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
This module will provide ProDy atom tools including selection and superposition.
"""
from pyworkflow.protocol import params

from pwem import *
from pwem.objects import AtomStruct, Transform, String
from pwem.protocols import EMProtocol

from pyworkflow.utils import *
from pyworkflow.protocol.params import PointerParam, StringParam

import prody

class ProDySelect(EMProtocol):
    """
    This protocol will perform atom selection
    """
    _label = 'Atom Selection'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label='ProDy Select')

        form.addParam('inputStructure', PointerParam, label="Input structure",
                      important=True,
                      pointerClass='AtomStruct',
                      help='The input structure can be an atomic model '
                           '(true PDB) or a pseudoatomic model\n'
                           '(an EM volume converted into pseudoatoms)')

        form.addParam('selection', StringParam, default="protein and name CA or nucleic and name P C4' C2",
                      label="selection string",
                      help='This determines which atoms are selected. '
                           'There is a rich selection engine with similarities to VMD '
                           '(see http://prody.csb.pitt.edu/tutorials/prody_tutorial/selection.html)')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        inputFn = self.inputStructure.get().getFileName()
        self._insertFunctionStep('selectionStep', inputFn)
        self._insertFunctionStep('createOutputStep')

    def selectionStep(self, inputFn):
        ag = prody.parsePDB(inputFn, alt='all')
        selection = ag.select(str(self.selection))

        self.pdbFileName = self._getPath('atoms.pdb')
        prody.writePDB(self.pdbFileName, selection)

    def createOutputStep(self):
        outputPdb = AtomStruct()
        outputPdb.setFileName(self.pdbFileName)
        self._defineOutputs(outputStructure=outputPdb)   


class ProDySuperpose(EMProtocol):
    """
    This protocol will perform atomic structure superposition
    """
    _label = 'Atom Superposition'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label='ProDy Superpose')

        form.addParam('mobStructure', PointerParam, label="Mobile structure",
                      important=True,
                      pointerClass='AtomStruct',
                      help='The structure to be moved can be an atomic model '
                           '(true PDB) or a pseudoatomic model\n'
                           '(an EM volume converted into pseudoatoms).'
                           'The two structures should have the same number of nodes.')

        form.addParam('tarStructure', PointerParam, label="Target structure",
                      important=True,
                      pointerClass='AtomStruct',
                      help='The target structure can be an atomic model '
                           '(true PDB) or a pseudoatomic model\n'
                           '(an EM volume converted into pseudoatoms)'
                           'The two structures should have the same number of nodes.')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        mobFn = self.mobStructure.get().getFileName()
        tarFn = self.tarStructure.get().getFileName()
        self._insertFunctionStep('superposeStep', mobFn, tarFn)
        self._insertFunctionStep('createOutputStep')

    def superposeStep(self, mobFn, tarFn):

        mob = prody.parsePDB(mobFn, alt='all')
        tar = prody.parsePDB(tarFn, alt='all')

        alg, self.T = prody.superpose(mob, tar)

        self.pdbFileName = self._getPath('atoms.pdb')
        prody.writePDB(self.pdbFileName, alg)

        self.matrixFileName = self._getPath('transformation.txt')
        prody.writeArray(self.matrixFileName, self.T.getMatrix())

    def createOutputStep(self):
        outputPdb = AtomStruct()
        outputPdb.setFileName(self.pdbFileName)

        outputTrans = Transform()
        outputTrans.setMatrix(self.T.getMatrix())

        self._defineOutputs(outputStructure=outputPdb,
                            outputTransformation=outputTrans)

