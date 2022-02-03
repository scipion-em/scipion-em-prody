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
This module will provide ProDy normal mode analysis using the anisotropic network model (ANM).
"""
from pyworkflow.protocol import params

from pwem import *
from pwem.objects import AtomStruct, String
from pwem.protocols import EMProtocol

from pyworkflow.utils import *
from pyworkflow.protocol.params import PointerParam, StringParam

import prody

class ProDySelect(EMProtocol):
    """
    This protocol will perform normal mode analysis using the anisotropic network model (ANM)
    """
    _label = 'Atom Selection'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label='ProDy ANM NMA')

        form.addParam('inputStructure', PointerParam, label="Input structure",
                      important=True,
                      pointerClass='AtomStruct',
                      help='The input structure can be an atomic model '
                           '(true PDB) or a pseudoatomic model\n'
                           '(an EM volume converted into pseudoatoms)')

        form.addParam('selection', StringParam, default="",
                      label="selection string",
                      help='This determines which atoms are used in the calculations. '
                           'If left empty, the default behaviour is to use '
                           '"protein and name CA or nucleic and name P C4\' C2" '
                           'for atomic structures and "all" for pseudoatoms.')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        inputFn = self.inputStructure.get().getFileName()
        self._insertFunctionStep('selectionStep', inputFn)
        self._insertFunctionStep('createOutputStep')

    def selectionStep(self, inputFn):
        if self.selection == "":
            self.selection = "protein and name CA or nucleic and name P C4' C2"

        ag = prody.parsePDB(inputFn, alt='all')
        selection = ag.select(str(self.selection))

        self.pdbFileName = self._getPath('atoms.pdb')
        prody.writePDB(self.pdbFileName, selection)

    def createOutputStep(self):
        outputPdb = AtomStruct()
        outputPdb.setFileName(self.pdbFileName)
        self._defineOutputs(outputStructure=outputPdb)   

