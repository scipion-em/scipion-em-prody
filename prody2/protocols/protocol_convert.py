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
from pyworkflow.protocol import Protocol, params, Integer

import os
from os.path import basename, exists, join
import math
import numpy as np

from pwem import *
from pwem.emlib import (MetaData, MDL_X, MDL_COUNT, MDL_NMA_MODEFILE, MDL_ORDER,
                        MDL_ENABLED, MDL_NMA_COLLECTIVITY, MDL_NMA_SCORE)
from pwem.objects import SetOfNormalModes, String
from pwem.protocols import EMProtocol

from pyworkflow.utils import *
from pyworkflow.utils.path import copyFile, createLink, makePath, cleanPath, moveFile
from pyworkflow.protocol.params import (PointerParam, IntParam, FloatParam, StringParam,
                                        BooleanParam, LEVEL_ADVANCED)

import prody

class ProDyWriteNMD(EMProtocol):
    """
    This protocol will convert a SetOfNormalModes from continuousflex NMA to ProDy NMD files for NMWiz
    """
    _label = 'ProDy writeNMD'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label='ProDy writeNMD')

        form.addParam('inputNMSet', PointerParam, label="Input SetOfNormalModes",
                      important=True,
                      pointerClass='SetOfNormalModes',
                      help='The input SetOfNormalModes can be from an atomic model '
                           '(true PDB) or a pseudoatomic model\n'
                           '(an EM volume converted into pseudoatoms)')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('convertModesStep')
        self._insertFunctionStep('createOutputStep')

    def convertModesStep(self):
        modes_path = os.path.split(self.inputNMSet.get().getFileName())[0]
        modes = prody.parseScipionModes(modes_path)

        pdb = glob(modes_path+"/*atoms.pdb")
        atoms = prody.parsePDB(pdb, altloc='all')

        prody.writeNMD(self._getPath('modes.nmd'), modes, atoms)

    def createOutputStep(self):
        outputNMSet = self._createSetOfNormalModes()
        _ = [outputNMSet.append(mode.clone()) for mode in self.inputNMSet.get().iterItems()]
        outputNMSet._nmdFileName = String(self._getPath('modes.nmd'))
        self._defineOutputs(outputModes=outputNMSet)
