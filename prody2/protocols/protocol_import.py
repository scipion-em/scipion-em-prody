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
This module will provide ProDy mode import tools.
"""
import os
import numpy as np

from pwem import *
from pwem.objects import AtomStruct, SetOfNormalModes, String
from pwem.protocols import EMProtocol, ProtImportFiles

from pyworkflow.utils import *
import pyworkflow.protocol.params as params

import prody
from prody.utilities import ZERO

NMD = 0
NPZ = 1
SCIPION = 2
GROMACS = 3

class ProDyImportModes(ProtImportFiles):
    """
    This protocol will import ProDy modes to a SetOfNormalModes object
    """
    _label = 'Import modes'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        importChoices = self._getImportChoices()
        filesCondition = self._getFilesCondition()

        form.addSection(label='Import modes')

        form.addParam('importFrom', params.EnumParam,
                      choices=importChoices, default=self._getDefaultChoice(),
                      label='Import from',
                      help='Select the type of import.')

        form.addParam('importType', params.EnumParam, choices=['nmd', 'npz', 'scipion', 'gromacs'],
                      default=NMD,
                      label='Type of modes to import',
                      help='ProDy can support import of modes in various file formats: \n'
                           'native nmd and npz files, a scipion-style directory for a SetOfNormalModes, '
                           'or a gromacs directory with ')

        form.addParam('filesPath', params.PathParam,
                      condition=filesCondition,
                      label="Files directory",
                      help="Directory with the files you want to import.\n\n"
                           "The path can also contain wildcards to select"
                           "from several folders. \n\n"
                           "Examples:\n"
                           "  ~/project/data/day??_files/\n"
                           "Each '?' represents one unknown character\n\n"
                           "  ~/project/data/day*_files/\n"
                           "'*' represents any number of unknown characters\n\n"
                           "  ~/project/data/day##_files/\n"
                           "'##' represents two digits that will be used as "
                           "file ID\n\n"
                           "NOTE: wildcard characters ('*', '?', '#') "
                           "cannot appear in the actual path.)")

        form.addParam('filesPattern', params.StringParam,
                      label='Pattern',
                      condition="importType!=%d" % SCIPION,
                      help="Pattern of the files to be imported.\n\n"
                           "The pattern can contain standard wildcards such as\n"
                           "*, ?, etc, or special ones like ### to mark some\n"
                           "digits in the filename as ID.\n\n"
                           "NOTE: wildcards and special characters "
                           "('*', '?', '#', ':', '%') cannot appear in the "
                           "actual path.\n\n"
                           "For gromacs modes, the first is for values and this is for vectors")

        form.addParam('filesPattern2', params.StringParam,
                      label='Pattern2',
                      condition="importType==%d" % GROMACS,
                      help="Pattern of the files to be imported.\n\n"
                           "The pattern can contain standard wildcards such as\n"
                           "*, ?, etc, or special ones like ### to mark some\n"
                           "digits in the filename as ID.\n\n"
                           "NOTE: wildcards and special characters "
                           "('*', '?', '#', ':', '%') cannot appear in the "
                           "actual path.\n\n"
                           "For gromacs modes, the first is for values and this is for vectors")

        form.addParam('inputStructure', params.PointerParam, label="Input structure",
                      pointerClass='AtomStruct',
                      help='The input structure can be an atomic model '
                           '(true PDB) or a pseudoatomic model\n'
                           '(an EM volume converted into pseudoatoms)')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('importModesStep')
        self._insertFunctionStep('createOutputStep')

    def importModesStep(self):
        files_paths = self.getMatchFiles()
        folder_path = os.path.split(files_paths[0])[0]
        self.pattern1 = os.path.split(files_paths[0])[1]

        pdb_filename = self.inputStructure.get().getFileName()
        
        if self.importType == NMD:
            if not self.pattern1.endswith('.nmd'):
                self.pattern1 += '.nmd'
            self.outModes, _ = prody.parseNMD(os.path.join(folder_path, self.pattern1))

        elif self.importType == NPZ:
            if not self.pattern1.endswith('.npz'):
                self.pattern1 += '.npz'
            self.outModes = prody.loadModel(os.path.join(folder_path, self.pattern1))

        elif self.importType == SCIPION:
            self.outModes = prody.parseScipionModes(folder_path)

        else:
            pattern2 = self.filesPattern2.get()
            self.outModes = prody.parseGromacsModes(folder_path, eigval_fname=pattern1,
                                                    eigvec_fname=pattern2, average_pdb=pdb_filename)

        prody.writeScipionModes(self._getPath(), self.outModes, write_star=True)
        
        if self.importType != NMD:
            atoms = prody.parsePDB(pdb_filename)
            self.nmdFileName = self._getPath('modes.nmd')
            prody.writeNMD(self.nmdFileName, self.outModes, atoms)
        else:
            self.nmdFileName = self.pattern1


    def createOutputStep(self):
        fnSqlite = self._getPath('modes.sqlite')
        nmSet = SetOfNormalModes(filename=fnSqlite)
        nmSet._nmdFileName = String(self.nmdFileName)

        inputPdb = self.inputStructure.get()
        nmSet.setPdb(inputPdb)

        self._defineOutputs(outputModes=nmSet)
        self._defineSourceRelation(self.inputStructure, nmSet)
