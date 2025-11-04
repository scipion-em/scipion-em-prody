# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     James Krieger (jamesmkrieger@gmail.com)
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
import os
from prody2 import Plugin

from pwem.objects import AtomStruct, EMFile, String
from pwem.protocols import EMProtocol

from pyworkflow.utils import glob
from pyworkflow.protocol.params import PointerParam, EnumParam, BooleanParam

NMA_METRIC_OVERLAP = 0
NMA_METRIC_COV_OVERLAP = 1
NMA_METRIC_RWSIP = 2

class ProDyCompare(EMProtocol):
    """
    This protocol will compare two SetOfNormalModes objects
    """
    _label = 'Compare modes'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label='ProDy compare')

        form.addParam('modes1', PointerParam, label="Input modes set 1",
                      important=True,
                      pointerClass='SetOfNormalModes',
                      help='The input modes can be a SetOfNormalModes '
                           'from an atomic model (true PDB) or a pseudoatomic model '
                           '(an EM volume compared into pseudoatoms) '
                           'or a SetOfPrincipalComponents.\n'
                           'The two sets should have the same number of nodes '
                           'unless one of them has exactly 1 mode in it.')

        form.addParam('modes2', PointerParam, label="Input modes set 2",
                      important=True,
                      pointerClass='SetOfNormalModes',
                      help='The input modes can be a SetOfNormalModes '
                           'from an atomic model (true PDB) or a pseudoatomic model '
                           '(an EM volume compared into pseudoatoms) '
                           'or a SetOfPrincipalComponents.\n'
                           'The two sets should have the same number of nodes '
                           'unless one of them has exactly 1 mode in it.')

        form.addParam('metric', EnumParam, choices=['Overlap', 'Covariance Overlap', 'RWSIP'],
                      default=NMA_METRIC_OVERLAP,
                      label='Comparison metric',
                      help='Modes can be compared pairwise using correlation cosine overlaps (inner products), '
                      'or in sets using either covariance overlap (Hess, Phys Rev E 2002; aka spectral overlap) '
                      'or the root weighted square inner product (RWSIP; Carnevale et al., J Phys Condens Matter 2007). \n'
                      'Covariance overlaps and RWSIPs are calculated over growing mode sets.\n'
                      'Zero eigenvalue modes are excluded from all calculations.')

        form.addParam('diag', BooleanParam, default=False, 
                      condition='metric==%d' % NMA_METRIC_OVERLAP,
                      label='Calculate diagonal values only',
                      help='Elect whether to calculate diagonal values only.')      

        form.addParam('match', BooleanParam, default=False, 
                      condition='metric!=%d' % NMA_METRIC_RWSIP,
                      label='Match modes',
                      help='Elect whether to match modes.')     

        form.addParam('norm', BooleanParam, default=True, 
                      condition='metric==%d' % NMA_METRIC_OVERLAP,
                      label='Normalise overlaps',
                      help='Elect whether to normalise vectors for overlaps or calculate raw dot products.')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('compareModesStep')
        self._insertFunctionStep('createOutputStep')

    def compareModesStep(self):
        modesPath1 = os.path.dirname(os.path.dirname(
            self.modes1.get()._getMapper().selectFirst().getModeFile()))

        pdb1 = glob(modesPath1+"/*atoms.pdb")
        if len(pdb1) != 0:
            pdb1 = pdb1[0]
        elif self.modes1.get().getPdb().getFileName() is not None:
            pdb1 = self.modes1.get().getPdb().getFileName()
        else:
            pdb1 = "None"

        modesFn1 = self.modes1.get().getFileName()

        modesPath2 = os.path.dirname(os.path.dirname(
            self.modes2.get()._getMapper().selectFirst().getModeFile()))
            
        pdb2 = glob(modesPath2+"/*atoms.pdb")
        if len(pdb2) != 0:
            pdb2 = pdb2[0]
        elif self.modes2.get().getPdb().getFileName() is not None:
            pdb2 = self.modes2.get().getPdb().getFileName()
        else:
            pdb2 = "None"

        modesFn2 = self.modes2.get().getFileName()

        args = '--inputPdbFns "{0}" --inputModesFns "{1}" --folder {2} '.format(
            ' '.join([pdb1, pdb2]),
            ' '.join([modesFn1, modesFn2]),
            self._getPath()
        )

        args += '--metric {0} '.format(self.metric.get())

        if self.match:
            args += '--match True '

        if self.diag:
            args += '--diag True '

        if self.norm:
            args += '--norm True '

        self.runJob(Plugin.getProgram('compare_modes.py', script=True), args)

    def createOutputStep(self):
        outputMatrix = EMFile(filename=self._getPath('matrix.txt'))

        if self.match:
            fnSqlite = self._getPath('modes.sqlite')
            inputClass = type(self.modes1.get())
            nmSet = inputClass(filename=fnSqlite)
            nmSet._nmdFileName = String(self.getNmdFileName)

            outputPdb = AtomStruct()
            outputPdb.setFileName(self._getPath('atoms.pdb'))
            nmSet.setPdb(outputPdb.get())

            outputMatch = EMFile(filename=self.getMatchIndsFn)
            nmSet._indsFileName = String(self.getMatchIndsFn)

            self._defineOutputs(matrixFile=outputMatrix,
                                matchFile=outputMatch,
                                outputModes=nmSet)
        else:
            self._defineOutputs(matrixFile=outputMatrix)

    def getMatchIndsFn(self):
        return self._getPath('matchInds.txt')

    def getNmdFileName(self):
        return glob(self._getPath() + '/matched_modes.*.nmd')
