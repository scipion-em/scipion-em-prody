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
import os
import numpy as np

from pwem import *
from pwem.objects import SetOfNormalModes, AtomStruct, EMFile, String
from pwem.protocols import EMProtocol

from pyworkflow.utils import *
from pyworkflow.protocol.params import PointerParam, EnumParam, BooleanParam

import prody

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

        form.addParam('modes1', PointerParam, label="Input SetOfNormalModes 1",
                      important=True,
                      pointerClass='SetOfNormalModes',
                      help='The input SetOfNormalModes can be from an atomic model '
                           '(true PDB) or a pseudoatomic model '
                           '(an EM volume compared into pseudoatoms).\n'
                           'The two sets should have the same number of nodes '
                           'unless one of them has exactly 1 mode in it.')

        form.addParam('modes2', PointerParam, label="Input SetOfNormalModes 2",
                      important=True,
                      pointerClass='SetOfNormalModes',
                      help='The input SetOfNormalModes can be from an atomic model '
                           '(true PDB) or a pseudoatomic model '
                           '(an EM volume compared into pseudoatoms).\n'
                           'The two sets should have the same number of nodes '
                           'unless one of them has exactly 1 mode in it.')

        form.addParam('metric', EnumParam, choices=['Overlap', 'Covariance Overlap', 'RWSIP'],
                      default=NMA_METRIC_OVERLAP,
                      label='Comparison metric',
                      help='Modes can be compared pairwise using correlation cosine overlaps (inner products), '
                      'or in sets using either covariance overlap (Hess, Phys Rev E 2002; aka spectral overlap) '
                      'or the root weighted square inner product (RWSIP; Carnevale et al., J Phys Condens Matter 2007). \n'
                      'Covariance overlaps and RWSIPs are calculated over growing mode sets.\n'
                      'The six zero modes are excluded from covariance overlap and RWSIP calculations.')

        form.addParam('diag', BooleanParam, default=False, 
                      condition='metric==%d' % NMA_METRIC_OVERLAP,
                      label='Calculate diagonal values only',
                      help='Elect whether to calculate diagonal values only.')      

        form.addParam('match', BooleanParam, default=False, 
                      condition='metric!=%d' % NMA_METRIC_RWSIP,
                      label='Match modes',
                      help='Elect whether to match modes.')     

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('compareModesStep')
        self._insertFunctionStep('createOutputStep')

    def compareModesStep(self):

        modes1_path = os.path.dirname(os.path.dirname(self.modes1.get()[1].getModeFile()))
        modes1 = prody.parseScipionModes(modes1_path)

        modes2_path = os.path.dirname(os.path.dirname(self.modes2.get()[1].getModeFile()))
        modes2 = prody.parseScipionModes(modes2_path)

        n_modes = np.max([modes1.numModes(), modes2.numModes()])
        min_n_modes = np.min([modes1.numModes(), modes2.numModes()])
        if modes1.numModes() != modes2.numModes() and min_n_modes != 1:
            raise ValueError('The two sets should have the same number of nodes '
                           'unless one of them has exactly 1 mode in it.')

        if min_n_modes != 1:
            mode_ens = prody.ModeEnsemble()
            mode_ens.addModeSet(modes1)
            mode_ens.addModeSet(modes2)
            if self.match:
                mode_ens.match()

                match_inds = prody.matchModes(modes1, modes2, index=True)

                prody.writeArray(self._getExtraPath('match_inds.txt'),
                                 np.array(match_inds, dtype=int)[1]+1,
                                 format='%3d')

                atoms = prody.parsePDB(self.modes1.get().getPdb().getFileName())
                prody.writeNMD(self._getExtraPath('matched_modes.nmd'), mode_ens[1], atoms)
                prody.writeScipionModes(self._getPath(), mode_ens[1], write_star=True)
        else:
            mode_ens = [modes1, modes2]
        
        if self.metric == NMA_METRIC_OVERLAP:
            self.matrix = prody.calcOverlap(mode_ens[0], mode_ens[1], diag=self.diag)
            if self.matrix.ndim == 1:
                self.matrix.reshape(-1, 1)

        else:
            self.matrix = np.empty((n_modes-6, 1))
            for i in range(6, n_modes):
                if self.metric == NMA_METRIC_COV_OVERLAP:
                    self.matrix[i-6, 0] = prody.calcEnsembleSpectralOverlaps(mode_ens[:, 6:i+1])[0, 1]
                else:
                    self.matrix[i-6, 0] = prody.calcRWSIP(mode_ens[0, 6:i+1], mode_ens[1, 6:i+1])

        prody.writeArray(self._getExtraPath('matrix.txt'), self.matrix)

    def createOutputStep(self):
        outputMatrix = EMFile(filename=self._getExtraPath('matrix.txt'))

        if self.match:
            fnSqlite = self._getPath('modes.sqlite')
            nmSet = SetOfNormalModes(filename=fnSqlite)
            nmSet._nmdFileName = String(self._getPath('modes.nmd'))

            outputPdb = AtomStruct()
            outputPdb.setFileName(self._getPath('atoms.pdb'))
            nmSet.setPdb(outputPdb.get())

            outputMatch = EMFile(filename=self._getExtraPath('match_inds.txt'))
            nmSet._indsFileName = String(self._getExtraPath('match_inds.txt'))

            self._defineOutputs(matrixFile=outputMatrix, matchFile=outputMatch, outputModes=nmSet)
        else:
            self._defineOutputs(matrixFile=outputMatrix)
