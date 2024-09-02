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
This module will provide ProDy mode algebra tools (linear combinations).
"""
import os
import numpy as np

from pwem.objects import SetOfNormalModes, String, Integer, CsvList

from pyworkflow.utils import glob, logger
from pyworkflow.protocol.params import (PointerParam, EnumParam, IntParam,
                                        StringParam, LEVEL_ADVANCED)

import prody
from prody2.protocols.protocol_modes_base import ProDyModesBase
from prody2 import Plugin

COEFF_POINTER = 0
COEFF_STRING = 1

class ProDyAlgebra(ProDyModesBase):
    """
    This protocol will add together components from a SetOfNormalModes object 
    with coefficients based on overlaps or user input
    """
    _label = 'Vector alegebra'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form, besidesAnimation=False):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label='ProDy linear algebra')
        form.addParam('modes', PointerParam, label='Input set of modes',
                      pointerClass='SetOfNormalModes',
                      help='The input modes can be a SetOfNormalModes '
                           'from an atomic model (true PDB) or a pseudoatomic model '
                           '(an EM volume compared into pseudoatoms)'
                           'or a SetOfPrincipalComponents.')

        form.addParam('coeffSource', EnumParam, choices=['Pointer', 'String'],
                    default=COEFF_STRING,
                    label='Type of edit',
                    help='Modes will be added together with these coefficients. If there are more modes '
                         'than coefficients then the remaining modes will be ignored.')

        form.addParam('coeffPointer', PointerParam,
                      label='Coefficients',
                      condition='coeffSource==%d' % COEFF_POINTER,
                      pointerClass='EMFile',
                      help='Atoms or pseudoatoms to use as new nodes.')   

        form.addParam('coeffString', StringParam, default='', 
                      condition='coeffSource==%d' % COEFF_STRING,
                      label='Coefficients',
                      help='Modes will be added together with these coefficients. If there are more modes '
                           'than coefficients then the remaining modes will be ignored.')

        form.addParam('numCoeffs', IntParam, default='-1',
                      label='Number of components',
                      help='This number of modes will be added together with coefficients. '
                           'The remaining modes will be ignored.')

        ProDyModesBase._defineParams(self, form, besidesAnimation=besidesAnimation)

    # --------------------------- STEPS functions ------------------------------
    # This is inherited from modes base protocol
    def _insertAllSteps(self, n=1, nzeros=0):
        self.nzero = nzeros
        super(ProDyAlgebra, self)._insertAllSteps(n=n, nzeros=nzeros)

    def computeModesStep(self):
        if self.coeffSource == COEFF_STRING:
            sep = ''
            coeffsString = self.coeffString.get()
            if ',' in coeffsString:
                sep += ','
            if ' ' in coeffsString:
                sep += ' '

            coeffs = np.array(coeffsString.split(sep), dtype=float)
            coeffsFn = self._getExtraPath('coeffs.txt')
            np.savetxt(coeffsFn, coeffs)
        else:
            coeffsFn = self.coeffPointer.get().getFileName()
            coeffs = np.loadtxt(coeffsFn)

        self.modesFn = self.modes.get().getFileName()
        self.atomsFn = self.modes.get().getPdb().getFileName()

        numCoeffs = min(self.numCoeffs.get(), len(coeffs))
        if numCoeffs == -1:
            numCoeffs = len(coeffs)

        args = '--modesFn {0} --atomsFn {1} --coeffsFn {2} --numCoeffs {3} ' \
            '--folder {4} --nmdFileName {5} --npzFileName {6}'.format(
                self.modesFn, self.atomsFn, coeffsFn, numCoeffs,
                self._getPath(), self.getNmdFileName(), self.getNpzFileName())

        self.runJob(Plugin.getProgram('algebra.py', script=True), args)

        self.outModesFn = self._getPath("modes.sqlite")
        self.coeffs = CsvList()
        self.coeffs._convertValue(["{:18.15f}".format(x) for x in coeffs[:numCoeffs]])

    def createOutputStep(self):
        fnSqlite = self._getPath('modes.sqlite')

        nmSet = SetOfNormalModes(filename=fnSqlite)
        nmSet._nmdFileName = String(self.getNmdFileName())

        pdb = self.modes.get().getPdb()
        nmSet.setPdb(pdb)
        os.symlink(os.path.abspath(pdb.getFileName()), 
                   self._getPath('atoms.pdb'))

        self._defineOutputs(outputModes=nmSet, coeffs=self.coeffs)
        self._defineSourceRelation(pdb, nmSet)

    def _summary(self):
        if not hasattr(self, 'outputModes'):
            summ = ['Output modes not ready yet']
        elif len(self.coeffs) == 1:
            summ = ['*1* mode scaled with coefficient *{0}*'.format(self.coeffs)]
        else:
            summ = ['*{0}* modes added with coefficients *{1}*'.format(
                    len(self.coeffs), self.coeffs)]
        return summ

    def getNmdFileName(self):
        return self._getPath('modes.nmd')
    
    def getNpzFileName(self):
        return self._getPath('modes.nma.npz')
