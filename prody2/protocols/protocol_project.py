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
This module will provide ProDy projection of structural ensembles on principal component or normal modes
"""
import numpy as np

from pwem.objects import SetOfAtomStructs
from pwem.protocols import EMProtocol

import pyworkflow.object as pwobj
from pyworkflow.protocol.params import (PointerParam, EnumParam, BooleanParam,
                                        MultiPointerParam, NumericRangeParam)
from pyworkflow.utils import getListFromRangeString

import prody
from prody2.constants import PROJ_COEFFS
from prody2 import fixVerbositySecondary, restoreVerbositySecondary

ONE = 0
TWO = 1
THREE = 2

class ProDyProject(EMProtocol):
    """
    This module will provide ProDy projection of structural ensembles on principal component or normal modes
    """
    _label = 'Projection'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        form.addSection(label='ProDy Projection')
        form.addParam('inputEnsemble', MultiPointerParam, label="Input ensemble(s)",
                      important=True,
                      pointerClass='SetOfAtomStructs,ProDyNpzEnsemble',
                      help='The input ensemble should be SetOfAtomStructs or ProDyNpzEnsemble '
                      'objects where all structures have the same number of atoms.')

        form.addParam('inputModes', PointerParam, label="Input set of modes",
                      important=True,
                      pointerClass='SetOfNormalModes,SetOfPrincipalComponents',
                      help='The input modes can come from Continuous-Flex NMA, ProDy NMA, or ProDy PCA.\n'
                           'The first modes from this set will be used. To use other modes, make a subset.')

        form.addParam('modeList', NumericRangeParam,
                      label="Modes selection", allowsNull=True, default="",
                      help='Select the normal modes that will be used for analysis.\n'
                           'If you leave this field empty, all the computed modes will be selected from.\n'
                           'If you only enter one number, all the computed modes will be selected from starting with that one.\n'
                           'You have several ways to specify the modes.\n'
                           '   Examples:\n'
                           ' "7,8-10" -> [7,8,9,10]\n'
                           ' "8, 10, 12" -> [8,10,12]\n'
                           ' "8 9, 10-12" -> [8,9,10,11,12])\n')
        
        form.addParam('numModes', EnumParam, choices=['1', '2', '3'],
                      label='Number of modes', default=TWO,
                      help='1, 2 or 3 modes can be used for projection')

        form.addParam('norm', BooleanParam, label="Normalize?", default=False,
                      help='Select whether to normalise projections.')

        form.addParam('rmsd', BooleanParam, label="RMSD scale?", default=True,
                      help='Select whether to scale projections to RMSDs.')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('computeStep')
        self._insertFunctionStep('createOutputStep')

    def computeStep(self):
        fixVerbositySecondary(self)

        inputModes = self.inputModes.get()
        modesPath = inputModes.getFileName()
        modes = prody.parseScipionModes(modesPath)

        if not self.modeList.empty():
            modeSelection = list(np.array(getListFromRangeString(self.modeList.get())) - 1)
            if len(modeSelection) == 1:
                modes = modes[modeSelection[0]:]
            else:
                modes = modes[modeSelection]

        prody.writeScipionModes(self._getPath(), modes, write_star=True)
        fnSqlite = self._getPath('modes.sqlite')
        inputClass = type(inputModes)
        self.outputModes = inputClass(filename=fnSqlite)

        self.proj = []
        for i, inputEnsemble in enumerate(self.inputEnsemble):
            ensGot = inputEnsemble.get()
            idSet = ensGot.getIdSet()
            if isinstance(ensGot, SetOfAtomStructs):
                ags = prody.parsePDB([tarStructure.getFileName() for tarStructure in ensGot])
                ens = prody.buildPDBEnsemble(ags, match_func=prody.sameChainPos, seqid=0., overlap=0., superpose=False, mapping=None)
                # the ensemble gets built exactly as the input is setup and nothing gets rejected
            else:
                ens = ensGot.loadEnsemble()

            projection = prody.calcProjection(ens, modes[:self.numModes.get()+1], rmsd=self.rmsd.get(),
                                              norm=self.norm.get())
            projDict = dict()
            for j, idx in enumerate(idSet):
                proj = projection[j]
                if isinstance(proj, float):
                    proj = [proj]

                projDict[idx] = proj
            self.proj.append(projDict)
            prody.writeArray(self._getPath('projection_{0}.csv'.format(i+1)), projection, 
                             format='%8.5f', delimiter=',')

        restoreVerbositySecondary(self)

    def createOutputStep(self):

        args = {}
        for self.ensId, inputEnsemble in enumerate(self.inputEnsemble): 
            ensGot = inputEnsemble.get()

            suffix = str(self.ensId+1)

            inputClass = type(ensGot)
            outSet = inputClass().create(self._getExtraPath(), suffix=suffix)
            outSet.copyItems(ensGot, updateItemCallback=self._setCoeffs)

            name = "outputEns" + suffix
            
            args[name] = outSet

        args["outputModes"] = self.outputModes

        self._defineOutputs(**args)

    # --------------------------- UTILS functions --------------------------------------------
    def _setCoeffs(self, item, row=None):
        # We provide data directly so don't need a row
        vector = pwobj.CsvList()
        vector._convertValue(["{:18.15f}".format(x) for x in (self.proj[self.ensId][item.getObjId()])])
        setattr(item, PROJ_COEFFS, vector)

    def _summary(self):
        if not hasattr(self, 'outputEns1'):
            summ = ['Projection not ready yet']
        else:
            summ = ['Projected structures onto *{0}* components'.format(self.numModes.get()+1)]
        return summ
        