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

from pwem.objects import SetOfPrincipalComponents, SetOfAtomStructs, EMFile
from pwem.protocols import EMProtocol

import pyworkflow.object as pwobj
from pyworkflow.protocol.params import PointerParam, EnumParam, MultiPointerParam

import prody
from prody2.constants import PROJ_COEFFS
from prody2.objects import ProDyNpzEnsemble

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

        form.addParam('numModes', EnumParam, choices=['1', '2', '3'],
                      label='Number of modes', default=TWO,
                      help='1, 2 or 3 modes can be used for projection')


    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('computeStep')
        self._insertFunctionStep('createOutputStep')

    def computeStep(self):
        # configure ProDy to automatically handle secondary structure information and verbosity
        oldSecondary = prody.confProDy("auto_secondary")
        oldVerbosity = prody.confProDy("verbosity")
        from pyworkflow import Config
        prodyVerbosity =  'none' if not Config.debugOn() else 'debug'
        prody.confProDy(auto_secondary=True, verbosity='{0}'.format(prodyVerbosity))

        modesPath = self.inputModes.get().getFileName()
        modes = prody.parseScipionModes(modesPath)

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

            projection = prody.calcProjection(ens, modes[:self.numModes.get()+1])
            projDict = dict()
            for i, idx in enumerate(idSet):
                projDict[idx] = projection[i]
            self.proj.append(projDict)
            prody.writeArray(self._getPath('projection_{0}.csv'.format(i+1)), projection, 
                             format='%8.5f', delimiter=',')

        # configure ProDy to restore secondary structure information and verbosity
        prody.confProDy(auto_secondary=oldSecondary, verbosity='{0}'.format(oldVerbosity))

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

        self._defineOutputs(**args)

    # --------------------------- UTILS functions --------------------------------------------
    def _setCoeffs(self, item, row=None):
        # We provide data directly so don't need a row
        vector = pwobj.CsvList()
        vector._convertValue(["{:18.15f}".format(x) for x in (self.proj[self.ensId][item.getObjId()])])
        setattr(item, PROJ_COEFFS, vector)

    def _summary(self):
        if not hasattr(self, 'outputStructures'):
            summ = ['Projection not ready yet']
        else:
            summ = ['Projected *{0}* structures onto *{1}* components'.format(
                   len(self.outputStructures), self.numModes.get()+1)]
        return summ
        