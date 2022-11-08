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
from pyworkflow.protocol import params

from os.path import basename, exists, join
import math
from multiprocessing import cpu_count

from pwem import *
import pwem.emlib.metadata as md
from pwem.objects import (SetOfPrincipalComponents, SetOfAtomStructs,
                          String, EMFile)
from pwem.protocols import EMProtocol

from pyworkflow.utils import *
from pyworkflow.utils.path import makePath
from pyworkflow.protocol.params import (PointerParam, IntParam, FloatParam, StringParam,
                                        BooleanParam, EnumParam, LEVEL_ADVANCED)

import prody

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
        form.addParam('inputEnsemble', PointerParam, label="Input ensemble",
                      important=True,
                      pointerClass='SetOfAtomStructs',
                      help='The input ensemble should be a SetOfAtomStructs '
                      'where all structures have the same number of atoms.')

        form.addParam('inputModes', PointerParam, label="Input set of modes",
                      important=True,
                      pointerClass='SetOfNormalModes',
                      help='The input modes can come from Continuous-Flex NMA, ProDy ANM NMA, or ProDy PCA.\n'
                           'The first modes from this set will be used. To use other modes, make a subset.')

        form.addParam('numModes', EnumParam, choices=['1', '2', '3'],
                      label='Number of modes',
                      help='1, 2 or 3 modes can be used for projection')
                 

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('computeStep')
        self._insertFunctionStep('createOutputStep')

    def computeStep(self):
        # configure ProDy to automatically handle secondary structure information and verbosity
        old_secondary = prody.confProDy("auto_secondary")
        old_verbosity = prody.confProDy("verbosity")
        from pyworkflow import Config
        prodyVerbosity =  'none' if not Config.debugOn() else 'debug'
        prody.confProDy(auto_secondary=True, verbosity='{0}'.format(prodyVerbosity))

        ags = prody.parsePDB([tarStructure.getFileName() for tarStructure in self.inputEnsemble.get()])
        ens = prody.buildPDBEnsemble(ags, match_func=prody.sameChainPos, seqid=0., overlap=0.)
        # the ensemble gets built exactly as the input is setup and nothing gets rejected

        modes_path = self.inputModes.get().getFileName()
        modes = prody.parseScipionModes(modes_path)

        self.proj = prody.calcProjection(ens, modes[:self.numModes.get()+1])
        prody.writeArray(self._getExtraPath('projection.txt'), self.proj)

        # configure ProDy to restore secondary structure information and verbosity
        prody.confProDy(auto_secondary=old_secondary, verbosity='{0}'.format(old_verbosity))

    def createOutputStep(self):
        inputEnsemble = self.inputEnsemble.get()
        coeffs = list(prody.parseArray(self._getExtraPath('projection.txt')))
        
        outSetAS = SetOfAtomStructs().create(self._getExtraPath())
        outSetAS.copyItems(inputEnsemble,
                          updateItemCallback=self._updateStructure,
                          itemDataIterator=iter(coeffs))
        
        outputProjection = EMFile(filename=self._getExtraPath('projection.txt'))
        self._defineOutputs(outputProjection=outputProjection,
                            outputStructures=outSetAS)


    # --------------------------- UTILS functions --------------------------------------------
    def _updateStructure(self, structure, value):
        setattr(structure, "_prody_pcaCoefficients", value)
