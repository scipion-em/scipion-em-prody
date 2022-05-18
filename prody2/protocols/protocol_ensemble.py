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
This module will provide ProDy ensemble tools.
"""
from pyworkflow.protocol import params

from pwem import *
from pwem.objects import AtomStruct, Transform, String, EMFile
from pwem.protocols import EMProtocol

from pyworkflow.utils import *
from pyworkflow.protocol.params import (PointerParam, StringParam, IntParam, FloatParam,
                                        EnumParam, LEVEL_ADVANCED)

import prody

STRUCTURE = 0
INDEX = 1

BEST_MATCH = 0
SAME_CHID = 1

class ProDyBuildPDBEnsemble(EMProtocol):
    """
    This protocol will build a PDBEnsemble
    """
    _label = 'buildPDBEnsemble'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label='ProDy buildPDBEnsemble')

        form.addParam('structures', PointerParam, label="Set of structures",
                      important=True,
                      pointerClass='SetOfAtomStructs',
                      help='The structures to be aligned must be atomic models.')

        form.addParam('refType', EnumParam, choices=['structure', 'index'], 
                      default=STRUCTURE,
                      label="Reference structure selection type",
                      help='The reference structure can be a separate structure or indexed from the set')

        form.addParam('refStructure', PointerParam, label="Reference structure",
                      condition="refType == %d" % STRUCTURE,
                      pointerClass='AtomStruct',
                      help='Select an atomic model as the reference structure')

        form.addParam('refIndex', IntParam, label="Reference structure index",
                      condition="refType == %d" % INDEX,
                      help='Select the index of the reference structure in the set, starting from 1')

        form.addParam('seqid', FloatParam, default=0.,
                      expertLevel=LEVEL_ADVANCED,
                      label="Sequence Identity Cut-off (%)",
                      help='Alignment mapping with lower sequence identity will not be accepted.\n'
                           'This can be a number between 0 and 100 or a decimal between 0 and 1')

        form.addParam('overlap', FloatParam, default=0.,
                      expertLevel=LEVEL_ADVANCED,
                      label="Overlap Cut-off (%)",
                      help='Alignment mapping with lower sequence coverage will not be accepted.\n'
                           'This can be a number between 0 and 100 or a decimal between 0 and 1') 

        form.addParam('matchFunc', EnumParam, choices=['bestMatch', 'sameChid'], 
                      default=SAME_CHID,
                      expertLevel=LEVEL_ADVANCED,
                      label="Chain matching function",
                      help='See http://prody.csb.pitt.edu/manual/release/v1.11_series.html for more details.\n'
                           'Custom match functions will be added soon.')    

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # handle inputs
        if self.refType.get() == STRUCTURE:
            ref = prody.parsePDB(self.refStructure.get().getFileName(), alt='all')
        else:
            ref = self.refIndex.get() - 1 # convert from Scipion (sqlite) to Prody (python) nomenclature

        tars = prody.parsePDB([tarStructure.getFileName() for tarStructure in self.structures.get()], alt='all')

        # actual steps
        self._insertFunctionStep('alignStep', ref, tars)
        self._insertFunctionStep('createOutputStep')

    def alignStep(self, ref, tars):
        """This step includes alignment mapping and superposition"""

        if self.matchFunc.get() == BEST_MATCH:
            match_func = prody.bestMatch
        else:
            match_func = prody.sameChid

        ens = prody.buildPDBEnsemble(tars, ref=ref,
                                     seqid=self.seqid.get(),
                                     overlap=self.overlap.get(),
                                     match_func=match_func)

        ens = prody.trimPDBEnsemble(ens, occupancy=1)

        self.T = [T.getMatrix() for T in ens.getTransformations()]

        self.pdbFileName = self._getPath('ensemble.pdb')
        prody.writePDB(self.pdbFileName, ens)

        self.dcdFileName = self._getPath('ensemble.dcd')
        prody.writeDCD(self.dcdFileName, ens)

        self.npzFileName = self._getPath('ensemble.ens.npz')
        prody.saveEnsemble(ens, self.npzFileName)

        self.matrixFileName = self._getPath('transformation_%05d.txt')
        [prody.writeArray(self.matrixFileName % i, T) for i, T in enumerate(self.T)]

    def createOutputStep(self):
        outputStructure = AtomStruct(filename=self.pdbFileName)

        outputDcd = EMFile(filename=self.dcdFileName)
        outputNpz = EMFile(filename=self.npzFileName)
        #outputTrans = [Transform(matrix=T) for T in self.T]

        self._defineOutputs(outputStructure=outputStructure,
                            outputDcd=outputDcd,
                            outputNpz=outputNpz)
                            #outputTransformations=outputTrans)

