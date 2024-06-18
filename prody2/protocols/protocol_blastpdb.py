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
This module will provide the ProDy interface for running an NCBI BLAST
search against the PDB
"""
from pwem.protocols import EMProtocol
from pyworkflow.protocol import params

import prody

from prody2.objects import SetOfBlastHits, createSetOfBlastResults
from prody2.constants import UNITE_CHAINS_LABEL, UNITE_CHAINS_HELP


class ProDyBlastPDB(EMProtocol):
    """
    This module will provide the ProDy interface for running an 
    NCBI BLAST search against the PDB
    """
    _label = 'BlastPDB'
    _possibleOutputs = {'outputResults': SetOfBlastHits}

    IMPORT_FROM_STRUCT = 0
    IMPORT_FROM_SEQ = 1
    USE_TEXT = 2

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        form.addSection(label='ProDy BlastPDB parsing')

        form.addParam('inputSeqData', params.EnumParam, choices=['AtomStruct', 'Sequence', 'text'],
                      label="Import sequence from",
                      default=self.USE_TEXT,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Extract sequence from AtomStruct, Sequence or text')

        form.addParam('inputStructure', params.PointerParam, label="Input structure",
                      condition='inputSeqData == IMPORT_FROM_STRUCT',
                      pointerClass='AtomStruct',
                      help='Provide any atomic structure with a sequence')

        form.addParam('uniteChains', params.BooleanParam, default=False,
                      condition='inputSeqData == IMPORT_FROM_STRUCT',
                      label=UNITE_CHAINS_LABEL,
                      help=UNITE_CHAINS_HELP)

        form.addParam('inputSequence', params.PointerParam, label="Input sequence",
                      condition='inputSeqData == IMPORT_FROM_SEQ',
                      pointerClass='Sequence',
                      help='Provide any Sequence object')

        form.addParam('inputSeqText', params.TextParam, width=50,
                       condition='inputSeqData == USE_TEXT', default="",
                       label='Input sequence',
                       help='Defined order of chains from custom matching.')       

        form.addParam('seqid', params.FloatParam, default=0.,
                      label="Sequence identity percentage cutoff",
                      help='Hits with lower percent sequence identity will not be included.\n'
                           'This can be a number between 0 and 100')

        form.addParam('overlap', params.FloatParam, default=0.,
                      label="Overlap percentage cutoff",
                      help='Hits with lower percent sequence coverage will not be included.\n'
                           'This can be a number between 0 and 100')

        form.addParam('separateChains', params.BooleanParam, default=False,
                      label="Whether to separate chains",
                      help="If true, a separate entry will be provided for each chain")

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('computeStep')
        self._insertFunctionStep('createOutputStep')

    def computeStep(self):
        """Run blastPDB and save the results to an xml file"""

        if self.inputSeqData.get() == self.IMPORT_FROM_STRUCT:
            ag = prody.parsePDB(self.inputStructure.get().getFileName(),
                                unite_chains=self.uniteChains.get())
            sequence = ag.ca.getSequence()
        
        elif self.inputSeqData.get() == self.IMPORT_FROM_SEQ:
            sequence = self.inputSequence.get().getSequence()

        else:
            sequence = self.inputSeqText.get()
    
        self.blastRec = prody.blastPDB(sequence=sequence, 
                                       filename=self._getExtraPath('results.xml'), 
                                       timeout=1e10) # ensure completion

    def createOutputStep(self):
        outputHitsSet = createSetOfBlastResults(self.blastRec, self,
                                                self.seqid.get(),
                                                self.overlap.get(),
                                                self.separateChains.get())
        self._defineOutputs(outputResults=outputHitsSet)

    def _summary(self):
        if not hasattr(self, 'outputResults'):
            summ = ['Output results are not ready yet']
        else:
            summ = ['blastPDB retrieved *{0}* hits with the input criteria'.format(
                len(self.outputResults))]
        return summ
