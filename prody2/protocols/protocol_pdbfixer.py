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
This module will provide the ProDy wrapper for OpenMM PDBFixer
"""
from pyworkflow.protocol import params

from os.path import basename, splitext


from pwem import *
from pwem.objects import AtomStruct
from pwem.protocols import EMProtocol

from pyworkflow.utils import *
from pyworkflow.protocol.params import PointerParam, FloatParam, LEVEL_ADVANCED

import prody
from prody2 import Plugin

class ProDyPDBFixer(EMProtocol):
    """
    This module will provide the ProDy wrapper for OpenMM PDBFixer
    """
    _label = 'PDBFixer'
    _possibleOutputs = {'outputStructure': AtomStruct}

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        form.addSection(label='ProDy PDBFixer modelling')
        
        form.addParam('inputStructure', PointerParam, label="Input structure",
                      important=True,
                      pointerClass='AtomStruct',
                      help='The input structure should be an AtomStruct. '
                      'Any AtomStruct from ProDy can be used for adding hydrogens '
                      'but one corresponding to an original PDB or mmCIF file with '
                      'SEQRES header data is needed to add missing residues.')

        form.addParam('pH', FloatParam, label="pH", default=7.4,
                      important=True, expertLevel=LEVEL_ADVANCED,
                      help='The pH to use for adding hydrogens.')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('computeStep')
        self._insertFunctionStep('createOutputStep')

    def computeStep(self):
        inputFn = self.inputStructure.get().getFileName()
        self.outputFn = self._getPath(splitext(basename(inputFn))[0] + '_fixed.pdb')

        args = '--inputFn {0} --pH {1} --outputFn {2}'.format(inputFn, self.pH.get(), self.outputFn)
        self.runJob(Plugin.getProgram('fixer.py', script=True), args)

    def createOutputStep(self):
        outAS = AtomStruct(self.outputFn)
        self._defineOutputs(outputStructure=outAS)

    def _summary(self):
        if not hasattr(self, 'outputStructure'):
            summ = ['Output structure not ready yet']
        else:
            input_ag = prody.parsePDB(self.inputStructure.get().getFileName())
            output_ag = prody.parsePDB(self.outputStructure.getFileName())
            summ = ['The new structure has *{0}* atoms from original *{1}* atoms'.format(
                   output_ag.numAtoms(), input_ag.numAtoms())]
            summ.append('The new structure has *{0}* protein residues '
                        'from original *{1}* protein residues'.format(
                        output_ag.ca.numAtoms(), input_ag.ca.numAtoms()))
        return summ
