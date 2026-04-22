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
This module will provide the ProDy wrapper for OpenMM PDBFixer
"""

from os.path import basename, splitext

from pwem.objects import AtomStruct, Integer
from pwem.protocols import EMProtocol

from pyworkflow.protocol.params import PointerParam, FloatParam, LEVEL_ADVANCED
from pyworkflow.utils import exists

from prody2.constants import N_ATOMS, N_RESIDUES, N_CHAINS
from prody2 import Plugin

SUMMARY_NO_OUTPUT = 'Output structure not ready yet'

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

        args = '--inputFn {0} --pH {1} --outputFn {2} --folder {3}'.format(
            inputFn, self.pH.get(), self.outputFn, self._getPath())
        self.runJob(Plugin.getProgram('fixer.py', script=True), args)

    def createOutputStep(self):
        with open(self._getPath('pdb_data.txt'), 'r') as fi:
            line = fi.readlines()[0]

        self.pdbFileName, numAtoms, numResidues, numChains = line.split('\t')
        if exists(self.pdbFileName):
            outputPdb = AtomStruct()
            setattr(outputPdb, N_ATOMS, Integer(numAtoms))
            setattr(outputPdb, N_RESIDUES, Integer(numResidues))
            setattr(outputPdb, N_CHAINS, Integer(numChains))
            outputPdb.setFileName(self.pdbFileName)
            self._defineOutputs(outputStructure=outputPdb)

    def _summary(self):
        if not hasattr(self, 'outputStructure'):
            summ = [SUMMARY_NO_OUTPUT]
        else:
            summ = ['The new structure has *{0}* residues '
                     'and *{1}* atoms in *{2}* chains'.format(
                     self.outputStructure.getAttributeValue(N_RESIDUES),
                     self.outputStructure.getAttributeValue(N_ATOMS),
                     self.outputStructure.getAttributeValue(N_CHAINS))]
        return summ
