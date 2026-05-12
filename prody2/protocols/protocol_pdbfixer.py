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

from os.path import basename, splitext

from pwem.objects import AtomStruct
from pwem.protocols import EMProtocol

from pyworkflow.protocol.params import PointerParam, FloatParam, LEVEL_ADVANCED

import prody
from prody2 import Plugin

class ProDyPDBFixer(EMProtocol):
    """
    Repairs and prepares atomic structures using the ProDy wrapper
    around OpenMM PDBFixer.

    The protocol is mainly used to add missing atoms, missing residues,
    and hydrogens to structural models.

    AI Generated:

    ProDy PDBFixer (ProDyPDBFixer) — User Manual
        Overview

        The ProDyPDBFixer protocol provides a structural preparation
        utility based on OpenMM PDBFixer.

        Its main purpose is to repair incomplete atomic structures so
        they can be used more reliably in downstream structural analysis,
        molecular simulation, or visualization workflows.

        In structural biology, experimentally determined models often
        contain incomplete regions such as:

            - missing side-chain atoms
            - missing backbone atoms
            - missing hydrogens
            - unresolved residues

        This protocol helps convert such incomplete models into more
        chemically consistent structural representations.

        Input Structure

        The protocol requires one input atomic structure.

        The input must be an:

            - AtomStruct

        Any atomic structure supported by ProDy can be used when the goal
        is simply to add hydrogens.

        However, if the user wants missing residues to be reconstructed,
        the input should preferably come from an original:

            - PDB file
            - mmCIF file

        containing SEQRES header information.

        This sequence information is important because residue recovery
        requires knowledge of the intended full protein sequence.

        pH Parameter

        The main user-controlled parameter is:

            pH

        This value determines the protonation conditions used when
        hydrogens are added.

        By default, the protocol uses:

            pH = 7.4

        This corresponds approximately to physiological conditions.

        Biologically, protonation state can strongly influence:

            - hydrogen bonding
            - electrostatic interactions
            - charge distribution
            - active-site chemistry

        For most structural analyses, the default is appropriate, but
        systems involving unusual environments may require a different
        pH.

        Computational Workflow

        The protocol performs two main steps.

        Structure Repair

        First, the input structure is passed to the PDBFixer backend.

        During this stage, the protocol can repair missing structural
        information and generate a corrected model.

        The repaired structure is written as a new PDB file whose name is
        derived from the original filename.

        Output Generation

        After repair, the corrected file is loaded into a new
        AtomStruct object.

        This repaired atomic structure becomes the protocol output.

        Output

        The protocol produces:

            - outputStructure

        This output is a new structural model containing the repaired
        atomic coordinates.

        It can be used directly in downstream workflows such as:

            - normal mode analysis
            - PCA
            - geometric measurements
            - molecular simulations
            - structural visualization

        Summary Information

        Once execution is complete, the protocol reports a summary
        comparing the original and repaired structures.

        The summary includes:

            - total number of atoms before and after fixing
            - total number of protein residues before and after fixing

        This gives the user a quick estimate of how much structural
        information was added.

        Biological Interpretation

        The biological relevance of this protocol lies in improving the
        chemical completeness of structural models.

        Missing atoms or incomplete residues can strongly affect
        downstream calculations.

        For example, missing atoms may distort:

            - steric contacts
            - residue packing
            - hydrogen bond networks
            - electrostatic interactions

        Adding hydrogens is especially important when structural analysis
        depends on realistic local chemistry.

        Practical Recommendations

        This protocol is particularly useful before performing analyses
        that are sensitive to structural completeness.

        Typical situations include:

            - preparing structures for molecular dynamics
            - refining structural models
            - comparing atomic-level geometries
            - generating chemically meaningful input for further ProDy
              analysis

        Users should keep in mind that missing residues added by PDBFixer
        are modeled rather than experimentally observed.

        Therefore, reconstructed regions should be interpreted with
        caution, especially when they correspond to flexible loops or
        poorly resolved segments.

        Final Perspective

        ProDyPDBFixer is best understood as a structural preparation
        protocol.

        Rather than asking:

            "What motions or measurements can be extracted?"

        it asks:

            "Is the structure chemically complete enough for reliable
            downstream analysis?"

        This makes it an important preprocessing step in many structural
        biology workflows.
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
            inputAg = prody.parsePDB(self.inputStructure.get().getFileName())
            outputAg = prody.parsePDB(self.outputStructure.getFileName())
            summ = ['The new structure has *{0}* atoms from original *{1}* atoms'.format(
                   outputAg.numAtoms(), inputAg.numAtoms())]
            summ.append('The new structure has *{0}* protein residues '
                        'from original *{1}* protein residues'.format(
                        outputAg.ca.numAtoms(), inputAg.ca.numAtoms()))
        return summ
