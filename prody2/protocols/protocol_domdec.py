# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     James Krieger (jmkrieger@cnb.csic.es)
# *              Ricardo Serrano Gutiérrez (rserranogut@hotmail.com)                 
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
This module will provide ProDy Dynamical Domain Decomposition using the Gaussian Network Modeling (GNM).
"""

import os

from pwem.objects import AtomStruct, EMFile
from pwem.protocols import EMProtocol

from pyworkflow.utils import glob
from pyworkflow.protocol.params import PointerParam, IntParam

import prody

class ProDyDomainDecomp(EMProtocol):
    """
    This protocol will perform dynamical domain decomposition.

    AI Generated:

    Domain Decomposition (ProDyDomainDecomp) — User Manual
        Overview

        The Domain Decomposition protocol identifies dynamical domains
        in a macromolecular structure using Gaussian Network Model (GNM)
        normal modes. Its goal is to partition the structure into groups
        of residues or pseudoatoms that move collectively as coherent
        dynamical units.

        In structural biology, this analysis is especially useful when
        studying domain organization, hinge motions, cooperative
        rearrangements, and long-range allosteric communication.
        Rather than focusing on individual fluctuations, the protocol
        reveals larger structural regions that behave as quasi-rigid
        blocks during collective motion.

        Inputs and General Workflow

        The protocol requires a SetOfNormalModes object generated from
        a Gaussian Network Model analysis.

        The input modes may originate from:

        - An atomic model (true PDB structure)
        - A pseudoatomic model derived from an EM volume

        The input must specifically contain GNM modes, since the
        decomposition method is designed for Gaussian Network Model
        dynamics.

        The user must also choose how many modes will be included in
        the decomposition.

        During execution, the protocol:

        1. Loads the selected GNM normal modes.
        2. Extracts the specified number of low-frequency modes.
        3. Computes dynamical domains based on correlated collective
           motion.
        4. Assigns a domain label to each atom or pseudoatom.
        5. Writes a PDB file where domain assignments are stored in
           the beta-factor column.

        Biological Meaning of Dynamical Domains

        Low-frequency GNM modes usually capture the largest-scale,
        biologically relevant collective motions.

        By combining these modes, the protocol identifies regions of
        the structure that move together. These coherent groups often
        correspond to biologically meaningful structural domains.

        In practice, such domains may represent:

        - Compact rigid-body units
        - Flexible domains connected by hinge regions
        - Substructures involved in cooperative functional motion
        - Regions participating in allosteric signal transmission

        This makes domain decomposition especially valuable when trying
        to understand how intrinsic dynamics relate to biological
        function.

        Choice of Number of Modes

        The number of selected modes strongly affects the decomposition.

        Using only a few low-frequency modes usually highlights the
        largest and most biologically relevant domain organization.

        Including more modes introduces finer dynamical detail, which
        can sometimes reveal smaller subdomains but may also produce
        more fragmented results.

        For most biological applications, beginning with only the first
        few non-trivial low-frequency modes is generally the most
        interpretable strategy.

        Output Representation

        The protocol produces a PDB file in which domain assignments
        are encoded in the beta-factor column.

        This allows immediate visualization using standard molecular
        graphics software.

        A VMD script is also generated automatically.

        In the VMD representation:

        - Atoms are colored according to domain identity
        - A bead representation is used for rapid visual inspection

        This makes it straightforward to identify coherent moving
        regions, domain boundaries, and hinge-like interfaces.

        Biological Interpretation

        The resulting decomposition helps answer biologically relevant
        structural questions such as:

        - Which regions move as cooperative rigid units?
        - Where are the boundaries between moving domains?
        - Which interfaces may behave as hinges?
        - Are known conformational transitions compatible with the
          intrinsic dynamical organization?

        Importantly, the protocol does not directly predict a specific
        conformational transition. Instead, it reveals how the
        structure is intrinsically organized into dynamical units.

        Practical Recommendations

        For meaningful interpretation:

        - Use only GNM-derived normal modes.
        - Focus first on low-frequency modes, since these usually
          capture the most biologically relevant collective motions.
        - Avoid using too many modes initially, as this may lead to
          over-fragmentation.
        - Visually inspect the output and compare the resulting
          dynamical domains with known structural domains, flexible
          linkers, or functional interfaces.

        If the requested number of modes exceeds the available number
        of modes in the input set, the protocol will report an invalid
        input condition.

        Outputs and Their Interpretation

        After execution, the protocol produces:

        - A PDB structure containing domain labels encoded in beta
          values
        - A VMD script for immediate visualization

        Together, these outputs allow both quantitative and visual
        interpretation of the dynamical partitioning.

        Final Perspective

        For structural biologists, domain decomposition provides a
        practical bridge between normal mode analysis and structural
        interpretation.

        Instead of inspecting abstract eigenvectors individually,
        the protocol summarizes intrinsic collective dynamics into
        biologically meaningful cooperative regions.

        This often gives a much clearer understanding of large-scale
        flexibility, especially in multi-domain proteins, assemblies,
        and systems undergoing collective functional motion.
    """
    _label = 'Domain Decomposition'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need params to belong to a section:
        form.addSection(label='ProDy DomainDecomp')

        form.addParam('modesGNM', PointerParam, label="Input SetOfNormalModes GNM",
                      important=True,
                      pointerClass='SetOfNormalModes',
                      help='The input SetOfNormalModes can be from an atomic model '
                           '(true PDB) or a pseudoatomic model '
                           '(an EM volume compared into pseudoatoms).\n'
                           'The set must be GNM modes')
        form.addParam('modeNumber', IntParam, default=2,
                label='Number of modes')
    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('computeDecompStep')
        self._insertFunctionStep('createOutputStep')

    def computeDecompStep(self):
        modesPath = os.path.dirname(os.path.dirname(self.modesGNM.get()[1].getModeFile()))
        modes = prody.parseScipionModes(self.modesGNM.get().getFileName(),
                                            pdb=glob(modesPath+"/*atoms.pdb"))

        numModes = self.modeNumber.get()
        
        try:
            mode = modes[:numModes] 
        except IndexError:
            return [self.errorMessage("Invalid number of modes *%d*\n"
                                      "Display the output Normal Modes to see "
                                      "the availables ones." % numModes,
                                      title="Invalid input")] 
        atoms = prody.parsePDB(glob(modesPath+"/*atoms.pdb"))

        domains = prody.calcGNMDomains(mode)

        self.pdbFilename = self._getPath("atoms.pdb")
        prody.writePDB(self.pdbFilename, atoms, beta=domains)

    def createOutputStep(self):
        fhCmd=open(self._getPath("domains.vmd"),'w')
        fhCmd.write("mol new %s\n" % self.pdbFilename)
        fhCmd.write("mol modcolor 0 0 Beta\n")
        fhCmd.write("mol modstyle 0 0 Beads\n")
        fhCmd.close()

        outputPdb = AtomStruct()
        outputPdb.setFileName(self.pdbFilename)

        outputvmd = EMFile()
        outputvmd.setFileName(self._getPath("domains.vmd"))
        
        self._defineOutputs(outputStructure=outputPdb, outputvmd=outputvmd)
