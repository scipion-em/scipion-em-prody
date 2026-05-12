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
    This module will provide the ProDy wrapper for OpenMM PDBFixer
    """
    """
    This protocol serves as a wrapper for the OpenMM PDBFixer tool, providing
    automated structural repair and preparation for macromolecular models.

    AI Generated:

    PDBFixer Modeling (ProDyPDBFixer) — User Manual
        Overview

        The PDBFixer protocol is designed to address common structural deficiencies in
        Protein Data Bank (PDB) files, such as missing atoms, incomplete residues,
        and the absence of hydrogen atoms. Its main purpose is to prepare "simulation-ready"
        models by ensuring structural integrity and proper protonation states. In
        typical structural biology workflows, this step is essential before performing
        energy minimization, molecular dynamics, or advanced flexibility analysis,
        where physical consistency is paramount.

        For a biological user, the most common applications include repairing
        loops that were not resolved in an experiment, standardizing a model
        for molecular simulations, or adding hydrogens to explore hydrogen-bonding
        networks. The protocol streamlines the transition from raw experimental
        coordinates to refined models suitable for biophysical computation.

        Inputs and General Workflow

        The protocol requires an input atomic structure (AtomStruct). While any
        standard model can be processed for hydrogen addition, the protocol
        reaches its full potential when provided with files containing original
        SEQRES header data (from PDB or mmCIF sources). This header data is
        biologically vital as it defines the complete sequence, allowing the
        software to identify and reconstruct residues that were missing from the
        experimental density.

        The workflow is automated and direct: the user provides the structure
        and specifies the environmental conditions. The protocol then invokes
        the OpenMM PDBFixer engine to analyze the structure, identify gaps,
        add heavy atoms, and finally protonate the molecule.

        Environmental Context: pH Sensitivity

        Protonation is a biologically critical process because the presence and
        position of hydrogen atoms depend heavily on the chemical environment.
        The protocol allows the user to specify a pH value, which defaults to
        physiological conditions (7.4).

        From a biological perspective, adjusting the pH is necessary when
        studying proteins that function in specific organelles (like the acidic
        lumen of a lysosome) or under non-standard experimental conditions.
        Correctly assigning hydrogen atoms based on pH ensures that the
        electrostatic properties and ionization states of amino acid side
        chains—such as Histidine, Lysine, and Glutamate—are accurately
        represented.

        Structural Reconstruction and Modeling

        Beyond simple atom addition, PDBFixer acts as a modeling tool. By
        reconciling the observed coordinates with the expected sequence, it
        can fill in missing side chains and loops. This reconstruction
        eliminates "gaps" in the protein backbone that would otherwise
        cause instabilities in downstream physics-based analyses.

        The protocol produces a repaired model that maintains the original
        spatial orientation while significantly increasing the total atom count.
        This "fixed" structure serves as the definitive reference for
        subsequent dynamic or static modeling steps.

        Outputs and Their Interpretation

        After execution, the protocol produces a new AtomStruct labeled as "fixed."
        This output is a complete PDB file containing all repaired atoms and
        newly added hydrogens.

        The protocol summary provides an immediate comparison of the atom and
        residue counts before and after the process. Biologically, an increase
        in the protein residue count (tracked via Alpha Carbons) confirms that
        missing segments have been successfully modeled, while the total atom
        count change reflects the addition of heavy atoms and the complete
        protonation of the system.

        Practical Recommendations

        In routine practice, it is advisable to inspect the resulting "fixed" PDB
        visually, particularly in regions where large loops were missing. If
        the reconstructed regions appear too extended or energetically
        unfavorable, subsequent energy minimization is highly recommended.

        When working with high-resolution structures that already contain some
        hydrogens, the protocol will intelligently handle the existing data to
        ensure a consistent protonation state across the entire complex.

        Final Perspective

        For most Scipion users, PDBFixer is a necessary bridge between experimental
        maps and theoretical models. Careful attention to the input pH and the
        inclusion of sequence metadata are the key elements for ensuring that
        the final model is not just mathematically complete, but biologically
        accurate for further scientific investigation.
    """