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
This module will provide ProDy atom tools including selection and superposition.
"""
from collections import OrderedDict
from os.path import basename, splitext, abspath

from pwem.objects import AtomStruct, SetOfAtomStructs, Transform, CsvList
from pwem.protocols import EMProtocol

from pyworkflow.utils import exists
from pyworkflow.protocol.params import (PointerParam, StringParam, FloatParam,
                                        BooleanParam, EnumParam, TextParam, IntParam,
                                        PathParam, MultiPointerParam, LEVEL_ADVANCED)

import prody
from pyworkflow.utils import logger

from prody2 import Plugin
from prody2.objects import Atom, SetOfAtoms
from prody2.constants import (NOTHING, PWALIGN, CEALIGN, DEFAULT,  # residue mapping methods
                              BEST_MATCH, SAME_CHID, SAME_POS, CUSTOM) # chain matching

def notFoundException(inputFn):
    return Exception("Atomic structure not found at *%s*" % inputFn)

UNITE_CHAINS_LABEL = "Unite chains in mmCIF segments"
UNITE_CHAINS_HELP = ('Elect whether to unite chains in mmCIF segments for each structure like ChimeraX. '
                     'Default is **False**, which means the smaller unit IDs are used for chains like PyMOL.')

IMPORT_FROM_ID_CONDITION = 'inputPdbData == IMPORT_FROM_ID'
SUMMARY_NO_OUTPUT = 'Output structure not ready yet'
NOT_DUMMY_SELSTR = "not dummy"


class ProDyAtomicBase(EMProtocol):
    """
    This protocol will perform atom selection
    """
    """
    Downloads and parses molecular dynamics trajectories from the
    BioExcel COVID-19 database, generating topology and trajectory
    objects that can be used for structural dynamics analysis.

    AI Generated:

    BioExcel CV19 Parser (ProDyBioExcelCV19) — User Manual
        Overview

        The BioExcel CV19 protocol provides an interface for accessing
        and processing molecular dynamics simulations deposited in the
        BioExcel COVID-19 database. Its primary purpose is to retrieve
        atomistic simulation trajectories together with their associated
        topology information and integrate them into structural biology
        and cryo-EM analysis workflows.

        From a biological perspective, this protocol allows researchers
        to study protein flexibility, conformational transitions, residue
        mobility, and large-scale molecular motions derived from molecular
        dynamics simulations. Such analyses are especially valuable when
        investigating viral proteins, membrane complexes, ligand-binding
        mechanisms, or dynamic conformational states that cannot be fully
        captured by static experimental structures alone.

        The protocol is particularly useful for integrating simulation
        data into hybrid structural biology workflows, where experimental
        cryo-EM maps and computational molecular dynamics trajectories
        complement one another to provide a more complete understanding
        of biomolecular behavior.

        Inputs and General Workflow

        The protocol requires a simulation accession identifier that
        corresponds to an entry in the BioExcel COVID-19 database.
        Using this accession, the workflow automatically retrieves
        the associated topology and trajectory files, including
        coordinate information in PDB format, structural topology
        definitions in PSF format, and trajectory frames stored in
        DCD format.

        During execution, the protocol launches an external BioExcel
        retrieval utility that downloads and prepares all required
        simulation files inside the workflow environment. Once the
        download is completed, the generated files are organized and
        converted into trajectory-compatible objects that can be used
        in downstream analyses.

        This automated retrieval process greatly simplifies access to
        large molecular dynamics datasets, allowing biological users
        to focus on structural interpretation rather than file handling
        or trajectory conversion.

        Atom Selection and Structural Reduction

        One of the most important biological options provided by this
        protocol is atom selection. Molecular dynamics trajectories
        can contain extremely large numbers of atoms, particularly
        for membrane proteins, viral assemblies, or solvated systems.
        Processing the full trajectory may therefore become expensive
        both computationally and analytically.

        To address this, the protocol allows the user to select only
        specific subsets of atoms. Carbon alpha selections provide a
        coarse-grained representation of protein backbone motion and
        are especially useful for large-scale conformational analyses,
        principal component calculations, or elastic network modeling.

        Backbone selections preserve the main structural scaffold of
        the protein while excluding side chains and solvent atoms,
        offering a balance between structural detail and computational
        efficiency. The combined backbone and carbon-alpha selection
        further reduces complexity while still retaining biologically
        meaningful information about protein architecture.

        From a biological interpretation standpoint, reduced atom
        selections are often sufficient for studying domain motions,
        global flexibility, and conformational transitions. Full
        atomistic representations become more important when analyzing
        local interactions, ligand binding, side-chain rearrangements,
        or solvent-mediated effects.

        Frame Selection and Trajectory Sampling

        The protocol also supports selective extraction of trajectory
        frames. Instead of processing an entire simulation, users may
        specify ranges or intervals of frames using compact selection
        expressions.

        Biologically, frame selection becomes especially important
        when studying long simulations that contain multiple temporal
        states or conformational transitions. Selecting only specific
        time windows allows researchers to focus on biologically
        relevant regions of the trajectory while reducing storage
        and computational requirements.

        For example, early frames may represent equilibration stages,
        whereas later frames may correspond to stabilized conformational
        states. Similarly, sparse frame sampling may be sufficient for
        exploratory analyses, while dense sampling becomes more useful
        when studying rapid structural fluctuations or detailed kinetic
        transitions.

        Output Trajectory Representation

        After processing, the protocol generates a molecular dynamics
        trajectory representation that includes the topology and the
        associated trajectory coordinates.

        When chemical support is enabled, the output is stored as a
        ProDyMDSystem object containing direct references to topology
        and trajectory files. This representation is optimized for
        downstream molecular dynamics analyses and structural dynamics
        workflows.

        Alternatively, the protocol can generate a SetOfTrajFrames
        object in which each frame is individually represented and
        assigned a statistical weight. This mode becomes useful in
        workflows requiring explicit frame manipulation, trajectory
        sampling, or integration into ensemble-based analyses.

        From a structural biology perspective, these trajectory objects
        provide the basis for downstream calculations such as RMSD,
        RMSF, normal mode analysis, principal component analysis,
        conformational clustering, and flexible fitting into cryo-EM
        density maps.

        Outputs and Biological Interpretation

        The final outputs include the processed topology structure,
        trajectory coordinates, and frame organization required for
        subsequent analyses. The protocol additionally reports the
        number of atoms and protein residues present in the generated
        trajectory.

        This information is biologically useful because it allows users
        to verify that the selected atom subset correctly represents
        the intended molecular system. In many practical workflows,
        validating atom counts and residue coverage is essential before
        beginning downstream structural interpretation.

        Practical Recommendations

        In routine molecular dynamics workflows, carbon-alpha or
        backbone selections are often sufficient for exploratory
        analyses and significantly reduce computational cost. Full
        atomistic trajectories should generally be reserved for cases
        where detailed interaction analysis is required.

        When working with very large simulations, selecting restricted
        frame ranges can greatly improve performance without losing
        biologically meaningful information. Researchers should also
        carefully consider whether trajectory averaging or sparse
        sampling may obscure transient but functionally important
        conformational states.

        For integrative cryo-EM studies, reduced trajectories are
        frequently advantageous because they simplify fitting and
        flexibility analyses while preserving the dominant structural
        motions relevant to experimental density interpretation.

        Final Perspective

        The BioExcel CV19 protocol provides a bridge between molecular
        dynamics simulations and integrative structural biology
        workflows. Rather than treating simulations as isolated
        computational datasets, the protocol enables their direct
        incorporation into biological interpretation pipelines.

        By facilitating controlled atom selection, frame sampling,
        and trajectory organization, the protocol allows researchers
        to efficiently explore conformational landscapes and dynamic
        structural behavior in biologically meaningful ways. This
        integration of simulation data with structural analysis
        workflows is increasingly important for understanding complex
        biomolecular systems, especially in modern cryo-EM and
        computational structural biology environments.
    """