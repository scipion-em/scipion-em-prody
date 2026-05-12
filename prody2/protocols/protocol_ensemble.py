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
from collections import OrderedDict
import numpy as np

from pwem.objects import (AtomStruct, SetOfAtomStructs, SetOfSequences,
                          EMFile)
from pwem.protocols import EMProtocol

from pyworkflow.utils import logger, getListFromRangeString, redStr
from pyworkflow.protocol.params import (PointerParam, MultiPointerParam,
                                        StringParam, IntParam, FloatParam,
                                        EnumParam, TextParam, NumericRangeParam,
                                        BooleanParam, LEVEL_ADVANCED)
from pyworkflow.object import Float

import prody
from prody2.objects import ProDyNpzEnsemble, TrajFrame
from prody2.constants import (NOTHING, PWALIGN, CEALIGN, DEFAULT,  # residue mapping methods
                              BEST_MATCH, SAME_CHID, SAME_POS, CUSTOM, # chain matching
                              ENSEMBLE_WEIGHTS)
from prody2 import parseMatchDict

import time

STRUCTURE = 0
INDEX = 1

BLAST = 0
DALI = 1

ENS_FILENAME = 'ensemble.dcd'

from prody2.objects import HAVE_CHEM
if HAVE_CHEM:
    from prody2.objects import DcdMDSystem

class ProDyBuildPDBEnsemble(EMProtocol):
    """
    This protocol will use ProDy's buildPDBEnsemble method to align atomic structures
    """
class ProDyBuildPDBEnsemble(EMProtocol):
    """
    Build PDB Ensemble (ProDyBuildPDBEnsemble) — User Manual

    Overview

    The Build PDB Ensemble protocol aligns multiple atomic structures into a shared
    coordinate system using ProDy's buildPDBEnsemble method. It creates a structural
    ensemble for comparative analysis, conformational studies, and downstream
    molecular interpretation.

    This protocol supports both user-provided structures and automatic retrieval of
    homologous structures via DALI searches, allowing flexible exploratory and
    in-depth structural analyses.

    Inputs and General Workflow

    The protocol accepts either explicit atomic structures or a PDB ID with a chain
    for automated DALI searches. Input structures can include all coordinate sets
    or only the active conformation. DALI-based retrieval supports filtering by
    RMSD, sequence identity, Z-score, and alignment length.

    Reference Structure Selection

    A reference structure defines the ensemble coordinate frame and can be provided
    explicitly or selected from the input set. The reference may optionally be removed
    after alignment to avoid biasing downstream analyses.

    Structural Alignment and Chain Matching

    Multiple chain matching strategies are available, including automatic best match,
    same-chain-ID, chain-position, or fully custom user-defined mapping. Custom chain
    matching is useful for multi-chain assemblies or structures with inconsistent
    chain labels.

    Residue Mapping and Structural Correspondence

    Residue correspondences are established using sequence alignment (Biopython pwalign),
    structural alignment (CE), or an automatic method. Thresholds for sequence identity
    and coverage ensure biologically meaningful mappings, especially for distant homologs.

    Atom Selection and Ensemble Construction

    Users define which atoms participate in the ensemble via selection strings
    (commonly alpha carbons). Unmapped structures are reported, and dummy atoms may
    be introduced for missing residues. Optional trimming removes poorly occupied atoms.

    Degeneracy and Multiple Coordinate Sets

    Supports NMR ensembles or multiple conformations per structure. Users can include
    all sets or restrict to the first active conformation. Including all sets captures
    structural heterogeneity.

    Output Generation and Ensemble Representation

    Outputs include a ProDy NPZ ensemble file, a multiple sequence alignment in FASTA,
    optionally aligned PDBs, and DCD trajectory files. Reference topologies are included
    when generating trajectories.

    Weight Management and Ensemble Metadata

    Structure or coordinate set weights are maintained throughout the workflow, allowing
    weighted ensemble analyses or integrative studies.

    Practical Recommendations

    Alpha-carbon selection with automatic chain matching is robust for homologous
    structures. Custom chain mapping should be used when chain identities differ.
    Trimming improves consistency for incomplete or flexible structures. DALI filters
    should be adjusted to prevent inclusion of distant homologs.

    Final Perspective

    The protocol provides a framework for organizing heterogeneous atomic models into
    biologically interpretable ensembles. Proper configuration of reference, chain
    matching, residue mapping, and trimming parameters is critical for meaningful
    structural analysis and conformational studies.
    """