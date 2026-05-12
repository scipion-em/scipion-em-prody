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
This module will provide ProDy principal component analysis (PCA) using atomic structures
"""

import numpy as np

from pwem.objects import SetOfAtomStructs, AtomStruct
from pwem.protocols import EMProtocol

from pyworkflow.protocol.params import (PointerParam, EnumParam, FloatParam,
                                        StringParam, BooleanParam, IntParam)
import pyworkflow.object as pwobj
import pyworkflow.utils as pwutils

from prody2.objects import (ProDyNpzEnsemble, TrajFrame, 
                            SetOfClassesTraj, ClassTraj)
from prody2.constants import ENSEMBLE_WEIGHTS
from prody2 import Plugin

import prody
import matplotlib.pyplot as plt


class ProDyRmsd(EMProtocol):
    """
    This protocol will perform ProDy principal component analysis (PCA) using atomic structures
    """
    """
    The ProDy RMSD protocol performs principal component analysis (PCA) and clustering
    on atomic structure ensembles using RMSD metrics. Its primary goal is to analyze
    structural variability, detect representative conformations, and optionally reorder
    or cluster the ensemble for downstream analyses. This protocol is particularly
    useful when handling multiple structures derived from simulations, experimental
    conditions, or modeling pipelines, as it helps extract biologically meaningful
    conformational differences.

    The protocol requires an input ensemble, which can be either a SetOfAtomStructs
    or a ProDyNpzEnsemble. All structures must contain the same number of atoms to
    ensure meaningful RMSD computations. Users can choose to cluster the ensemble
    using hierarchical clustering or k-medoids. Hierarchical clustering optionally
    reorders structures based on RMSD trees and allows defining subgroups via an RMSD
    threshold, while k-medoids partitions the ensemble into a fixed number of clusters
    and identifies representative medoids. Cluster weights are assigned proportionally
    based on subgroup sizes, reflecting the relative importance of each conformer.

    During execution, the protocol calculates pairwise RMSDs, generates clustering trees
    or medoid assignments, computes representative structures, and can reorder the
    ensemble according to similarity. Users may also generate representative PDB files
    for individual structures. Outputs include sets of clustered structures, each
    represented by a ClassTraj object, and optionally the reordered ensemble and
    PDB structures. Representative structures are weighted to preserve ensemble statistics,
    enabling downstream PCA or comparative analysis.

    For practical use, hierarchical clustering is recommended to explore relationships
    between conformers, while k-medoids is suited for fixed-size clustering. Generating
    representative PDB files facilitates visualization and integration with other
    modeling workflows. Proper selection of RMSD thresholds, cluster numbers, and tree
    methods ensures biologically meaningful grouping without losing structural detail.

    ProDy RMSD is not just a computational protocol but a tool to extract biologically
    relevant conformational information from atomic ensembles. Thoughtful application
    of clustering, weighting, and optional reordering provides insight into molecular
    flexibility, dominant states, and representative structures suitable for further
    analysis or publication-quality figures.
    """