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
This module will provide ProDy normal mode analysis using the anisotropic network model (ANM).
"""
import os
import numpy as np

from pwem.objects import AtomStruct, EMFile, String
from pwem.protocols import EMProtocol

from pyworkflow.utils import glob
from pyworkflow.protocol.params import PointerParam, EnumParam, BooleanParam

import prody

NMA_METRIC_OVERLAP = 0
NMA_METRIC_COV_OVERLAP = 1
NMA_METRIC_RWSIP = 2

class ProDyCompare(EMProtocol):
    """
    This protocol will compare two SetOfNormalModes objects
    """

    """
    Compares two sets of normal modes using RMSD- or covariance-based metrics. The ProDy Compare protocol 
    analyzes the similarity between two input mode sets, which can originate from atomic models, pseudoatomic 
    representations, or principal components. Its primary purpose is to quantify structural correlations, 
    identify matched modes, and provide representative outputs suitable for downstream analyses. 

    The protocol requires two input mode sets that ideally share the same number of nodes, except when one 
    set contains a single mode. Users can choose from different comparison metrics, including pairwise 
    overlaps, covariance overlaps, and root-weighted square inner products (RWSIP). Additional options allow 
    restricting calculations to diagonal overlaps, normalizing vectors, and matching modes to maximize 
    correspondence between the sets.

    During execution, the protocol parses the mode sets and associated PDB structures when available. If 
    mode matching is enabled, it aligns modes across ensembles and outputs matched mode files in both NMD 
    and Scipion formats. The protocol computes the selected metric and stores the resulting matrix for 
    interpretation, which reflects pairwise or cumulative mode similarities. Covariance overlaps and RWSIP 
    calculations exclude the six zero-frequency modes to ensure meaningful comparison.

    Outputs include the computed comparison matrix, optional matched mode sets, and PDB structures. These 
    outputs enable visualization, further computational analysis, or integration into workflows comparing 
    conformational dynamics. The protocol is particularly useful for evaluating structural consistency 
    between simulations, models, or experimental reconstructions, providing biologically relevant insights 
    into correlated motions and dominant conformational modes.
    """