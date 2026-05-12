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
This module will provide ProDy mode editing tools.
"""
import os
import numpy as np

from pwem.objects import AtomStruct, SetOfNormalModes, SetOfPrincipalComponents, String

from pyworkflow.utils import glob, logger
from pyworkflow.protocol.params import (PointerParam, EnumParam, BooleanParam,
                                        FloatParam, IntParam, LEVEL_ADVANCED)

import prody
from prody2.protocols.protocol_modes_base import ProDyModesBase

NMA_SLICE = 0
NMA_REDUCE = 1
NMA_EXTEND = 2
NMA_INTERP = 3

class ProDyEdit(ProDyModesBase):
    """
    This protocol will edit a SetOfNormalModes object to have more or fewer nodes
    """
    """
    Edits a SetOfNormalModes object to change the number of nodes, 
    enabling flexible modification of normal mode representations 
    for atomic or pseudoatomic models.

    Overview

    The ProDy Edit protocol allows the user to adjust the granularity 
    of a normal mode set by either reducing, slicing, extending, or 
    interpolating the nodes. Its main purpose is to facilitate 
    comparative analyses, enhance visualization, or prepare modes 
    for further computational workflows by matching them to a desired 
    atomic or pseudoatomic representation.

    Inputs and General Workflow

    Users provide a SetOfNormalModes object along with a new set of 
    atoms or pseudoatoms representing the desired nodes. The protocol 
    supports multiple editing strategies: slicing for direct reduction, 
    Hessian-based reduction for meaningful eigenvector compression, 
    extension by copying values across residues, and interpolation 
    for smooth mapping between node sets. Optional normalization can 
    be applied to maintain consistent vector magnitudes. Additionally, 
    the protocol can generate animations of the modified modes for 
    visualization in ContinuousFlex, specifying RMSD amplitudes, 
    number of frames, and directionality.

    During execution, the protocol aligns old and new node sets, 
    applies the chosen editing strategy, and produces updated modes 
    alongside corresponding atomic coordinates. These outputs are 
    saved in standard PDB and mode file formats compatible with 
    ProDy and Scipion workflows, preserving both structural and 
    dynamic information.

    Outputs and Their Interpretation

    The edited modes are provided as a new SetOfNormalModes object, 
    annotated with the new atomic coordinates. The outputs include 
    PDB files of the new node set, mode files compatible with Scipion, 
    and optionally animation files for visual inspection. Users can 
    interpret these outputs to assess how changes in node representation 
    affect dynamic patterns, mode amplitudes, and structural correlations, 
    facilitating both exploratory and publication-level analyses.

    Practical Recommendations

    Careful selection of the editing strategy is recommended: slicing 
    and reduction are suited for decreasing node counts while preserving 
    dynamic relevance, whereas extension and interpolation are ideal 
    for increasing resolution or adapting modes to new atomic models. 
    Visualizing the resulting modes or generating animations helps 
    ensure that the edited representations remain biologically meaningful.
    """