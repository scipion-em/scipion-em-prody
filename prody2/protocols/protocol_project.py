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
This module will provide ProDy projection of structural ensembles on principal component or normal modes
"""
import numpy as np
import os

from pwem.objects import SetOfAtomStructs
from pwem.protocols import EMProtocol

import pyworkflow.object as pwobj
from pyworkflow.protocol.params import (PointerParam, EnumParam, BooleanParam,
                                        MultiPointerParam, NumericRangeParam)
from pyworkflow.utils import getListFromRangeString, glob

import prody
from prody2.constants import PROJ_COEFFS

ONE = 0
TWO = 1
THREE = 2

class ProDyProject(EMProtocol):
    """
    This module will provide ProDy projection of structural ensembles on principal component or normal modes
    """


class ProDyProject(EMProtocol):
    """
    AI Generated Summary:

    Structural Projection (ProDyProject) — User Manual

    OVERVIEW
    The Projection protocol maps structural ensembles onto a reduced coordinate
    system defined by Normal Mode Analysis (NMA) or Principal Component
    Analysis (PCA). Its primary goal is to simplify high-dimensional molecular
    movements into discrete coefficients that describe a structure's position
    along specific biological pathways or functional transitions.

    INPUTS AND COMPATIBILITY
    The protocol requires two essential components:
    - Structural Ensembles: Atomic sets (PDBs) or ProDy NPZ ensembles.
    - Reference Modes: Principal Components or Normal Modes.
    A fundamental constraint is that all input structures must share the
    identical atom count as the model used to generate the reference modes
    to ensure mathematical validity during vector overlap calculations.

    MODE SELECTION AND SUBSPACE
    Users can define the dimensionality of the projection space:
    - Selective Indexing: Specific modes can be chosen using range strings
      (e.g., "1, 2, 5-10").
    - Multi-Mode Projection: The protocol allows projecting onto 1, 2, or 3
      modes simultaneously, which is ideal for creating 2D or 3D
      conformational landscapes to identify structural clusters.

    SCALING AND NORMALIZATION
    To enhance biological interpretability, the protocol provides:
    - RMSD Scaling: Expresses the projection in Angstroms, allowing
      researchers to quantify how far a structure has moved along a mode.
    - Normalization: Adjusts coefficients relative to the mode magnitude,
      facilitating comparison between different structural models.

    OUTPUTS AND DATA EXPORT
    The results are integrated into the output ensembles as new attributes
    (projection coefficients). Additionally, the protocol exports:
    - CSV Files: Containing raw coefficients and structural weights for
      external statistical analysis.
    - NMD Files: For visual inspection of the projection subspace in
      standard NMA viewers.
    """
