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
This module will provide ProDy normal mode analysis (NMA) using the the rotation and translation of blocks (RTB) framework.
"""
from pwem.emlib import (MetaData, MDL_NMA_MODEFILE, MDL_ORDER,
                        MDL_ENABLED, MDL_NMA_COLLECTIVITY, MDL_NMA_SCORE, 
                        MDL_NMA_EIGENVAL)
from pwem.objects import SetOfNormalModes, String

from pyworkflow.utils import glob, redStr
from pyworkflow.protocol.params import (PointerParam, IntParam, FloatParam,
                                        BooleanParam, EnumParam, LEVEL_ADVANCED)

import prody
from prody2.protocols.protocol_modes_base import ProDyModesBase

BLOCKS_FROM_RES = 0
BLOCKS_FROM_SECSTR = 1

class ProDyRTB(ProDyModesBase):
    """
    This protocol will perform normal mode analysis (NMA) using the rotation and translation of blocks (RTB) framework
    """

    """
    AI Generated Summary:

    RTB Normal Mode Analysis (ProDyRTB) — User Manual

    OVERVIEW
    The RTB (Rotation-Translation of Blocks) protocol is a high-performance 
    implementation of Normal Mode Analysis. It is specifically designed to 
    handle massive macromolecular complexes by partitioning the structure into 
    rigid blocks, thereby reducing the degrees of freedom and making the 
    calculation of global motions computationally feasible for systems that 
    would otherwise exhaust system memory.

    BLOCK DEFINITION AND COARSE-GRAINING
    The efficiency of the method relies on how the molecule is divided:
    - Residue-Based: Groups a specific number of residues into single blocks.
    - Secondary Structure: Groups atoms based on biological units like alpha 
      helices and beta sheets, preserving essential structural domains.
    - Optimization: Automatically splits long blocks and merges short ones 
      to maintain numerical stability and structural relevance.

    PHYSICAL MODELING (ENM)
    The protocol builds a Hessian matrix based on an Elastic Network Model:
    - Cut-off Distance: Controls the interaction range between blocks (default 
      15 Å for C-alpha).
    - Spring Constant (Gamma): Defines the strength of the virtual springs 
      connecting the blocks.
    - Hessian Matrix: Describes the potential energy landscape of the 
      partitioned system.

    PERFORMANCE AND VALIDATION
    - Turbo Mode: An optimized, memory-intensive algorithm for rapid 
      eigenvector decomposition.
    - Sparse Matrices: Automatically utilized if the system detects potential 
      MemoryErrors, ensuring robustness for extremely large assemblies.
    - Collectivity Filtering: Identifies collective functional motions vs. 
      local fluctuations, excluding the first six rigid-body modes.

    OUTPUTS AND VISUALIZATION
    The protocol generates a 'SetOfNormalModes' and NMD files. It also 
    produces VMD-compatible animations where the user can adjust RMSD 
    amplitude to visually inspect how the protein 'breathes' or 'twists' 
    along the calculated vectors, providing a direct link between 
    mathematical modes and biological function.
    """
