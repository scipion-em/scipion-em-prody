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

class  ProDyDomainDecomp(EMProtocol):
    """
    This protocol will perform dynamical domain decomposition
    """
    """
    Performs dynamical domain decomposition of a set of GNM normal modes 
    to identify structurally coherent domains within a protein or pseudoatomic model.

    Overview

    The ProDy Domain Decomposition protocol analyzes input normal modes 
    derived from Gaussian Network Models (GNM) to partition a structure 
    into dynamically correlated domains. Its main purpose is to reveal 
    functionally relevant substructures and flexible regions, helping 
    to interpret collective motions and structural modularity.

    Inputs and General Workflow

    Users provide a SetOfNormalModes object containing GNM modes, along 
    with the number of modes to analyze. The protocol parses the mode 
    files and the corresponding atomic coordinates, selecting the 
    specified number of modes. It then calculates domains based on 
    correlated motions using ProDy's GNM decomposition algorithms.

    The output structure is annotated with domain information encoded 
    in the B-factor column of the PDB, allowing immediate visualization 
    of dynamic regions. The protocol also generates a VMD script for 
    convenient graphical inspection, coloring domains and rendering 
    them as beads to illustrate the modular organization.

    Outputs and Their Interpretation

    The primary output is a PDB file with domain annotations, complemented 
    by a VMD script for visual exploration. The domains highlight regions 
    of coordinated motion, enabling users to understand which parts of 
    the structure move together and how flexibility is distributed. 
    These insights are valuable for interpreting allosteric effects, 
    conformational transitions, or potential sites for functional 
    regulation.

    Practical Recommendations

    Choosing the appropriate number of modes is critical: too few modes 
    may miss relevant correlations, while too many may produce noisy 
    partitions. Users are encouraged to inspect the output PDB and 
    VMD visualization to validate the biological relevance of the 
    decomposed domains and adjust the mode selection as necessary.
    """