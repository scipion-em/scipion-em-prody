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
This module will provide ProDy deformation vector analysis.
"""
from pwem.emlib import MetaData, MDL_NMA_MODEFILE, MDL_NMA_ATOMSHIFT
from pwem.objects import AtomStruct, SetOfNormalModes, String
from pwem.protocols import EMProtocol

from pyworkflow.utils import join, makePath
from pyworkflow.protocol.params import (PointerParam, FloatParam, IntParam, 
                                        BooleanParam, LEVEL_ADVANCED)

import prody
import math

class ProDyDefvec(EMProtocol):
    """
    This protocol will perform deformation vector analysis
    """
    """
    Performs deformation vector analysis to explore structural differences 
    between two atomic or pseudoatomic models. The ProDy Deformation 
    protocol calculates the displacement vectors required to move a 
    mobile structure toward a target structure, enabling visualization, 
    quantification, and animation of conformational changes.

    Overview

    The protocol requires two input structures: a mobile structure, which 
    will be transformed, and a target structure, which defines the 
    desired conformation. Both structures should have matching numbers 
    of nodes. The main purpose is to provide insight into the structural 
    transitions between conformations and to generate animations of 
    atomic displacements along computed deformation vectors.

    Inputs and General Workflow

    Users provide the mobile and target structures, along with optional 
    parameters for animation, such as RMSD amplitude, number of frames, 
    and whether to include positive and negative directions. The protocol 
    first parses the input structures, computes RMSD if not provided, 
    and calculates the deformation vector representing the displacement 
    of each atom from mobile to target.

    The protocol generates a mode object representing the deformation 
    vector, writes it in Scipion mode and NMD formats, and optionally 
    animates the transformation using the specified RMSD and frame 
    settings. Visualization scripts compatible with VMD are also 
    produced to facilitate immediate inspection.

    Atom Shift Calculation

    Beyond animation, the protocol computes per-atom displacement 
    magnitudes and stores them as metadata. These profiles allow 
    quantitative analysis of the most mobile regions and identification 
    of atoms undergoing significant structural shifts.

    Outputs and Their Interpretation

    The primary output is a SetOfNormalModes object encapsulating the 
    deformation vector, associated mode files, and reference to the 
    mobile structure. Users can inspect animations, RMSD-based 
    displacement profiles, and transformed coordinates to interpret 
    conformational differences. The protocol is particularly useful 
    for comparing alternative structural states, exploring flexibility, 
    and generating illustrative animations for biological interpretation.

    Practical Recommendations

    In typical workflows, RMSD amplitudes can be left at zero to use 
    the actual structural deviation. The number of frames should be 
    chosen to balance smooth animations and computational cost. 
    Including both positive and negative directions provides a complete 
    view of the conformational landscape, while selective inclusion 
    may focus on biologically relevant transitions.
    """