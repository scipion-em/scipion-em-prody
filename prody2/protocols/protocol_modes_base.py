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
This module will provide ProDy normal mode analysis (NMA) using the anisotropic network model (ANM).
"""
import math
from os.path import exists, join

from pwem.emlib import (MetaData, MDL_NMA_MODEFILE, MDL_ORDER,
                        MDL_ENABLED, MDL_NMA_COLLECTIVITY, MDL_NMA_SCORE, 
                        MDL_NMA_ATOMSHIFT, MDL_NMA_EIGENVAL)
from pwem.objects import AtomStruct, SetOfNormalModes, String
from pwem.protocols import EMProtocol

from pyworkflow.utils import glob, redStr
from pyworkflow.utils.path import makePath
from pyworkflow.protocol.params import (PointerParam, IntParam, FloatParam, StringParam,
                                        BooleanParam, LEVEL_ADVANCED)

import prody


class ProDyModesBase(EMProtocol):
    """
    This protocol acts as a base class for various kinds of mode analysis,
    providing easier access to the qualify and animate steps.
    Currently, the only child class is ProDyEdit.
    """

    """
    AI Generated Summary:

    ProDy Modes Base (ProDyModesBase) — User Manual

    OVERVIEW
    The ProDy Modes Base protocol serves as the fundamental architectural layer
    for Elastic Network Model (ENM) analysis. It centralizes the core logic
    required to compute and validate macromolecular dynamics, ensuring that
    child protocols—such as ProDyEdit—follow a consistent biophysical workflow.

    PHYSICAL PARAMETERS AND MODELING
    The protocol transforms structural data (PDB files or pseudoatoms from EM
    volumes) into a spring-node network.
    - Cut-off Distance: Defines the reach of atomic interactions (e.g., 15 Å
      for C-alpha).
    - Spring Constant (Gamma): Determines the stiffness of the connections.
    These parameters dictate the vibrational frequencies that describe how
    the protein naturally flexes and moves.

    MODE QUALIFICATION AND SCORING
    Not every mathematical mode is biologically relevant. This protocol
    implements a 'qualification' step:
    - Collectivity: Measures if a motion is global (affecting the whole protein)
      or local. Modes below a 0.15 threshold are typically filtered out.
    - Rigid-Body Filtering: The first six modes (translations and rotations)
      are deselected to focus strictly on internal structural deformations.
    - Scoring: Modes are ranked and assigned scores to help users prioritize
      which movements to use in downstream image analysis.

    ANIMATION AND DISPLACEMENT
    To provide a tangible understanding of protein dynamics, the protocol
    includes an animation engine:
    - RMSD Amplitude: Controls the scale of the visualized motion.
    - VMD Integration: Automatically generates scripts for external high-quality
      rendering.
    - Atom Shifts: Quantifies the displacement of each node, identifying
      flexible hinges and stable structural cores.

    OUTPUTS AND INTEGRATION
    The final output is a 'SetOfNormalModes', a standardized metadata package
    containing vectors, eigenvalues, and collectivity scores. This output
    is fully compatible with the ContinuousFlex framework, enabling
    advanced studies of structural heterogeneity and functional transitions.
    """
