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
from multiprocessing import cpu_count

import prody
from prody2 import Plugin
from prody2.protocols.protocol_modes_base import ProDyModesBase

from pwem.emlib import (MetaData, MDL_NMA_MODEFILE, MDL_ORDER,
                        MDL_ENABLED, MDL_NMA_COLLECTIVITY, MDL_NMA_SCORE, 
                        MDL_NMA_ATOMSHIFT, MDL_NMA_EIGENVAL)
from pwem.objects import SetOfNormalModes, String

from pyworkflow.utils import glob, redStr
from pyworkflow.utils.path import makePath
from pyworkflow.protocol.params import (PointerParam, IntParam, FloatParam, StringParam,
                                        BooleanParam, LEVEL_ADVANCED)

vecStr = "vec.%d"

class ProDyANM(ProDyModesBase):
    """
    This protocol will perform normal mode analysis (NMA) using the anisotropic network model (ANM)
    """
    """
    Performs normal mode analysis (NMA) using the anisotropic network model (ANM) to explore
    collective motions of macromolecular structures.

    AI Generated:

    ANM Normal Mode Analysis (ProDyANM) — User Manual

        Overview

        The ProDy ANM protocol computes normal modes for atomic or pseudoatomic models
        using the anisotropic network model. This approach captures intrinsic motions
        and flexibility within biomolecules, allowing the study of functional dynamics,
        conformational changes, or principal component analysis of structural ensembles.

        Inputs and General Workflow

        The protocol requires a single input structure, which can be a PDB atomic model
        or a pseudoatomic model derived from an EM volume. Users specify the number of
        modes to calculate, cutoff distances for atomic interactions, and spring constants
        governing the network. Optional advanced parameters include sparse matrices,
        KDTree optimization, explicit membrane modeling, and turbo computation mode
        for enhanced performance.

        ANM modes are computed through sequential steps. The protocol first calculates
        the Hessian matrix and derives mode vectors. Modes are then evaluated for
        collectivity, with less collective modes optionally deselected based on a user-defined
        threshold. Animations of selected modes can be generated with configurable RMSD,
        number of frames, and positive or negative direction motion.

        Outputs and Interpretation

        Upon completion, the protocol produces a set of normal modes saved in Scipion
        format alongside a metadata file containing eigenvalues, collectivity scores,
        and mode enablement. Each mode captures a distinct collective motion of the
        input structure, providing insights into flexibility and conformational transitions.

        Practical Recommendations

        For routine use, it is recommended to verify input parameters such as the number
        of modes and cutoff distance. Visual inspection of mode animations can help
        validate biological relevance. Advanced options, including membrane modeling
        or sparse computation, may be leveraged for specialized cases. Users should
        ensure that computed modes are consistent with expected physical behavior
        and biological function.

        Final Perspective

        The ProDy ANM protocol integrates rigorous computational methods with user-friendly
        configuration to facilitate structural dynamics analysis. It provides a robust
        framework for exploring macromolecular flexibility, interpreting functional
        motions, and generating biologically meaningful insights that support downstream
        modeling, simulation, or comparative studies.
    """