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
This module will provide ProDy linear discriminant analysis (LRA) using atomic structures
"""
from collections import OrderedDict
import numpy as np

from pwem.objects import Float, String
from pyworkflow.utils import getListFromRangeString
from pyworkflow.protocol.params import (MultiPointerParam, IntParam, FloatParam,
                                        BooleanParam, StringParam, TextParam, 
                                        NumericRangeParam, 
                                        LEVEL_ADVANCED, Float)

from prody2.protocols.protocol_modes_base import ProDyModesBase
from prody2.objects import SetOfLogisticModes, loadAndWriteEnsemble
from prody2.constants import PRODY_FRACT_VARS
from prody2 import parseMatchDict

import prody


class ProDyLRA(ProDyModesBase):
    """
    This protocol will perform ProDy logistic regression analysis (LRA) using atomic structures
    """
    """
    This protocol performs Logistic Regression Analysis (LRA) using ProDy to identify 
    structural components that best discriminate between two different functional 
    or conformational states.

    AI Generated:

    Logistic Regression Analysis (ProDyLRA) — User Manual
        Overview

        The Logistic Regression Analysis protocol is designed to identify the collective 
        structural motions that characterize the difference between two distinct groups 
        of structures. By applying a supervised learning approach to structural 
        ensembles, LRA finds "modes" that are optimized to separate defined biological 
        classes, such as "open" vs "closed" states or "ligand-bound" vs "apo" forms.

        For a biological researcher, this tool is invaluable for moving beyond simple 
        Principal Component Analysis (PCA). While PCA finds directions of maximum 
        variance, LRA finds the specific directions that are most biologically 
        relevant to the functional transition being studied.

        Inputs and General Workflow

        The protocol requires one or more input ensembles, which can be provided as 
        atomic structures, ProDy NPZ files, or DCD trajectories. A critical 
        requirement is that the user must define exactly two classes for the 
        analysis. If more or fewer than two classes are detected in the labeling 
        dictionary, the protocol will signal a validation error.

        The workflow involves labeling each structure in the ensemble with a class 
        identifier. The algorithm then calculates the logistic regression components 
        that maximize the separation between these classes. To ensure the results are 
        statistically robust and not due to random variation, the protocol includes 
        a shuffling parameter to assess the significance of the findings.

        Class Labeling and Customization

        Labeling is the core biological driver of this protocol. Users can provide 
        a custom class label dictionary to group their structural data. This allows 
        the user to explicitly tell the software which structures belong to which 
        biological state.

        Advanced users can manipulate the label dictionary by inserting custom 
        labels at specific indices or recovering specific label numbers. This 
        flexibility is essential when dealing with complex datasets where structures 
        may be interleaved or require manual sorting into functional categories.

        Atoms Selection and Degeneracy

        To focus the analysis on relevant structural features, a selection string 
        is provided. By default, "name CA" (Alpha Carbons) is used, which is 
        sufficient for most protein systems to capture global conformational 
        changes while reducing computational noise.

        The degeneracy option allows users to decide whether to use only the 
        first conformation of each structure or all available coordinate sets. 
        This is particularly useful when importing multi-model PDB files or 
        trajectories where some frames may be redundant or represent the 
        same equilibrium state.

        Animation and Visualization

        Once the LRA modes are calculated, the protocol allows for their 
        visualization through animations. These animations depict how the 
        structure moves along the discriminative components.

        The user can control the RMSD Amplitude to define how far the atoms 
        move in the animation, and set the number of frames to ensure a smooth 
        transition. Options to include both positive and negative directions 
        allow for a full view of the structural transition between the two 
        defined biological states.

        Outputs and Their Interpretation

        After execution, the protocol produces a SetOfLogisticModes. This 
        output includes the calculated components, which can be explored in 
        standard NMD viewers. The summary provides a clear count of how many 
        LRA components were calculated relative to the number of structures 
        and atoms analyzed.

        The protocol also outputs the processed ensemble and the reference 
        structure used for the calculation. This ensures that the results 
        are perfectly mapped back to the physical model, allowing for a 
        direct biological interpretation of which residues are key to 
        the conformational switch.

        Practical Recommendations

        For optimal results, users should ensure their ensembles are well-aligned 
        before running LRA. It is highly recommended to start with a clean CA-only 
        selection to identify the main hinges and domains involved in the 
        transition. If the "Number of Shuffles" shows that the components are 
        easily replicated by random labeling, the user should re-evaluate the 
        biological consistency of their class assignments.

        Final Perspective

        Logistic Regression Analysis turns structural data into a diagnostic 
        tool. It doesn't just ask "how does this molecule move?", but rather 
        "which specific movements define the difference between these two states?". 
        By providing a clear mathematical bridge between classification and 
        molecular dynamics, LRA helps researchers pinpoint the physical 
        basis of biological regulation.
    """