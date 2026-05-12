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
This module will provide ProDy mode algebra tools (linear combinations).
"""
import os
import numpy as np

from pwem.objects import SetOfNormalModes, String, Integer, CsvList

from pyworkflow.utils import glob, logger
from pyworkflow.protocol.params import (PointerParam, EnumParam, IntParam,
                                        StringParam, LEVEL_ADVANCED)

import prody
from prody2.protocols.protocol_modes_base import ProDyModesBase

COEFF_POINTER = 0
COEFF_STRING = 1

class ProDyAlgebra(ProDyModesBase):
    """
    This protocol will add together components from a SetOfNormalModes object 
    with coefficients based on overlaps or user input
    """

    """
    Performs linear algebra operations on a set of normal modes to generate
    combined mode vectors based on user-provided or calculated coefficients.

    AI Generated:

    ProDy Algebra (ProDyAlgebra) — User Manual

        Overview

        The ProDy Algebra protocol allows users to combine normal modes
        obtained from atomic or pseudoatomic models into new composite
        vectors. Each input mode is scaled by a corresponding coefficient,
        which can be provided either as a numeric string or as an external
        file. This process is biologically meaningful when exploring
        collective motions of macromolecules, analyzing principal
        components of structural variability, or testing hypotheses about
        conformational changes.

        Inputs and General Workflow

        The protocol requires a SetOfNormalModes as input. This set may
        originate from atomic PDB models, pseudoatomic representations
        derived from EM volumes, or principal component analysis of
        atomic ensembles. Each mode in the set represents a directional
        deformation of the structure.

        Users must specify coefficients that determine how each mode
        contributes to the final vector. Coefficients can be entered
        directly as a comma- or space-separated string, or imported
        from a file. The number of coefficients applied may be limited
        to a subset of the total modes, allowing selective combination
        of dominant motions.

        The protocol parses the input modes and associated atomic
        coordinates, scales each mode by its coefficient, and sums them
        to generate a composite vector. This vector is then converted
        into a new Normal Mode Analysis object suitable for further
        analysis or visualization.

        Outputs and Interpretation

        After execution, the protocol produces a set of combined modes
        saved in Scipion format and a corresponding NMD file containing
        the atomic coordinates and mode vectors. The computed coefficients
        are stored alongside the output to ensure reproducibility and
        traceability.

        Biologically, the resulting vectors can be used to explore
        functional motions, simulate conformational transitions, or
        serve as input for downstream modeling and analysis protocols.
        The integrity of coefficients and alignment with the atomic
        structure are critical for meaningful interpretations.

        Practical Recommendations

        For routine use, it is advisable to verify the correspondence
        between input modes and coefficients and to ensure that the
        number of coefficients matches the number of modes intended
        for combination. Visual inspection of the resulting modes in
        molecular visualization tools is recommended to confirm that
        the generated vectors reflect biologically plausible motions.

        When coefficients are derived from experimental data, care
        should be taken to normalize or scale values appropriately
        to prevent overrepresentation of certain modes. For exploratory
        analyses, selecting a subset of dominant modes can simplify
        interpretation and highlight major structural transitions.

        Final Perspective

        The ProDy Algebra protocol provides a flexible and powerful
        framework for mode combination in structural biology. By
        enabling custom linear combinations of normal modes, it allows
        researchers to probe conformational landscapes, generate
        biologically relevant motions, and prepare inputs for
        advanced modeling and analysis pipelines.
    """