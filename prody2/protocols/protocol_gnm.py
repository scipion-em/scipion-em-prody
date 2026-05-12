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
This module will provide ProDy normal mode analysis (NMA) using the Gaussian network model (GNM).
"""
from os.path import exists, join
import math

from pwem.emlib import (MetaData, MDL_NMA_MODEFILE, MDL_ORDER,
                        MDL_ENABLED, MDL_NMA_COLLECTIVITY, MDL_NMA_SCORE, 
                        MDL_NMA_ATOMSHIFT, MDL_NMA_EIGENVAL)
from pwem.objects import String, EMFile
from pwem.protocols import EMProtocol

from pyworkflow.utils import glob, redStr
from pyworkflow.utils.path import makePath
from pyworkflow.protocol.params import (PointerParam, IntParam, FloatParam, StringParam,
                                        BooleanParam, LEVEL_ADVANCED)

import prody
from prody2.objects import SetOfGnmModes
from prody2 import Plugin

class ProDyGNM(EMProtocol):
    """
    This protocol will perform normal mode analysis (NMA) using the Gaussian network model (GNM)
    """
class ProDyGNM(EMProtocol):
    """
    Performs normal mode analysis (NMA) on atomic or pseudoatomic structures using
the Gaussian Network Model (GNM).

    AI Generated:

    GNM Analysis (ProDyGNM) — User Manual

        Overview

        The ProDyGNM protocol performs normal mode analysis on molecular structures
using the Gaussian Network Model. Its main goal is to analyze intrinsic collective
motions of proteins and assemblies, providing insights into flexibility, conformational
changes, and functional dynamics. The protocol works with either atomic models or
pseudoatomic representations derived from EM volumes, supporting both exploratory
and publication-level structural studies.

        Inputs and General Workflow

        Users provide an input structure (atomic or pseudoatomic) and select the number
of modes to compute. Optional parameters include a cut-off distance for interactions,
spring constant for the network, explicit membrane modeling, collectivity threshold
for mode selection, and inclusion of zero eigenvalues. The protocol is organized
into sequential steps: computation of normal modes, qualification of modes based on
collectivity, calculation of atomic displacements, and generation of outputs.

        Mode Computation

        Normal modes are computed using the ProDy GNM or ExGNM model. Input structures
are parsed and written to temporary PDB files. Depending on the membrane setting,
modes are labeled accordingly. Users can include or exclude zero eigenvalue modes.
The resulting mode covariance and cross-correlation matrices are saved for further
analysis.

        Mode Qualification

        Modes are evaluated based on collectivity and eigenvalues. Low-collectivity or
zero-eigenvalue modes are deselected. Scores are assigned to modes according to
their collectivity, allowing prioritization of biologically meaningful movements.
Metadata files store mode properties for downstream use.

        Atomic Shifts and Displacements

        The protocol calculates maximum atomic displacements across modes, recording
which mode contributes most to each atomic shift. Distance profiles are saved for
each mode, facilitating detailed structural interpretation and identification of
dynamic hotspots.

        Outputs

        Final outputs include a set of GNM modes (SQLite database), covariance and
cross-correlation matrices, and metadata describing mode properties and atomic shifts.
The protocol links the outputs to the original input structure, ensuring traceability
and reproducibility in downstream analyses.

        Practical Recommendations

        For routine analysis, selecting alpha-carbon nodes with default cut-off and
spring parameters is robust. Membrane modeling is recommended when analyzing
transmembrane proteins. Collectivity thresholds help filter non-meaningful modes.
Careful inspection of modes and atomic shifts supports biologically relevant
interpretation of protein dynamics.

        Final Perspective

        ProDyGNM integrates computational NMA with metadata management to provide
a comprehensive view of protein flexibility. Proper configuration of input
parameters, mode selection, and post-analysis interpretation is essential for
extracting meaningful insights into conformational dynamics and functional motions.
    """
