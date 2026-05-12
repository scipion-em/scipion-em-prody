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
This module will provide ProDy principal component analysis (PCA) using atomic structures
"""

from multiprocessing import cpu_count

from pwem.emlib import (MetaData, MDL_NMA_MODEFILE, MDL_ORDER,
                        MDL_ENABLED, MDL_NMA_COLLECTIVITY, MDL_NMA_SCORE, 
                        MDL_NMA_EIGENVAL)
from pwem.objects import SetOfPrincipalComponents, String, AtomStruct, EMFile

from pyworkflow.utils import glob, redStr, copyFile
from pyworkflow.protocol.params import (MultiPointerParam, IntParam, FloatParam,
                                        BooleanParam, StringParam,
                                        LEVEL_ADVANCED)
from pyworkflow.object import Float

from prody2.protocols.protocol_modes_base import ProDyModesBase
from prody2.objects import replaceCoordsets, loadAndWriteEnsemble
from prody2.constants import PRODY_FRACT_VARS
from prody2 import Plugin

import prody

from prody2.objects import HAVE_CHEM, DcdMDSystem
if HAVE_CHEM:
    POINTER_CLASSES = 'SetOfAtomStructs, ProDyNpzEnsemble, DcdMDSystem'
else:
    POINTER_CLASSES = 'SetOfAtomStructs, ProDyNpzEnsemble'

class ProDyPCA(ProDyModesBase):
    """
    This protocol will perform ProDy principal component analysis (PCA) using atomic structures
    """

    """
    AI Generated Summary:

    Principal Component Analysis (ProDyPCA) — User Manual

    OVERVIEW
    The Principal Component Analysis protocol is the primary tool for reducing 
    the dimensionality of structural ensembles. It identifies the dominant 
    directions of motion—Principal Components—that account for the largest 
    variations in a set of structures, effectively filtering biological signal 
    from structural noise.

    INPUTS AND PRE-PROCESSING
    The protocol supports multi-source ensembles including atomic structures, 
    DCD trajectories, and ProDy ensembles. To capture internal dynamics rather 
    than global rotation, the system can either keep existing alignment or 
    perform an iterative superposition. Furthermore, the degeneracy parameter 
    allows the choice between using all coordinate sets or just the 
    representative first conformation of each structure.

    CORE ANALYSIS AND METRICS
    The protocol computes the covariance matrix of atomic positions to extract 
    fractional variance, which quantifies the percentage of total movement 
    explained by each mode. It also analyzes cross-correlation to understand 
    the coordinated movement between different residues or domains. Parallel 
    execution is optimized for performance, using multi-threaded processing 
    during the intensive covariance calculations.

    BIOLOGICAL FILTERING AND VALIDATION
    Collectivity serves as a key metric where high collectivity modes represent 
    coordinated domain movements and low collectivity indicates local 
    fluctuations. The protocol generates visual animations in both positive 
    and negative directions to confirm that mathematical components align with 
    plausible biological transitions. Selection strings are typically focused 
    on the protein backbone (C-alpha) to analyze global fold changes.

    OUTPUTS
    The results are packaged into a 'SetOfPrincipalComponents' object complete 
    with eigenvalues and fractional variance metadata. Additionally, the 
    protocol exports raw covariance and correlation matrices as external files 
    to facilitate advanced statistical validation and further research.
    """
