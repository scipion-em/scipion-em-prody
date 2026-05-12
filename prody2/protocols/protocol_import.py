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
This module will provide ProDy mode import tools.
"""
import os
import numpy as np

from pwem.objects import (String, AtomStruct, SetOfAtomStructs, EMFile,
                          SetOfNormalModes, SetOfPrincipalComponents)
from pwem.protocols import ProtImportFiles

from prody2.objects import (ProDyNpzEnsemble, TrajFrame,
                            SetOfGnmModes, SetOfLogisticModes)
from prody2.constants import ENSEMBLE_WEIGHTS

import pyworkflow.object as pwobj
import pyworkflow.protocol.params as params
from pyworkflow.utils import logger

import prody
from prody.dynamics.gnm import ZERO

NMD = 0
MODES_NPZ = 1
SCIPION = 2
GROMACS = 3

PDB_FILENAME = 'atoms.pdb'
PSF_FILENAME = 'atoms.psf'
DCD_FILENAME = 'ensemble.dcd'

filesPatternHelp = """Pattern of the files to be imported.\n\n
The pattern can contain standard wildcards such as\n
*, ?, etc, or special ones like ### to mark some\n
digits in the filename as ID.\n\n
NOTE: wildcards and special characters 
('*', '?', '#', ':', '%') cannot appear in the actual path.\n\n
For gromacs modes, the first is for values and this is for vectors"""

class ProDyImportModes(ProtImportFiles):
    """
    This protocol will import ProDy modes to a SetOfNormalModes object
    """
    """
    AI Generated: ProDy Import Manual

    OVERVIEW
    This protocol serves as a bridge for importing structural dynamics data 
    into Scipion. It handles both Normal Modes (vibrational analysis) and 
    Ensembles (conformational sets), allowing for a seamless transition 
    from static structures to dynamic biological insights.

    IMPORTING MODES
    The protocol transforms mathematical vectors from ProDy or Gromacs 
    into a 'SetOfNormalModes'. This is essential for studying protein 
    flexibility. It supports several formats:
    - Native NMD and NPZ files.
    - Gromacs directories (eigenvalues and eigenvectors).
    - Scipion-style directories.
    Because modes describe movement, they must be linked to an 'Input 
    Structure' (PDB or pseudoatoms) to map the flexibility onto the 
    physical biological model.

    IMPORTING ENSEMBLES
    The ensemble protocol manages collections of protein structures, such 
    as those from Molecular Dynamics or NMR.
    - Standardization: It can focus on specific regions using 'Selection 
      Strings' (e.g., 'protein' or 'name CA') to eliminate noise from 
      solvent or flexible loops.
    - Spatial Alignment: It offers 'Superpose' and 'Iterpose' options 
      to remove global rotation and translation, ensuring that observed 
      differences represent true internal conformational changes.
    - Flexibility: Data can be loaded via file patterns or by pointing 
      directly to existing Scipion objects like trajectories or atom sets.

    OUTPUTS AND BEST PRACTICES
    The system generates standardized Scipion objects (PCA, GNM, LRA, or 
    NMA modes) and high-performance NPZ ensembles. For best results, 
    users should ensure that the reference structure matches the atom 
    count of the imported ensemble. Using C-alpha selections is 
    recommended for large datasets to maintain efficiency without 
    losing the primary biological signal.
    """
