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
This module will provide ProDy distance and angle measurement for structural ensembles
"""
import numpy as np

from pwem.objects import SetOfAtomStructs
from pwem.protocols import EMProtocol

import pyworkflow.object as pwobj
from pyworkflow.protocol import params

import prody
from prody2.constants import MEASURES

DISTANCE = 0
ANGLE = 1
DIHEDRAL = 2

selstrHelp = '''The distance, angle or dihedral will be calculated between the centers of 2, 3 or 4 selections.
There is a rich selection engine with similarities to VMD. 
See http://http://www.bahargroup.org/prody/tutorials/prody_tutorial/selection.html'''

defaultSelstr = "protein and name CA or nucleic and name P C4' C2"

class ProDyMeasure(EMProtocol):
    """
    This module will provide ProDy distance and angle measurement for structural ensembles
    """

"""
Calculates geometric properties such as distances, angles, and dihedral angles across
structural ensembles. This protocol allows for the tracking of structural changes
and conformational dynamics by measuring specific relationships between defined atom selections.

```
AI Generated:

Structural Measurement (ProDyMeasure) — User Manual
    Overview

    The ProDy Measure protocol provides a robust way to quantify geometric relationships 
    within one or more structural ensembles. Its main purpose is to transform complex 3D 
    conformational changes into discrete numerical data—distances, angles, or dihedrals—so 
    they can be compared, plotted, or used for statistical analysis. In cryo-EM and 
    structural biology workflows, this step is vital for characterizing the range of 
    motion in flexible complexes or verifying specific functional states.

    For a biological user, the most common applications include monitoring the distance 
    between two domains during a catalytic cycle, measuring the hinge angle of a 
    molecular motor, or tracking the dihedral rotation of a specific side chain across 
    an MD trajectory or a set of reconstructed maps converted to pseudoatoms.

    Inputs and General Workflow

    The protocol accepts multiple input ensembles, which can be provided as sets of atomic 
    structures or ProDy NPZ files. A fundamental requirement is that all structures 
    within the input ensembles must share the same number of atoms to ensure consistency 
    during the measurement loop.

    The workflow centers on defining "selections" that represent the points of interest. 
    The protocol calculates the geometric center of each selection for every frame in 
    the ensemble. Once these centers are established, the protocol applies the 
    appropriate geometric formula based on the selected measurement type.

    Measurement Types and Selections

    The protocol offers three distinct modes of geometric analysis, each requiring a 
    different number of point selections:

    Distance measurement is the most straightforward, requiring two selection strings. 
    It calculates the Euclidean distance between the centers of mass of selection 1 
    and selection 2. This is ideal for tracking the "opening" or "closing" of a pocket.

    Angle measurement requires three selections. It calculates the angle formed at the 
    vertex (selection 2) by the lines connecting to selections 1 and 3. This is useful 
    for quantifying hinge-like movements in multi-domain proteins.

    Dihedral measurement is the most complex, requiring four selections. It calculates 
    the torsion angle between the planes defined by points (1,2,3) and (2,3,4). This is 
    the standard way to measure twisting motions or bond rotations.

    Atom Selection Logic

    Selection strings are powered by the ProDy selection language, allowing users to 
    target specific residues, atom names, or chains. From a biological perspective, 
    selecting a robust set of atoms (like C-alpha atoms of a stable helix) to define a 
    center is generally preferred over selecting a single atom, as it provides a more 
    statistically stable representation of a domain's position.

    If the ensemble contains structural discrepancies, the protocol includes an 
    automated "trimming" step to ensure that the atom selections remain valid across 
    all coordinate sets, preventing crashes during center-of-mass calculations.

    Outputs and Data Interpretation

    Upon completion, the protocol generates new ensembles that are copies of the 
    inputs but augmented with a new attribute: the calculated measure for each 
    individual structure. These values are also exported to CSV files, formatted 
    for easy import into external plotting or spreadsheet software.

    Biologically, these outputs allow the user to correlate structural identifiers 
    with geometric values. For instance, one can identify which subset of an 
    ensemble exhibits a "closed" conformation by filtering based on the distance 
    attribute now attached to the output structures.

    Practical Recommendations

    In routine practice, it is best to verify selection strings on a single 
    representative PDB file before running the protocol on a large ensemble. 
    If the resulting distances or angles appear erratic, consider broadening 
    the selection string to include more stable backbone atoms, which reduces 
    noise caused by local side-chain fluctuations.

    When measuring complex motions, it is often helpful to run multiple instances 
    of the protocol—for example, combining distance and angle measurements—to 
    obtain a multi-dimensional view of the protein's conformational landscape.

    Final Perspective

    For most researchers, measuring a distance or an angle is the first step toward 
    turning a visual observation into a publication-quality statistic. By automating 
    this across entire ensembles, the ProDy Measure protocol bridges the gap between 
    qualitative structural inspection and quantitative biophysical characterization.
"""

