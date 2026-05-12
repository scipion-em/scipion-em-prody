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
This module will provide the ProDy interface for parsing files from 
the BioExcel CV19 database
"""
from pwem.protocols import EMProtocol
from pyworkflow.protocol import params
from pyworkflow.object import Float

import prody

from prody2 import Plugin
from prody2.constants import ENSEMBLE_WEIGHTS
from prody2.objects import SetOfTrajFrames, TrajFrame

from prody2.objects import HAVE_CHEM
if HAVE_CHEM:
    from prody2.objects import DcdMDSystem as ProDyMDSystem
else:
    ProDyMDSystem = SetOfTrajFrames


class ProDyBioExcelCV19(EMProtocol):
    """
    This module will provide the ProDy interface for parsing files from 
    the BioExcel CV19 database
    """
    """
    Downloads and parses molecular dynamics trajectories from the BioExcel COVID-19 database.
    The protocol provides an automated interface between ProDy and BioExcel simulation repositories,
    allowing users to retrieve trajectory datasets together with their associated topology and
    structural information for downstream structural biology and molecular dynamics analysis.

    AI Generated:

    BioExcel CV19 Parser (ProDyBioExcelCV19) — User Manual

        Overview

        The ProDyBioExcelCV19 protocol is designed to retrieve and prepare molecular dynamics
        simulations hosted in the BioExcel COVID-19 database. Its primary objective is to simplify
        the access to curated MD trajectories associated with biologically relevant systems,
        particularly viral proteins and biomolecular complexes studied during the COVID-19 pandemic.

        In practical workflows, this protocol enables researchers to import simulation trajectories
        directly into Scipion and ProDy-based environments without manually downloading,
        organizing, or converting trajectory files. The protocol automatically generates the
        structural topology, coordinate files, and trajectory representations needed for
        downstream analysis such as conformational exploration, flexibility studies,
        principal component analysis, or ensemble characterization.

        For biological users, this protocol is especially useful when studying protein dynamics,
        ligand interactions, conformational transitions, or large-scale structural fluctuations
        extracted from publicly available simulation datasets. It provides a reproducible and
        automated mechanism for integrating BioExcel resources into structural biology workflows.

        Inputs and Simulation Retrieval

        The protocol requires a simulation accession identifier corresponding to an entry in the
        BioExcel COVID-19 database. This accession uniquely defines the simulation dataset to
        retrieve and determines which trajectory, topology, and coordinate files will be parsed.

        During execution, the protocol launches the BioExcel parsing utility and stores all
        retrieved files inside the protocol working directory. The generated files include a PDB
        structure file, a PSF topology file, and a DCD trajectory file. Together, these files
        describe both the static molecular architecture and the temporal evolution of the system.

        From a biological perspective, the accession number represents a complete simulation
        experiment. Different accessions may correspond to distinct proteins, variants,
        ligand-binding conditions, solvent environments, or simulation protocols. Careful
        selection of the dataset is therefore essential for biologically meaningful analysis.

        Atom Selection and Structural Reduction

        One of the most important features of this protocol is the ability to reduce the trajectory
        to specific atom subsets before importing it into the workflow environment. This greatly
        decreases computational cost and simplifies downstream analyses.

        The protocol supports several atom selection modes, including carbon-only selections,
        backbone atoms, or backbone carbon atoms. These reduced representations are commonly used
        in large-scale conformational studies where the global protein motion is more relevant than
        atomistic side-chain detail.

        Selecting only backbone atoms is particularly useful for coarse-grained analyses,
        normal mode calculations, trajectory clustering, or essential dynamics studies.
        Carbon-alpha selections are frequently employed when comparing conformational landscapes
        across multiple trajectories because they provide a compact representation of protein
        geometry while preserving the major structural transitions.

        From a computational perspective, reduced atom selections significantly improve memory
        efficiency and analysis speed, especially when working with long trajectories or large
        molecular assemblies.

        Frame Selection and Temporal Sampling

        The protocol also allows selective extraction of trajectory frames. Instead of loading
        the entire simulation, users may specify frame ranges or intervals using compact range
        expressions such as "1-5" or "10:20:2".

        This functionality is biologically important because molecular dynamics trajectories often
        contain thousands or millions of frames, many of which may be redundant for a particular
        analysis. Temporal subsampling enables efficient exploration of conformational states
        while reducing storage and processing requirements.

        In practical applications, users commonly extract representative frame intervals for
        visualization, ensemble averaging, conformational clustering, or flexible fitting.
        Sparse sampling may also be useful for identifying slow collective motions or monitoring
        large conformational transitions over time.

        Care should nevertheless be taken when selecting only a subset of frames, since excessive
        reduction may remove transient intermediate states or rare biologically relevant events.

        Output Formats and Trajectory Representation

        After processing, the protocol generates a trajectory object that can be represented in
        two different ways depending on the selected configuration.

        When MDSystem output is enabled, the protocol creates a ProDyMDSystem object that stores
        the structural coordinates, topology, and trajectory in a unified representation. This
        format is particularly suitable for downstream molecular dynamics analysis pipelines,
        including flexibility studies, ensemble calculations, and structural motion analysis.

        Alternatively, the protocol can generate a SetOfTrajFrames object, where each frame of
        the trajectory is represented independently. This representation is useful for workflows
        focused on frame-by-frame processing, classification, scoring, or visualization.

        Each generated frame is associated with a normalized ensemble weight, allowing the
        trajectory to be interpreted as a statistical ensemble of conformational states.

        Biological Interpretation of the Outputs

        The resulting trajectories provide a dynamic description of biomolecular behavior rather
        than a single static structure. Unlike crystallographic or cryo-EM models, molecular
        dynamics trajectories capture fluctuations, transitions, and transient interactions that
        may be essential for biological function.

        Protein flexibility, domain rearrangements, loop dynamics, and ligand-binding events can
        often be studied more effectively through trajectory analysis than through static models
        alone. For viral systems, these motions may reveal mechanisms related to infectivity,
        immune escape, or drug recognition.

        The protocol summary reports the total number of atoms together with the number of protein
        residues contained in the parsed trajectory. These values provide a quick estimate of the
        system size and complexity before downstream processing begins.

        Practical Recommendations

        In routine biological workflows, backbone or carbon-alpha selections are generally the
        best starting point for exploratory analyses because they substantially reduce
        computational cost while preserving the dominant structural motions.

        Full-atom trajectories are more appropriate when studying detailed intermolecular
        interactions, ligand binding, hydrogen-bond networks, or local conformational effects.
        However, these trajectories require significantly larger computational resources.

        Frame subsampling should be selected carefully depending on the biological question.
        Dense frame sampling is preferable for kinetic or transition-state analyses, whereas
        sparse sampling is often sufficient for visualization or global conformational studies.

        When working with large MD datasets, it is usually advisable to begin with reduced atom
        selections and limited frame ranges before scaling to full trajectories.

        Final Perspective

        Molecular dynamics simulations provide a time-resolved view of biomolecular systems that
        complements experimental structural biology techniques. The ProDyBioExcelCV19 protocol
        simplifies access to these datasets by integrating BioExcel simulation repositories
        directly into ProDy and Scipion workflows.

        For most structural biology users, the protocol serves as a bridge between public
        molecular simulation databases and advanced conformational analysis tools, enabling
        efficient exploration of protein flexibility, structural heterogeneity, and dynamic
        biological mechanisms.
    """