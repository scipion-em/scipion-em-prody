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
This module will provide the ClustENM(D) hybrid simulation method from ProDy, combining clustering, ENM NMA and MD.
"""

from multiprocessing import cpu_count
import os

from pwem.objects import AtomStruct, SetOfAtomStructs
from pwem.protocols import EMProtocol

import pyworkflow.object as pwobj
from pyworkflow.protocol.params import (IntParam, FloatParam, StringParam, BooleanParam,
                                        EnumParam, MultiPointerParam, LEVEL_ADVANCED)

import prody
from prody2.constants import ENSEMBLE_WEIGHTS
from prody2.objects import ProDyNpzEnsemble, TrajFrame
from prody2 import Plugin

IMP = 0
EXP = 1

from pyworkflow.utils import logger

class ProDyClustENM(EMProtocol):
    """
    This protocol will provide the ClustENM and ClustENMD hybrid simulation methods from ProDy, combining clustering, ENM NMA, minimisation and MD.
    """
    """
    Performs conformational sampling of biomolecular structures using 
    the ClustENM and ClustENMD methodologies implemented in ProDy. 
    The protocol combines Elastic Network Model (ENM) normal mode 
    analysis, clustering, energy minimization, and optional molecular 
    dynamics simulations to explore biologically relevant conformational 
    landscapes of proteins and macromolecular assemblies.

    AI Generated:

    ClustENM(D) (ProDyClustENM) — User Manual

        Overview

        The ClustENM(D) protocol is designed to generate and refine 
        alternative conformations of biomolecular structures by combining 
        coarse-grained normal mode analysis with iterative structural 
        sampling and molecular dynamics refinement. The method provides 
        an efficient way to explore large-scale collective motions that 
        are often associated with biological function, such as domain 
        rearrangements, hinge bending, loop movements, or allosteric 
        transitions.

        In structural biology workflows, this protocol is particularly 
        useful for studying conformational heterogeneity, generating 
        structural ensembles for flexible fitting, investigating protein 
        dynamics, or preparing candidate conformations for downstream 
        cryo-EM or molecular docking analyses. Unlike conventional MD 
        simulations alone, ClustENM(D) focuses sampling along collective 
        low-frequency motions predicted by ENM normal mode analysis, 
        allowing efficient exploration of biologically meaningful states 
        with reduced computational cost.

        Inputs and General Workflow

        The protocol requires one or more atomic structures as input. 
        These structures are typically protein models in PDB format and 
        act as the starting conformations for conformational sampling. 
        Each input structure is independently processed through iterative 
        generations of normal mode analysis, conformer generation, 
        clustering, and refinement.

        During each generation, the protocol computes low-frequency ENM 
        modes and generates new conformers by perturbing the structure 
        along collective motions. The resulting conformers are then 
        clustered to reduce redundancy and identify representative states. 
        Optionally, the protocol performs energy minimization and short 
        molecular dynamics simulations to refine generated conformations 
        and improve their physical realism.

        This iterative workflow allows the protocol to progressively 
        expand the conformational landscape while maintaining biologically 
        meaningful structural diversity.

        Conformational Sampling and Generational Expansion

        One of the central concepts of ClustENM(D) is the generation-based 
        exploration of conformational space. Each generation begins from 
        the representative conformers retained from the previous cycle. 
        New conformations are produced by displacing structures along 
        selected normal modes, with the average RMSD controlling the 
        magnitude of structural perturbation.

        The number of modes determines how many collective motions are 
        sampled during conformer generation. In most biological systems, 
        low-frequency modes capture the largest functional motions, and 
        using a small number of modes typically provides stable and 
        interpretable results. Excessively large numbers of modes may 
        introduce unrealistic local distortions rather than biologically 
        relevant collective movements.

        The number of conformers generated per structure controls the 
        breadth of sampling. Larger values increase structural diversity 
        but also raise computational cost. Similarly, increasing the 
        number of generations expands conformational exploration but may 
        lead to progressively less physically relevant states if sampling 
        becomes too aggressive.

        Clustering and Structural Diversity

        Clustering plays an essential role in controlling redundancy and 
        preserving representative conformational states. After each 
        generation, conformers are grouped according to structural 
        similarity, allowing the protocol to retain representative 
        structures while discarding highly redundant conformations.

        The protocol supports clustering either by specifying a maximum 
        number of clusters or by defining an RMSD threshold. The 
        maxclust strategy is generally more efficient for large datasets 
        or many generations, while RMSD threshold clustering provides 
        more direct structural control over ensemble diversity.

        From a biological perspective, clustering prevents oversampling 
        of nearly identical conformations and helps maintain a balanced 
        representation of distinct structural states. Choosing overly 
        permissive thresholds may retain excessive redundancy, whereas 
        very strict clustering can remove meaningful intermediate states.

        Elastic Network Model Parameters

        The ENM component of the protocol models the structure as a 
        network of interacting nodes connected by springs. The spring 
        constant determines the strength of interactions between residues, 
        while the cutoff distance defines which atoms are considered 
        connected in the elastic network.

        The default cutoff values generally perform well for globular 
        proteins and standard biomolecular systems. However, highly 
        elongated structures, membrane proteins, or flexible assemblies 
        may require parameter adjustment to better capture collective 
        dynamics.

        Advanced computational options such as sparse Hessian matrices, 
        KDTree construction, and turbo mode allow optimization of memory 
        usage and computational speed. These settings are particularly 
        relevant for large macromolecular complexes or high-throughput 
        ensemble generation workflows.

        Molecular Dynamics Refinement

        The protocol optionally integrates short molecular dynamics 
        simulations after conformer generation. This refinement stage 
        improves structural realism by relaxing steric clashes and 
        allowing local adaptation of the perturbed conformations.

        Simulations may be performed using either implicit or explicit 
        solvent environments. Implicit solvent simulations are generally 
        faster and computationally lighter, making them suitable for 
        exploratory ensemble generation. Explicit solvent simulations 
        provide more physically realistic environments but require 
        substantially greater computational resources.

        Additional parameters such as temperature, ionic strength, 
        minimization tolerance, and simulation length allow fine control 
        over the refinement process. In biological applications, moderate 
        simulation lengths are often sufficient to stabilize generated 
        conformers without drifting excessively from the intended ENM 
        perturbations.

        Outlier Detection and Ensemble Quality

        For implicit solvent simulations, the protocol can automatically 
        detect and exclude energetic outliers using modified z-score 
        analysis. This filtering step removes conformations with unusually 
        unfavorable potential energies that may correspond to unstable or 
        nonphysical structural states.

        Biologically, outlier filtering improves the quality and 
        interpretability of the final ensemble by reducing the presence 
        of unrealistic conformations. However, users should exercise 
        caution because highly flexible or partially unfolded states may 
        occasionally appear as energetic outliers despite having potential 
        biological relevance.

        Flexible Fitting to Experimental Volumes

        An important optional feature of the protocol is the ability to 
        perform filtering against experimental density maps. When enabled, 
        generated conformations are evaluated according to their agreement 
        with input volumes, similarly to approaches used in flexible 
        fitting workflows.

        This functionality is particularly useful in cryo-EM studies where 
        conformational ensembles need to remain compatible with 
        experimental density data. Structures that significantly reduce 
        cross-correlation with the target map can be filtered out, helping 
        guide sampling toward experimentally supported conformations.

        The protocol also supports iterative replacement of filtered 
        conformers, allowing the ensemble size to remain approximately 
        constant while enforcing map consistency.

        Outputs and Their Interpretation

        After execution, the protocol produces sets of refined atomic 
        structures corresponding to the generated conformational ensemble. 
        Each output structure is associated with ensemble weights derived 
        from the clustering and sampling procedure.

        In addition, the protocol generates compressed ProDy ensemble 
        files containing coordinate sets and metadata for all sampled 
        conformers. These ensembles can be reused for downstream 
        structural analysis, visualization, dimensionality reduction, or 
        comparison with experimental data.

        Biologically, the resulting ensemble should be interpreted as a 
        representation of accessible collective motions rather than an 
        exact thermodynamic distribution. The generated conformations are 
        most valuable for exploring plausible functional transitions and 
        identifying structurally meaningful dynamic states.

        Practical Recommendations

        For most biological applications, beginning with a small number 
        of modes and moderate RMSD perturbations provides the most stable 
        and interpretable results. Excessive perturbation amplitudes or 
        too many generations may generate unrealistic conformations that 
        deviate from experimentally plausible structures.

        Implicit solvent refinement is generally recommended for rapid 
        exploratory studies, while explicit solvent simulations are more 
        suitable for detailed structural refinement or publication-level 
        analyses. When studying large conformational transitions, combining 
        several generations with careful clustering often provides a good 
        balance between diversity and structural realism.

        When fitting against cryo-EM maps, users should carefully choose 
        map thresholds and simulated volume resolutions to avoid 
        overfitting noise or introducing artificial structural bias.

        Final Perspective

        ClustENM(D) provides an efficient hybrid framework for exploring 
        biomolecular flexibility by combining coarse-grained collective 
        motion analysis with atomistic refinement techniques. Rather than 
        relying solely on long-timescale molecular dynamics simulations, 
        the protocol guides sampling toward biologically relevant motions 
        predicted by elastic network theory.

        For structural biologists, this approach offers a practical way 
        to investigate conformational variability, generate flexible 
        structural ensembles, and bridge computational modeling with 
        experimental structural data such as cryo-EM density maps. 
        Careful parameter selection, balanced conformational sampling, 
        and biologically informed interpretation remain essential for 
        obtaining meaningful and reliable results.
    """