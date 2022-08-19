=======================
Scipion ProDy plugin
=======================

This plugin provide a wrapper around `ProDy <https://github.com/prody/prody>`_ software: A Python Package for Protein Dynamics Analysis

Installation
-------------

You will need to use 3.0+ version of Scipion to be able to run these protocols. To install the plugin, you have one option so far:

Developer's version

   * download repository

    .. code-block::

        git clone https://github.com/scipion-em/scipion-em-prody.git

   * install

    .. code-block::

       scipion3 installp -p ./scipion-em-prody --devel

ProDy software will be installed automatically with the plugin but you can also use an another version 
by installing that one in your scipion3 environment.

**Important:** you need to have conda (miniconda3 or anaconda3) pre-installed to use this program.

Configuration variables
-----------------------
*CONDA_ACTIVATION_CMD*: If undefined, it will rely on conda command being in the
PATH (not recommended), which can lead to execution problems mixing scipion
python with conda ones. One example of this could can be seen below but
depending on your conda version and shell you will need something different:
CONDA_ACTIVATION_CMD = eval "$(/extra/miniconda3/bin/conda shell.bash hook)"


Protocols
----------

* ProDyANM for anisotropic network model (ANM) normal mode analysis
* ProDyRTB for rotation and translation blocks (RTB) normal mode analysis
* ProDyDefvec for calculating deformation vectors between two structures

* ProDyEdit to edit normal modes to different nodes/atoms, including slice and extend
* ProDyCompare for comparing motions between different calculations including from continuousflex

* ProDySelect for atom selection
* ProDyAlign for alignment two structures, including atom matching and superposition

* ProDyImportModes for importing modes calculated outside Scipion, including from Gromacs

* ProDyBuildPDBEnsemble for building ensembles from heterogeneous atomic structures
* ProDyImportEnsemble for importing ensembles calculated outside Scipion

* ProDyPCA for principal component analysis based on ensembles of atomic structures
* ProDyProject for projecting ensembles of structures onto principal components or normal modes

* ProDyGNM for Gaussian network model (GNM) analysis
* ProDyDomainDecomp for dynamical domain decomposition based on GNM modes

Viewers
----------

* ProDyComparisonsViewer for viewing matrices and bar graphs quantifying motion similarities
* ProDyModeViewer for viewing modes of motion in the VMD plugin NMWiz

* ProDyGNMViewer for viewing GNM mode shapes and covariance or cross-correlation matrices
* ProDyDomainViewer for viewing dynamical domains

* ProDyProjectionsViewer for viewing conformational landscapes from projecting ensembles onto modes

References
-----------

1. Zhang S, Krieger JM, Zhang Y, Kaya C, Kaynak B, Mikulska-Ruminska K, Doruker P, Li H, Bahar I (2021). ProDy 2.0: Increased scale and scope after 10 years of protein dynamics modelling with Python. Bioinformatics, btab187.
2. Bakan A, Meireles LM, Bahar I ProDy: Protein Dynamics Inferred from Theory and Experiments 2011 Bioinformatics 27(11):1575-1577
3. Bakan A, Dutta A, Mao W, Liu Y, Chennubhotla C, Lezon TR, Bahar I Evol and ProDy for Bridging Protein Sequence Evolution and Structural Dynamics 2014 Bioinformatics 30(18):2681-2683
