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

       scipion installp -p path_to_scipion-em-prody --devel

ProDy software will be installed automatically with the plugin but you can also use an existing installation by providing *PRODY_ENV_ACTIVATION* (see below).

**Important:** you need to have conda (miniconda3 or anaconda3) pre-installed to use this program.

Configuration variables
-----------------------
*CONDA_ACTIVATION_CMD*: If undefined, it will rely on conda command being in the
PATH (not recommended), which can lead to execution problems mixing scipion
python with conda ones. One example of this could can be seen below but
depending on your conda version and shell you will need something different:
CONDA_ACTIVATION_CMD = eval "$(/extra/miniconda3/bin/conda shell.bash hook)"

*PRODY_ENV_ACTIVATION* (default = conda activate prody-github):
Command to activate the ProDy environment.


Protocols
----------

* ProDy ANM for anisotropic network model normal mode analysis

References
-----------

1. Zhang S, Krieger JM, Zhang Y, Kaya C, Kaynak B, Mikulska-Ruminska K, Doruker P, Li H, Bahar I (2021). ProDy 2.0: Increased scale and scope after 10 years of protein dynamics modelling with Python. Bioinformatics, btab187.
2. Bakan A, Meireles LM, Bahar I ProDy: Protein Dynamics Inferred from Theory and Experiments 2011 Bioinformatics 27(11):1575-1577
3. Bakan A, Dutta A, Mao W, Liu Y, Chennubhotla C, Lezon TR, Bahar I Evol and ProDy for Bridging Protein Sequence Evolution and Structural Dynamics 2014 Bioinformatics 30(18):2681-2683
