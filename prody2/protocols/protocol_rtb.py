# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     James Krieger (jamesmkrieger@gmail.com)
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
This module will provide ProDy normal mode analysis (NMA) using the the rotation and translation of blocks (RTB) framework.
"""
from pwem.objects import SetOfNormalModes, String, AtomStruct
from pyworkflow.protocol.params import (PointerParam, IntParam, FloatParam,
                                        BooleanParam, EnumParam, LEVEL_ADVANCED)

from prody2.protocols.protocol_modes_base import ProDyModesBase
from prody2 import Plugin, copyConvertPDB

BLOCKS_FROM_RES = 0
BLOCKS_FROM_SECSTR = 1

class ProDyRTB(ProDyModesBase):
    """
    Performs normal mode analysis using the Rotation Translation of
    Blocks (RTB) framework.

    The protocol approximates collective structural motions by grouping
    atoms into blocks and computing their coupled rigid-body motions.

    AI Generated:

    ProDy RTB Normal Mode Analysis (ProDyRTB) — User Manual
        Overview

        The ProDyRTB protocol performs normal mode analysis using the
        Rotation Translation of Blocks (RTB) formalism.

        Its main purpose is to reduce the computational cost of normal
        mode analysis while preserving the biologically meaningful
        collective motions of large macromolecular systems.

        Instead of treating every atom independently, RTB groups atoms
        into blocks and models the motion of those blocks as rigid-body
        translations and rotations.

        This makes the protocol especially useful for:

            - large proteins
            - multi-domain complexes
            - coarse-grained structural models
            - pseudoatomic EM models

        Biological Motivation

        In many biological systems, large-scale functional motions often
        involve coordinated movement of structural regions rather than
        isolated atomic fluctuations.

        RTB exploits this idea by representing groups of residues as
        collective moving units.

        This allows efficient approximation of motions such as:

            - domain rearrangements
            - hinge bending
            - subunit displacement
            - collective conformational transitions

        Input Structure

        The protocol requires one input atomic structure.

        The input can be:

            - a standard atomic model
            - a pseudoatomic model

        The structure is loaded with secondary structure information
        whenever available.

        Block Definition

        A key feature of RTB is how the structural blocks are defined.

        The protocol supports two strategies.

        Residue-Based Blocks

        Blocks can be defined by assigning a fixed number of residues to
        each block.

        The parameter:

            res_per_block

        controls the target number of residues per block.

        This is a simple and robust option when the user wants uniform
        coarse-graining.

        Secondary Structure Blocks

        Blocks can also be defined using secondary structure
        information.

        In this case, blocks tend to follow biologically meaningful
        structural elements such as:

            - helices
            - beta strands
            - compact structural segments

        This can often provide more interpretable collective motions.

        Block Refinement

        Several parameters control block refinement.

        Shortest Block

        Very small blocks can be unstable.

        Blocks shorter than the selected threshold are merged with the
        previous block.

        Longest Block

        Very large blocks may be overly rigid.

        Blocks longer than the threshold are split.

        Distance-Based Splitting

        The parameter:

            min_dist_cutoff

        allows block subdivision based on internal structural distance.

        Residues that are too far apart are not forced into the same
        block.

        Biologically, this helps avoid grouping disconnected structural
        regions into one rigid body.

        Elastic Network Parameters

        Cutoff Distance

        The cutoff defines which block interactions are included in the
        elastic network.

        For most alpha-carbon models, the default value of 15 Å is often
        appropriate.

        Shorter values may be preferable for denser atomic models.

        Spring Constant

        The spring constant controls the strength of inter-block
        coupling.

        This defines the stiffness of the RTB elastic network.

        Computational Workflow

        The protocol performs the following main steps.

        Structure Preparation

        The input structure is loaded and converted into a block-mapped
        representation.

        The resulting block assignment is stored together with an atom
        mapping structure.

        Hessian Construction

        The RTB Hessian matrix is constructed using the selected blocks
        and elastic network parameters.

        If memory limitations occur, the protocol automatically switches
        to sparse matrix representation.

        Mode Calculation

        Normal modes are computed from the RTB Hessian.

        The user specifies the desired number of modes.

        Two additional options control this stage.

        Zero Eigenvalues

        The user may decide whether zero-eigenvalue modes are retained.

        When enabled, the first six rigid-body modes are preserved.

        Turbo Mode

        Turbo mode uses a faster but more memory-intensive matrix
        decomposition.

        If memory becomes limiting, the protocol automatically falls back
        to a non-turbo calculation.

        Output Files

        The protocol exports the computed RTB modes in several formats:

            - Scipion mode files
            - NMD visualization file
            - NPZ model file including matrices

        These outputs allow visualization and downstream structural
        analysis.

        Mode Qualification

        Each computed mode is evaluated according to:

            - collectivity
            - eigenvalue
            - ranking score
            - enabled/disabled state

        If zero modes are included, the first six rigid-body modes are
        automatically disabled for interpretation.

        Modes below the selected collectivity threshold are also
        deselected.

        This filtering helps focus on biologically meaningful collective
        motions.

        Animation

        The protocol automatically generates animations of the computed
        RTB modes.

        Animation parameters include:

            - RMSD amplitude
            - number of frames
            - positive direction
            - negative direction

        These animations are especially useful because RTB often captures
        large-amplitude domain-scale rearrangements that are easy to
        interpret visually.

        Atom Shift Profiles

        The protocol also computes atom displacement profiles for the
        selected modes.

        This identifies which regions of the structure undergo the
        largest displacements.

        Biologically, this helps detect:

            - flexible hinges
            - moving domains
            - collective deformation hotspots

        Output

        The final output is a:

            - SetOfNormalModes

        The output modes remain linked to the original input structure.

        This allows downstream interpretation in structural and
        conformational analyses.

        Biological Interpretation

        RTB should be understood as a coarse-grained approximation of
        normal mode analysis.

        It is particularly useful when the biological question concerns
        large-scale collective motion rather than local atomic detail.

        Compared with fully atomistic normal mode analysis, RTB often
        provides:

            - faster computation
            - better scalability
            - clearer interpretation of domain-level motion

        Practical Recommendations

        Residue-based blocks are usually a good starting point for
        exploratory analyses.

        Secondary-structure-based blocks may be preferable when the user
        wants motions that align more closely with biologically defined
        structural elements.

        For very large systems, RTB is often much more practical than
        full atomistic normal mode analysis.

        Final Perspective

        ProDyRTB is best understood as a scalable collective-motion
        approximation.

        Rather than asking:

            "How does every atom fluctuate?"

        it asks:

            "How do coherent structural blocks move relative to one
            another?"

        This makes it especially powerful for studying large
        macromolecular rearrangements.
    """
    _label = 'RTB NMA'
    _possibleOutputs = {'outputModes': SetOfNormalModes}

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form, besidesAnimation=False):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label='ProDy RTB NMA')

        form.addParam('inputStructure', PointerParam, label="Input structure",
                      important=True,
                      pointerClass='AtomStruct',
                      help='The input structure can be an atomic model '
                           '(true PDB) or a pseudoatomic model\n'
                           '(an EM volume converted into pseudoatoms)')

        form.addParam('numberOfModes', IntParam, default=20,
                      label='Number of modes',
                      help='The maximum number of modes allowed by the method for '
                           'atomic normal mode analysis is 3 times the '
                           'number of nodes (Calpha atoms or pseudoatoms).')

        form.addParam('blockDef', EnumParam, choices=['res', 'secstr'],
                      label="Block definition type",
                      default=BLOCKS_FROM_RES,
                      display=EnumParam.DISPLAY_HLIST,
                      help='Define blocks using either a number of residues or secondary structure information')

        form.addParam('res_per_block', IntParam, default=10,
                      condition='blockDef==%d' % BLOCKS_FROM_RES,
                      label="Number of residues per block",
                      help='All blocks will have this number of residues except the last one')

        form.addParam('shortest_block', IntParam, default=4,
                      expertLevel=LEVEL_ADVANCED,
                      label='Number of residues in shortest block',
                      help='Blocks with fewer residues will be combined into the previous block. '
                           'Fewer than 4 can be problematic.')

        form.addParam('longest_block', IntParam, default=20,
                      expertLevel=LEVEL_ADVANCED,
                      label='Number of residues in longest block',
                      help='Blocks with more residues will be split in half')

        form.addParam('min_dist_cutoff', FloatParam, default=20.,
                      expertLevel=LEVEL_ADVANCED,
                      label='Distance cutoff for splitting blocks',
                      help='Distance of a residue from others beyond which '
                           'it is not included in the same block based on a distance tree. '
                           'This is calculated using ProDy function findSubgroups.')

        form.addParam('cutoff', FloatParam, default=15.,
                      label="Cut-off distance (A)",
                      help='Atoms or pseudoatoms beyond this distance will not interact.\n'
                           'For Calpha atoms, the default distance of 15 A works well in the majority of cases. '
                           'For all atoms, a shorter distance such as 5 or 7 A is recommended.\n'
                           'For fewer atoms or pseudoatoms, set this according to the level of coarse-graining '
                           '(see Doruker et al., J Comput Chem 2002 though values may differ for RTB).')

        form.addParam('gamma', FloatParam, default=1.,
                      expertLevel=LEVEL_ADVANCED,
                      label="Spring constant",
                      help='This number or function determines the strength of the springs.\n'
                           'More sophisticated options are available within the ProDy API and '
                           'the resulting modes can be imported back into Scipion.\n'
                           'See http://http://www.bahargroup.org/prody/tutorials/enm_analysis/gamma.html')

        form.addParam('collectivityThreshold', FloatParam, default=0.15,
                      expertLevel=LEVEL_ADVANCED,
                      label='Threshold on collectivity',
                      help='Collectivity degree is related to the number of atoms or pseudoatoms that are affected by '
                      'the mode, and it is normalized between 0 and 1. Modes below this threshold are deselected in '
                      'the modes metadata file as these modes are much less collective. \n'
                      'For no deselection, this parameter should be set to 0 . \n'
                      'Zero modes 1-6 are always deselected as they are related to rigid-body movements. \n'
                      'The modes metadata file can be used to see which modes are more collective '
                      'in order to decide which modes to use at the image analysis step.')

        form.addParam('zeros', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Include zero eigvals",
                      help='Elect whether modes with zero eigenvalues will be kept.')

        form.addParam('turbo', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Use turbo mode",
                      help='Elect whether to use a memory intensive, but faster way to calculate modes.')

        form.addParam('sparse', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Use sparse matrices?",
                      help='This saves memory at the expense of computational time.')
        form.addParam('kdtree', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Use KDTree for building Hessian matrix?",
                      help='This takes more computational time.')
        form.addParam('membrane', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Use explicit membrane model?",
                      help='An explicit lattice elastic network is used to model the membrane. '
                      'This option requires a protein oriented with opm or ppm.')

        form.addSection(label='Animation')        
        form.addParam('rmsd', FloatParam, default=5,
                      label='RMSD Amplitude (A)',
                      help='Used only for animations of computed normal modes. '
                      'This is the maximal amplitude with which atoms or pseudoatoms are moved '
                      'along normal modes in the animations. \n')
        form.addParam('n_steps', IntParam, default=10,
                      expertLevel=LEVEL_ADVANCED,
                      label='Number of frames',
                      help='Number of frames used in each direction of animations.')
        form.addParam('pos', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Include positive direction",
                      help='Elect whether to animate in the positive mode direction.')
        form.addParam('neg', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Include negative direction",
                      help='Elect whether to animate in the negative mode direction.')

        form.addParam('registerAnimations', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Register animation pdbs as outputs",
                      help='Elect whether to register multi-state pdbs from animations as outputs.')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps

        # Link the input
        inputFn = self.inputStructure.get().getFileName()
        numModes = self.numberOfModes.get()
        self.gnm = False

        self.nzeros = 6 if self.zeros.get() else 0

        self._insertFunctionStep('computeModesStep', inputFn, numModes)
        self._insertFunctionStep('animateModesStep', self.rmsd.get(), self.n_steps.get(),
                                 self.neg.get(), self.pos.get(), self.nzeros)
        self._insertFunctionStep('qualifyModesStep', numModes,
                                 self.collectivityThreshold.get())
        self._insertFunctionStep('computeAtomShiftsStep', numModes, self.nzeros)
        self._insertFunctionStep('createOutputStep')

    def computeModesStep(self, inputFn='', n=20):       
        self.atomsFn = self._getPath('atoms.pdb')
        copyConvertPDB(inputFn, self.atomsFn)

        args = '{0} -s "all" --altloc "all"  --hessian --export-scipion --npzmatrices ' \
            '--npz -o {1} -p {2} -n {3} -g {4} -c "{5}" -P {6}'.format(self.atomsFn,
                self._getPath(), self.getPrefix(), n, self.gamma.get(),
                self.cutoff.get(), self.numberOfThreads.get())

        if self.sparse.get():
            args += ' --sparse-hessian'

        if self.kdtree.get():
            args += ' --use-kdtree'

        if self.zeros.get():
            args += ' --zero-modes'

        if self.turbo.get():
            args += ' --turbo'

        if self.blockDef.get() == BLOCKS_FROM_RES:
            args += ' --block-input-type 1 --res-per-block {0}'.format(self.res_per_block.get())
        else:
            args += ' --block-input-type 2'

        args += f' --res-per-block {self.res_per_block.get()} --shortest-block {self.shortest_block.get()}'
        args += f' --longest-block {self.longest_block.get()} --min-block-dist-cutoff {self.min_dist_cutoff.get()}'

        self.runJob(Plugin.getProgram('rtb'), args)

    def createOutputStep(self):
        fnSqlite = self._getPath('modes.sqlite')
        nmSet = SetOfNormalModes(filename=fnSqlite)
        nmSet._nmdFileName = String(self._getPath('modes.rtb.nmd'))

        inputPdb = self.inputStructure.get()
        nmSet.setPdb(inputPdb)

        self._defineOutputs(outputModes=nmSet)
        self._defineSourceRelation(self.inputStructure, nmSet)

        if self.registerAnimations.get():
            args = {}
            for i in range(len(nmSet)):
                name = "animation" + str(i+1)
                args[name] = AtomStruct(self.getAnimationPdbPath(i))

            self._defineOutputs(**args)

    def getPrefix(self):
        return 'modes.rtb'
