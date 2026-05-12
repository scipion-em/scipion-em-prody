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
This module will provide ProDy normal mode analysis (NMA) using the anisotropic network model (ANM).
"""
from multiprocessing import cpu_count

import prody
from prody2 import Plugin
from prody2.protocols.protocol_modes_base import ProDyModesBase

from pwem.emlib import (MetaData, MDL_NMA_MODEFILE, MDL_ORDER,
                        MDL_ENABLED, MDL_NMA_COLLECTIVITY, MDL_NMA_SCORE, 
                        MDL_NMA_ATOMSHIFT, MDL_NMA_EIGENVAL)
from pwem.objects import SetOfNormalModes, String

from pyworkflow.utils import glob, redStr
from pyworkflow.utils.path import makePath
from pyworkflow.protocol.params import (PointerParam, IntParam, FloatParam, StringParam,
                                        BooleanParam, LEVEL_ADVANCED)

vecStr = "vec.%d"

class ProDyANM(ProDyModesBase):
    """
    This protocol performs normal mode analysis (NMA) using the
    anisotropic network model (ANM).

    AI Generated:

    ANM NMA (ProDyANM) — User Manual
        Overview

        The ANM NMA protocol computes collective motions of a molecular
        structure using the Anisotropic Network Model (ANM). ANM is one of the
        most widely used coarse-grained approaches for studying large-scale
        structural dynamics in proteins, nucleic acids, and macromolecular
        assemblies.

        Rather than simulating atomic trajectories over time, ANM estimates the
        intrinsic directions in which a structure can move most easily around
        its equilibrium conformation. These motions often correspond to
        biologically meaningful conformational changes such as domain closure,
        hinge bending, subunit rearrangements, breathing motions, or ligand
        gating.

        For cryo-EM and structural biology users, ANM is especially useful when
        exploring functional flexibility, interpreting structural variability,
        generating candidate motions for flexible fitting, or selecting
        collective deformation coordinates for downstream analysis.

        Inputs and General Workflow

        The protocol requires a single input structure.

        This structure can be a conventional atomic model (for example a PDB
        file) or a pseudoatomic model derived from an EM density map.

        The protocol constructs an elastic network where nodes correspond to
        atoms or pseudoatoms, and springs connect nearby nodes. From this
        network, the Hessian matrix is built and diagonalized to obtain the
        normal modes.

        The resulting modes describe preferred directions of collective motion.

        Biological Interpretation of ANM

        ANM is best viewed as a model of intrinsic structural mechanics.

        Low-frequency non-zero modes often describe collective motions that are
        most relevant biologically. These are typically the modes associated
        with functional conformational changes.

        High-frequency modes generally represent more localized fluctuations and
        are often less informative when studying large-scale biological motion.

        Six zero modes correspond to rigid-body translations and rotations.
        These do not describe internal flexibility.

        In practice, biological interpretation usually focuses on the first few
        low-frequency non-zero modes.

        Number of Modes

        The Number of modes parameter controls how many modes are computed.

        For exploratory analysis, values around 10 to 20 are often sufficient.

        If the goal is to characterize broader conformational variability or to
        provide a richer basis for downstream flexible fitting, larger values
        may be useful.

        The theoretical maximum number of modes is three times the number of
        nodes.

        In most biological applications, computing very large numbers of modes
        rarely provides major practical benefit unless downstream dimensionality
        reduction or clustering is planned.

        Cutoff Distance

        The cutoff distance is one of the most important parameters in ANM.

        It determines which nodes interact through springs.

        Biologically, this defines the effective mechanical connectivity of the
        structure.

        For C-alpha models, the default value of 15 Å usually works well.

        Slightly larger values such as 18 Å may improve robustness in some
        proteins, especially elongated or multi-domain systems.

        For all-atom models, much smaller cutoffs such as 5–7 Å are generally
        more appropriate.

        For pseudoatomic models, the optimal value depends on the level of
        coarse-graining and particle density.

        If the cutoff is too small, the network may become poorly connected and
        modes may become unstable or fragmented.

        If the cutoff is too large, the model becomes overly rigid and may lose
        biologically meaningful flexibility.

        A practical biological strategy is to begin with the default and adjust
        only if the computed modes appear unphysical or overly localized.

        Spring Constant (Gamma)

        Gamma controls the stiffness of the elastic springs.

        In many biological applications, the default value of 1 is entirely
        sufficient because relative mode shapes matter more than absolute
        frequencies.

        More advanced users may introduce structure-dependent gamma functions
        when modeling specific physical hypotheses, but this is usually not
        necessary for standard exploratory structural analysis.

        Collectivity Threshold

        Collectivity is a particularly useful biological descriptor.

        It measures how broadly distributed a motion is across the structure.

        Modes with high collectivity involve large fractions of the molecule and
        often correspond to biologically relevant collective rearrangements.

        Modes with low collectivity tend to be more localized and may reflect
        local flexibility rather than global conformational change.

        The collectivity threshold allows automatic deselection of poorly
        collective modes.

        For many biological analyses, the default value provides a useful first
        filter.

        Setting the threshold to zero disables deselection entirely.

        Zero Modes

        The protocol can optionally retain zero eigenvalue modes.

        These correspond to rigid-body motions and generally do not provide
        information about internal structural flexibility.

        In most biological analyses, these are not of primary interest.

        However, keeping them may be useful for technical completeness or
        specialized downstream workflows.

        Sparse, KDTree, and Turbo Options

        These parameters mainly affect computational performance rather than
        biological interpretation.

        Sparse matrices reduce memory usage at the cost of speed.

        KDTree changes how neighbors are identified during network
        construction.

        Turbo mode uses a faster but more memory-intensive diagonalization
        strategy.

        For most users, the default settings are appropriate.

        Explicit Membrane Model

        For membrane proteins, an explicit membrane elastic network can be
        included.

        This is particularly relevant when the mechanical environment of the
        lipid bilayer strongly influences the dominant motions.

        Biologically, this can improve interpretation of channels,
        transporters, and membrane-associated assemblies.

        This option should only be used when the structure has already been
        properly oriented relative to the membrane.

        Animation and Visual Interpretation

        The protocol automatically generates animations of the computed modes.

        These animations are extremely useful for biological interpretation.

        They help reveal whether a mode corresponds to hinge closure, domain
        rotation, interface breathing, gate opening, or other collective
        rearrangements.

        RMSD amplitude controls the visual excursion along the mode.

        Larger amplitudes make motions easier to inspect but can exaggerate
        structural changes beyond realistic physical scales.

        Number of frames determines smoothness of the animation.

        Positive and negative directions simply explore both directions along
        the same mode vector.

        Outputs and Their Interpretation

        The main output is a SetOfNormalModes object.

        Each mode includes:

        - an eigenvector describing the direction of motion
        - an eigenvalue related to stiffness
        - collectivity information
        - metadata indicating whether the mode passed collectivity filtering

        The protocol also produces visualization files compatible with ProDy
        and ContinuousFlex viewers.

        These outputs can be used directly in downstream analyses such as mode
        comparison, deformation fitting, image analysis, or structural
        interpretation.

        Practical Recommendations

        For most biological systems, a good starting point is:

        - 10 to 20 modes
        - cutoff near 15 Å for C-alpha models
        - default collectivity filtering

        If modes appear fragmented or excessively localized, increasing the
        cutoff is often the first parameter worth testing.

        If the structure is a membrane protein, consider the membrane option
        only if the orientation is biologically meaningful.

        In practice, visual inspection of the first few non-zero collective
        modes usually provides the most biologically useful information.

        Final Perspective

        ANM does not attempt to reproduce exact physical trajectories.

        Instead, it identifies the easiest collective deformations allowed by
        the architecture of the structure.

        For structural biology users, this makes ANM especially powerful for
        understanding how molecular architecture constrains biological motion.

        The most reliable biological conclusions usually come from combining
        ANM with structural knowledge, biochemical context, and direct visual
        inspection of the dominant collective modes.
    """
    _label = 'ANM NMA'
    _possibleOutputs = {'outputModes': SetOfNormalModes}

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form, besidesAnimation=False):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        cpus = cpu_count()//2 # don't use everything
        form.addParallelSection(threads=cpus, mpi=0)

        # You need a params to belong to a section:
        form.addSection(label='ProDy ANM NMA')

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

        form.addParam('cutoff', StringParam, default=15.,
                      label="Cut-off distance (A)",
                      help='Atoms or pseudoatoms beyond this distance will not interact.\n'
                           'For Calpha atoms, the default distance of 15 A works well in the majority of cases although '
                           '18 A may sometimes be better, see Eyal et al., Bioinformatics 2006.\n'
                           'For all atoms, a shorter distance such as 5 or 7 A is recommended, see Tirion et al., Phys Rev Lett 1996.\n'
                           'For other levels of coarse-graining including pseudoatoms, see Doruker et al., J Comput Chem 2002.\n'
                           'It is also possible to use other functions for the cutoff e.g. 2.9 * math.log(numResidues) - 2.9 for ed-ENM, '
                           'replacing numResidues with the actual number of residues')
        form.addParam('gamma', StringParam, default=1.,
                      expertLevel=LEVEL_ADVANCED,
                      label="Spring constant",
                      help='This number or function determines the strength of the springs.\n'
                           'Besides pre-defined Gamma functions such as GammaStructureBased from Lezon et al., PLoS Comput Biol 2010 '
                           'and GammaED from Orellana et al., J Chem Theory Comput 2010, '
                           'more sophisticated options are available within the ProDy API and '
                           'the resulting modes can be imported back into Scipion.\n'
                           'See http://http://www.bahargroup.org/prody/tutorials/enm_analysis/gamma.html')
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
                      label="Include zero eigvals",
                      help='Elect whether modes with zero eigenvalues will be kept.')
        form.addParam('turbo', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Use turbo mode",
                      help='Elect whether to use a memory intensive, but faster way to calculate modes.')

        form.addSection(label='Animation')        
        form.addParam('rmsd', FloatParam, default=5,
                      label='RMSD Amplitude (A)',
                      help='Used only for animations of computed normal modes. '
                      'This is the maximal amplitude with which atoms or pseudoatoms are moved '
                      'along normal modes in the animations. \n')
        form.addParam('numSteps', IntParam, default=10,
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

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self, n=20, nzeros=6):
        # Insert processing steps

        # Link the input
        inputFn = self.inputStructure.get().getFileName()
        numModes = self.numberOfModes.get()

        self.gnm = False
        self.nzeros = 6 if self.zeros.get() else 0

        self._insertFunctionStep('computeModesStep', inputFn, numModes)
        self._insertFunctionStep('qualifyModesStep', numModes,
                                 self.collectivityThreshold.get())
        self._insertFunctionStep('animateModesStep', self.rmsd.get(), self.numSteps.get(),
                                 self.neg.get(), self.pos.get(), self.nzeros)
        self._insertFunctionStep('computeAtomShiftsStep', numModes, self.nzeros)
        self._insertFunctionStep('createOutputStep')

    def computeModesStep(self, inputFn='', n=20):
        """Compute ANM normal modes"""

        self.pdbFileName = self._getPath('atoms.pdb')
        self.atoms = prody.parsePDB(inputFn, alt='all')
        prody.writePDB(self.pdbFileName, self.atoms)

        if self.membrane.get():
            self.prefix = 'modes.exanm'
        else:
            self.prefix = 'modes.anm'
        filename = self.prefix + '.npz'

        args = '{0} -s "all" --altloc "all"  --hessian --export-scipion --npzmatrices ' \
            '--npz -o {1} -p {2} -n {3} -g {4} -c "{5}" -P {6}'.format(self.pdbFileName,
                                                                         self._getPath(),
                                                                         self.prefix, n,
                                                                         self.gamma.get(),
                                                                         self.cutoff.get(),
                                                                         self.numberOfThreads.get())

        if self.sparse.get():
            args += ' --sparse-hessian'

        if self.kdtree.get():
            args += ' --use-kdtree'

        if self.zeros.get():
            args += ' --zero-modes'
            self.startMode = 6
        else:
            self.startMode = 0

        if self.turbo.get():
            args += ' --turbo'

        if self.membrane.get():
            args += ' --membrane'

        self.runJob(Plugin.getProgram('anm'), args)
        self.outModes = prody.loadModel(self._getPath(filename))

    def qualifyModesStep(self, numberOfModes, collectivityThreshold=0.15, suffix=''):
        self._enterWorkingDir()

        fnVec = glob("modes/vec.*")

        if len(fnVec) < numberOfModes:
            msg = "There are only %d modes instead of %d. "
            msg += "Check the number of modes you asked to compute and/or consider increasing cut-off distance. "
            msg += "The maximum number of modes allowed by the method for ANM normal mode analysis is "
            msg += "3 times the number of nodes (atoms or pseudoatoms; %d). "
            self.warning(redStr(msg % (len(fnVec), numberOfModes, self.atoms.numAtoms()*3)))

        mdOut = MetaData()
        collectivityList = list(prody.calcCollectivity(self.outModes))
        eigvals = self.outModes.getEigvals()

        for n in range(len(fnVec)):
            collectivity = collectivityList[n]

            objId = mdOut.addObject()
            modefile = self._getPath("modes", vecStr % (n + 1))
            mdOut.setValue(MDL_NMA_MODEFILE, modefile, objId)
            mdOut.setValue(MDL_ORDER, int(n + 1), objId)

            eigval = eigvals[n]
            mdOut.setValue(MDL_NMA_EIGENVAL, eigval, objId)

            if eigval > prody.utilities.ZERO:
                mdOut.setValue(MDL_ENABLED, 1, objId)
            else:
                mdOut.setValue(MDL_ENABLED, -1, objId)

            mdOut.setValue(MDL_NMA_COLLECTIVITY, collectivity, objId)
            
            if collectivity < collectivityThreshold:
                mdOut.setValue(MDL_ENABLED, -1, objId)

        idxSorted = [i[0] for i in sorted(enumerate(collectivityList), key=lambda x: x[1], reverse=True)]

        score = []
        for _ in range(len(fnVec)):
            score.append(0)

        modeNum = []
        l = 0
        for k in range(len(fnVec)):
            modeNum.append(k)
            l += 1

        for i in range(len(fnVec)):
            score[idxSorted[i]] = idxSorted[i] + modeNum[i] + 2
        i = 0
        for objId in mdOut:
            score[i] = float(score[i]) / (2.0 * l)
            mdOut.setValue(MDL_NMA_SCORE, score[i], objId)
            i += 1
        mdOut.write("modes%s.xmd" % suffix)

        self._leaveWorkingDir()
        
        prody.writeScipionModes(self._getPath(), self.outModes, scores=score, only_sqlite=True,
                                collectivityThreshold=collectivityThreshold)

    def createOutputStep(self):
        fnSqlite = self._getPath('modes.sqlite')
        nmSet = SetOfNormalModes(filename=fnSqlite)
        nmSet._nmdFileName = String(self._getPath(self.prefix + '.nmd'))

        inputPdb = self.inputStructure.get()
        nmSet.setPdb(inputPdb)

        self._defineOutputs(outputModes=nmSet)
        self._defineSourceRelation(self.inputStructure, nmSet)

    def _summary(self):
        if not hasattr(self, 'outputModes'):
            summ = ['Output modes not ready yet']
        else:
            modes = prody.parseScipionModes(self.outputModes.getFileName())

            summ = ['*{0}* ANM modes calculated for *{1}* nodes'.format(
                    modes.numModes(), modes.numAtoms())]
        return summ

