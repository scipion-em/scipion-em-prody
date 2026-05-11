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
    Performs Principal Component Analysis (PCA) on structural ensembles
    using ProDy.

    The protocol identifies the dominant collective structural variations
    present across an ensemble of conformations.

    AI Generated:

    ProDy Principal Component Analysis (ProDyPCA) — User Manual
        Overview

        The ProDyPCA protocol performs principal component analysis on
        structural ensembles.

        Its main purpose is to identify the dominant directions of
        structural variability present in a collection of conformations.

        In structural biology, PCA is commonly used to detect collective
        motions that naturally emerge from experimental ensembles,
        molecular dynamics trajectories, or collections of related
        structures.

        Unlike supervised methods such as logistic regression, PCA does
        not use predefined classes.

        Instead, it asks:

            "What are the dominant structural fluctuations sampled by
            the ensemble?"

        Input Data

        The protocol accepts structural ensembles provided as:

            - SetOfAtomStructs
            - ProDy ensemble objects
            - molecular dynamics trajectory systems (DcdMDSystem)

        All conformations must represent the same molecular system and
        contain equivalent atoms.

        If the input is a molecular dynamics trajectory, the protocol
        extracts:

            - the trajectory file
            - the associated reference structure

        If the input is a structural ensemble, it is first converted into
        the internal ProDy representation.

        Atom Selection

        The parameter:

            selstr

        determines which atoms are included in the PCA calculation.

        Recommended common choices are:

            - "all"
            - "name CA"

        Using alpha carbons often provides a robust description of
        large-scale collective protein motions while reducing noise.

        Number of Components

        The user specifies how many principal components to compute.

        The theoretical maximum is:

            number of conformations - 1

        This differs from normal mode analysis, where the upper limit is
        determined by the number of structural nodes.

        In practice, the first few principal components usually capture
        most biologically meaningful structural variability.

        Structural Alignment

        PCA is highly sensitive to structural alignment.

        The protocol provides two options:

            - keep the input alignment
            - realign the conformations before analysis

        Keep Alignment

        When alignment is preserved, the protocol assumes the input
        conformations are already in a common structural frame.

        This is appropriate when the ensemble has been carefully prepared
        beforehand.

        Realignment

        If alignment is not preserved, the protocol realigns the
        structures after trajectory generation.

        This is particularly important when translational or rotational
        differences would otherwise dominate the covariance matrix.

        Biologically, proper alignment is critical because PCA should
        capture internal conformational variability rather than rigid-body
        displacement.

        Covariance Analysis

        PCA is based on the covariance matrix of atomic displacements.

        During execution, the protocol computes this covariance matrix and
        derives the principal components from it.

        Each component represents an independent direction of structural
        variance.

        The associated eigenvalues quantify the amount of variance
        captured by each component.

        Fractional Variance

        For every principal component, the protocol computes the
        fractional variance.

        This indicates how much of the total structural variance is
        explained by each mode.

        Biologically, this helps identify which components dominate the
        ensemble dynamics.

        A few large fractional variances often indicate a relatively
        simple collective motion landscape.

        Cross-Correlation Matrix

        The protocol also computes a cross-correlation matrix between
        atomic displacements.

        This matrix describes how atomic motions are correlated across
        the ensemble.

        Positive correlations indicate atoms moving together.

        Negative correlations indicate atoms moving in opposite
        directions.

        This information is particularly useful for studying:

            - long-range coupling
            - domain communication
            - allosteric behavior

        Workflow

        The protocol performs the following steps:

            1. Load the ensemble or trajectory.
            2. Select the requested atoms.
            3. Compute the PCA model.
            4. Parse the generated principal components.
            5. Compute fractional variances.
            6. Compute cross-correlation matrices.
            7. Rank and qualify the resulting modes.
            8. Generate animations.
            9. Export outputs.

        Mode Qualification

        Each principal component is evaluated using:

            - collectivity
            - eigenvalue
            - ranking score
            - enable/disable flag

        Unlike elastic network normal mode analysis, PCA does not
        automatically exclude rigid-body modes.

        Components may optionally be filtered according to a
        collectivity threshold.

        By default, the threshold is zero because biologically important
        PCA modes may not always be highly collective.

        Animation

        The protocol automatically generates animations for the computed
        principal components.

        Animation parameters include:

            - RMSD amplitude
            - number of frames
            - positive direction
            - negative direction

        These animations provide a visual representation of the structural
        displacement associated with each principal component.

        This is often the most intuitive way to interpret the biological
        meaning of a component.

        Output Data

        The protocol generates a set of principal components together
        with several associated outputs.

        Main outputs include:

            - outputModes
            - optional aligned outputEnsemble
            - covariance matrix file
            - cross-correlation matrix file

        Each principal component is linked to:

            - eigenvalue
            - collectivity
            - score
            - fractional variance

        The components are also associated with the average reference
        structure.

        Biological Interpretation

        The biological meaning of PCA is fundamentally different from
        energy-based normal mode analysis.

        PCA does not predict possible motions.

        Instead, it describes motions that are actually sampled in the
        structural ensemble.

        This makes PCA especially powerful for studying:

            - experimentally observed heterogeneity
            - molecular dynamics trajectories
            - conformational continua
            - dominant collective fluctuations

        Practical Recommendations

        In most structural biology applications, the first few principal
        components contain the most interpretable motions.

        It is generally useful to inspect together:

            - fractional variance
            - cross-correlation
            - animations

        This combined interpretation often reveals whether structural
        variability reflects:

            - domain motion
            - hinge bending
            - flexible loops
            - collective rearrangements

        Final Perspective

        ProDyPCA is best understood as an unsupervised structural
        dimensionality reduction method.

        Rather than asking:

            "Which motions are theoretically accessible?"

        it asks:

            "Which motions are actually sampled by the ensemble?"

        This makes it especially useful when the scientific goal is to
        characterize experimentally observed or simulated structural
        variability.
    """
    _label = 'PCA'
    _possibleOutputs = {'outputModes': SetOfPrincipalComponents}
    _nmdFileName = 'modes.pca.nmd'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form, besidesAnimation=False):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        cpus = cpu_count()//2 # don't use everything
        form.addParallelSection(threads=cpus, mpi=0)

        form.addSection(label='ProDy PCA')
        form.addParam('inputEnsemble', MultiPointerParam, label="Input ensemble",
                      important=True,
                      pointerClass=POINTER_CLASSES,
                      help='The input ensemble should be a SetOfAtomStructs '
                      'where all structures have the same number of atoms or a ProDy ensemble.')
        form.addParam('degeneracy', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Take only first conformation from each structure/set",
                      help='Elect whether only the active coordinate set (**True**) or all the coordinate sets '
                           '(**False**) of each structure should be added to the ensemble. Default is **True**.')
        form.addParam('numberOfModes', IntParam, default=5,
                      label='Number of modes',
                      help='The maximum number of modes allowed by the method for '
                           'atomic normal mode analysis is 3 times the '
                           'number of nodes (Calpha atoms or pseudoatoms).')
        form.addParam('collectivityThreshold', FloatParam, default=0, # important modes may well not be collective
                      expertLevel=LEVEL_ADVANCED,
                      label='Threshold on collectivity',
                      help='Collectivity degree is related to the number of atoms or pseudoatoms that are affected by '
                      'the mode, and it is normalized between 0 and 1. Modes below this threshold are deselected in '
                      'the modes metadata file as these modes are much less collective. \n'
                      'For no deselection, this parameter should be set to 0 . \n')
        form.addParam('selstr', StringParam, default="all",
                      label="Selection string",
                      help='Selection string for atoms to include in the calculation.\n'
                           'It is recommended to use "all" (default) or "name CA"')
        form.addParam('keepAlignment', BooleanParam, default=True,
                      label="Keep alignment", help="The alternative is to realign the structures")

        form.addSection(label='Animation')        
        form.addParam('rmsd', FloatParam, default=2,
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

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        n = self.numberOfModes.get()

        self.gnm = False
        nzeros = 0

        self._insertFunctionStep('computeModesStep', n)
        self._insertFunctionStep('qualifyModesStep', n,
                                 self.collectivityThreshold.get())
        self._insertFunctionStep('computeAtomShiftsStep', n, nzeros)
        self._insertFunctionStep('animateModesStep', self.rmsd.get(), self.n_steps.get(),
                                 self.neg.get(), self.pos.get(), 0)
        self._insertFunctionStep('createOutputStep')

    def computeModesStep(self, n=5):
        if (len(self.inputEnsemble)==1 and
            isinstance(self.inputEnsemble[0].get(), DcdMDSystem)):
                self.npz = None
                system = self.inputEnsemble[0].get()
                self.dcdFileName = system.getTrajectoryFile()

                self.pdbFileName = self._getPath('atoms.pdb')
                copyFile(system.getSystemFile(), self.pdbFileName)
                self.averageStructure = AtomStruct()
                self.averageStructure.setFileName(self.pdbFileName)
        else:
            loadAndWriteEnsemble(self) # creates self.npz, self.dcdFileName, self.pdbFileName and others

        args = '{0} --pdb {1} -s "{2}" ' \
               '--covariance --export-scipion --npz --npzmatrices' \
               ' -o {3} -p modes.pca -n {4} -P {5}'.format(self.dcdFileName,
                                                           self.pdbFileName,
                                                           self.selstr.get(),
                                                           self._getPath(), n,
                                                           self.numberOfThreads.get())
        if self.keepAlignment.get():
            args += " --aligned"

        self.runJob(Plugin.getProgram('pca'), args)
        
        self.outModes, self.atoms = prody.parseNMD(self._getPath(self._nmdFileName),
                                                   type=prody.PCA)
        
        crossCorr = prody.calcCrossCorr(self.outModes)
        prody.writeArray(self._getPath('modes.pca_crossCorr.txt'), crossCorr)

        if not self.keepAlignment.get():
            dcdEnsemble = prody.parseDCD(self._getPath('ensemble.dcd'))
            dcdEnsemble.iterpose()

            if self.npz is not None:
                self.npz2 = replaceCoordsets(self.npz, dcdEnsemble.getCoordsets(),
                                            suffix='_aligned', iterpose=False,
                                            coords=dcdEnsemble.getCoords())
            else:
                self.npz2 = None
        else:
            self.npz2 = self.npz
        
        self.fract_vars = prody.calcFractVariance(self.outModes)
        prody.writeArray(self._getPath('pca_fract_vars.txt'), self.fract_vars)

    def qualifyModesStep(self, numberOfModes, collectivityThreshold=0, suffix=None):
        self._enterWorkingDir()

        fnVec = glob("modes/vec.*")

        if len(fnVec) < numberOfModes:
            msg = "There are only %d modes instead of %d. "
            msg += "Check the number of modes you asked to compute and/or consider increasing cut-off distance. "
            msg += "The maximum number of modes allowed by the method for atomic principal component analysis is "
            msg += "the number of structures - 1 (%d). "
            self.warning(redStr(msg % (len(fnVec), numberOfModes, self.ens.numConfs())))

        mdOut = MetaData()
        collectivityList = list(prody.calcCollectivity(self.outModes))
        eigvals = self.outModes.getEigvals()

        for n in range(len(fnVec)):
            collectivity = collectivityList[n]

            objId = mdOut.addObject()
            modefile = self._getPath("modes", "vec.%d" % (n + 1))
            mdOut.setValue(MDL_NMA_MODEFILE, modefile, objId)
            mdOut.setValue(MDL_ORDER, int(n + 1), objId)

            mdOut.setValue(MDL_ENABLED, 1, objId)
            mdOut.setValue(MDL_NMA_COLLECTIVITY, collectivity, objId)
            mdOut.setValue(MDL_NMA_EIGENVAL, eigvals[n], objId)

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
        mdOut.write("modes.xmd")

        self._leaveWorkingDir()
        
        prody.writeScipionModes(self._getPath(), self.outModes, scores=score, only_sqlite=True,
                                collectivityThreshold=collectivityThreshold)


    def createOutputStep(self):
        fnSqlite = self._getPath('modes.sqlite')
        nmSet = SetOfPrincipalComponents(filename=fnSqlite)
        nmSet._nmdFileName = String(self._getPath(self._nmdFileName))

        self.fractVarsDict = {}
        for i, item in enumerate(nmSet):
            self.fractVarsDict[item.getObjId()] = self.fract_vars[i]

        outSet = SetOfPrincipalComponents().create(self._getPath())
        outSet.copyItems(nmSet, updateItemCallback=self._setFractVars)
        outSet._nmdFileName = String(self._getPath(self._nmdFileName))

        inputPdb = self.averageStructure
        self._defineOutputs(refPdb=inputPdb)
        outSet.setPdb(inputPdb)

        self._defineOutputs(outputModes=outSet)
        self._defineSourceRelation(inputPdb, outSet)

        if self.npz2 is not None:
            self._defineOutputs(outputEnsemble=self.npz2)

        outputMatrixCov = EMFile(filename=self._getExtraPath('modes.pca_covariance.txt'))
        outputMatrixCrosCor = EMFile(filename=self._getExtraPath('modes.pca_crossCorr.txt'))
        self._defineOutputs(matrixFileCC=outputMatrixCrosCor,
                            matrixFileCV=outputMatrixCov)

    def _summary(self):
        if not hasattr(self, 'outputModes'):
            summ = ['Output modes not ready yet']
        else:
            modes = prody.parseScipionModes(self.outputModes.getFileName())

            if hasattr(self, 'outputEnsemble'):
                ens = self.outputEnsemble.loadEnsemble()

                summ = ['*{0}* principal components calculated from *{1}* structures of *{2}* atoms'.format(
                        modes.numModes(), ens.numConfs(), ens.numAtoms())]
            else:
                summ = ['*{0}* principal components calculated'.format(modes.numModes())]
        return summ

    def _setFractVars(self, item, row=None):
        # We provide data directly so don't need a row
        fractVar = Float(self.fractVarsDict[item.getObjId()])
        setattr(item, PRODY_FRACT_VARS, fractVar)
