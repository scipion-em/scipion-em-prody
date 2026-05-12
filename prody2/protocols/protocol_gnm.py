# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     James Krieger (jmkrieger@cnb.csic.es)
# *              Ricardo Serrano Gutiérrez (rserranogut@hotmail.com)                 
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
This module will provide ProDy normal mode analysis (NMA) using the Gaussian network model (GNM).
"""
from os.path import exists, join
import math

from pwem.emlib import (MetaData, MDL_NMA_MODEFILE, MDL_ORDER,
                        MDL_ENABLED, MDL_NMA_COLLECTIVITY, MDL_NMA_SCORE, 
                        MDL_NMA_ATOMSHIFT, MDL_NMA_EIGENVAL)
from pwem.objects import String, EMFile
from pwem.protocols import EMProtocol

from pyworkflow.utils import glob, redStr
from pyworkflow.utils.path import makePath
from pyworkflow.protocol.params import (PointerParam, IntParam, FloatParam, StringParam,
                                        BooleanParam, LEVEL_ADVANCED)

import prody
from prody2.objects import SetOfGnmModes
from prody2 import Plugin

class ProDyGNM(EMProtocol):
    """
    Performs Gaussian Network Model (GNM) normal mode analysis on an atomic
    structure or pseudoatomic model.

    AI Generated:

    GNM Analysis (ProDyGNM) — User Manual
        Overview

        The GNM Analysis protocol performs normal mode analysis using the
        Gaussian Network Model (GNM), a coarse-grained elastic network
        approach widely used to characterize collective motions in proteins
        and other macromolecular assemblies.

        In structural biology, GNM is commonly used to identify intrinsic
        flexibility patterns encoded by the native structure. Rather than
        simulating explicit time evolution, the method analyzes how the
        topology of interatomic contacts gives rise to preferred collective
        fluctuations.

        For biological users, this protocol is particularly useful when
        studying domain motions, flexible regions, hinge behavior, or
        identifying collective motions that may be functionally relevant.

        Inputs and General Workflow

        The protocol requires a single input structure.

        The input may be either:

            - a standard atomic structure (for example, a PDB model)
            - a pseudoatomic representation derived from volumetric data

        The workflow follows these steps:

            1. Read the input structure.
            2. Build the GNM elastic network.
            3. Compute normal modes.
            4. Evaluate collectivity and eigenvalues.
            5. Generate covariance and cross-correlation matrices.
            6. Estimate atom-wise displacement profiles.

        This produces both mode metadata and auxiliary outputs useful for
        interpretation and downstream analysis.

        Number of Modes

        The number of modes defines how many normal modes are computed.

        In GNM, low-frequency modes usually capture large-scale collective
        motions that are often biologically meaningful.

        In practical biological interpretation:

            - low-order modes often correspond to global collective motions
            - higher-order modes often describe more localized fluctuations

        Computing too many modes is not always necessary. In many practical
        analyses, a moderate number of modes is sufficient to characterize
        dominant structural flexibility.

        Cutoff Distance

        The cutoff defines which atoms or pseudoatoms interact in the elastic
        network.

        This parameter is biologically important because it controls network
        connectivity.

            - shorter cutoffs produce more local interactions
            - larger cutoffs produce more global connectivity

        For Cα-based protein models, values around the default range usually
        work well. For pseudoatomic models or sparse systems, somewhat larger
        cutoffs may be needed to maintain meaningful connectivity.

        If the cutoff is too small, the network may become fragmented or may
        produce fewer usable modes than expected.

        Spring Constant

        The spring constant controls interaction strength between connected
        nodes.

        In most biological applications, the absolute value is less important
        than the relative fluctuation patterns between residues.

        Therefore, the default value is generally sufficient for exploratory
        analyses unless a specific calibrated elastic network model is being
        used.

        Membrane-Aware GNM

        The protocol optionally supports an explicit membrane elastic network.

        This mode is intended for membrane proteins that have already been
        oriented consistently relative to the membrane, for example using OPM
        or PPM orientations.

        Biologically, this can improve the realism of the fluctuation model
        because membrane constraints often strongly affect collective motions
        of transmembrane assemblies.

        This option should generally be used only when membrane orientation
        is structurally meaningful.

        Zero Eigenvalue Modes

        The protocol can optionally keep modes with zero eigenvalues.

        In elastic network analysis, zero modes usually correspond to
        trivial rigid-body motions rather than internal conformational
        flexibility.

        For most biological interpretation, users are typically more
        interested in non-zero internal modes.

        However, retaining zero modes may still be useful for technical
        inspection or advanced downstream analyses.

        Mode Collectivity

        One of the most biologically useful outputs is mode collectivity.

        Collectivity measures how broadly motion is distributed across the
        structure.

            - high collectivity:
              many atoms participate in the motion

            - low collectivity:
              motion is localized to fewer atoms

        The collectivity threshold allows automatic deselection of highly
        localized modes.

        In practice, this is useful because biologically relevant global
        motions are often more collective than highly localized fluctuations.

        The protocol records collectivity values for all modes and uses them
        to rank and annotate the resulting mode metadata.

        Covariance and Cross-Correlation Matrices

        After mode calculation, the protocol computes:

            - covariance matrix
            - normalized cross-correlation matrix

        These matrices are extremely valuable for biological interpretation.

        Covariance reflects the magnitude of coupled fluctuations.

        Cross-correlation reveals whether regions move:

            - together (positive correlation)
            - oppositely (negative correlation)
            - independently (near zero correlation)

        In proteins, correlated motions often help identify dynamic domains,
        communication pathways, or long-range allosteric coupling.

        Atom Shift Profiles

        The protocol also computes atom-wise displacement amplitudes across
        modes.

        For each atom or pseudoatom, it records:

            - the largest observed displacement
            - the mode where that displacement occurs

        Biologically, this provides a simple way to identify:

            - highly mobile regions
            - flexible loops
            - hinge zones
            - localized hotspots of structural motion

        Outputs and Their Interpretation

        The protocol generates several outputs.

        outputModes

            A structured set of GNM normal modes including metadata such as:

                - eigenvalues
                - collectivity
                - ranking score
                - enable/disable flags

        matrixFileCC

            Cross-correlation matrix between nodes.

        matrixFileCV

            Covariance matrix of structural fluctuations.

        Additional metadata files are also generated for atom shift
        distributions and per-mode displacement profiles.

        Practical Recommendations

        For routine protein flexibility analysis:

            - start with moderate mode numbers
            - use default spring constant
            - choose a reasonable cutoff based on model granularity

        If too few modes are obtained, increasing the cutoff is often the
        most useful first adjustment.

        For membrane proteins, use the membrane option only when the input
        structure has biologically meaningful membrane orientation.

        When selecting modes for downstream interpretation, collectivity is
        often one of the most informative criteria.

        Final Perspective

        GNM does not simulate atomistic trajectories.

        Instead, it provides a physically intuitive description of intrinsic
        structural flexibility encoded by the contact topology.

        For structural biologists, its main strength lies in rapidly
        identifying collective motions that may underlie biological
        function, conformational change, or long-range communication
        within macromolecular assemblies.
    """
    _label = 'GNM analysis'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label='ProDy GNM NMA')

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

        form.addParam('cutoff', FloatParam, default=10,
                      expertLevel=LEVEL_ADVANCED,
                      label="Cut-off distance (A)",
                      help='Atoms or pseudoatoms beyond this distance will not interact. \n'
                           'For Calpha atoms, the default distance of 7.5 A works well in the majority of cases. \n'
                           'For all atoms, a shorter distance is recommended.'
                           'For fewer atoms or pseudoatoms, a longer distance is recommended.')

        form.addParam('gamma', StringParam, default=1.,
                      expertLevel=LEVEL_ADVANCED,
                      label="Spring constant",
                      help='This number or function determines the strength of the springs.\n'
                           'More sophisticated options are available within the ProDy API and '
                           'the resulting modes can be imported back into Scipion.\n'
                           'See http://http://www.bahargroup.org/prody/tutorials/enm_analysis/gamma.html')

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
                      'Modes 1-6 are always deselected as they are related to rigid-body movements. \n'
                      'The modes metadata file can be used to see which modes are more collective '
                      'in order to decide which modes to use at the image analysis step.')

        form.addParam('zeros', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Include zero eigvals",
                      help='Elect whether modes with zero eigenvalues will be kept.')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps

        # Link the input
        inputFn = self.inputStructure.get().getFileName()
        self.structureEM = self.inputStructure.get().getPseudoAtoms()
        n = self.numberOfModes.get()

        self._insertFunctionStep('computeModesStep', inputFn, n)
        self._insertFunctionStep('qualifyModesStep', n,
                                 self.collectivityThreshold.get(),
                                 self.structureEM)
        self._insertFunctionStep('computeAtomShiftsStep', n)
        self._insertFunctionStep('createOutputStep')

    def computeModesStep(self, inputFn, n):

        if self.structureEM:
            self.pdbFileName = self._getPath('pseudoatoms.pdb')
        else:
            self.pdbFileName = self._getPath('atoms.pdb')

        ag = prody.parsePDB(inputFn, alt='all')
        prody.writePDB(self.pdbFileName, ag)

        if self.membrane.get():
            self.prefix = 'modes.exgnm'
        else:
            self.prefix = 'modes.gnm'
        filename = self.prefix + '.npz'

        args = '{0} -s "all" --altloc "all" --kirchhoff --export-scipion --npz --npzmatrices ' \
               '-o {1} -p {2} -n {3} -g {4} -c {5} -P {6}'.format(self.pdbFileName,
                                                              self._getPath(),
                                                              self.prefix, n,
                                                              self.gamma.get(),
                                                              self.cutoff.get(),
                                                              self.numberOfThreads.get())

        if self.zeros.get():
            args += ' --zero-modes'
            self.startMode = 1
        else:
            self.startMode = 0
        
        if self.membrane.get():
            args += ' --membrane'

        self.runJob(Plugin.getProgram('gnm'), args)

        self.gnm = prody.loadModel(self._getPath(filename))
        covariances = prody.calcCrossCorr(self.gnm[self.startMode:], norm=False)
        prody.writeArray(self._getExtraPath('modes_covariance.txt'), covariances)

        crossCorr = prody.calcCrossCorr(self.gnm[self.startMode:])
        prody.writeArray(self._getExtraPath('modes_crossCorr.txt'), crossCorr)

    def qualifyModesStep(self, numberOfModes, collectivityThreshold, structureEM, suffix=''):
        self._enterWorkingDir()

        fnVec = glob("modes/vec.*")

        if len(fnVec) < numberOfModes:
            msg = "There are only %d modes instead of %d. "
            msg += "Check the number of modes you asked to compute and/or consider increasing cut-off distance. "
            msg += "The maximum number of modes allowed by the method for GNM normal mode analysis is "
            msg += "1 times the number of nodes (atoms or pseudoatoms; %d). "
            self.warning(redStr(msg % (len(fnVec), numberOfModes, self.atoms.numAtoms())))

        mdOut = MetaData()
        collectivityList = list(prody.calcCollectivity(self.gnm))
        eigvals = self.gnm.getEigvals()

        vecStr = "vec.%d"

        for n in range(len(fnVec)):
            collectivity = collectivityList[n]

            objId = mdOut.addObject()
            modefile = self._getPath("modes", vecStr % (n + 1))
            mdOut.setValue(MDL_NMA_MODEFILE, modefile, objId)
            mdOut.setValue(MDL_ORDER, int(n + 1), objId)

            mdOut.setValue(MDL_NMA_COLLECTIVITY, collectivity, objId)

            eigval = eigvals[n]
            mdOut.setValue(MDL_NMA_EIGENVAL, eigval, objId)

            if eigval > prody.utilities.ZERO:
                mdOut.setValue(MDL_ENABLED, 1, objId)
            else:
                mdOut.setValue(MDL_ENABLED, -1, objId)

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
        
        prody.writeScipionModes(self._getPath(), self.gnm, scores=score, only_sqlite=True,
                                collectivityThreshold=collectivityThreshold)

    def computeAtomShiftsStep(self, numberOfModes):
        fnOutDir = self._getExtraPath("distanceProfiles")
        makePath(fnOutDir)
        maxShift=[]
        maxShiftMode=[]
        vecStr = "vec.%d"
        for n in range(self.startMode+1, numberOfModes+1):
            fnVec = self._getPath("modes", vecStr % n)
            if exists(fnVec):
                fhIn = open(fnVec)
                md = MetaData()
                atomCounter = 0
                for line in fhIn:
                    d = abs(float(line))
                    if n==self.startMode+1:
                        maxShift.append(d)
                        maxShiftMode.append(self.startMode+1)
                    else:
                        if d>maxShift[atomCounter]:
                            maxShift[atomCounter]=d
                            maxShiftMode[atomCounter]=n
                    atomCounter+=1
                    md.setValue(MDL_NMA_ATOMSHIFT,d,md.addObject())
                md.write(join(fnOutDir,"vec%d.xmd" % n))
                fhIn.close()
                
        md = MetaData()
        for i, _ in enumerate(maxShift):
            fnVec = self._getPath("modes", vecStr % (maxShiftMode[i]+1))
            if exists(fnVec):
                objId = md.addObject()
                md.setValue(MDL_NMA_ATOMSHIFT, maxShift[i],objId)
                md.setValue(MDL_NMA_MODEFILE, fnVec, objId)
        md.write(self._getExtraPath('maxAtomShifts.xmd'))

    def createOutputStep(self):
        outputMatrixCov = EMFile(filename=self._getExtraPath('modes_covariance.txt'))
        outputMatrixCrosCor = EMFile(filename=self._getExtraPath('modes_crossCorr.txt'))

        fnSqlite = self._getPath('modes.sqlite')
        nmSet = SetOfGnmModes(filename=fnSqlite)
        nmSet._nmdFileName = String(self._getPath(self.prefix + '.nmd'))

        inputPdb = self.inputStructure.get()
        nmSet.setPdb(inputPdb)

        self._defineOutputs(outputModes=nmSet,
                            matrixFileCC=outputMatrixCrosCor,
                            matrixFileCV=outputMatrixCov)
        self._defineSourceRelation(self.inputStructure, nmSet)

