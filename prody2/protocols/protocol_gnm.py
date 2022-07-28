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
from os.path import basename, exists, join
import math

from pwem import *
from pwem.emlib import (MetaData, MDL_NMA_MODEFILE, MDL_ORDER,
                        MDL_ENABLED, MDL_NMA_COLLECTIVITY, MDL_NMA_SCORE, 
                        MDL_NMA_ATOMSHIFT, MDL_NMA_EIGENVAL)
from pwem.objects import AtomStruct, SetOfNormalModes, String, EMFile
from pwem.protocols import EMProtocol

from pyworkflow.utils import *
from pyworkflow.utils.path import makePath
from pyworkflow.protocol.params import (PointerParam, IntParam, FloatParam, StringParam,
                                        BooleanParam, LEVEL_ADVANCED)

import prody

class ProDyGNM(EMProtocol):
    """
    This protocol will perform normal mode analysis (NMA) using the Gaussian network model (GNM)
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
                           'For pseudoatoms, set this according to the level of coarse-graining '
                           '(see Doruker et al., J Comput Chem 2002). \n'
                           'For all atoms, a shorter distance is recommended.')

        form.addParam('gamma', FloatParam, default=1.,
                      expertLevel=LEVEL_ADVANCED,
                      label="Spring constant",
                      help='This number or function determines the strength of the springs.\n'
                           'More sophisticated options are available within the ProDy API and '
                           'the resulting modes can be imported back into Scipion.\n'
                           'See http://prody.csb.pitt.edu/tutorials/enm_analysis/gamma.html')

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

        self.model_type = 'gnm'
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

        if self.zeros.get():
            args = ('gnm {0} -s "all" --altloc "all" --zero-modes --kirchhoff '
                        '--export-scipion --npz -o {1} -p modes -n {2} -g {3} -c {4}'.format(self.pdbFileName,
                                                                                             self._getPath(), n,
                                                                                             self.gamma.get(),
                                                                                             self.cutoff.get()))
            self.startMode = 1
        else:
            args = ('gnm {0} -s "all" --altloc "all" --kirchhoff '
                        '--export-scipion --npz -o {1} -p modes -n {2} -g {3} -c {4}'.format(self.pdbFileName,
                                                                                             self._getPath(), n,
                                                                                             self.gamma.get(),
                                                                                             self.cutoff.get()))
            self.startMode = 0

        if self.membrane.get():
            args += ' --membrane'

        self.runJob('prody', args)
        
        self.gnm = prody.loadModel(self._getPath('modes.gnm.npz'))
        
        eigvecs = self.gnm.getEigvecs()
        eigvals = self.gnm.getEigvals()
        kirchhoff = prody.parseArray(self._getPath('modes_kirchhoff.txt'))

        self.gnm.setKirchhoff(kirchhoff)
        self.gnm.setEigens(eigvecs, eigvals)
        prody.saveModel(self.gnm, self._getPath('modes.gnm.npz'), matrices=True)

        covariances = prody.calcCrossCorr(self.gnm[self.startMode:], norm=False)
        prody.writeArray(self._getExtraPath('modes_covariance.txt'), covariances)

        crossCorr = prody.calcCrossCorr(self.gnm[self.startMode:])
        prody.writeArray(self._getExtraPath('modes_crossCorr.txt'), crossCorr)

    def qualifyModesStep(self, numberOfModes, collectivityThreshold, structureEM, suffix=''):
        self._enterWorkingDir()

        fnVec = glob("modes/vec.*")

        if len(fnVec) < numberOfModes:
            msg = "There are only %d modes instead of %d. "
            msg += "Check the number of modes you asked to compute and/or consider increasing cut-off distance."
            msg += "The maximum number of modes allowed by the method for atomic normal mode analysis is "
            msg += "3 times the number of nodes (pseudoatoms or Calphas). "
            self._printWarnings(redStr(msg % (len(fnVec), numberOfModes)))

        mdOut = MetaData()
        collectivityList = list(prody.calcCollectivity(self.gnm))
        eigvals = self.gnm.getEigvals()

        for n in range(len(fnVec)):
            collectivity = collectivityList[n]

            objId = mdOut.addObject()
            modefile = self._getPath("modes", "vec.%d" % (n + 1))
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
        for j in range(len(fnVec)):
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
        
        for n in range(self.startMode+1, numberOfModes+1):
            fnVec = self._getPath("modes", "vec.%d" % n)
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
            fnVec = self._getPath("modes", "vec.%d" % (maxShiftMode[i]+1))
            if exists(fnVec):
                objId = md.addObject()
                md.setValue(MDL_NMA_ATOMSHIFT, maxShift[i],objId)
                md.setValue(MDL_NMA_MODEFILE, fnVec, objId)
        md.write(self._getExtraPath('maxAtomShifts.xmd'))

    def createOutputStep(self):
        outputMatrixCov = EMFile(filename=self._getExtraPath('modes_covariance.txt'))
        outputMatrixCrosCor = EMFile(filename=self._getExtraPath('modes_crossCorr.txt'))

        fnSqlite = self._getPath('modes.sqlite')
        nmSet = SetOfNormalModes(filename=fnSqlite)
        nmSet._nmdFileName = String(self._getPath('modes.nmd'))

        inputPdb = self.inputStructure.get()
        nmSet.setPdb(inputPdb)

        self._defineOutputs(outputModes=nmSet, matrixFileCC=outputMatrixCrosCor, matrixFileCV=outputMatrixCov)
        self._defineSourceRelation(self.inputStructure, nmSet)

