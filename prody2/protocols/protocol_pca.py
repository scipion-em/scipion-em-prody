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
from pyworkflow.protocol import params

from os.path import basename, exists, join
import math
from multiprocessing import cpu_count

from pwem import *
from pwem.emlib import (MetaData, MDL_NMA_MODEFILE, MDL_ORDER,
                        MDL_ENABLED, MDL_NMA_COLLECTIVITY, MDL_NMA_SCORE, 
                        MDL_NMA_ATOMSHIFT, MDL_NMA_EIGENVAL)
from pwem.objects import SetOfAtomStructs, SetOfPrincipalComponents, String
from pwem.protocols import EMProtocol

from pyworkflow.utils import *
from pyworkflow.utils.path import makePath
from pyworkflow.protocol.params import (PointerParam, IntParam, FloatParam, StringParam,
                                        BooleanParam, LEVEL_ADVANCED)

from prody2.protocols.protocol_modes_base import ProDyModesBase

import prody
import matplotlib.pyplot as plt


class ProDyPCA(ProDyModesBase):
    """
    This protocol will perform ProDy principal component analysis (PCA) using atomic structures
    """
    _label = 'PCA'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        cpus = cpu_count()//2 # don't use everything
        form.addParallelSection(threads=cpus, mpi=0)

        form.addSection(label='ProDy PCA')
        form.addParam('inputEnsemble', PointerParam, label="Input ensemble",
                      important=True,
                      pointerClass='SetOfAtomStructs',
                      help='The input ensemble should be a SetOfAtomStructs '
                      'where all structures have the same number of atoms.')
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

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps

        # Link the input
        ags = prody.parsePDB([tarStructure.getFileName() for tarStructure in self.inputEnsemble.get()])
        ens = prody.buildPDBEnsemble(ags, match_func=prody.sameChainPos, seqid=0., overlap=0., superpose=False)
        # the ensemble gets built exactly as the input is setup and nothing gets rejected
        
        self.dcdFileName = self._getPath('ensemble.dcd')
        prody.writeDCD(self.dcdFileName, ens)

        self.model_type = 'pca'
        n = self.numberOfModes.get()

        self._insertFunctionStep('computeModesStep', ens, n)
        self._insertFunctionStep('qualifyModesStep', n,
                                 self.collectivityThreshold.get())
        self._insertFunctionStep('computeAtomShiftsStep', n)
        self._insertFunctionStep('animateModesStep', n,
                                 self.rmsd.get(), self.n_steps.get(),
                                 self.neg.get(), self.pos.get(), 0)
        self._insertFunctionStep('createOutputStep')

    def computeModesStep(self, ens, n):
        # configure ProDy to automatically handle secondary structure information and verbosity
        old_secondary = prody.confProDy("auto_secondary")
        old_verbosity = prody.confProDy("verbosity")
        
        from pyworkflow import Config
        prodyVerbosity =  'none' if not Config.debugOn() else 'debug'
        prody.confProDy(auto_secondary=True, verbosity='{0}'.format(prodyVerbosity))

        self.pdbFileName = self._getPath('atoms.pdb')
        prody.writePDB(self.pdbFileName, ens.getAtoms())

        self.runJob('prody', 'pca {0} --pdb {1} -s "all" --covariance --export-scipion'
                    ' -o {2} -p modes -n {3} -P {4} --aligned'.format(self.dcdFileName,
                                                                      self.pdbFileName,
                                                                      self._getPath(), n,
                                                                      self.numberOfThreads.get()))
        
        self.outModes, self.atoms = prody.parseNMD(self._getPath('modes.nmd'), type=prody.PCA)
        
        eigvecs = self.outModes.getEigvecs()
        eigvals = self.outModes.getEigvals()
        cov = prody.parseArray(self._getPath('modes_covariance.txt'))

        self.outModes.setCovariance(cov)
        self.outModes.setEigens(eigvecs, eigvals)
        prody.saveModel(self.outModes, self._getPath('modes.pca.npz'), matrices=True)
        
        plt.figure()
        prody.showFractVars(self.outModes)
        plt.savefig(self._getPath('pca_fract_vars.png'))
        
        self.fract_vars = prody.calcFractVariance(self.outModes)
        prody.writeArray(self._getExtraPath('fract_vars.txt'), self.fract_vars)
        
        # configure ProDy to restore secondary structure information and verbosity
        prody.confProDy(auto_secondary=old_secondary, verbosity='{0}'.format(old_verbosity))

    def qualifyModesStep(self, numberOfModes, collectivityThreshold, suffix=''):
        self._enterWorkingDir()

        fnVec = glob("modes/vec.*")

        if len(fnVec) < numberOfModes:
            msg = "There are only %d modes instead of %d. "
            msg += "Check the number of modes you asked to compute and/or consider increasing cut-off distance."
            msg += "The maximum number of modes allowed by the method for atomic normal mode analysis is "
            msg += "3 times the number of nodes (pseudoatoms or Calphas). "
            self._printWarnings(redStr(msg % (len(fnVec), numberOfModes)))

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
        
        prody.writeScipionModes(self._getPath(), self.outModes, scores=score, only_sqlite=True,
                                collectivityThreshold=collectivityThreshold)

    def computeAtomShiftsStep(self, numberOfModes):
        fnOutDir = self._getExtraPath("distanceProfiles")
        makePath(fnOutDir)
        maxShift=[]
        maxShiftMode=[]
        
        for n in range(7, numberOfModes+1):
            fnVec = self._getPath("modes", "vec.%d" % n)
            if exists(fnVec):
                fhIn = open(fnVec)
                md = MetaData()
                atomCounter = 0
                for line in fhIn:
                    x, y, z = map(float, line.split())
                    d = math.sqrt(x*x+y*y+z*z)
                    if n==7:
                        maxShift.append(d)
                        maxShiftMode.append(7)
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
        fnSqlite = self._getPath('modes.sqlite')
        nmSet = SetOfPrincipalComponents(filename=fnSqlite)
        nmSet._nmdFileName = String(self._getPath('modes.nmd'))
        
        outputFractVars = EMFile(filename=self._getExtraPath('fract_vars.txt'))

        self._defineOutputs(outFractVars=outputFractVars,
                            outputModes=nmSet)

