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
from pwem.objects import SetOfAtomStructs, SetOfPrincipalComponents, String, AtomStruct

from pyworkflow.utils import glob, redStr
from pyworkflow.protocol.params import (PointerParam, IntParam, FloatParam,
                                        BooleanParam, StringParam,
                                        LEVEL_ADVANCED, Float)

from prody2.protocols.protocol_modes_base import ProDyModesBase
from prody2.objects import ProDyNpzEnsemble, TrajFrame
from prody2.constants import PCA_FRACT_VARS
from prody2 import Plugin

import prody
import matplotlib.pyplot as plt


class ProDyPCA(ProDyModesBase):
    """
    This protocol will perform ProDy principal component analysis (PCA) using atomic structures
    """
    _label = 'PCA'
    _possibleOutputs = {'outputModes': SetOfPrincipalComponents}

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
                      pointerClass='SetOfAtomStructs, ProDyNpzEnsemble',
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
        self._insertFunctionStep('animateModesStep', n,
                                 self.rmsd.get(), self.n_steps.get(),
                                 self.neg.get(), self.pos.get(), 0)
        self._insertFunctionStep('createOutputStep')

    def computeModesStep(self, n=5):
        # configure ProDy to automatically handle secondary structure information and verbosity
        self.oldSecondary = prody.confProDy("auto_secondary")
        self.oldVerbosity = prody.confProDy("verbosity")
        
        from pyworkflow import Config
        prodyVerbosity =  'none' if not Config.debugOn() else 'debug'
        prody.confProDy(auto_secondary=True, verbosity='{0}'.format(prodyVerbosity))

        loadAndWriteEnsemble(self) # creates self.npz and others
        
        # configure ProDy to restore secondary structure information and verbosity
        prody.confProDy(auto_secondary=self.oldSecondary, verbosity='{0}'.format(self.oldVerbosity))

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
        
        self.outModes, self.atoms = prody.parseNMD(self._getPath('modes.pca.nmd'), type=prody.PCA)
        if not self.keepAlignment.get():
            dcdEnsemble = prody.parseDCD(self._getPath('ensemble.dcd'))
            dcdEnsemble.iterpose()
            self.npz.replaceCoordsets(dcdEnsemble.getCoordsets(), suffix='_aligned')
        
        plt.figure()
        prody.showFractVars(self.outModes)
        prody.showCumulFractVars(self.outModes, 'r')
        plt.savefig(self._getPath('pca_fract_vars.png'))
        
        self.fract_vars = prody.calcFractVariance(self.outModes)
        prody.writeArray(self._getPath('pca_fract_vars.txt'), self.fract_vars)

    def qualifyModesStep(self, numberOfModes, collectivityThreshold):
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
        nmSet._nmdFileName = String(self._getPath('modes.pca.nmd'))

        self.fractVarsDict = {}
        for i, item in enumerate(nmSet):
            self.fractVarsDict[item.getObjId()] = self.fract_vars[i]

        outSet = SetOfPrincipalComponents().create(self._getPath())
        outSet.copyItems(nmSet, updateItemCallback=self._setFractVars)

        inputPdb = self.averageStructure
        self._defineOutputs(refPdb=inputPdb)
        outSet.setPdb(inputPdb)

        self._defineOutputs(outputModes=outSet, outputEnsemble=self.npz)
        self._defineSourceRelation(inputPdb, outSet)

    def _summary(self):
        if not hasattr(self, 'outputModes'):
            summ = ['Output modes not ready yet']
        else:
            modes = prody.parseScipionModes(self.outputModes.getFileName())
            ens = self.outputEnsemble.loadEnsemble()

            summ = ['*{0}* principal components calculated from *{1}* structures of *{2}* atoms'.format(
                    modes.numModes(), ens.numConfs(), ens.numAtoms())]
        return summ

    def _setFractVars(self, item, row=None):
        # We provide data directly so don't need a row
        fractVar = Float(self.fractVarsDict[item.getObjId()])
        setattr(item, PCA_FRACT_VARS, fractVar)

def loadAndWriteEnsemble(cls):
    """Handle inputs to load ensemble into ProDy and write outputs"""

    inputEnsemble = cls.inputEnsemble.get()

    if isinstance(inputEnsemble, SetOfAtomStructs):
        ags = prody.parsePDB([tarStructure.getFileName() for tarStructure in inputEnsemble])
        cls.ens = prody.buildPDBEnsemble(ags, match_func=prody.sameChainPos, seqid=0., 
                                            overlap=0., superpose=False, degeneracy=cls.degeneracy.get())
        # the ensemble gets built exactly as the input is setup and nothing gets rejected
    else:
        cls.ens = inputEnsemble.loadEnsemble()

    cls.ens.select(cls.selstr.get())

    avgStruct = cls.ens.getAtoms()
    avgStruct.setCoords(cls.ens.getCoords())

    cls.pdbFileName = cls._getPath('atoms.pdb')
    prody.writePDB(cls.pdbFileName, avgStruct)
    cls.averageStructure = AtomStruct()
    cls.averageStructure.setFileName(cls.pdbFileName)

    cls.dcdFileName = cls._getPath('ensemble.dcd')
    prody.writeDCD(cls.dcdFileName, cls.ens)

    cls.npzFileName = cls._getPath('ensemble.ens.npz')
    prody.saveEnsemble(cls.ens, cls.npzFileName)
    cls.npz = ProDyNpzEnsemble().create(cls._getExtraPath())
    for j in range(cls.ens.numConfs()):
        frame = TrajFrame((j+1, cls.npzFileName), objLabel=cls.ens.getLabels()[j])
        cls.npz.append(frame)
