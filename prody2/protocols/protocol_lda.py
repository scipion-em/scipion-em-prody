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
This module will provide ProDy linear discriminant analysis (LDA) using atomic structures
"""
from collections import OrderedDict

from pwem import *
from pwem.emlib import (MetaData, MDL_NMA_MODEFILE, MDL_ORDER,
                        MDL_ENABLED, MDL_NMA_COLLECTIVITY, MDL_NMA_SCORE, 
                        MDL_NMA_EIGENVAL)
from pwem.objects import SetOfAtomStructs, SetOfPrincipalComponents, String, AtomStruct

from pyworkflow.utils import *
from pyworkflow.protocol.params import (PointerParam, IntParam, FloatParam,
                                        BooleanParam, StringParam, TextParam, 
                                        NumericRangeParam, 
                                        LEVEL_ADVANCED, Float)

from prody2.protocols.protocol_modes_base import ProDyModesBase
from prody2.objects import ProDyNpzEnsemble, TrajFrame, SetOfLdaModes
from prody2.constants import LDA_FRACT_VARS

import prody
import matplotlib.pyplot as plt


class ProDyLDA(ProDyModesBase):
    """
    This protocol will perform ProDy linear discriminant analysis (LDA) using atomic structures
    """
    _label = 'LDA'
    _possibleOutputs = {'outputModes': SetOfLdaModes}

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:

        form.addSection(label='ProDy LDA')
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
        form.addParam('numberOfModes', IntParam, default=1,
                      label='Number of modes',
                      help='The maximum number of modes allowed by the method for '
                           'atomic linear discriminant analysis is the number of classes - 1.')
        form.addParam('numberOfShuffles', IntParam, default=10,
                      label='Number of random shuffles',
                      help='The class labels will be shuffled this many times for LDA to '
                           'assess random variation.')
        form.addParam('selstr', StringParam, default="name CA",
                      label="Selection string",
                      help='Selection string for atoms to include in the calculation.\n'
                           'It is recommended to use "name CA" (default)')
        
        group = form.addGroup('Class labels')
        group.addParam('chainOrders', TextParam, width=60, default='{}',
                       label='Custom class label dictionary',
                       help='Defined labels for classes. These can be any string including numbers')
        group.addParam('insertOrder', NumericRangeParam, default='1',
                       label='Insert label index',
                       help='Insert the class label with the specified index into the label dict.\n'
                            'The default (when empty) is the last position.')
        group.addParam('customOrder', StringParam, default='1',
                       label='Custom label to insert at the specified index',
                       help='Enter the desired label here.\n'
                            'The default (when empty) is the number 1.')
        group.addParam('label', StringParam, default='',
                       label='Ensemble label for item with the specified number for recovering custom class labels',
                       help='This cannot be changed by the user and is for display only.')
        group.addParam('recoverOrder', StringParam, default='1',
                       label='Recover custom label number',
                       help='Enter the desired class label index here.\n'
                            'Recover the class label with the specified index from the label dict.')

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

        self.model_type = 'lda'
        n = self.numberOfModes.get()

        self.gnm = False
        nzeros = 0

        self._insertFunctionStep('computeModesStep', n)
        self._insertFunctionStep('qualifyModesStep', n)
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

        inputEnsemble = self.inputEnsemble.get()
        if isinstance(inputEnsemble, SetOfAtomStructs):
            ags = prody.parsePDB([tarStructure.getFileName() for tarStructure in inputEnsemble])
            self.ens = prody.buildPDBEnsemble(ags, match_func=prody.sameChainPos, seqid=0., 
                                              overlap=0., superpose=False,
                                              degeneracy=self.degeneracy.get())
            # the ensemble gets built exactly as the input is setup and nothing gets rejected
        else:
            self.ens = inputEnsemble.loadEnsemble()

        self.ens.select(self.selstr.get())

        self.atoms = self.ens.getAtoms()
        self.atoms.setCoords(self.ens.getCoords())

        self.pdbFileName = self._getPath('atoms.pdb')
        prody.writePDB(self.pdbFileName, self.atoms)
        self.averageStructure = AtomStruct()
        self.averageStructure.setFileName(self.pdbFileName)

        self.dcdFileName = self._getPath('ensemble.dcd')
        prody.writeDCD(self.dcdFileName, self.ens)

        self.npzFileName = self._getPath('ensemble.ens.npz')
        prody.saveEnsemble(self.ens, self.npzFileName)
        self.npz = ProDyNpzEnsemble().create(self._getExtraPath())
        for j in range(self.ens.numConfs()):
            frame = TrajFrame((j+1, self.npzFileName), objLabel=self.ens.getLabels()[j])
            self.npz.append(frame)

        # configure ProDy to restore secondary structure information and verbosity
        prody.confProDy(auto_secondary=self.oldSecondary, verbosity='{0}'.format(self.oldVerbosity))

        labelsMap = eval(self.chainOrders.get())
        self.classes = list(labelsMap.values())

        self.outModes = prody.LDA()
        self.outModes.calcModes(self.ens, self.classes,
                                self.numberOfModes.get(),
                                n_shuffles=self.numberOfShuffles.get())
        
        plt.figure()
        prody.showFractVars(self.outModes)
        prody.showCumulFractVars(self.outModes, 'r')
        plt.savefig(self._getPath('lda_fract_vars.png'))
        
        self.fract_vars = prody.calcFractVariance(self.outModes)
        prody.writeArray(self._getPath('lda_fract_vars.txt'), self.fract_vars)

        prody.writeScipionModes(self._getPath(), self.outModes)
        prody.writeNMD(self._getPath('modes.nmd'), self.outModes, self.atoms)
        prody.saveModel(self.outModes, self._getPath('modes.lda.npz'), matrices=True)

    def qualifyModesStep(self, numberOfModes):
        fnVec = glob(self._getPath("modes/vec.*"))

        if len(fnVec) < numberOfModes:
            msg = "There are only %d modes instead of %d. "
            msg += "Check the number of modes you asked to compute. "
            msg += "The maximum number of components allowed by the method for linear discriminant analysis is "
            msg += "the number of classes - 1 (%d)."
            self.warning(redStr(msg % (len(fnVec), numberOfModes, len(set(self.classes))-1)))

        mdOut = MetaData()
        collectivity = prody.calcCollectivity(self.outModes)
        if isinstance(collectivity, float):
            collectivityList = [collectivity]
        else:
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

        idxSorted = [i[0] for i in sorted(enumerate(collectivityList),
                                          key=lambda x: x[1], reverse=True)]

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
        mdOut.write(self._getPath("modes.xmd"))
        
        prody.writeScipionModes(self._getPath(), self.outModes, scores=score,
                                only_sqlite=True)

    def createOutputStep(self):
        fnSqlite = self._getPath('modes.sqlite')
        nmSet = SetOfLdaModes(filename=fnSqlite)
        nmSet._nmdFileName = String(self._getPath('modes.nmd'))


        self.fractVarsDict = {}
        for i, item in enumerate(nmSet):
            self.fractVarsDict[item.getObjId()] = self.fract_vars[i]

        outSet = SetOfLdaModes().create(self._getPath())
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

            summ = ['*{0}* LDA components calculated from *{1}* structures of *{2}* atoms'.format(
                    modes.numModes(), ens.numConfs(), ens.numAtoms())]
        return summ

    def _setFractVars(self, item, row=None):
        # We provide data directly so don't need a row
        fractVar = Float(self.fractVarsDict[item.getObjId()])
        setattr(item, LDA_FRACT_VARS, fractVar)

    def createMatchDic(self, index, label=None):

        # configure ProDy to automatically handle secondary structure information and verbosity
        oldSecondary = prody.confProDy("auto_secondary")
        oldVerbosity = prody.confProDy("verbosity")
        
        from pyworkflow import Config
        prodyVerbosity =  'none' if not Config.debugOn() else 'debug'
        prody.confProDy(auto_secondary=False, verbosity='{0}'.format(prodyVerbosity))
        
        self.matchDic = eval(self.chainOrders.get())
        self.labels = list(self.matchDic.keys())
        self.classes = list(self.matchDic.values())

        # reinitialise to update with new keys
        # that are still ordered correctly
        self.matchDic = OrderedDict()

        if self.labels == []:
            inputEnsemble = self.inputEnsemble.get()
            if isinstance(inputEnsemble, SetOfAtomStructs):
                ags = prody.parsePDB([tarStructure.getFileName() for tarStructure in inputEnsemble])
                self.ens = prody.buildPDBEnsemble(ags, match_func=prody.sameChainPos, seqid=0., 
                                                overlap=0., superpose=False, degeneracy=self.degeneracy.get())
                # the ensemble gets built exactly as the input is setup and nothing gets rejected
            else:
                self.ens = inputEnsemble.loadEnsemble()

            self.labels = self.ens.getLabels()
            self.classes = list(np.ones(len(self.labels), dtype=str))

        if not isinstance(self.labels[0], tuple):
            self.labels = [(i+1, label) for i, label in enumerate(self.labels)]

        inds = [item-1 for item in getListFromRangeString(index)]
        for idx in inds:
            self.classes[idx] = self.customOrder.get()
        
        # configure ProDy to restore secondary structure information and verbosity
        prody.confProDy(auto_secondary=oldSecondary, verbosity='{0}'.format(oldVerbosity))

        self.matchDic.update(zip(self.labels, self.classes))
        return self.matchDic
