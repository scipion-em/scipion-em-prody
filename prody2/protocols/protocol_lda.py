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

from pwem.emlib import (MetaData, MDL_NMA_MODEFILE, MDL_ORDER,
                        MDL_ENABLED, MDL_NMA_COLLECTIVITY, MDL_NMA_SCORE, 
                        MDL_NMA_EIGENVAL)
from pwem.objects import SetOfAtomStructs, String, AtomStruct

from pyworkflow.utils import glob, redStr, getListFromRangeString
from pyworkflow.protocol.params import (MultiPointerParam, IntParam, FloatParam,
                                        BooleanParam, StringParam, TextParam, 
                                        NumericRangeParam, 
                                        LEVEL_ADVANCED, Float)

from prody2.protocols.protocol_modes_base import ProDyModesBase
from prody2.protocols.protocol_pca import loadAndWriteEnsemble
from prody2.objects import ProDyNpzEnsemble, TrajFrame, SetOfLdaModes
from prody2.constants import PRODY_FRACT_VARS

import prody
import matplotlib.pyplot as plt
import numpy as np


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
        form.addParam('inputEnsemble', MultiPointerParam, label="Input ensemble(s)",
                      important=True,
                      pointerClass='SetOfAtomStructs, ProDyNpzEnsemble',
                      help='Each input ensemble should be a SetOfAtomStructs or a ProDy NPZ ensemble.')
        form.addParam('degeneracy', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      condition='isinstance(inputEnsemble, SetOfAtomStructs)',
                      label="Take only first conformation from each structure/set",
                      help='Elect whether only the active coordinate set (**True**) or all the coordinate sets '
                           '(**False**) of each structure should be added to the ensemble. Default is **True**.')
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
    def _insertAllSteps(self, n=1, nzeros=0):
        # Insert processing steps
        labelsMap = self.createMatchDic(self.insertOrder.get())
        self.classes = list(labelsMap.values())
        numModes = len(set(self.classes)) - 1
        self.gnm = False
        self.nzero = nzeros

        self._insertFunctionStep('computeModesStep', numModes)
        self._insertFunctionStep('qualifyModesStep', numModes, 0.)
        self._insertFunctionStep('computeAtomShiftsStep', numModes, nzeros)
        self._insertFunctionStep('animateModesStep', numModes,
                                 self.rmsd.get(), self.n_steps.get(),
                                 self.neg.get(), self.pos.get(), 0)
        self._insertFunctionStep('createOutputStep')

    def computeModesStep(self, n=1):
        # configure ProDy to automatically handle secondary structure information and verbosity
        self.oldSecondary = prody.confProDy("auto_secondary")
        self.oldVerbosity = prody.confProDy("verbosity")
        
        from pyworkflow import Config
        prodyVerbosity =  'none' if not Config.debugOn() else 'debug'
        prody.confProDy(auto_secondary=True, verbosity='{0}'.format(prodyVerbosity))

        loadAndWriteEnsemble(self)
        self.atoms = self.ens.getAtoms()

        # configure ProDy to restore secondary structure information and verbosity
        prody.confProDy(auto_secondary=self.oldSecondary, verbosity='{0}'.format(self.oldVerbosity))

        self.outModes = prody.LDA()
        self.outModes.calcModes(self.ens, self.classes, n,
                                n_shuffles=self.numberOfShuffles.get())
        
        plt.figure()
        prody.showFractVars(self.outModes)
        prody.showCumulFractVars(self.outModes, 'r')
        plt.savefig(self._getPath('lda_fract_vars.png'))
        
        self.fract_vars = prody.calcFractVariance(self.outModes)
        prody.writeArray(self._getPath('lda_fract_vars.txt'), self.fract_vars)

        prody.writeScipionModes(self._getPath(), self.outModes)
        self._nmdFileName = String(self._getPath('modes.lda.nmd'))
        prody.writeNMD(self._nmdFileName.get(), self.outModes, self.atoms)
        prody.saveModel(self.outModes, self._getPath('modes.lda.npz'), matrices=True)

    def createOutputStep(self):
        fnSqlite = self._getPath('modes.sqlite')
        nmSet = SetOfLdaModes(filename=fnSqlite)
        nmSet._nmdFileName = self._nmdFileName

        self.fractVarsDict = {}
        for i, item in enumerate(nmSet):
            self.fractVarsDict[item.getObjId()] = self.fract_vars[i]

        outSet = SetOfLdaModes().create(self._getPath())
        outSet.copyItems(nmSet, updateItemCallback=self._setFractVars)
        outSet._nmdFileName = self._nmdFileName

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
        setattr(item, PRODY_FRACT_VARS, fractVar)

    def createMatchDic(self, index, label=None):

        # configure ProDy to automatically handle secondary structure information and verbosity
        oldSecondary = prody.confProDy("auto_secondary")
        oldVerbosity = prody.confProDy("verbosity")
        
        from pyworkflow import Config
        prodyVerbosity =  'none' if not Config.debugOn() else 'debug'
        prody.confProDy(auto_secondary=False, verbosity='{0}'.format(prodyVerbosity))

        parseMatchDict(self)
        self.classes = list(self.matchDic.values())

        # reinitialise to update with new keys
        # that are still ordered correctly
        self.matchDic = OrderedDict()

        if self.labels == []:
            loadAndWriteEnsemble(self)
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

def parseMatchDict(cls):
    if cls.chainOrders.get() != "":
        cls.matchDic = eval(cls.chainOrders.get())
    else:
        cls.matchDic = OrderedDict()

    if not isinstance(cls.matchDic, OrderedDict):
        cls.matchDic = OrderedDict()

    cls.labels = list(cls.matchDic.keys())
