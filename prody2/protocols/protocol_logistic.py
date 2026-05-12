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
This module will provide ProDy linear discriminant analysis (LRA) using atomic structures
"""
from collections import OrderedDict
import numpy as np

from pwem.objects import Float, String
from pyworkflow.utils import getListFromRangeString
from pyworkflow.protocol.params import (MultiPointerParam, IntParam, FloatParam,
                                        BooleanParam, StringParam, TextParam, 
                                        NumericRangeParam, 
                                        LEVEL_ADVANCED, Float)

from prody2.protocols.protocol_modes_base import ProDyModesBase
from prody2.objects import SetOfLogisticModes, loadAndWriteEnsemble
from prody2.constants import PRODY_FRACT_VARS
from prody2 import parseMatchDict

import prody


class ProDyLRA(ProDyModesBase):
    """
    Performs ProDy Logistic Regression Analysis (LRA) on ensembles of atomic
    structures in order to identify structural variations that best separate
    two predefined classes.

    AI Generated:

    ProDy Logistic Regression Analysis (ProDyLRA) — User Manual
        Overview

        The ProDyLRA protocol applies logistic regression analysis (LRA)
        to a structural ensemble. Its main purpose is to detect collective
        structural changes that discriminate between two classes of
        conformations.

        In practical structural biology workflows, this protocol is useful
        when comparing two functional states of the same macromolecule,
        such as open versus closed conformations, ligand-bound versus
        ligand-free states, or wild-type versus mutant ensembles.

        Unlike classical normal mode analysis, LRA does not describe the
        dominant fluctuations of the ensemble alone. Instead, it finds the
        directions in conformational space that maximize separation between
        two user-defined groups.

        Input Data

        The protocol accepts structural ensembles provided as:

            - SetOfAtomStructs
            - ProDy NPZ ensembles
            - DCD molecular dynamics systems

        These ensembles must represent comparable conformations of the same
        molecular system.

        If atomic structures contain multiple coordinate sets, the user may
        optionally keep only the first conformation from each structure.
        This can be useful when each structure should contribute only one
        representative conformation.

        Atom Selection

        A selection string defines which atoms are included in the analysis.

        By default, the protocol uses:

            "name CA"

        This selects alpha carbons only, which is generally recommended for
        protein structural analyses because it reduces noise while preserving
        large-scale collective motions.

        More detailed selections are possible, but overly large selections
        may increase noise and computational cost.

        Class Labels

        Logistic regression requires exactly two classes.

        The protocol allows the user to assign custom class labels to the
        ensemble entries. Labels can represent any biologically meaningful
        grouping, for example:

            - state A vs state B
            - bound vs unbound
            - mutant vs wild type

        Labels are internally stored in an ordered dictionary that maps each
        ensemble element to a class.

        The protocol validates this step before execution. If the number of
        unique class labels is not exactly two, execution stops.

        This restriction is important because the underlying implementation
        performs binary logistic regression.

        Random Shuffling

        The parameter:

            numberOfShuffles

        controls how many random permutations of the class labels are
        generated.

        These shuffles estimate how much class separation could arise by
        chance alone.

        Biologically, this provides a simple way to assess whether the
        observed discriminative mode reflects meaningful structural
        differences rather than random variation in the ensemble.

        Workflow

        The protocol follows these main steps:

            1. Load and preprocess the structural ensemble.
            2. Generate the class label mapping.
            3. Compute the number of logistic modes.

        The number of modes is determined as:

            number of unique classes - 1

        Since LRA requires two classes, this normally produces one
        discriminative mode.

        Mode Computation

        During execution, the protocol creates a ProDy LRA object and
        computes logistic regression modes from the ensemble and the
        associated class labels.

        The analysis uses the specified number of shuffled label trials.

        The resulting model is written to disk in several formats:

            - Scipion modes format
            - NMD format for visualization
            - NPZ model file including matrices

        These outputs allow later visualization and downstream analysis.

        Output Generation

        After mode calculation, the protocol creates a
        SetOfLogisticModes object.

        Each generated mode is assigned a fractional variance value.

        In this implementation, every mode receives:

            fractional variance = 1

        This reflects that logistic regression modes are discriminative
        directions rather than conventional variance-explaining PCA modes.

        The output includes:

            - outputModes
            - outputEnsemble

        The resulting mode set is linked to the reference average structure,
        allowing direct structural interpretation.

        Animation

        The protocol also supports mode animation.

        Animation parameters include:

            - RMSD amplitude
            - number of frames
            - positive direction
            - negative direction

        These animations provide a visual representation of the structural
        displacement associated with the discriminative logistic mode.

        From a biological perspective, animation helps interpret which
        regions of the molecule contribute most strongly to class
        separation.

        Summary Information

        Once finished, the protocol reports a summary describing:

            - number of LRA modes
            - number of conformations analyzed
            - number of atoms included

        This provides a quick overview of the scale of the calculation.

        Practical Interpretation

        The most important biological meaning of ProDyLRA is that it
        identifies motions associated with class discrimination rather than
        simply structural variability.

        Therefore, large-amplitude motions found by LRA may not be the most
        frequent motions in the ensemble. Instead, they are the motions most
        strongly associated with the biological difference encoded in the
        labels.

        This makes the protocol particularly useful for studying:

            - conformational transitions
            - functional state changes
            - mutation-induced structural shifts
            - ligand-dependent rearrangements

        Practical Recommendations

        For most protein applications, using alpha carbons only is usually
        sufficient.

        The biological relevance of the results depends strongly on the
        quality of the class definition. Poorly defined classes may produce
        discriminative modes that are mathematically valid but biologically
        difficult to interpret.

        It is also important that both classes contain representative and
        sufficiently sampled conformations.

        Final Perspective

        ProDyLRA is best understood as a supervised structural analysis
        method.

        Rather than asking:

            "What motions dominate the ensemble?"

        it asks:

            "What motions best distinguish the two biological states?"

        This makes it especially powerful when the scientific question is
        focused on structural determinants of functional differences.
    """
    _label = 'LRA'
    _possibleOutputs = {'outputModes': SetOfLogisticModes}
    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form, besidesAnimation=False):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:

        form.addSection(label='ProDy LRA')
        form.addParam('inputEnsemble', MultiPointerParam, label="Input ensemble(s)",
                      important=True,
                      pointerClass='SetOfAtomStructs, ProDyNpzEnsemble, DcdMDSystem',
                      help='Each input ensemble should be a SetOfAtomStructs or a ProDy NPZ ensemble.')
        form.addParam('degeneracy', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      condition='isinstance(inputEnsemble, SetOfAtomStructs)',
                      label="Take only first conformation from each structure/set",
                      help='Elect whether only the active coordinate set (**True**) or all the coordinate sets '
                           '(**False**) of each structure should be added to the ensemble. Default is **True**.')
        form.addParam('numberOfShuffles', IntParam, default=10,
                      label='Number of random shuffles',
                      help='The class labels will be shuffled this many times for LRA to '
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
        self._insertFunctionStep('animateModesStep', self.rmsd.get(), self.n_steps.get(),
                                 self.neg.get(), self.pos.get(), 0)
        self._insertFunctionStep('createOutputStep')

    def computeModesStep(self, n=1):
        loadAndWriteEnsemble(self)
        self.atoms = self.ens.getAtoms()

        self.outModes = prody.LRA()
        self.outModes.calcModes(self.ens, self.classes,
                                n_shuffles=self.numberOfShuffles.get())

        prody.writeScipionModes(self._getPath(), self.outModes)
        self._nmdFileName = String(self._getPath('modes.logreg.nmd'))
        prody.writeNMD(self._nmdFileName.get(), self.outModes, self.atoms)
        prody.saveModel(self.outModes, self._getPath('modes.logreg.npz'), matrices=True)

    def createOutputStep(self):
        fnSqlite = self._getPath('modes.sqlite')
        nmSet = SetOfLogisticModes(filename=fnSqlite)
        nmSet._nmdFileName = self._nmdFileName

        self.fractVarsDict = {}
        for _, item in enumerate(nmSet):
            self.fractVarsDict[item.getObjId()] = 1

        outSet = SetOfLogisticModes().create(self._getPath())
        outSet.copyItems(nmSet, updateItemCallback=self._setFractVars)
        outSet._nmdFileName = self._nmdFileName

        inputPdb = self.averageStructure
        self._defineOutputs(refPdb=inputPdb)
        outSet.setPdb(inputPdb)

        self._defineOutputs(outputModes=outSet, outputEnsemble=self.npz)
        self._defineSourceRelation(inputPdb, outSet)

    def _validate(self):
        errors = []
        labelsMap = self.createMatchDic(self.insertOrder.get())
        numClasses = len(set(list(labelsMap.values())))
        if numClasses != 2:
            errors.append('The number of class labels should be 2')

        return errors

    def _summary(self):
        if not hasattr(self, 'outputModes'):
            summ = ['Output modes not ready yet']
        else:
            modes = prody.parseScipionModes(self.outputModes.getFileName())
            ens = self.outputEnsemble.loadEnsemble()

            summ = ['*{0}* LRA components calculated from *{1}* structures of *{2}* atoms'.format(
                    modes.numModes(), ens.numConfs(), ens.numAtoms())]
        return summ

    def _setFractVars(self, item, row=None):
        # We provide data directly so don't need a row
        fractVar = Float(self.fractVarsDict[item.getObjId()])
        setattr(item, PRODY_FRACT_VARS, fractVar)

    def createMatchDic(self, index, label=None):
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

        self.matchDic.update(zip(self.labels, self.classes))
        return self.matchDic
