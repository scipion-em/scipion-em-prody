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

import numpy as np

from pwem.objects import SetOfAtomStructs, AtomStruct
from pwem.protocols import EMProtocol

from pyworkflow.protocol.params import (PointerParam, EnumParam, FloatParam,
                                        StringParam, BooleanParam)
import pyworkflow.object as pwobj

from prody2.objects import ProDyNpzEnsemble, TrajFrame
from prody2.constants import ENSEMBLE_WEIGHTS

import prody
import matplotlib.pyplot as plt


class ProDyRmsd(EMProtocol):
    """
    This protocol will perform ProDy principal component analysis (PCA) using atomic structures
    """
    _label = 'RMSD Reorder and Cluster'
    _possibleOutputs = {'outputEnsemble': ProDyNpzEnsemble}

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """

        form.addSection(label='ProDy RMSD')
        form.addParam('inputEnsemble', PointerParam, label="Input ensemble",
                      important=True,
                      pointerClass='SetOfAtomStructs, ProDyNpzEnsemble',
                      help='The input ensemble should be a SetOfAtomStructs '
                      'where all structures have the same number of atoms.')
        
        form.addParam('treeMethod', EnumParam, choices=['upgma', 'nj',
                                                        'single', 'average',
                                                        'ward', 'other'],
                      label="RMSD tree method", default=0,
                      display=EnumParam.DISPLAY_HLIST,
                      help='Method for calculating RMSD tree.\n'
                           'Acceptable options are ``"upgma"``, ``"nj"``, or methods '
                           'supported by :func:`~scipy.cluster.hierarchy.linkage` such as '
                           '``"single"``, ``"average"``, ``"ward"``, etc.')
        
        form.addParam('otherMethod', StringParam, default="None",
                      condition="treeMethod==5",
                      label="Other tree method",
                      help='You can type another tree method here')
        
        form.addParam('doReorder', BooleanParam, default=True,
                      label="Reorder ensemble?",
                      help='Whether to reorder ensemble')
        
        form.addParam('doCluster', BooleanParam, default=True,
                      label="Cluster ensemble?",
                      help='Whether to cluster ensemble')
        
        form.addParam('rmsdThreshold', FloatParam, default=1.,
                      condition="doCluster==True",
                      label='RMSD threshold',
                      help='RMSD threshold for clustering the tree')
        
        form.addParam('writePDBFiles', BooleanParam, default=False,
                      label="Whether to write lots of PDB files",
                      help='These will be registered as output too')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps

        # Link the input
        inputEnsemble = self.inputEnsemble.get()
        if isinstance(inputEnsemble, SetOfAtomStructs):
            ags = prody.parsePDB([tarStructure.getFileName() for tarStructure in inputEnsemble])
            self.ens = prody.buildPDBEnsemble(ags, match_func=prody.sameChainPos, seqid=0., overlap=0., superpose=False)
            # the ensemble gets built exactly as the input is setup and nothing gets rejected
        else:
            self.ens = inputEnsemble.loadEnsemble()

        self._insertFunctionStep('ensembleModificationStep')
        self._insertFunctionStep('createOutputStep')

    def ensembleModificationStep(self):

        matrix = self.ens.getRMSDs(pairwise=True)
        labels = self.ens.getLabels()

        if len(labels) > 50:
            allticks = False
        else:
            allticks = True

        plt.figure()
        prody.showMatrix(matrix, allticks=allticks)
        plt.tight_layout()
        plt.savefig(self._getExtraPath('rmsd_matrix'))

        tree = prody.calcTree(labels, matrix)

        plt.figure()
        prody.showTree(tree)
        plt.tight_layout()
        plt.axis('off')
        plt.savefig(self._getExtraPath('rmsd_tree'))

        reordRMSDs, reordIndices = prody.reorderMatrix(labels,
                                            matrix,
                                            tree)
        
        reordLabels = list(np.array(labels)[reordIndices])

        plt.figure()
        prody.showMatrix(reordRMSDs, allticks=allticks)
        plt.tight_layout()
        plt.savefig(self._getExtraPath('reordered_matrix'))

        if self.doReorder.get():
            self.ens = self.ens[reordIndices]
            self.ensBaseName = self._getExtraPath('rmsd_reord_ensemble')
            labels = reordLabels
        else:
            self.ensBaseName = self._getExtraPath('ensemble')

        if not self.doCluster.get():
            idx = range(self.ens.numConfs())
        else:
            subgroups = prody.findSubgroups(tree, self.rmsdThreshold.get())
            self.weights = np.zeros(len(subgroups), dtype=int)
            allWeights = np.zeros(self.ens.numCoordsets(), dtype=int)
            idx = np.zeros(len(subgroups), dtype=int)
            for i, sg in enumerate(subgroups):
                sgIdx = [labels.index(label) for label in sg]
                submatrix = matrix[sgIdx, :][:, sgIdx]
                idx[i] = sgIdx[np.argmin(np.mean(submatrix, axis=0))]

                weight = len(sg)
                self.weights[i] = weight
                allWeights[idx[i]] = weight

        prody.writePDB(self.ensBaseName, self.ens)

        self.ens.setData('size', allWeights)
        prody.saveEnsemble(self.ens, self.ensBaseName)

        self.npz = ProDyNpzEnsemble().create(self._getExtraPath())

        if self.writePDBFiles.get():
            ag = self.ens.getAtoms().copy()
            self.pdbs = SetOfAtomStructs().create(self._getExtraPath())

        for j in sorted(idx):
            label = self.ens.getLabels()[j]
            frame = TrajFrame((j+1, self.ensBaseName+'.ens.npz'), objLabel=label)
            self.npz.append(frame)
            
            if self.writePDBFiles.get():
                ag.setCoords(self.ens.getCoordsets()[j])
                filename = self._getExtraPath('{:06d}_{:s}.pdb'.format(j+1, label))
                prody.writePDB(filename, ag)
                pdb = AtomStruct(filename)
                self.pdbs.append(pdb)

    def createOutputStep(self):

        args = {}

        outNpz = ProDyNpzEnsemble().create(self._getPath())
        outNpz.copyItems(self.npz, updateItemCallback=self._setWeights)
        args["outputNpz"] = outNpz

        if self.writePDBFiles.get():
            outSetAS = SetOfAtomStructs().create(self._getPath())
            outSetAS.copyItems(self.pdbs, updateItemCallback=self._setWeights)
            args["outputStructures"] = outSetAS

        self._defineOutputs(**args)
        
    def _summary(self):
        if not hasattr(self, 'outputNpz'):
            summ = ['Output ensemble not ready yet']
        else:
            ens = self.outputNpz.loadEnsemble()
            summ = ['Output ensemble has *{0}* structures of *{1}* atoms'.format(
                   ens.numConfs(), ens.numAtoms())]
        return summ

    def _setWeights(self, item, row=None):
            weight = pwobj.Integer(self.weights[item.getObjId()-1])
            setattr(item, ENSEMBLE_WEIGHTS, weight)
