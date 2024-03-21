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
                                        StringParam, BooleanParam, IntParam)
import pyworkflow.object as pwobj
import pyworkflow.utils as pwutils

from prody2.objects import (ProDyNpzEnsemble, TrajFrame, 
                            SetOfClassesTraj, ClassTraj)
from prody2.constants import ENSEMBLE_WEIGHTS
from prody2 import Plugin

import prody
import matplotlib.pyplot as plt


class ProDyRmsd(EMProtocol):
    """
    This protocol will perform ProDy principal component analysis (PCA) using atomic structures
    """
    _label = 'RMSD Cluster'
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

        form.addParam('doCluster', BooleanParam, default=True,
                      label="Cluster ensemble?",
                      help='Whether to cluster ensemble')

        form.addParam('clusteringMethod', EnumParam, choices=['hierarchical', 'kmedoids'],
                      default=0,
                      display=EnumParam.DISPLAY_HLIST,
                      label="Clustering method",
                      condition='doCluster',
                      help='Method used for clustering using RMSD')
                
        form.addParam('doReorder', BooleanParam, default=False,
                      label="Reorder ensemble?",
                      condition='clusteringMethod==0',
                      help='Whether to reorder ensemble based on RMSD tree')
        
        form.addParam('treeMethod', EnumParam, choices=['upgma', 'nj',
                                                        'single', 'average',
                                                        'ward', 'other'],
                      condition='clusteringMethod==0',
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

        form.addParam('rmsdThreshold', FloatParam, default=1.,
                      condition="doCluster==True and clusteringMethod==0",
                      label='RMSD threshold',
                      help='RMSD threshold for clustering the tree')
        
        form.addParam('nClusters', IntParam, default=2,
                      condition="doCluster==True and clusteringMethod==1",
                      label='Number of clusters',
                      help='Kmedoids will create this many clusters')
        
        form.addParam('writePDBFiles', BooleanParam, default=False,
                      label="Write representative PDB files",
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

        self.ensBaseName = self._getExtraPath('ensemble')
        ensFn = prody.saveEnsemble(self.ens, self.ensBaseName)

        allWeights = self.ens.getData('size')
        if allWeights is None:
            allWeights = np.ones(self.ens.numConfs(), dtype=float)
        
        if not self.doCluster.get():
            repIdx = range(self.ens.numConfs())
        else:
            if self.clusteringMethod.get() == 0:

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
                plt.close()

                tree = prody.calcTree(labels, matrix)

                plt.figure()
                prody.showTree(tree)
                plt.tight_layout()
                plt.axis('off')
                plt.savefig(self._getExtraPath('rmsd_tree'))
                plt.close()

                reordRMSDs, reordIndices = prody.reorderMatrix(labels,
                                                               matrix,
                                                               tree)

                plt.figure()
                prody.showMatrix(reordRMSDs, allticks=allticks)
                plt.tight_layout()
                plt.savefig(self._getExtraPath('reordered_matrix'))
                plt.close()

                classLabels = np.zeros(self.ens.numCoordsets(), dtype=int)
                subgroups = prody.findSubgroups(tree, self.rmsdThreshold.get())
                self.weights = np.zeros(len(subgroups), dtype=float)
                repIdx = np.zeros(len(subgroups), dtype=int)
                sgIdx = []
                for i, sg in enumerate(subgroups):
                    sgIdx.append([labels.index(label) for label in sg])
                    submatrix = matrix[sgIdx[i], :][:, sgIdx[i]]
                    repIdx[i] = sgIdx[i][np.argmin(np.mean(submatrix, axis=0))]

                    weight = len(sg)
                    self.weights[i] = allWeights[repIdx[i]] * weight
                    allWeights[sgIdx[i]] *= weight
                    classLabels[sgIdx[i]] = i

                if self.doReorder.get():
                    self.ens = self.ens[reordIndices]
            else:
                args = '--inputEns {0} --nClusters {1} --outputDir {2}'.format(ensFn, self.nClusters.get(), 
                                                                               self._getExtraPath())
                self.runJob(Plugin.getProgram('kmedoids.py', script=True), args)

                classLabels = np.loadtxt(self._getExtraPath("cluster_labels.txt"))
                repIdx = np.loadtxt(self._getExtraPath("cluster_medoids.txt"), dtype=int)
                weights = np.loadtxt(self._getExtraPath("cluster_counts.txt"))
                
                sgIdx = [np.nonzero(classLabels==label)[0] for label in np.unique(classLabels)]
                self.weights = np.zeros(len(weights), dtype=float)
                for i, weight in enumerate(weights):
                    weight /= weights.sum()
                    allWeights[sgIdx[i]] *= weight
                    self.weights[i] = allWeights[repIdx[i]]

        prody.writePDB(self.ensBaseName, self.ens)
        self.ens.setData('size', allWeights)
        prody.saveEnsemble(self.ens, self.ensBaseName)

        if self.writePDBFiles.get():
            ag = self.ens.getAtoms().copy()
            self.pdbs = SetOfAtomStructs().create(self._getExtraPath())

        self.npzClasses = SetOfClassesTraj().create(self._getExtraPath())
        frames = ProDyNpzEnsemble().create(self._getExtraPath())
        for i, label in enumerate(self.ens.getLabels()):
            frames.append(TrajFrame((i+1, self.ensBaseName+'.ens.npz'), 
                                    objLabel=label, weight=allWeights[i]))

        for i, sg in enumerate(sgIdx):
            newClass = ClassTraj().create(self._getExtraPath(), suffix=i+1)
            newClass.setRef(frames[repIdx[i]])
            self.npzClasses.append(newClass)
            for j in sg:
                newClass.append(frames[int(j+1)])
            self.npzClasses.update(newClass)
            
            if self.writePDBFiles.get():
                ag.setCoords(self.ens.getCoordsets()[i])
                filename = self._getExtraPath('{:06d}_{:s}.pdb'.format(i+1, 
                                                                       self.ens.getLabels()[i]))
                prody.writePDB(filename, ag)
                pdb = AtomStruct(filename)
                self.pdbs.append(pdb)

        self.npzClasses.write()

    def createOutputStep(self):
        args = {}
        args["outputClasses"] = self.npzClasses

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
            weight = pwobj.Float(self.weights[item.getObjId()-1])
            setattr(item, ENSEMBLE_WEIGHTS, weight)

    def _createSetOfClassesTraj(self, frameSet, suffix=''):
        classes = self.__createSet(SetOfClassesTraj,
                                   self._getExtraPath('classesTraj%s.sqlite'), 
                                   suffix)
        classes.setImages(frameSet)
        return classes

    def __createSet(self, SetClass, template, suffix, **kwargs):
        """ Create a set and set the filename using the suffix.
        If the file exists, it will be delete. """
        setFn = self._getPath(template % suffix)
        # Close the connection to the database if
        # it is open before deleting the file
        pwutils.cleanPath(setFn)
        setObj = SetClass(filename=setFn, **kwargs)
        return setObj
