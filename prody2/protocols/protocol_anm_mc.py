# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     James Krieger (jamesmkrieger@gmail.com)
# *
# * Centro Nacional de Biotecnologia, CSIC
# * Francis Crick Institute, London, UK
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
This module will provide ProDy normal mode random walks using the anisotropic network model Monte Carlo algorithm (ANM MC).
"""
from multiprocessing import cpu_count
import numpy as np
import os
import prody

from prody2 import Plugin
from prody2.objects import DcdMDSystem

from pwem.objects import AtomStruct, SetOfAtomStructs
from pwem.protocols import EMProtocol

from pyworkflow.protocol import params

class ProDyANMMC(EMProtocol):
    """
    This protocol will generate random walks in normal mode space using ANM Monte Carlo
    """
    _label = 'ANM MC walks'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        cpus = cpu_count()//2 # don't use everything
        form.addParam('binThreads', params.IntParam,
                      label='threads',
                      default=cpus,
                      help='Number of threads used by ProDy each time it is called in the protocol execution. For '
                           'example, if 3 Scipion threads and 6 ProDy threads are set, 2 random walks will run '
                           'at the same time with each call of ProDy using 6 threads for normal mode analysis, so '
                           '12 threads will be used at the same time. This should be balanced for memory and efficiency.')
        form.addParallelSection(threads=1, mpi=0)
        # You need a params to belong to a section:
        form.addSection(label='ProDy ANM MC')

        form.addParam('numberOfWalks', params.IntParam, default=2,
                      label='Number of random walks',
                      help='This protocol can run multiple random walks in parallel or serial depending on threads')

        form.addParam('startingStructure', params.PointerParam, label="Starting structure",
                      important=True,
                      pointerClass='AtomStruct',
                      help='The starting structure should have just representative atoms '
                            '(e.g. CA) for efficient normal mode analysis')

        form.addParam('targetStructure', params.PointerParam, label="Target structure (optional)",
                      allowsNull=True,
                      pointerClass='AtomStruct',
                      help='The target structure, if provided, should have matching atoms '
                            'to the starting structure')
        


    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        stepIds = []
        for i in range(self.numberOfWalks.get()):
            comp = self._insertFunctionStep('computeStep', i,
                                            prerequisites=[],
                                            needsGPU=False)
            outputs = self._insertFunctionStep('createOutputsStep', i,
                                               prerequisites=comp, 
                                               needsGPU=False)
            
            stepIds.append(outputs)

        self._insertFunctionStep('createOutputStep',
                                 prerequisites=stepIds, needsGPU=False)

    def computeStep(self):

        args = f"{self.startingStructure.get().getFileName()}"
        self.runJob(Plugin.getProgram(
            'anmmc.py',
            location=os.path.join(prody.__path__, 'dynamics/comd.py')), args)

    def createOutputsStep(self, i):

        suffix = str(i+1)
        direc = self._getDirectory(i)

        structs = SetOfAtomStructs.create(self._getExtraPath(), suffix=suffix)
        for filename in sorted(os.listdir(os.path.join(direc, 'pdbs'))):
            pdb = AtomStruct(os.path.join(direc, 'pdbs', filename))
            structs.append(pdb)

        if not os.path.exists(os.path.join(direc, 'weights.txt')):
            args = '--path {0} --filename {1}'.format(direc, 'pdbs.ens.npz')
            self.runJob(Plugin.getProgram('ensemble_weights.py', script=True), args)

        self.weights = np.loadtxt(os.path.join(direc, 'weights.txt'))
        if self.weights.ndim == 0:
            self.weights = self.weights.reshape(-1)

        if self.doFitting.get():
            self.ccs = np.loadtxt(os.path.join(direc, 'pdbs_cc.txt'))
            if self.ccs.ndim == 0:
                self.ccs = self.ccs.reshape(-1)

        self.labels = np.loadtxt(os.path.join(direc, 'labels.txt'), dtype=str)
        if self.labels.ndim == 0:
            self.labels = self.labels.reshape(-1)
        if len(self.labels) == 0:
            self.labels = np.arange(len(self.weights))

        outSetAS = SetOfAtomStructs().create(self._getPath(), suffix=suffix)
        outSetAS.copyItems(structs, updateItemCallback=self._setWeightsCCs)

        self.args["outputStructures" + suffix] = outSetAS

        self.ensBaseName = os.path.join(direc, 'pdbs')
        npz = ProDyNpzEnsemble().create(self._getExtraPath(), suffix=suffix)
        for j in range(len(self.weights)):
            if self.doFitting.get():
                frame = TrajFrame((j+1, self.ensBaseName+'.ens.npz'),
                                  objLabel=self.labels[j],
                                  weight=self.weights[j],
                                  cc=self.ccs[j])
            else:
                frame = TrajFrame((j+1, self.ensBaseName+'.ens.npz'),
                                  objLabel=self.labels[j],
                                  weight=self.weights[j])
            npz.append(frame)

        self.args["outputNpz" + suffix] = npz

    def _setWeightsCCs(self, item, row=None):
        weight = pwobj.Float(self.weights[item.getObjId()-1])
        setattr(item, ENSEMBLE_WEIGHTS, weight)

        if self.doFitting:
            cc = pwobj.Float(self.ccs[item.getObjId()-1])
            setattr(item, ENSEMBLE_CCS, cc)

    def createOutputStep(self):
        self._defineOutputs(**self.args)