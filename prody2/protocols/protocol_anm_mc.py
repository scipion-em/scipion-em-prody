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
import os

from prody2 import Plugin, ENSEMBLE_WEIGHTS

from pwem import Config
from pwem.objects import AtomStruct, SetOfAtomStructs, Float
from pwem.protocols import EMProtocol

from pyworkflow.protocol import params

ALTERNATING = 0
ONEWAY = 1
SERIAL = 2

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
        form.addParallelSection(threads=2, mpi=0)

        # You need a params to belong to a section:
        form.addSection(label='ProDy ANM MC')

        form.addParam('numberOfWalks', params.IntParam, default=2,
                      label='Number of random walks',
                      help='This protocol can run multiple random walks in parallel or serial '
                        'depending on threads')

        form.addParam('startingStructure', params.PointerParam, label="Starting structure",
                      pointerClass='AtomStruct',
                      help='The starting structure should have just representative atoms '
                            '(e.g. CA) for efficient normal mode analysis')

        form.addParam('useTarget', params.BooleanParam, default=False,
                      label='Whether to use a target structure.',
                      help='If using a target, steps will be accepted depending on approaching it')

        form.addParam('targetStructure', params.PointerParam,
                      label="Target structure",
                      condition='useTarget==True',
                      pointerClass='AtomStruct',
                      help='The target structure should have matching atoms '
                            'to the starting structure. Steps will be accepted or rejected '
                            'with a certain probability based on an energy from '
                            'contacts agreeing with the target (depending on the acceptance ratio)')

        form.addParam('anmmcSteps', params.IntParam,
                      label="Number of ANM MC steps",
                      default=1000000,
                      help='This is a very large number of steps and should probably be reduced when combining with CoMD')

        form.addParam('useCoMD', params.BooleanParam, default=False,
                      label='Whether to use the collective MD.',
                      help='This algorithm adds targeted MD and energy minimisation for all-atom relaxation. '
                            'This could be much slower and may reduce conformational changes')

        form.addParam('comdGens', params.IntParam,
                      label="Number of CoMD generations",
                      default=6,
                      condition='useCoMD==True',
                      help='Each generation runs many steps of ANM MC and then runs targeted MD towards the '
                           'resulting structure and then minimises the output')

        form.addParam('comdDirectionMode', params.EnumParam,
                      choices=['Alternating', 'One way', 'Serial'],
                      label="CoMD direction mode",
                      default=ALTERNATING,
                      condition='useCoMD==True and useTarget==True',
                      help="Direction scheduling mode for runs starting from initial and target structures. "
                            "Alternating is classical in CoMD and Adaptive ANM, but other modes should be good too.")

        form.addParam('devi', params.FloatParam,
                      label="Maximum deviation per step (A)",
                      default=0.5,
                      help='Each step is scaled by the mode frequency and this scale factor '
                            'such that steps along the slowest mode with largest amplitude have '
                            'this step size in Angstroms')

        form.addParam('stepcutoff', params.FloatParam,
                      label="Maximum total RMSD (A)",
                      default=2.,
                      help='The random walk is stopped when the RMSD exceeds this value. '
                            'Unreasonable deformations may occur if this is too high in one walk. '
                            'It may be better to use the output of one run and the input for another '
                            'so that the normal modes are recalculated')

        form.addParam('acceptance_ratio', params.FloatParam, label="Acceptance ratio",
                      default=0.9,
                      condition='useTarget==True',
                      help='This parameter scales the probability of accepting '
                            'moves in the wrong direction')
        
        form.addParam('cutoff', params.FloatParam, default=15,
                      expertLevel=params.LEVEL_ADVANCED,
                      label="ANM cut-off distance (A)",
                      help='Atoms beyond this distance will not interact')

        form.addParam('useAllAtoms', params.BooleanParam, default=True,
                      label='Whether to use all atoms.',
                      help='Otherwise, CA atoms are selected')


    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):

        self.args = {}

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

    def computeStep(self, i):
        direc = self._getDirectory(i)
        if not os.path.exists(direc):
            os.mkdir(direc)

        startingStructureFn = self.startingStructure.get().getFileName()
        if self.useTarget.get():
            targetStructureFn = self.targetStructure.get().getFileName()
        else:
            targetStructureFn = startingStructureFn

        args = f"{startingStructureFn} {targetStructureFn} "*2 # repeats make sense with CoMD code

        args += f"{i+1} {self.devi.get()} {self.stepcutoff.get()} {self.acceptance_ratio.get()} "
        args += f"{self.cutoff.get()} {self.anmmcSteps.get()} "

        args += f"{os.path.join(direc, f'run_{i+1}_final_structure.dcd')} "
        args += f"{int(self.useAllAtoms.get())} 1 1 "  # these numbers 1 are for saving all coordinate sets and writing pdbs
        args += f"{self.useCoMD.get()} "
        if self.useCoMD.get():
            args += f"{self.comdGens.get()} "
            if self.useTarget.get():
                args += f"{self.comdDirectionMode.get()}"
            else:
                args += f"{ONEWAY}"

        self.runJob(
            Plugin.getProgram(
                'comd.py',
                location=os.path.join(
                    Config.EM_ROOT,
                    "prody-github",
                    "ProDy",
                    "prody",
                    "dynamics"
                )
            ),
            args
        )

    def createOutputsStep(self, i):

        suffix = str(i+1)
        direc = self._getDirectory(i)

        structs = SetOfAtomStructs.create(self._getExtraPath(), suffix=suffix)
        for filename in sorted(os.listdir(os.path.join(direc))):
            if filename.endswith(".pdb"):
                pdb = AtomStruct(os.path.join(direc, filename))
                structs.append(pdb)
        self.args["outputStructures" + suffix] = structs

    def createOutputStep(self):
        outSetAS = SetOfAtomStructs.create(self._getExtraPath())
        for key, output in self.args.items():
            if key.startswith("outputStructures"):
                outSetAS.copyItems(output, updateItemCallback=self._cleanIds)
        self.args["outputStructures"] = outSetAS

        self._defineOutputs(**self.args)

    def _getDirectory(self, i):
        suffix = str(i+1)
        return self._getPath('walk_{0}'.format(suffix))

    def _cleanIds(self, item, row=None):
        item.cleanObjId()
        setattr(item, ENSEMBLE_WEIGHTS, Float(1))

    def _validate(self):
        errors = []
        if not (self.tarStructure.hasValue() or (self.use_trans.get() 
                                                 and self.transformation.hasValue())):
            errors.append('A target structure or transformation matrix must be provided')

        return errors
