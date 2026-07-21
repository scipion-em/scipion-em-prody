# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     James Krieger (jamesmkrieger@gmail.com)
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
This module will provide the ClustENM(D) hybrid simulation method from ProDy, combining clustering, ENM NMA and MD.
"""

from multiprocessing import cpu_count
import numpy as np
import os

from pwem.objects import AtomStruct, SetOfAtomStructs
from pwem.protocols import EMProtocol

import pyworkflow.object as pwobj
from pyworkflow.protocol.params import (IntParam, FloatParam, StringParam, BooleanParam,
                                        EnumParam, MultiPointerParam, LEVEL_ADVANCED,
                                        USE_GPU, GPU_LIST)
from pyworkflow.protocol import STEPS_PARALLEL, STEPS_SERIAL

from prody2.constants import ENSEMBLE_WEIGHTS, ENSEMBLE_CCS
from prody2.objects import ProDyNpzEnsemble, TrajFrame
from prody2 import Plugin

IMP = 0
EXP = 1

from pyworkflow.utils import logger

class ProDyClustENM(EMProtocol):
    """
    Performs hybrid conformational sampling using the ProDy ClustENM/ClustENMD framework.

    AI Generated:

    ClustENM(D) (ProDyClustENM) — User Manual
        Overview

        The ClustENM(D) protocol explores conformational variability of one or more
        atomic structures by combining normal mode analysis (ENM), conformer generation,
        clustering, energy minimization, and optional molecular dynamics refinement.

        In practical structural biology terms, this protocol is designed to generate
        plausible alternative conformations starting from one or more input atomic
        models. It is particularly useful when studying intrinsic flexibility,
        conformational transitions, domain motions, or preparing structural ensembles
        for downstream fitting, analysis, or interpretation.

        Unlike a single minimization or a standard MD trajectory, ClustENM(D)
        generates multiple branches of structural alternatives across successive
        generations. This makes it especially valuable when the biological goal is to
        explore possible motions rather than only refine a single structure.

        Inputs and General Workflow

        The protocol requires one or more input atomic structures. Each structure is
        processed independently.

        For every input structure, the workflow follows this general scheme:

            1. Perform elastic network normal mode analysis.
            2. Generate new conformers along selected low-frequency modes.
            3. Cluster the generated conformers.
            4. Refine representative conformers by minimization.
            5. Optionally perform short molecular dynamics simulations.
            6. Repeat this process for the desired number of generations.

        The result is a structured ensemble of conformations that represent accessible
        motions around the starting structure.

        Number of Modes and Conformer Generation

        The number of normal modes determines how many collective motions are used to
        generate structural perturbations.

        In biological applications, low-frequency modes often correspond to large-scale
        collective motions such as hinge bending, domain rearrangements, or breathing
        motions. Because of this, a small number of modes (commonly 3–5) is usually
        sufficient for exploratory analyses.

        The number of conformers controls how broadly each conformational branch is
        sampled. Higher values increase diversity but also increase computational cost.

        A practical strategy is:

            - small exploratory runs: few conformers and few generations
            - broader sampling: more conformers and more generations

        RMSD Sampling Amplitude

        The RMSD parameter defines the average displacement of newly generated
        conformers relative to their parent structure.

        Biologically, this controls how far the exploration moves away from the
        current conformation.

            - small RMSD values favor local exploration
            - larger RMSD values allow broader conformational searches

        Different RMSD values can be assigned to successive generations. This is often
        useful when beginning with broader exploration and gradually refining later
        generations.

        Clustering Strategy

        After conformer generation, structures are clustered to remove redundancy and
        retain representative states.

        Two clustering strategies are available:

            - maxclust:
              limits the maximum number of clusters. This is generally more efficient
              for large searches.

            - threshold:
              groups structures according to an RMSD cutoff. This can be useful when
              structural similarity has a clear biological meaning.

        In most practical workflows, maxclust is often easier to control when many
        conformers or multiple generations are used.

        Normal Mode Analysis Parameters

        Several advanced parameters control ENM behavior:

            - gamma:
              spring constant controlling interaction strength

            - cutoff:
              distance threshold defining which Cα atoms interact

            - sparse:
              saves memory at the cost of longer computation

            - kdtree:
              alternative Hessian construction strategy

            - turbo:
              faster but more memory-demanding mode calculation

        For most biological applications, the default cutoff and gamma values are
        suitable starting points unless very unusual systems are being studied.

        Simulation and Refinement

        After conformer generation, structures can be refined using energy
        minimization and optional molecular dynamics.

        If simulation is enabled, the protocol performs:

            - minimization
            - optional heating
            - short MD sampling

        This helps remove unrealistic distortions introduced by mode-based
        perturbations and improves physical plausibility.

        Solvent Models

        Two solvent models are available:

            - implicit solvent:
              computationally cheaper and generally suitable for exploratory
              conformational sampling

            - explicit solvent:
              more realistic but significantly more expensive

        For most routine conformational exploration, implicit solvent is usually the
        preferred starting choice.

        Explicit solvent becomes more relevant when physical detail is especially
        important, for example when local side-chain packing or solvent-mediated
        effects may matter.

        Outlier Filtering

        In implicit solvent mode, conformers with unusually unfavorable energies can
        be filtered automatically using modified z-score statistics.

        From a biological perspective, this helps remove highly distorted or unstable
        structures that are less likely to represent meaningful conformational states.

        The outlier threshold should usually be kept conservative unless aggressive
        sampling is intentionally being performed.

        Optional Fitting to Experimental Volumes

        The protocol can optionally filter generated conformers against one or more
        experimental volumes.

        In this mode, simulated density maps are generated from candidate conformers
        and compared against target maps.

        This becomes particularly useful in cryo-EM workflows when one wants to
        retain only conformers compatible with experimental density.

        Practical biological applications include:

            - exploring flexible fitting candidates
            - selecting conformers consistent with low-resolution maps
            - filtering out conformers that deviate too strongly from experimental data

        If enabled, the protocol can optionally resample conformers to replace those
        rejected during the fitting stage.

        Parallel Execution

        Conformer generation can be parallelized across CPU threads.

        This primarily accelerates ENM/NMA-based sampling and is especially useful
        when processing multiple structures or larger conformational searches.

        Outputs and Their Interpretation

        For each input structure, the protocol produces:

            - outputStructuresN:
              a set of sampled atomic conformers

            - outputNpzN:
              an ensemble representation of the same conformers with associated weights

        Each conformer receives a weight derived from the ensemble statistics.

        Biologically, these weights can be interpreted as relative representation
        within the sampled ensemble, although they should not automatically be treated
        as rigorous thermodynamic populations.

        Practical Recommendations

        For exploratory biological studies, a good starting strategy is:

            - 2 generations
            - 3 to 5 modes
            - moderate RMSD (~1 Å)
            - implicit solvent
            - modest clustering

        For broader conformational searches:

            - increase number of generations
            - increase number of conformers
            - gradually tune RMSD and clustering thresholds

        When fitting to cryo-EM maps, careful attention should be paid to map
        resolution and threshold selection, since overly aggressive filtering may
        discard biologically relevant alternatives.

        Final Perspective

        ClustENM(D) is best understood not as a conventional molecular dynamics
        protocol, but as a structured conformational exploration framework.

        For structural biologists, its main value lies in efficiently sampling
        physically plausible collective motions that may correspond to biologically
        meaningful functional transitions.

        When used carefully, it provides an effective bridge between coarse-grained
        normal mode analysis and more detailed atomistic refinement.
    """
    _label = 'ClustENM(D)'
    _possibleOutputs = {'outputStructures1': SetOfAtomStructs,
                        'outputNpz1': ProDyNpzEnsemble}
    stepsExecutionMode = STEPS_PARALLEL

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        cpus = cpu_count()//2 # don't use everything
        form.addParam('binThreads', IntParam,
                      label='threads',
                      default=cpus,
                      help='Number of threads used by ProDy each time it is called in the protocol execution. For '
                           'example, if 3 Scipion threads and 6 ProDy threads are set, the structures will be '
                           'processed in groups of 2 at the same time with a call of ProDy with 6 threads each, so '
                           '12 threads will be used at the same time. Beware the memory of your machine has '
                           'memory enough to load together the number of structures specified by Scipion threads.')
        form.addParallelSection(threads=1, mpi=0)

        form.addSection(label='ClustENM(D)')
        form.addParam('inputStructures', MultiPointerParam, label="Input structures",
                      important=True,
                      pointerClass='AtomStruct',
                      help='Each input structures should be an atomic model')
        form.addParam('mergeInputs', EnumParam,
                      choices=['Separate run per structure', 'Merge into one multi-start run'],
                      default=0, display=EnumParam.DISPLAY_HLIST,
                      label="Handling of multiple input structures",
                      help='Separate: run ClustENM(D) independently on each input structure, giving one '
                           'output ensemble per structure (the original behaviour).\n'
                           'Merge: seed a single ClustENM(D) run with all the input structures at once '
                           '(multi-start) -- they form the initial population together and the generations '
                           'then proceed jointly -- giving one combined output ensemble. The inputs must '
                           'share the same topology (same atoms after fixing) to be merged.')
        form.addParam('numberOfModes', IntParam, default=3,
                      label='Number of modes',
                      help='The maximum number of modes allowed by the method for '
                           'atomic normal mode analysis is 3 times the '
                           'number of nodes (Calpha atoms), but we suggest 3 to 5.')
        form.addParam('n_gens', IntParam, default=2,
                      label='Number of generations',
                      help='Number of generations of NMA, clustering and refinement')     
        form.addParam('n_confs', IntParam, default=10,
                      label='Number of conformers from each existing conformer',
                      help='Number of new conformers to be generated based on any conformer '
                           'from the previous generation')    
        form.addParam('sim', BooleanParam, default=True,
                      label="Whether to run a short MD simulation as well as minimisation",
                      help='This includes a heating-up phase until the desired temperature is reached plus '
                           'the numbers of steps set below.')
        form.addParam('parallel', BooleanParam, default=False,
                      label='Whether to use parallel threads for conformer generation.',
                      help='This will only affect the ENM NMA steps')  
        form.addParam('rmsd', StringParam, default="1.",
                      label="Average RMSD (A) of the new conformers from source conformer",
                      help='Average RMSD of the new conformers with respect to the conformer'
                           'from which they are generated \n'
                           'A tuple of floats can be given, e.g. (1.0, 1.5, 1.5) for subsequent generations.')
        form.addParam('clusterMode', EnumParam, choices=['maxclust', 'threshold'], default=0,
                      label="Method for clustering for each generation",
                      help='Either maxclust or RMSD threshold should be given! For large number of '
                           'generations and/or structures, specifying maxclust is more efficient.')
        form.addParam('maxclust', StringParam,
                      condition='clusterMode==0',
                      default='None',
                      label="Maximum number of clusters for each generation",
                      help='A tuple of floats can be given, e.g. (10, 30, 50) for subsequent generations.')
        form.addParam('threshold', StringParam, condition='clusterMode==1',
                      default='None',
                      label="RMSD threshold (A) to apply when forming clusters",
                      help='A tuple of floats can be given, e.g. (1.0, 1.5, 1.5) for subsequent generations.\n'
                           'This parameter has been used in ClustENMv1, setting it to 75%% of the maximum RMSD for sampling. '
                           'For the current version (v2), this should be chosen carefully for efficiency')                          

        form.addSection(label='NMA')
        form.addParam('gamma', FloatParam, default=1.,
                      expertLevel=LEVEL_ADVANCED,
                      label="Spring constant",
                      help='This number or function determines the strength of the springs.\n'
                           'More sophisticated options are available within the ProDy API and '
                           'the resulting modes can be imported back into Scipion.\n'
                           'See http://http://www.bahargroup.org/prody/tutorials/enm_analysis/gamma.html')
        form.addParam('cutoff', FloatParam, default=15.,
                      expertLevel=LEVEL_ADVANCED,
                      label="Cut-off distance (A)",
                      help='Calpha atoms beyond this distance will not interact. \n'
                           'The default distance of 15 A works well in the majority of cases.')
        form.addParam('sparse', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Use sparse matrices?",
                      help='This saves memory at the expense of computational time.')
        form.addParam('kdtree', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Use KDTree for building Hessian matrix?",
                      help='This takes more computational time.')
        form.addParam('turbo', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Use turbo mode",
                      help='Elect whether to use a memory intensive, but faster way to calculate modes.')

        form.addSection(label='Simulation')
        form.addParam('solvent', EnumParam, choices=['implicit', 'explicit'],
                      label="Solvent type", default=IMP,
                      display=EnumParam.DISPLAY_HLIST,
                      help='Choose whether to use implicit or explicit solvent')
        form.addParam('padding', FloatParam, default=1.,
                      expertLevel=LEVEL_ADVANCED,
                      label="Padding distance to use for solvation",
                      help='Padding distance to use for the solvent box')
        form.addParam('ionicStrength', FloatParam, default=0.,
                      condition="solvent==%d" % EXP,
                      label="Total concentration of ions (both positive and negative) to add in mol/L",
                      help='This does not include ions that are added to neutralize the system.')
        form.addParam('force_field', StringParam, default="None",
                      expertLevel=LEVEL_ADVANCED,
                      label="Force field to use",
                      help="If set to None, Implicit solvent force field is ('amber99sbildn.xml', 'amber99_obc.xml') "
                           "and Explicit solvent force field is ('amber14-all.xml', 'amber14/tip3pfb.xml').\n"
                           "Any other tuple of protein and water force fields implemented in OpenMM should work.")

        form.addParam('tolerance', FloatParam, default=10.,
                      expertLevel=LEVEL_ADVANCED,
                      label="Energy tolerance to which the system should be minimized in kJ/mole",
                      help='Energy tolerance for stopping energy minimisation')
        form.addParam('maxIterations', IntParam, default=0,
                      expertLevel=LEVEL_ADVANCED,
                      label="Maximum number of iterations to perform during energy minimization",
                      help='If this is 0 (default), minimization is continued until the results converge without '
                           'regard to how many iterations it takes')

        form.addParam('parallelSim', IntParam, default=1,
                      condition="mergeInputs==1",
                      label="Parallel simulation workers",
                      help='Only for a merged (multi-start) run. Number of worker processes to run the '
                           'per-conformer energy minimisation/MD in parallel (1 = off, the default; 0 = as '
                           'many as CPUs). This is independent of the conformer-generation "parallel" option. '
                           'When GPUs are selected the workers are spread one-per-GPU across the chosen GPUs. '
                           'Most useful for MD or large merged ensembles; for short minimisation the '
                           'per-worker start-up (a fresh OpenMM/CUDA context) can outweigh the gain.')

        simTrue = "sim==True"
        form.addParam('temp', FloatParam, default=303.15,
                      expertLevel=LEVEL_ADVANCED,
                      condition=simTrue,
                      label="Temperature (K)",
                      help='Temperature (K) at which the simulations are conducted')                           
        form.addParam('t_steps_i', IntParam, default=1000,
                      expertLevel=LEVEL_ADVANCED,
                      condition=simTrue,
                      label="Number of 2 fs MD time steps for the initial starting structure",
                      help='Default value is good for reducing possible drift from the starting structure') 
        form.addParam('t_steps_g', StringParam, default="7500",
                      expertLevel=LEVEL_ADVANCED,
                      condition=simTrue,
                      label="Number of 2 fs MD time steps for each conformer from each generation",
                      help="A tuple of integers can be given, e.g. (3000, 5000, 7000) for subsequent generations.")

        form.addParam('outlier', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      condition="solvent==%d" % IMP,
                      label="Exclude conformers detected as outliers in each generation?",
                      help="Outliers, if any, are detected by the modified z-scores of the conformers' "
                           "potential energies over a generation.\n"
                           "Note: It is automatically set to False when explicit solvent is being used")                      
        form.addParam('mzscore', FloatParam, default=3.5,
                      expertLevel=LEVEL_ADVANCED,
                      condition="outlier==True",
                      label="Modified z-score threshold to label conformers as outliers",
                      help='Modified z-score threshold to label conformers as outliers')   

        form.addSection(label='Fitting')
        form.addParam('doFitting', BooleanParam, default=False,
                      label="Whether to do fitting to volumes like MDeNM-EMFit?",
                      help="If selected, this will filter structures to those that do not reduce the cross-correlation much")
        fittingCondition = 'doFitting==True'
        form.addParam('inputVolumes', MultiPointerParam, label="Target volumes",
                      important=True, allowsNull=True, condition=fittingCondition,
                      pointerClass='Volume',
                      help='If fitting, there should be the same number of volumes as models or just one for all of them')
        form.addParam('fitResolution', FloatParam, default=5.,
                      expertLevel=LEVEL_ADVANCED,
                      condition=fittingCondition,
                      label="Resolution for simulated volumes (A)",
                      help='Resolution (A) for simulating volumes to compare against the target')
        form.addParam('replaceFiltered', BooleanParam, default=False, condition=fittingCondition,
                      label="Whether to sample again to replace filtered conformations?",
                      help="If selected, this will sample and filter structures repeatedly until the selected number are kept")
        form.addParam('mapCutoff', FloatParam, default=0.1,
                      expertLevel=LEVEL_ADVANCED,
                      condition=fittingCondition,
                      label="Intensity threshold for target maps",
                      help='Minimum intensity cutoff for reading target maps to avoid noise')

        form.addHidden(USE_GPU, BooleanParam, default=False,
                       label="Use GPU for execution",
                       help="This protocol has both CPU and GPU implementation.\
                       Select the one you want to use. Be aware that the GPU program is new and may have problems")

        form.addHidden(GPU_LIST, StringParam, default='0',
                       expertLevel=LEVEL_ADVANCED,
                       label="Choose GPU IDs",
                       help="Add a list of GPU devices that can be used")

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):

        self.args = {}

        # Insert processing steps
        pdbs = [struct.get().getFileName() for struct in self.inputStructures]
        merge = self.mergeInputs.get() == 1

        if self.doFitting.get():
            if self.inputVolumes is not None:
                self.volumes = [vol.get().getFileName() for vol in self.inputVolumes]
            else:
                self.volumes = []

            if len(self.volumes) < len(pdbs) and len(self.volumes) != 1:
                if len(self.volumes) != 0:
                    logger.warning("Ignoring volumes as the number of them does not match structures.")
                self.volumes = None

            if not merge and len(pdbs) == 1 and self.volumes and len(self.volumes) > 1:
                pdbs = [pdbs[0] for _ in self.volumes]

        # merged multi-start: a single run seeded with all structures (comma-separated, parsed by the
        # clustenm app); otherwise one run per structure (the original behaviour).
        if merge:
            runInputs = [",".join(pdbs)]
            self.stepsExecutionMode = STEPS_SERIAL
        else:
            runInputs = pdbs
            if len(runInputs) == 1:
                self.stepsExecutionMode = STEPS_SERIAL

        stepIds = []
        for i, pdb in enumerate(runInputs):
            comp = self._insertFunctionStep('computeStep', i, pdb,
                                            prerequisites=[],
                                            needsGPU=self.useGpu)
            outputs = self._insertFunctionStep('createOutputsStep', i,
                                               prerequisites=comp,
                                               needsGPU=False)

            stepIds.append(outputs)

        self._insertFunctionStep('createOutputStep',
                                 prerequisites=stepIds, needsGPU=False)

    def computeStep(self, i, pdb):
        direc = self._getDirectory(i)
        if not os.path.exists(direc):
            os.mkdir(direc)

        args = '"{0}" --ngens {1} --number-of-modes {2} --nconfs {3} --rmsd {4} -c {5} -g {6} ' \
               '--solvent {7} --force_field {8} --ionicStrength {9} --padding {10} --temp {11} ' \
                '--t_steps_i {12} --t_steps_g {13} --select all ' \
               '--tolerance {14} --maxIterations {15} -o {16} --file-prefix pdbs --multiple -P {17}'.format(
                   pdb, self.n_gens.get(), self.numberOfModes.get(),
                   self.n_confs.get(), self.rmsd.get(), self.cutoff.get(), self.gamma.get(),
                   self._getSolvent(), self.force_field.get(), self.ionicStrength.get(), self.padding.get(),
                   self.temp.get(), self.t_steps_i.get(), self.t_steps_g.get(),
                   self.tolerance.get(), self.maxIterations.get(), direc, self.numberOfThreads.get())

        if self.n_gens.get() > 0:
            args += ' --maxclust "{0}" --threshold "{1}"'.format(self.maxclust.get(), self.threshold.get())

        if self.sim.get() is False:
            args += ' --no-sim'

        if self.sparse.get():
            args += ' --sparse-hessian'

        if self.kdtree.get():
            args += ' --use-kdtree'

        if self.turbo.get():
            args += ' --turbo'

        if self.parallel.get():
            args += ' --parallel'

        if self.outlier.get():
            args += ' --mzscore {0}'.format(self.mzscore.get())
        else:
            args += ' --no-outlier'

        if self.doFitting.get():
            args += ' --fitmap {0} --fit_resolution {1} --map_cutoff {2}'.format(
                self.volumes[i], self.fitResolution.get(), self.mapCutoff.get())

            if self.replaceFiltered.get():
                args += ' --replace_filtered'

        if self.useGpu:
            gpuId = self._stepsExecutor.getGpuList()
            if isinstance(gpuId, int):
                gpuStr = str(gpuId)
            else:
                gpuStr = ','.join([str(g) for g in gpuId])

            args += ' --platform CUDA --gpu-ids {0}'.format(gpuStr)
        else:
            args += ' --platform CPU'

        # merged multi-start run: optionally parallelise the per-conformer minimisation/MD across
        # worker processes, spread one-per-GPU over the selected GPUs (0-based DeviceIndex within the
        # visible set that Scipion exposes via CUDA_VISIBLE_DEVICES).
        if self.mergeInputs.get() == 1 and self.parallelSim.get() != 1:
            args += ' --parallel_sim {0}'.format(self.parallelSim.get())
            if self.useGpu:
                ngpu = len(gpuStr.split(','))
                args += ' --sim_devices {0}'.format(','.join(str(k) for k in range(ngpu)))

        if not os.path.exists(os.path.join(direc, 'pdbs.ens.npz')):
            self.runJob('export OPENMM_CPU_THREADS={0} && '.format(
                self.numberOfThreads.get()
                ) + Plugin.getProgram('clustenm'), args)

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

    def _summary(self):
        if not hasattr(self, 'outputStructures1'):
            summ = ['Output not ready yet']
        elif self.mergeInputs.get() == 1:
            summ = ['ClustENM completed *{0}* generations for a single merged (multi-start) run of '
                    '*{1}* input structures'.format(self.n_gens.get(), len(self.inputStructures))]
        else:
            numStructs = (self.numberOfSteps-1)//2 # 1 collect output step, 2 steps per struct
            if numStructs == 1:
                summ = ['ClustENM completed *{0}* generations for *{1}* structure'.format(
                        self.n_gens.get(), numStructs)]
            else:
                summ = ['ClustENM completed *{0}* generations for *{1}* structures'.format(
                        self.n_gens.get(), numStructs)]
        return summ

    def _getDirectory(self, i):
        suffix = str(i+1)
        return self._getPath('clustenm_{0}'.format(suffix))

    def _getSolvent(self):
        if self.solvent.get() == IMP:
            return 'imp'
        
        return 'exp'
