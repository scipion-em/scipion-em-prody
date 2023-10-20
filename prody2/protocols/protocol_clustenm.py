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
This module will provide the ClustENM(D) hybrid simulation method from ProDy, combining clustering, ENM NMA and MD.
"""

from multiprocessing import cpu_count
import os

from pwem.objects import AtomStruct, SetOfAtomStructs
from pwem.protocols import EMProtocol

import pyworkflow.object as pwobj
from pyworkflow.protocol.params import (IntParam, FloatParam, StringParam, BooleanParam,
                                        EnumParam, MultiPointerParam, LEVEL_ADVANCED)

import prody
from prody2.constants import ENSEMBLE_WEIGHTS
from prody2.objects import ProDyNpzEnsemble, TrajFrame
from prody2 import Plugin

IMP = 0
EXP = 1

import logging
logger = logging.getLogger(__name__)

class ProDyClustENM(EMProtocol):
    """
    This protocol will provide the ClustENM and ClustENMD hybrid simulation methods from ProDy, combining clustering, ENM NMA, minimisation and MD.
    """
    _label = 'ClustENM(D)'
    _possibleOutputs = {'outputTraj1': SetOfAtomStructs}

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        cpus = cpu_count()//2 # don't use everything
        form.addParallelSection(threads=cpus, mpi=0)

        form.addSection(label='ClustENM(D)')
        form.addParam('inputStructures', MultiPointerParam, label="Input structures",
                      important=True,
                      pointerClass='AtomStruct',
                      help='Each input structures should be an atomic model')
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
                           'See http://prody.csb.pitt.edu/tutorials/enm_analysis/gamma.html')
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
        form.addParam('mapCutoff', FloatParam, default=0.1,
                      expertLevel=LEVEL_ADVANCED,
                      condition=fittingCondition,
                      label="Intensity threshold for target maps",
                      help='Minimum intensity cutoff for reading target maps to avoid noise')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):

        self.args = {}

        # Insert processing steps
        pdbs = [struct.get().getFileName() for struct in self.inputStructures]

        if self.inputVolumes is not None:
            self.volumes = [vol.get().getFileName() for vol in self.inputVolumes]
        else:
            self.volumes = []

        if len(self.volumes) < len(pdbs) and len(self.volumes) != 1:
            if len(self.volumes) != 0:
                logger.warning("Ignoring volumes as the number of them does not match structures.")
            self.volumes = None


        if self.solvent.get() == IMP:
            self.solvent = 'imp'
        else:
            self.solvent = 'exp'

        for i, pdb in enumerate(pdbs):
            self._insertFunctionStep('computeStep', i, pdb)

        self._insertFunctionStep('createOutputStep')

    def computeStep(self, i, pdb):

        suffix = str(i+1)
        direc = self._getPath('clustenm_{0}'.format(suffix))
        if not os.path.exists(direc):
            os.mkdir(direc)

        args = '{0} --ngens {1} --number-of-modes {2} --nconfs {3} --rmsd {4} -c {5} -g {6} --maxclust "{7}" --threshold "{8}" ' \
               '--solvent {9} --force_field {10} --ionicStrength {11} --padding {12} --temp {13} --t_steps_i {14} --t_steps_g {15} ' \
               '--tolerance {16} --maxIterations {17} -o {18} --file-prefix pdbs --multiple -P {19}'.format(pdb, self.n_gens.get(), self.numberOfModes.get(),
                    self.n_confs.get(), self.rmsd.get(), self.cutoff.get(), self.gamma.get(), 
                    self.maxclust.get(), self.threshold.get(),
                    self.solvent, self.force_field.get(), self.ionicStrength.get(), self.padding.get(),
                    self.temp.get(), self.t_steps_i.get(), self.t_steps_g.get(),
                    self.tolerance.get(), self.maxIterations.get(), direc, self.numberOfThreads.get())
        
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

        if self.volumes is not None:
            args += ' --fitmap {0} --fit_resolution {1} --map_cutoff {2}'.format(self.volumes[i],
                                                                                 self.fitResolution.get(),
                                                                                 self.mapCutoff.get())

        self.runJob(Plugin.getProgram('clustenm'), args)

        structs = SetOfAtomStructs.create(self._getExtraPath())
        for filename in os.listdir(os.path.join(direc, 'pdbs')):
            pdb = AtomStruct(os.path.join(direc, 'pdbs', filename))
            structs.append(pdb)

        ens = prody.loadEnsemble(os.path.join(direc, 'pdbs.ens.npz'))
        self.weights = ens.getSizes()

        outSetAS = SetOfAtomStructs().create(self._getPath(), suffix=suffix)
        outSetAS.copyItems(structs, updateItemCallback=self._setWeights)
        self.args["outputStructures" + suffix] = outSetAS

        self.ensBaseName = os.path.join(direc, 'pdbs')
        npz = ProDyNpzEnsemble().create(self._getExtraPath(), suffix=suffix)
        for j in range(ens.numCoordsets()):
            frame = TrajFrame((j+1, self.ensBaseName+'.ens.npz'), objLabel=ens.getLabels()[j])
            npz.append(frame)

        outNpz = ProDyNpzEnsemble().create(self._getPath(), suffix=suffix)
        outNpz.copyItems(npz, updateItemCallback=self._setWeights)
        self.args["outputNpz" + suffix] = outNpz

    def _setWeights(self, item, row=None):
            weight = pwobj.Integer(self.weights[item.getObjId()-1])
            setattr(item, ENSEMBLE_WEIGHTS, weight)

    def createOutputStep(self):
        self._defineOutputs(**self.args)

    def _summary(self):
        if not hasattr(self, 'outputStructures1'):
            summ = ['Output not ready yet']
        else:
            summ = ['ClustENM completed *{0}* generations for *{1}* structures'.format(
                    self.n_gens.get(), self.numberOfSteps-1)]
        return summ

