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
This module will provide the ClustENM hybrid simulation method from ProDy, combining clustering, ENM NMA and MD.
"""

from multiprocessing import cpu_count

from pwem import *
from pwem.objects import AtomStruct
from pwem.protocols import EMProtocol

from prody2.objects import ProDyNpzEnsemble, TrajFrame

from pyworkflow.utils import *
from pyworkflow.protocol.params import (PointerParam, IntParam, FloatParam, StringParam,
                                        BooleanParam, EnumParam, LEVEL_ADVANCED)

import prody

IMP = 0
EXP = 1

class ProDyClustENM(EMProtocol):
    """
    This protocol will provide the ClustENM and ClustENMD hybrid simulation methods from ProDy, combining clustering, ENM NMA, minimisation and MD.
    """
    _label = 'ClustENM(D)'
    _possibleOutputs = {'outputNpz': ProDyNpzEnsemble}

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        cpus = cpu_count()//2 # don't use everything
        form.addParallelSection(threads=cpus, mpi=0)

        form.addSection(label='ClustENM(D)')
        form.addParam('inputStructure', PointerParam, label="Input structure",
                      important=True,
                      pointerClass='AtomStruct',
                      help='The input structure can be an atomic model '
                           '(true PDB) or a pseudoatomic model\n'
                           '(an EM volume converted into pseudoatoms)')
        form.addParam('numberOfModes', IntParam, default=3,
                      label='Number of modes',
                      help='The maximum number of modes allowed by the method for '
                           'atomic normal mode analysis is 3 times the '
                           'number of nodes (Calpha atoms), but we suggest 3 to 5.')
        form.addParam('n_gens', IntParam, default=5,
                      label='Number of generations',
                      help='Number of generations of NMA, clustering and refinement')     
        form.addParam('n_confs', IntParam, default=50,
                      label='Number of conformers from each existing conformer',
                      help='Number of new conformers to be generated based on any conformer '
                           'from the previous generation')    
        form.addParam('sim', BooleanParam, default=True,
                      #expertLevel=LEVEL_ADVANCED,
                      label="Whether to run a short MD simulation as well as minimisation",
                      help='This includes a heating-up phase until the desired temperature is reached plus '
                           'the numbers of steps set below.')
        form.addParam('parallel', BooleanParam, default=False,
                      label='Whether to use parallel threads for conformer generation.',
                      help='This will only affect the ENM NMA steps')  
        form.addParam('rmsd', StringParam, default="1.",
                      #expertLevel=LEVEL_ADVANCED,
                      label="Average RMSD (A) of the new conformers from source conformer",
                      help='Average RMSD of the new conformers with respect to the conformer'
                           'from which they are generated \n'
                           'A tuple of floats can be given, e.g. (1.0, 1.5, 1.5) for subsequent generations.')
        form.addParam('clusterMode', EnumParam, choices=['maxclust', 'threshold'],
                      label="Method for clustering for each generation",
                      help='Either maxclust or RMSD threshold should be given! For large number of '
                           'generations and/or structures, specifying maxclust is more efficient.')
        form.addParam('maxclust', StringParam, default="None",
                      condition='clusterMode==0',
                      label="Maximum number of clusters for each generation",
                      help='A tuple of floats can be given, e.g. (10, 30, 50) for subsequent generations.')
        form.addParam('threshold', StringParam, default="None",
                      condition='clusterMode==1',
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
                      #expertLevel=LEVEL_ADVANCED,
                      display=EnumParam.DISPLAY_HLIST,
                      help='Choose whether to use implicit or explicit solvent')
        form.addParam('padding', FloatParam, default=1.,
                      expertLevel=LEVEL_ADVANCED,
                      label="Padding distance to use for solvation",
                      help='Padding distance to use for the solvent box')
        form.addParam('ionicStrength', FloatParam, default=0.,
                      condition="solvent==%d" % EXP,
                      #expertLevel=LEVEL_ADVANCED, 
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

        form.addParam('temp', FloatParam, default=303.15,
                      expertLevel=LEVEL_ADVANCED,
                      condition="sim==True",
                      label="Temperature (K)",
                      help='Temperature (K) at which the simulations are conducted')                           
        form.addParam('t_steps_i', IntParam, default=1000,
                      expertLevel=LEVEL_ADVANCED,
                      condition="sim==True",
                      label="Number of 2 fs MD time steps for the initial starting structure",
                      help='Default value is good for reducing possible drift from the starting structure') 
        form.addParam('t_steps_g', StringParam, default="7500",
                      expertLevel=LEVEL_ADVANCED,
                      condition="sim==True",
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



    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('computeStep')
        self._insertFunctionStep('createOutputStep')

    def computeStep(self):
        atoms = prody.parsePDB(self.inputStructure.get().getFileName())

        if self.solvent.get() == IMP:
            self.solvent = 'imp'
        else:
            self.solvent = 'exp'

        nproc = self.numberOfThreads.get()
        if nproc:
            try:
                from threadpoolctl import threadpool_limits
            except ImportError:
                raise ImportError('Please install threadpoolctl to control threads')

            with threadpool_limits(limits=nproc, user_api="blas"):
                ens = prody.ClustENM(str(atoms))
                ens.setAtoms(atoms)
                ens.run(n_gens=self.n_gens.get(), n_modes=self.numberOfModes.get(),
                        n_confs=self.n_confs.get(), rmsd=eval(self.rmsd.get()),
                        cutoff=self.cutoff.get(), gamma=self.gamma.get(), 
                        maxclust=eval(self.maxclust.get()), threshold=eval(self.threshold.get()),
                        solvent=self.solvent, force_field=eval(self.force_field.get()),  
                        sim=self.sim.get(), temp=self.temp.get(),
                        t_steps_i=self.t_steps_i.get(), t_steps_g=eval(self.t_steps_g.get()),
                        outlier=self.outlier.get(), mzscore=self.mzscore.get(),
                        sparse=self.sparse.get(), kdtree=self.kdtree.get(), turbo=self.turbo.get(), 
                        parallel=self.parallel.get())

        self.outFileName = self._getPath('clustenm')
        prody.saveEnsemble(ens, self.outFileName+'.ens.npz')
        prody.writePDB(self.outFileName+'.pdb', ens)
        prody.writeDCD(self.outFileName+'.dcd', ens)

    def createOutputStep(self):
        outputPdb = AtomStruct()
        outputPdb.setFileName(self.outFileName+'.pdb')
        self._defineOutputs(outputStructure=outputPdb,
                            outputNPZ=self.npz)

    def _summary(self):
        if not hasattr(self, 'outputNPZ'):
            sum = ['Output ensemble not ready yet']
        else:
            ens = self.outputNPZ.loadEnsemble()
            sum = ['ClustENM completed *{0}* generations and generated *{1}* structures of *{2}* atoms'.format(
                    self.n_gens.get(), ens.numConfs(), ens.numAtoms())]
        return sum

