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
This module will provide ProDy normal mode analysis (NMA) using the anisotropic network model (ANM).
"""
import math
from multiprocessing import cpu_count
from os.path import exists, join

import prody
from prody2 import Plugin
from prody2.protocols.protocol_modes_base import ProDyModesBase

from pwem.emlib import (MetaData, MDL_NMA_MODEFILE, MDL_ORDER,
                        MDL_ENABLED, MDL_NMA_COLLECTIVITY, MDL_NMA_SCORE, 
                        MDL_NMA_ATOMSHIFT, MDL_NMA_EIGENVAL)
from pwem.objects import SetOfNormalModes, String

from pyworkflow.utils import glob, redStr
from pyworkflow.utils.path import makePath
from pyworkflow.protocol.params import (PointerParam, IntParam, FloatParam, StringParam,
                                        BooleanParam, LEVEL_ADVANCED)

vecStr = "vec.%d"

class ProDyANM(ProDyModesBase):
    """
    This protocol will perform normal mode analysis (NMA) using the anisotropic network model (ANM)
    """
    _label = 'ANM NMA'
    _possibleOutputs = {'outputModes': SetOfNormalModes}

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        cpus = cpu_count()//2 # don't use everything
        form.addParallelSection(threads=cpus, mpi=0)

        # You need a params to belong to a section:
        form.addSection(label='ProDy ANM NMA')

        form.addParam('inputStructure', PointerParam, label="Input structure",
                      important=True,
                      pointerClass='AtomStruct',
                      help='The input structure can be an atomic model '
                           '(true PDB) or a pseudoatomic model\n'
                           '(an EM volume converted into pseudoatoms)')

        form.addParam('numberOfModes', IntParam, default=20,
                      label='Number of modes',
                      help='The maximum number of modes allowed by the method for '
                           'atomic normal mode analysis is 3 times the '
                           'number of nodes (Calpha atoms or pseudoatoms).')

        form.addParam('cutoff', StringParam, default=15.,
                      label="Cut-off distance (A)",
                      help='Atoms or pseudoatoms beyond this distance will not interact.\n'
                           'For Calpha atoms, the default distance of 15 A works well in the majority of cases although '
                           '18 A may sometimes be better, see Eyal et al., Bioinformatics 2006.\n'
                           'For all atoms, a shorter distance such as 5 or 7 A is recommended, see Tirion et al., Phys Rev Lett 1996.\n'
                           'For other levels of coarse-graining including pseudoatoms, see Doruker et al., J Comput Chem 2002.\n'
                           'It is also possible to use other functions for the cutoff e.g. 2.9 * math.log(numResidues) - 2.9 for ed-ENM, '
                           'replacing numResidues with the actual number of residues')
        form.addParam('gamma', StringParam, default=1.,
                      expertLevel=LEVEL_ADVANCED,
                      label="Spring constant",
                      help='This number or function determines the strength of the springs.\n'
                           'Besides pre-defined Gamma functions such as GammaStructureBased from Lezon et al., PLoS Comput Biol 2010 '
                           'and GammaED from Orellana et al., J Chem Theory Comput 2010, '
                           'more sophisticated options are available within the ProDy API and '
                           'the resulting modes can be imported back into Scipion.\n'
                           'See http://prody.csb.pitt.edu/tutorials/enm_analysis/gamma.html')
        form.addParam('sparse', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Use sparse matrices?",
                      help='This saves memory at the expense of computational time.')
        form.addParam('kdtree', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Use KDTree for building Hessian matrix?",
                      help='This takes more computational time.')

        form.addParam('membrane', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Use explicit membrane model?",
                      help='An explicit lattice elastic network is used to model the membrane. '
                      'This option requires a protein oriented with opm or ppm.')

        form.addParam('collectivityThreshold', FloatParam, default=0.15,
                      expertLevel=LEVEL_ADVANCED,
                      label='Threshold on collectivity',
                      help='Collectivity degree is related to the number of atoms or pseudoatoms that are affected by '
                      'the mode, and it is normalized between 0 and 1. Modes below this threshold are deselected in '
                      'the modes metadata file as these modes are much less collective. \n'
                      'For no deselection, this parameter should be set to 0 . \n'
                      'Zero modes 1-6 are always deselected as they are related to rigid-body movements. \n'
                      'The modes metadata file can be used to see which modes are more collective '
                      'in order to decide which modes to use at the image analysis step.')

        form.addParam('zeros', BooleanParam, default=True,
                      label="Include zero eigvals",
                      help='Elect whether modes with zero eigenvalues will be kept.')
        form.addParam('turbo', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Use turbo mode",
                      help='Elect whether to use a memory intensive, but faster way to calculate modes.')

        form.addSection(label='Animation')        
        form.addParam('rmsd', FloatParam, default=5,
                      label='RMSD Amplitude (A)',
                      help='Used only for animations of computed normal modes. '
                      'This is the maximal amplitude with which atoms or pseudoatoms are moved '
                      'along normal modes in the animations. \n')
        form.addParam('numSteps', IntParam, default=10,
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
    def _insertAllSteps(self, n=20, nzeros=6,
                        collecThreshold=0.15):
        # Insert processing steps

        # Link the input
        inputFn = self.inputStructure.get().getFileName()
        numModes = self.numberOfModes.get()

        self.gnm = False

        self._insertFunctionStep('computeModesStep', inputFn, numModes)
        self._insertFunctionStep('qualifyModesStep', numModes,
                                 self.collectivityThreshold.get())
        self._insertFunctionStep('animateModesStep', numModes,
                                 self.rmsd.get(), self.numSteps.get(),
                                 self.neg.get(), self.pos.get())
        self._insertFunctionStep('computeAtomShiftsStep', numModes)
        self._insertFunctionStep('createOutputStep')

    def computeModesStep(self, inputFn='', n=20):
        # configure ProDy to automatically handle secondary structure information and verbosity
        self.oldSecondary = prody.confProDy("auto_secondary")
        self.oldVerbosity = prody.confProDy("verbosity")

        self.pdbFileName = self._getPath('atoms.pdb')
        self.atoms = prody.parsePDB(inputFn, alt='all')
        prody.writePDB(self.pdbFileName, self.atoms)

        args = '{0} -s "all" --altloc "all"  --hessian --export-scipion --npzmatrices ' \
            '--npz -o {1} -p modes -n {2} -g {3} -c "{4}" -P {5}'.format(self.pdbFileName,
                                                                         self._getPath(), n,
                                                                         self.gamma.get(),
                                                                         self.cutoff.get(),
                                                                         self.numberOfThreads.get())

        if self.sparse.get():
            args += ' --sparse-hessian'

        if self.kdtree.get():
            args += ' --use-kdtree'

        if self.zeros.get():
            args += ' --zero-modes'
            self.startMode = 6
        else:
            self.startMode = 0

        if self.turbo.get():
            args += ' --turbo'

        if self.membrane.get():
            args += ' --membrane'
            filename = 'modes.exanm.npz'
        else:
            filename = 'modes.anm.npz'

        self.runJob(Plugin.getProgram('anm'), args)

        from pyworkflow import Config
        prodyVerbosity =  'none' if not Config.debugOn() else 'debug'
        prody.confProDy(auto_secondary=True, verbosity='{0}'.format(prodyVerbosity))
        
        self.outModes = prody.loadModel(self._getPath(filename))

    def qualifyModesStep(self, numberOfModes, collectivityThreshold=0.15, suffix=''):
        self._enterWorkingDir()

        fnVec = glob("modes/vec.*")

        if len(fnVec) < numberOfModes:
            msg = "There are only %d modes instead of %d. "
            msg += "Check the number of modes you asked to compute and/or consider increasing cut-off distance. "
            msg += "The maximum number of modes allowed by the method for ANM normal mode analysis is "
            msg += "3 times the number of nodes (atoms or pseudoatoms; %d). "
            self.warning(redStr(msg % (len(fnVec), numberOfModes, self.atoms.numAtoms()*3)))

        mdOut = MetaData()
        collectivityList = list(prody.calcCollectivity(self.outModes))
        eigvals = self.outModes.getEigvals()

        for n in range(len(fnVec)):
            collectivity = collectivityList[n]

            objId = mdOut.addObject()
            modefile = self._getPath("modes", vecStr % (n + 1))
            mdOut.setValue(MDL_NMA_MODEFILE, modefile, objId)
            mdOut.setValue(MDL_ORDER, int(n + 1), objId)

            eigval = eigvals[n]
            mdOut.setValue(MDL_NMA_EIGENVAL, eigval, objId)

            if eigval > prody.utilities.ZERO:
                mdOut.setValue(MDL_ENABLED, 1, objId)
            else:
                mdOut.setValue(MDL_ENABLED, -1, objId)

            mdOut.setValue(MDL_NMA_COLLECTIVITY, collectivity, objId)
            
            if collectivity < collectivityThreshold:
                mdOut.setValue(MDL_ENABLED, -1, objId)

        idxSorted = [i[0] for i in sorted(enumerate(collectivityList), key=lambda x: x[1], reverse=True)]

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
        mdOut.write("modes%s.xmd" % suffix)

        self._leaveWorkingDir()
        
        prody.writeScipionModes(self._getPath(), self.outModes, scores=score, only_sqlite=True,
                                collectivityThreshold=collectivityThreshold)

    def createOutputStep(self):
        fnSqlite = self._getPath('modes.sqlite')
        nmSet = SetOfNormalModes(filename=fnSqlite)
        nmSet._nmdFileName = String(self._getPath('modes.nmd'))

        inputPdb = self.inputStructure.get()
        nmSet.setPdb(inputPdb)

        self._defineOutputs(outputModes=nmSet)
        self._defineSourceRelation(self.inputStructure, nmSet)

    def _summary(self):
        if not hasattr(self, 'outputModes'):
            summ = ['Output modes not ready yet']
        else:
            modes = prody.parseScipionModes(self.outputModes.getFileName())

            summ = ['*{0}* ANM modes calculated for *{1}* nodes'.format(
                    modes.numModes(), modes.numAtoms())]
        return summ

