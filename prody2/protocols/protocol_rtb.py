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
This module will provide ProDy normal mode analysis (NMA) using the the rotation and translation of blocks (RTB) framework.
"""
from pyworkflow.protocol import params

from os.path import basename, exists, join
import math

from pwem import *
from pwem.emlib import (MetaData, MDL_NMA_MODEFILE, MDL_ORDER,
                        MDL_ENABLED, MDL_NMA_COLLECTIVITY, MDL_NMA_SCORE, 
                        MDL_NMA_ATOMSHIFT, MDL_NMA_EIGENVAL)
from pwem.objects import AtomStruct, SetOfNormalModes, String
from pwem.protocols import EMProtocol

from pyworkflow.utils import *
from pyworkflow.utils.path import makePath
from pyworkflow.protocol.params import (PointerParam, IntParam, FloatParam, StringParam,
                                        BooleanParam, EnumParam, LEVEL_ADVANCED)

import prody

BLOCKS_FROM_RES = 0
BLOCKS_FROM_SECSTR = 1

class ProDyRTB(EMProtocol):
    """
    This protocol will perform normal mode analysis (NMA) using the rotation and translation of blocks (RTB) framework
    """
    _label = 'RTB NMA'
    _possibleOutputs = {'outputModes': SetOfNormalModes}

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label='ProDy RTB NMA')

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

        form.addParam('blockDef', EnumParam, choices=['res', 'secstr'],
                      label="Block definition type",
                      default=BLOCKS_FROM_RES,
                      display=EnumParam.DISPLAY_HLIST,
                      help='Define blocks using either a number of residues or secondary structure information')

        form.addParam('res_per_block', IntParam, default=10,
                      condition='blockDef==%d' % BLOCKS_FROM_RES,
                      label="Number of residues per block",
                      help='All blocks will have this number of residues except the last one')

        form.addParam('shortest_block', IntParam, default=4,
                      label='Number of residues in shortest block',
                      help='Blocks with fewer residues will be combined into the previous block. '
                           'Fewer than 4 can be problematic.')

        form.addParam('longest_block', IntParam, default=20,
                      label='Number of residues in longest block',
                      help='Blocks with more residues will be split in half')

        form.addParam('min_dist_cutoff', FloatParam, default=20.,
                      expertLevel=LEVEL_ADVANCED,
                      label='Distance cutoff for splitting blocks',
                      help='Distance of a residue from others beyond which '
                           'it is not included in the same block based on a distance tree. '
                           'This is calculated using ProDy function findSubgroups.')

        form.addParam('cutoff', FloatParam, default=15.,
                      expertLevel=LEVEL_ADVANCED,
                      label="Cut-off distance (A)",
                      help='Atoms or pseudoatoms beyond this distance will not interact.\n'
                           'For Calpha atoms, the default distance of 15 A works well in the majority of cases. '
                           'For all atoms, a shorter distance such as 5 or 7 A is recommended.\n'
                           'For fewer atoms or pseudoatoms, set this according to the level of coarse-graining '
                           '(see Doruker et al., J Comput Chem 2002 though values may differ for RTB).')

        form.addParam('gamma', FloatParam, default=1.,
                      expertLevel=LEVEL_ADVANCED,
                      label="Spring constant",
                      help='This number or function determines the strength of the springs.\n'
                           'More sophisticated options are available within the ProDy API and '
                           'the resulting modes can be imported back into Scipion.\n'
                           'See http://prody.csb.pitt.edu/tutorials/enm_analysis/gamma.html')

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
                      expertLevel=LEVEL_ADVANCED,
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
    def _insertAllSteps(self):
        # Insert processing steps

        # Link the input
        inputFn = self.inputStructure.get().getFileName()
        self.structureEM = self.inputStructure.get().getPseudoAtoms()

        self.model_type = 'rtb'
        n = self.numberOfModes.get()

        self._insertFunctionStep('computeModesStep', inputFn, n)
        self._insertFunctionStep('animateModesStep', n,
                                 self.rmsd.get(), self.n_steps.get(),
                                 self.neg.get(), self.pos.get())
        self._insertFunctionStep('qualifyModesStep', n,
                                 self.collectivityThreshold.get(),
                                 self.structureEM)
        self._insertFunctionStep('computeAtomShiftsStep', n)
        self._insertFunctionStep('createOutputStep')

    def computeModesStep(self, inputFn, n):
        # configure ProDy to automatically handle secondary structure information and verbosity
        self.oldSecondary = prody.confProDy("auto_secondary")
        self.oldVerbosity = prody.confProDy("verbosity")
        from pyworkflow import Config
        prodyVerbosity =  'none' if not Config.debugOn() else 'debug'
        prody.confProDy(auto_secondary=True, verbosity='{0}'.format(prodyVerbosity))
        
        if self.structureEM:
            self.pdbFileName = self._getPath('pseudoatoms.pdb')
        else:
            self.pdbFileName = self._getPath('atoms.pdb')

        atoms = prody.parsePDB(inputFn, alt='all', secondary=True)

        if self.blockDef.get() == BLOCKS_FROM_RES:
            self.blocks, self.amap = prody.assignBlocks(atoms, res_per_block=self.res_per_block.get(),
                                                        shortest_block=self.shortest_block.get(),
                                                        longest_block=self.longest_block.get(),
                                                        min_dist_cutoff=self.min_dist_cutoff.get())
        else:
            self.blocks, self.amap = prody.assignBlocks(atoms, secstr=True,
                                                        shortest_block=self.shortest_block.get(),
                                                        longest_block=self.longest_block.get(),
                                                        min_dist_cutoff=self.min_dist_cutoff.get())

        prody.writePDB(self.pdbFileName, self.amap)

        self.rtb = prody.RTB()
        try:
            self.rtb.buildHessian(self.amap, self.blocks, cutoff=self.cutoff.get(),
                                gamma=self.gamma.get())
        except MemoryError as err:
            prody.LOGGER.warn("{0} so using sparse matrix".format(err))
            self.rtb.buildHessian(self.amap, self.blocks, cutoff=self.cutoff.get(),
                                    gamma=self.gamma.get(), sparse=True)

        try:
            self.rtb.calcModes(n, zeros=self.zeros.get(), turbo=self.turbo.get())
        except MemoryError as err:
            prody.LOGGER.warn("{0} so using not using turbo decomposition".format(err))
            self.rtb.calcModes(n, zeros=self.zeros.get(), turbo=False)

        if self.zeros.get():
            self.startMode = 6
        else:
            self.startMode = 0
        
        prody.writeScipionModes(self._getPath(), self.rtb)
        prody.writeNMD(self._getPath('modes.nmd'), self.rtb, self.amap)
        prody.saveModel(self.rtb, self._getPath('modes.rtb.npz'), matrices=True)

    def animateModesStep(self, numberOfModes, rmsd, n_steps, pos, neg):

        animations_dir = self._getExtraPath('animations')
        makePath(animations_dir)
        for i, mode in enumerate(self.rtb[self.startMode:]):
            modenum = i+self.startMode+1
            fnAnimation = join(animations_dir, "animated_mode_%03d"
                               % modenum)
             
            self.outAtoms = prody.traverseMode(mode, self.amap, rmsd=rmsd,
                                               n_steps=n_steps,
                                               pos=pos, neg=neg)
            prody.writePDB(fnAnimation+".pdb", self.outAtoms)

            fhCmd=open(fnAnimation+".vmd",'w')
            fhCmd.write("mol new %s.pdb\n" % fnAnimation)
            fhCmd.write("animate style Rock\n")
            fhCmd.write("display projection Orthographic\n")
            if self.structureEM:
                fhCmd.write("mol modcolor 0 0 Beta\n")
                fhCmd.write("mol modstyle 0 0 Beads 1.0 8.000000\n")
            else:
                fhCmd.write("mol modcolor 0 0 Index\n")

                if self.amap.select('name P') is not None:
                    num_p_atoms = self.amap.select('name P').numAtoms()
                else:
                    num_p_atoms = 0

                if self.amap.ca is not None:
                    num_ca_atoms = self.amap.ca.numAtoms()
                else:
                    num_ca_atoms = 0

                num_rep_atoms = num_ca_atoms + num_p_atoms
                if num_rep_atoms == self.amap.numAtoms():
                    fhCmd.write("mol modstyle 0 0 Beads 2.000000 8.000000\n")
                    # fhCmd.write("mol modstyle 0 0 Beads 1.800000 6.000000 "
                    #         "2.600000 0\n")
                else:
                    fhCmd.write("mol modstyle 0 0 NewRibbons 1.800000 6.000000 "
                            "2.600000 0\n")
            fhCmd.write("animate speed 0.5\n")
            fhCmd.write("animate forward\n")
            fhCmd.close()    

    def qualifyModesStep(self, numberOfModes, collectivityThreshold, structureEM, suffix=''):
        self._enterWorkingDir()

        fnVec = glob("modes/vec.*")

        if len(fnVec) < numberOfModes:
            msg = "There are only %d modes instead of %d. "
            msg += "Check the number of modes you asked to compute and/or consider increasing cut-off distance."
            msg += "The maximum number of modes allowed by the method for atomic normal mode analysis is "
            msg += "3 times the number of nodes (pseudoatoms or Calphas). "
            self._printWarnings(redStr(msg % (len(fnVec), numberOfModes)))

        mdOut = MetaData()
        collectivityList = list(prody.calcCollectivity(self.rtb))
        eigvals = self.rtb.getEigvals()

        for n in range(len(fnVec)):
            collectivity = collectivityList[n]

            objId = mdOut.addObject()
            modefile = self._getPath("modes", "vec.%d" % (n + 1))
            mdOut.setValue(MDL_NMA_MODEFILE, modefile, objId)
            mdOut.setValue(MDL_ORDER, int(n + 1), objId)

            if n >= self.startMode:
                mdOut.setValue(MDL_ENABLED, 1, objId)
            else:
                mdOut.setValue(MDL_ENABLED, -1, objId)

            mdOut.setValue(MDL_NMA_COLLECTIVITY, collectivity, objId)
            mdOut.setValue(MDL_NMA_EIGENVAL, eigvals[n] , objId)

            if collectivity < collectivityThreshold:
                mdOut.setValue(MDL_ENABLED, -1, objId)

        idxSorted = [i[0] for i in sorted(enumerate(collectivityList), key=lambda x: x[1], reverse=True)]

        score = []
        for j in range(len(fnVec)):
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
        
        prody.writeScipionModes(self._getPath(), self.rtb, scores=score, only_sqlite=True,
                                collectivityThreshold=collectivityThreshold)

    def computeAtomShiftsStep(self, numberOfModes):
        fnOutDir = self._getExtraPath("distanceProfiles")
        makePath(fnOutDir)
        maxShift=[]
        maxShiftMode=[]
        
        for n in range(self.startMode+1, numberOfModes+1):
            fnVec = self._getPath("modes", "vec.%d" % n)
            if exists(fnVec):
                fhIn = open(fnVec)
                md = MetaData()
                atomCounter = 0
                for line in fhIn:
                    x, y, z = map(float, line.split())
                    d = math.sqrt(x*x+y*y+z*z)
                    if n==self.startMode+1:
                        maxShift.append(d)
                        maxShiftMode.append(self.startMode+1)
                    else:
                        if d>maxShift[atomCounter]:
                            maxShift[atomCounter]=d
                            maxShiftMode[atomCounter]=n
                    atomCounter+=1
                    md.setValue(MDL_NMA_ATOMSHIFT,d,md.addObject())
                md.write(join(fnOutDir,"vec%d.xmd" % n))
                fhIn.close()
        md = MetaData()
        for i, _ in enumerate(maxShift):
            fnVec = self._getPath("modes", "vec.%d" % (maxShiftMode[i]+1))
            if exists(fnVec):
                objId = md.addObject()
                md.setValue(MDL_NMA_ATOMSHIFT, maxShift[i],objId)
                md.setValue(MDL_NMA_MODEFILE, fnVec, objId)
        md.write(self._getExtraPath('maxAtomShifts.xmd'))
        
        # configure ProDy to restore secondary structure information and verbosity
        prody.confProDy(auto_secondary=self.oldSecondary, 
                        verbosity='{0}'.format(self.oldVerbosity))

    def createOutputStep(self):
        fnSqlite = self._getPath('modes.sqlite')
        nmSet = SetOfNormalModes(filename=fnSqlite)
        nmSet._nmdFileName = String(self._getPath('modes.nmd'))

        inputPdb = self.inputStructure.get()
        nmSet.setPdb(inputPdb)

        self._defineOutputs(outputModes=nmSet)
        self._defineSourceRelation(self.inputStructure, nmSet)

