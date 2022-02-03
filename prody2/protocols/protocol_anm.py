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
from pyworkflow.protocol import params

from os.path import basename, exists, join
import math

from pwem import *
from pwem.emlib import (MetaData, MDL_NMA_MODEFILE, MDL_ORDER,
                        MDL_ENABLED, MDL_NMA_COLLECTIVITY, MDL_NMA_SCORE, 
                        MDL_NMA_ATOMSHIFT)
from pwem.objects import AtomStruct, SetOfNormalModes, String
from pwem.protocols import EMProtocol

from pyworkflow.utils import *
from pyworkflow.utils.path import makePath
from pyworkflow.protocol.params import (PointerParam, IntParam, FloatParam, StringParam,
                                        BooleanParam, LEVEL_ADVANCED)

import prody

class ProDyANM(EMProtocol):
    """
    This protocol will perform normal mode analysis (NMA) using the anisotropic network model (ANM)
    """
    _label = 'ANM analysis'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label='ProDy ANM NMA')

        form.addParam('inputStructure', PointerParam, label="Input structure",
                      important=True,
                      pointerClass='AtomStruct',
                      help='The input structure can be an atomic model '
                           '(true PDB) or a pseudoatomic model\n'
                           '(an EM volume converted into pseudoatoms)')

        form.addParam('selection', StringParam, default="",
                      label="selection string",
                      help='This determines which atoms are used in the calculations. '
                           'If left empty, the default behaviour is to use '
                           '"protein and name CA or nucleic and name P C4\' C2" '
                           'for atomic structures and "all" for pseudoatoms.')

        form.addParam('numberOfModes', IntParam, default=20,
                      label='Number of modes',
                      help='The maximum number of modes allowed by the method for '
                           'atomic normal mode analysis is 3 times the '
                           'number of nodes (Calpha atoms or pseudoatoms).')

        form.addParam('cutoff', FloatParam, default=15.,
                      expertLevel=LEVEL_ADVANCED,
                      label="Cut-off distance (A)",
                      help='Atoms or pseudoatoms beyond this distance will not interact. \n'
                           'For Calpha atoms, the default distance of 15 A works well in the majority of cases. \n'
                           'For pseudoatoms, set this according to the level of coarse-graining '
                           '(see Doruker et al., J Comput Chem 2002). \n'
                           'For all atoms, Tirion (Phys. Rev. Lett. 1996) found ~5 A works well '
                           '(2 vdW radii of ~1.5 A + ~2 A)')

        form.addParam('gamma', StringParam, default=1.,
                      expertLevel=LEVEL_ADVANCED,
                      label="Spring constant",
                      help='This determines the strength of the strings and can be a function. \n'
                           'This allows the use of ENMs based on pseudoatom, atom or residue properties including radii, '
                           'secondary structure, distance-dependent potentials, and other custom functions. \n'
                           'See http://prody.csb.pitt.edu/tutorials/enm_analysis/gamma.html')

        form.addParam('collectivityThreshold', FloatParam, default=0.15,
                      expertLevel=LEVEL_ADVANCED,
                      label='Threshold on collectivity',
                      help='Collectivity degree is related to the number of atoms or pseudoatoms that are affected by '
                      'the mode, and it is normalized between 0 and 1. Modes below this threshold are deselected in '
                      'the modes metadata file as these modes are much less collective. \n'
                      'For no deselection, this parameter should be set to 0 . \n'
                      'Modes 1-6 are always deselected as they are related to rigid-body movements. \n'
                      'The modes metadata file can be used to see which modes are more collective '
                      'in order to decide which modes to use at the image analysis step.')

        form.addParam('sparse', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Use sparse matrices",
                      help='Elect whether to use Scipy sparse matrices. \n'
                           'This allows efficient usage of memory at the cost of computation speed')

        form.addParam('kdtree', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Use KDTree for building Hessian matrix",
                      help='Elect whether to use KDTree for building Hessian matrix. This is usually slower.')

        form.addParam('zeros', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="include zero eigvals",
                      help='Elect whether modes with zero eigenvalues will be kept.')

        form.addParam('turbo', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Use turbo mode",
                      help='Elect whether to use a memory intensive, but faster way to calculate modes.')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps

        # Link the input
        inputFn = self.inputStructure.get().getFileName()
        self.structureEM = self.inputStructure.get().getPseudoAtoms()

        self.model_type = 'anm'
        n = self.numberOfModes.get()

        self._insertFunctionStep('computeModesStep', inputFn, n)
        self._insertFunctionStep('qualifyModesStep', n,
                                 self.collectivityThreshold.get(),
                                 self.structureEM)
        self._insertFunctionStep('computeAtomShiftsStep', n)
        self._insertFunctionStep('createOutputStep')

    def computeModesStep(self, inputFn, n):
        
        if self.structureEM:
            if self.selection == "":
                self.selection = "all"
            self.pdbFileName = self._getPath('pseudoatoms.pdb')
        else:
            if self.selection == "":
                self.selection = "protein and name CA or nucleic and name P C4' C2"
            self.pdbFileName = self._getPath('atoms.pdb')

        ag = prody.parsePDB(inputFn, alt='all')
        selection = ag.select(str(self.selection))
        prody.writePDB(self.pdbFileName, selection)

        self.runJob('prody', 'anm {0} -s "all" --altloc "all" --zero-modes '
                    '--export-scipion --npz -o {1} -p modes -n {2}'.format(self.pdbFileName,
                                                                           self._getPath(), n))
        self.anm = prody.loadModel(self._getPath('modes.anm.npz'))

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
        collectivityList = list(prody.calcCollectivity(self.anm))

        for n in range(len(fnVec)):
            collectivity = collectivityList[n]

            objId = mdOut.addObject()
            modefile = self._getPath("modes", "vec.%d" % (n + 1))
            mdOut.setValue(MDL_NMA_MODEFILE, modefile, objId)
            mdOut.setValue(MDL_ORDER, int(n + 1), objId)

            if n >= 6:
                mdOut.setValue(MDL_ENABLED, 1, objId)
            else:
                mdOut.setValue(MDL_ENABLED, -1, objId)
            mdOut.setValue(MDL_NMA_COLLECTIVITY, collectivity, objId)

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
        
        prody.writeScipionModes(self._getPath(), self.anm, scores=score, only_sqlite=True,
                              collectivityThreshold=collectivityThreshold)

    def computeAtomShiftsStep(self, numberOfModes):
        fnOutDir = self._getExtraPath("distanceProfiles")
        makePath(fnOutDir)
        maxShift=[]
        maxShiftMode=[]
        
        for n in range(7, numberOfModes+1):
            fnVec = self._getPath("modes", "vec.%d" % n)
            if exists(fnVec):
                fhIn = open(fnVec)
                md = MetaData()
                atomCounter = 0
                for line in fhIn:
                    x, y, z = map(float, line.split())
                    d = math.sqrt(x*x+y*y+z*z)
                    if n==7:
                        maxShift.append(d)
                        maxShiftMode.append(7)
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

    def createOutputStep(self):
        fnSqlite = self._getPath('modes.sqlite')
        nmSet = SetOfNormalModes(filename=fnSqlite)
        nmSet._nmdFileName = String(self._getPath('modes.nmd'))

        outputPdb = AtomStruct()
        outputPdb.setFileName(self.pdbFileName)
        nmSet.setPdb(outputPdb.get())

        self._defineOutputs(outputModes=nmSet, outputStructure=outputPdb)
        self._defineSourceRelation(outputPdb, nmSet)

