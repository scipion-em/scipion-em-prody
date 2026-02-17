# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     James Krieger (jamesmkrieger@gmail.com)
# *              Ricardo Serrano Gutiérrez (rserranogut@hotmail.com)                 
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
This module will provide ProDy normal mode analysis (NMA) using the Gaussian network model (GNM).
"""
from pwem.objects import String, EMFile, SetOfNormalModes

from pyworkflow.protocol.params import (PointerParam, IntParam, FloatParam, StringParam,
                                        BooleanParam, LEVEL_ADVANCED)

from prody2.protocols.protocol_modes_base import ProDyModesBase
from prody2 import Plugin, copyConvertPDB

class ProDyGNM(ProDyModesBase):
    """
    This protocol will perform normal mode analysis (NMA) using the Gaussian network model (GNM)
    """
    _label = 'GNM analysis'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label='ProDy GNM NMA')

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

        form.addParam('cutoff', FloatParam, default=10,
                      expertLevel=LEVEL_ADVANCED,
                      label="Cut-off distance (A)",
                      help='Atoms or pseudoatoms beyond this distance will not interact. \n'
                           'For Calpha atoms, the default distance of 7.5 A works well in the majority of cases. \n'
                           'For all atoms, a shorter distance is recommended.'
                           'For fewer atoms or pseudoatoms, a longer distance is recommended.')

        form.addParam('gamma', StringParam, default=1.,
                      expertLevel=LEVEL_ADVANCED,
                      label="Spring constant",
                      help='This number or function determines the strength of the springs.\n'
                           'More sophisticated options are available within the ProDy API and '
                           'the resulting modes can be imported back into Scipion.\n'
                           'See http://http://www.bahargroup.org/prody/tutorials/enm_analysis/gamma.html')

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
                      'Modes 1-6 are always deselected as they are related to rigid-body movements. \n'
                      'The modes metadata file can be used to see which modes are more collective '
                      'in order to decide which modes to use at the image analysis step.')

        form.addParam('zeros', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Include zero eigvals",
                      help='Elect whether modes with zero eigenvalues will be kept.')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps

        # Link the input
        inputFn = self.inputStructure.get().getFileName()
        self.structureEM = self.inputStructure.get().getPseudoAtoms()
        n = self.numberOfModes.get()
        nzeros = self.getNzero()

        self._insertFunctionStep(self.computeModesStep, inputFn, n)
        self._insertFunctionStep(self.qualifyModesStep, n,
                                 self.collectivityThreshold.get(),
                                 self.structureEM)
        self._insertFunctionStep(self.computeAtomShiftsStep, n, nzeros)
        self._insertFunctionStep(self.createOutputStep)

    def computeModesStep(self, inputFn, n):

        if self.structureEM:
            self.pdbFileName = self._getPath('pseudoatoms.pdb')
        else:
            self.pdbFileName = self._getPath('atoms.pdb')

        copyConvertPDB(inputFn, self.pdbFileName)

        args = '{0} -s "all" --altloc "all" --kirchhoff --export-scipion --npz --npzmatrices ' \
               '-o {1} -p {2} -n {3} -g {4} -c {5} -P {6}' \
                ' --covariance --cross-correlations'.format(
                    self.pdbFileName, self._getPath(), self.getPrefix(), n,
                    self.gamma.get(), self.cutoff.get(), self.numberOfThreads.get())

        if self.zeros.get():
            args += ' --zero-modes'
            self.startMode = 1
        else:
            self.startMode = 0
        
        if self.membrane.get():
            args += ' --membrane'

        self.runJob(Plugin.getProgram('gnm'), args)

    # def computeAtomShiftsStep(self, numberOfModes):
    #     fnOutDir = self._getExtraPath("distanceProfiles")
    #     makePath(fnOutDir)
    #     maxShift=[]
    #     maxShiftMode=[]
    #     vecStr = "vec.%d"
    #     for n in range(self.startMode+1, numberOfModes+1):
    #         fnVec = self._getPath("modes", vecStr % n)
    #         if exists(fnVec):
    #             fhIn = open(fnVec)
    #             md = MetaData()
    #             atomCounter = 0
    #             for line in fhIn:
    #                 d = abs(float(line))
    #                 if n==self.startMode+1:
    #                     maxShift.append(d)
    #                     maxShiftMode.append(self.startMode+1)
    #                 else:
    #                     if d>maxShift[atomCounter]:
    #                         maxShift[atomCounter]=d
    #                         maxShiftMode[atomCounter]=n
    #                 atomCounter+=1
    #                 md.setValue(MDL_NMA_ATOMSHIFT,d,md.addObject())
    #             md.write(join(fnOutDir,"vec%d.xmd" % n))
    #             fhIn.close()
                
    #     md = MetaData()
    #     for i, _ in enumerate(maxShift):
    #         fnVec = self._getPath("modes", vecStr % (maxShiftMode[i]+1))
    #         if exists(fnVec):
    #             objId = md.addObject()
    #             md.setValue(MDL_NMA_ATOMSHIFT, maxShift[i],objId)
    #             md.setValue(MDL_NMA_MODEFILE, fnVec, objId)
    #     md.write(self._getExtraPath('maxAtomShifts.xmd'))

    def createOutputStep(self):
        outputMatrixCov = EMFile(filename=self._getExtraPath('modes_covariance.txt'))
        outputMatrixCrosCor = EMFile(filename=self._getExtraPath('modes_cross-correlations.txt'))

        fnSqlite = self._getPath('modes.sqlite')
        nmSet = SetOfNormalModes(filename=fnSqlite)
        nmSet._nmdFileName = String(self._getPath(self.getPrefix() + '.nmd'))

        inputPdb = self.inputStructure.get()
        nmSet.setPdb(inputPdb)

        self._defineOutputs(outputModes=nmSet,
                            matrixFileCC=outputMatrixCrosCor,
                            matrixFileCV=outputMatrixCov)
        self._defineSourceRelation(self.inputStructure, nmSet)

    def getPrefix(self):
        if self.membrane.get():
            return 'modes.exgnm'
        else:
            return 'modes.gnm'

    def getNzero(self):
        if self.zeros.get():
            return 1
        else:
            return 0
