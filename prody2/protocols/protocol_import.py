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
This module will provide ProDy mode import tools.
"""
import os
import numpy as np

from pwem.objects import (String, AtomStruct, SetOfAtomStructs,
                          SetOfNormalModes, SetOfPrincipalComponents)
from pwem.protocols import ProtImportFiles

from prody2.objects import ProDyNpzEnsemble, TrajFrame

import pyworkflow.protocol.params as params

import prody
from prody.dynamics.gnm import ZERO

NMD = 0
modes_NPZ = 1
SCIPION = 2
GROMACS = 3

filesPatternHelp = """Pattern of the files to be imported.\n\n
The pattern can contain standard wildcards such as\n
*, ?, etc, or special ones like ### to mark some\n
digits in the filename as ID.\n\n
NOTE: wildcards and special characters 
('*', '?', '#', ':', '%') cannot appear in the actual path.\n\n
For gromacs modes, the first is for values and this is for vectors"""

class ProDyImportModes(ProtImportFiles):
    """
    This protocol will import ProDy modes to a SetOfNormalModes object
    """
    _label = 'Import modes'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        importChoices = self._getImportChoices()
        filesCondition = self._getFilesCondition()

        form.addSection(label='Import modes')

        form.addParam('importFrom', params.EnumParam,
                      choices=importChoices, default=self._getDefaultChoice(),
                      label='Import from',
                      help='Select the type of import.')

        form.addParam('importType', params.EnumParam,
                      choices=['prody nmd', 'prody npz', 'scipion', 'gromacs'],
                      default=NMD,
                      label='Type of modes to import',
                      help='ProDy can support import of modes in various file formats: \n'
                           'native nmd and npz files, a scipion-style directory for a SetOfNormalModes, '
                           'or a gromacs directory with files containing eigenvectors and eigenvalues')

        form.addParam('filesPath', params.PathParam,
                      condition=filesCondition,
                      label="Files directory",
                      help="Directory with the files you want to import.\n\n"
                           "The path can also contain wildcards to select"
                           "from several folders. \n\n"
                           "Examples:\n"
                           "  ~/project/data/day??_files/\n"
                           "Each '?' represents one unknown character\n\n"
                           "  ~/project/data/day*_files/\n"
                           "'*' represents any number of unknown characters\n\n"
                           "  ~/project/data/day##_files/\n"
                           "'##' represents two digits that will be used as "
                           "file ID\n\n"
                           "NOTE: wildcard characters ('*', '?', '#') "
                           "cannot appear in the actual path.)")

        form.addParam('filesPattern', params.StringParam,
                      label='File pattern',
                      condition="importType!=%d" % SCIPION,
                      help=filesPatternHelp)

        form.addParam('filesPattern2', params.StringParam,
                      label='File pattern for eigenvalues',
                      condition="importType==%d" % GROMACS,
                      help=filesPatternHelp)

        form.addParam('inputStructure', params.PointerParam, label="Input structure",
                      pointerClass='AtomStruct', condition="importType!=%d" % NMD,
                      help='The input structure can be an atomic model '
                           '(true PDB) or a pseudoatomic model\n'
                           '(an EM volume converted into pseudoatoms)')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('importModesStep')
        self._insertFunctionStep('createOutputStep')

    def importModesStep(self):
        # configure ProDy to automatically handle secondary structure information and verbosity
        from pyworkflow import Config
        global oldSecondary; oldSecondary = prody.confProDy("auto_secondary")
        global oldVerbosity; oldVerbosity = prody.confProDy("verbosity")
        prodyVerbosity =  'none' if not Config.debugOn() else 'debug'
        prody.confProDy(auto_secondary=True, verbosity='{0}'.format(prodyVerbosity))

        filesPaths = self.getMatchFiles()

        if self.importType == SCIPION:
            filePath = filesPaths[0]
        else:
            folderPath = os.path.split(filesPaths[0])[0]
            self.pattern1 = os.path.split(filesPaths[0])[1]

        pdbFilename = self.inputStructure.get().getFileName()
        
        if self.importType == NMD:
            if not self.pattern1.endswith('.nmd'):
                self.pattern1 += '.nmd'
            self.outModes, _ = prody.parseNMD(os.path.join(folderPath, self.pattern1))

        elif self.importType == modes_NPZ:
            if not self.pattern1.endswith('.npz'):
                self.pattern1 += '.npz'
            self.outModes = prody.loadModel(os.path.join(folderPath, self.pattern1))

        elif self.importType == SCIPION:
            self.outModes = prody.parseScipionModes(filePath, pdb=pdbFilename)

        else:
            pattern2 = self.filesPattern2.get()
            self.outModes = prody.parseGromacsModes(folderPath, eigval_fname=self.pattern1,
                                                    eigvec_fname=pattern2, average_pdb=pdbFilename)

        prody.writeScipionModes(self._getPath(), self.outModes, write_star=True)
        
        if self.importType != NMD:
            atoms = prody.parsePDB(pdbFilename)
            self.nmdFileName = self._getPath('modes.nmd')
            prody.writeNMD(self.nmdFileName, self.outModes, atoms)
        else:
            self.nmdFileName = self.pattern1
            
        # configure ProDy to restore secondary structure information and verbosity
        prody.confProDy(auto_secondary=oldSecondary, verbosity='{0}'.format(oldVerbosity))

    def createOutputStep(self):
        fnSqlite = self._getPath('modes.sqlite')

        if self.outModes.getEigvals()[0] <= self.outModes.getEigvals()[1] or self.outModes.getEigvals()[0] < ZERO:
            nmSet = SetOfNormalModes(filename=fnSqlite)
        else:
            nmSet = SetOfPrincipalComponents(filename=fnSqlite)

        nmSet._nmdFileName = String(self.nmdFileName)

        inputPdb = self.inputStructure.get()
        nmSet.setPdb(inputPdb)

        self._defineOutputs(outputModes=nmSet)
        self._defineSourceRelation(self.inputStructure, nmSet)

PDB = 0
DCD = 1
ens_NPZ = 2

NO_SUP = 0
YES_SUP = 1
ITERPOSE = 2

class ProDyImportEnsemble(ProtImportFiles):
    """
    This protocol will import and optionally trim a ProDy Ensemble
    """
    _label = 'Import/Trim Ensemble'
    _possibleOutputs = {'outputAtomStructs': SetOfAtomStructs,
                        'outputNpz': ProDyNpzEnsemble}

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """

        form.addSection(label='Import ensemble')

        importChoices = self._getImportChoices()
        filesCondition = self._getFilesCondition()

        form.addParam('importFrom', params.EnumParam,
                      choices=importChoices, default=self._getDefaultChoice(),
                      label='Import from',
                      help='Select the type of import.')

        form.addParam('importType', params.EnumParam, choices=['pdb', 'dcd', 'ens.npz'],
                      default=PDB, condition=filesCondition,
                      label='Type of ensemble file to import',
                      help='ProDy can support import of ensembles in various file formats: \n'
                           'pdb, dcd and the native ens.npz format')

        form.addParam('importPointer', params.PointerParam,
                      pointerClass='AtomStruct,SetOfAtomStructs,ProDyNpzEnsemble',
                      condition='importFrom!=0',
                      label='Ensemble file to import',
                      help='ProDy can support import of ensembles in various file formats: \n'
                           'pdb, dcd and the native ens.npz format')

        form.addParam('filesPath', params.PathParam,
                      label="Files directory", condition=filesCondition,
                      help="Directory with the files you want to import.\n\n"
                           "The path can also contain wildcards to select"
                           "from several folders. \n\n"
                           "Examples:\n"
                           "  ~/project/data/day??_files/\n"
                           "Each '?' represents one unknown character\n\n"
                           "  ~/project/data/day*_files/\n"
                           "'*' represents any number of unknown characters\n\n"
                           "  ~/project/data/day##_files/\n"
                           "'##' represents two digits that will be used as "
                           "file ID\n\n"
                           "NOTE: wildcard characters ('*', '?', '#') "
                           "cannot appear in the actual path.)")

        form.addParam('filesPattern', params.StringParam,
                      label='Pattern', condition=filesCondition,
                      help=filesPatternHelp)

        form.addParam('inputStructure', params.PointerParam, label="Input structure",
                      pointerClass='AtomStruct', condition="importType==%d or importFrom!=0" % DCD,
                      help='The input structure can be an atomic model '
                           '(true PDB) or a pseudoatomic model '
                           '(an EM volume converted into pseudoatoms) '
                           'The input structure should have the same number of atoms '
                           'as the original ensemble.')

        form.addParam('superpose', params.EnumParam, label="Superpose?",
                      choices=['No', 'Once', 'Iteratively'], default=NO_SUP,
                      help='Elect whether and how to superpose the ensemble')

        form.addParam('savePDBs', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Save PDB files?",
                      help='This takes up storage and time, but '
                           'may be helpful for interfacing with ContinuousFlex.')

        form.addParam('selstr', params.StringParam, default="all",
                      label="Selection string",
                      help='Selection string for atoms to include in the ensemble.\n'
                           'It is recommended to use "all" (default), "protein" or "name CA"')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('importEnsembleStep')
        self._insertFunctionStep('createOutputStep')

    def importEnsembleStep(self):
        # configure ProDy to automatically handle secondary structure information and verbosity
        from pyworkflow import Config
        global oldSecondary; oldSecondary = prody.confProDy("auto_secondary")
        global oldVerbosity; oldVerbosity = prody.confProDy("verbosity")
        prodyVerbosity =  'none' if not Config.debugOn() else 'debug'
        prody.confProDy(auto_secondary=True, verbosity='{0}'.format(prodyVerbosity))
        
        if self.importFrom.get() == 0:
            filesPaths = self.getMatchFiles()
            folderPath = os.path.split(filesPaths[0])[0]
            self.pattern1 = os.path.split(filesPaths[0])[1]
            
            if self.importType == PDB:
                if not (self.pattern1.endswith('.pdb') or self.pattern1.endswith('.cif')):
                    self.pattern1 += '.pdb'
                self.atoms = prody.parsePDB(os.path.join(folderPath, self.pattern1))
                self.outEns = prody.PDBEnsemble(self.atoms)

            elif self.importType == DCD:
                if not self.pattern1.endswith('.dcd'):
                    self.pattern1 += '.dcd'
                self.outEns = prody.PDBEnsemble(prody.parseDCD(os.path.join(folderPath, self.pattern1)))
                self.atoms = prody.parsePDB(self.inputStructure.get().getFileName())

            elif self.importType == ens_NPZ:
                if not self.pattern1.endswith('.ens.npz'):
                    self.pattern1 += '.ens.npz'
                self.outEns = prody.PDBEnsemble(prody.loadEnsemble(os.path.join(folderPath, self.pattern1)))
                self.atoms = self.outEns.getAtoms()
        else:
            point = self.importPointer.get()
            if isinstance(point, ProDyNpzEnsemble):
                self.outEns = point.loadEnsemble()
            elif isinstance(point, SetOfAtomStructs):
                ags = prody.parsePDB([struct.getFileName() for struct in point])
                self.outEns = prody.PDBEnsemble()
                self.outEns.setCoords(ags[0])
                for ag in ags:
                    self.outEns.addCoordset(ag)
            else:
                self.outEns = prody.PDBEnsemble(prody.parsePDB(point.getFileName()))

            self.atoms = prody.parsePDB(self.inputStructure.get().getFileName())

        self.outEns.setAtoms(self.atoms)

        selstr = self.selstr.get()
        self.outEns.setAtoms(self.atoms.select(selstr))
        self.outEns = prody.trimPDBEnsemble(self.outEns) # hard

        if self.superpose == YES_SUP:
            self.outEns.superpose()
        elif self.superpose == ITERPOSE:
            self.outEns.iterpose()
            
        if self.savePDBs.get():
            self.pdbs = SetOfAtomStructs().create(self._getExtraPath())
            for i, coordset in enumerate(self.outEns.getCoordsets()):
                atoms = self.atoms.select(selstr).copy()
                atoms.setCoords(coordset)
                filename = self._getExtraPath('{:s}_{:06d}.pdb'.format(atoms.getTitle(), i))
                prody.writePDB(filename, atoms)
                pdb = AtomStruct(filename)
                self.pdbs.append(pdb)

        self.filename = prody.saveEnsemble(self.outEns, self._getExtraPath('ensemble.ens.npz'))

        self.npz = ProDyNpzEnsemble().create(self._getExtraPath())
        for j in range(self.outEns.numConfs()):
            self.npz.append(TrajFrame((j+1, self.filename)))

        # configure ProDy to restore secondary structure information and verbosity
        prody.confProDy(auto_secondary=oldSecondary, verbosity='{0}'.format(oldVerbosity))

    def createOutputStep(self):
        if self.savePDBs.get():
            self._defineOutputs(outputNpz=self.npz, outputStructures=self.pdbs)
        else:
            self._defineOutputs(outputNpz=self.npz)

    def _summary(self):
        if not hasattr(self, 'outputNpz'):
            summ = ['Output ensemble not ready yet']
        elif not hasattr(self, 'sumOutput'):
            ens = self.outputNpz.loadEnsemble()
            summ = ['Ensemble imported with *{0}* structures of *{1}* atoms'.format(
                   ens.numConfs(), ens.numAtoms())]
            self.sumOutput = True
        return summ

    def _getImportChoices(self):
        return ['files', 'pointer']

