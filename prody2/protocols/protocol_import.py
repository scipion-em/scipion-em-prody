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

from pwem.objects import (String, AtomStruct, SetOfAtomStructs, EMFile,
                          SetOfNormalModes, SetOfPrincipalComponents)
from pwem.protocols import ProtImportFiles

from prody2.objects import (ProDyNpzEnsemble, TrajFrame,
                            SetOfGnmModes, SetOfLdaModes)
from prody2.constants import ENSEMBLE_WEIGHTS
from prody2 import fixVerbositySecondary, restoreVerbositySecondary

import pyworkflow.object as pwobj
import pyworkflow.protocol.params as params
from pyworkflow.utils import logger

import prody
from prody.dynamics.gnm import ZERO

NMD = 0
MODES_NPZ = 1
SCIPION = 2
GROMACS = 3

PDB_FILENAME = 'atoms.pdb'
PSF_FILENAME = 'atoms.psf'
ENS_FILENAME = 'ensemble.dcd'

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
        fixVerbositySecondary(self)

        filesPaths = self.getMatchFiles()

        if self.importType == SCIPION:
            filePath = filesPaths[0]
        else:
            folderPath = os.path.split(filesPaths[0])[0]
            self.pattern1 = os.path.split(filesPaths[0])[1]

        if self.inputStructure.get() is not None:
            pdbFilename = self.inputStructure.get().getFileName()
        
        if self.importType == NMD:
            if not self.pattern1.endswith('.nmd'):
                self.pattern1 += '.nmd'

            if self.pattern1.find('pca') != -1:
                prodyType = prody.PCA
            elif self.pattern1.find('lda') != -1:
                prodyType = prody.LDA
            elif self.pattern1.find('gnm') != -1:
                prodyType = prody.GNM
            else:
                prodyType = prody.NMA

            self.outModes, self.atoms = prody.parseNMD(os.path.join(folderPath, self.pattern1),
                                                       type=prodyType)

            if self.inputStructure.get() is None:
                pdbFilename = prody.writePDB(self._getExtraPath('atoms'), self.atoms)
                self.inputStructure = AtomStruct(filename=pdbFilename)

        elif self.importType == MODES_NPZ:
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
            self.atoms = prody.parsePDB(pdbFilename)
            if isinstance(self.outModes, prody.PCA):
                self.nmdFileName = self._getPath('modes.pca.nmd')
            elif isinstance(self.outModes, prody.GNM):
                self.nmdFileName = self._getPath('modes.gnm.nmd')
            else:
                self.nmdFileName = self._getPath('modes.nmd')
            prody.writeNMD(self.nmdFileName, self.outModes, self.atoms)
        else:
            self.nmdFileName = self.pattern1
            
        restoreVerbositySecondary(self)

    def createOutputStep(self):
        fnSqlite = self._getPath('modes.sqlite')

        if isinstance(self.outModes, prody.GNM) or self.outModes.numAtoms() < self.atoms.numAtoms():
            nmSet = SetOfGnmModes(filename=fnSqlite)
        elif (self.outModes.getEigvals()[0] <= self.outModes.getEigvals()[1]
            or self.outModes.getEigvals()[0] < ZERO):
            nmSet = SetOfNormalModes(filename=fnSqlite)
        elif isinstance(self.outModes, prody.LDA):
            nmSet = SetOfLdaModes(filename=fnSqlite)
        else:
            nmSet = SetOfPrincipalComponents(filename=fnSqlite)

        nmSet._nmdFileName = String(self.nmdFileName)

        inputPdb = self.inputStructure.get()
        nmSet.setPdb(inputPdb)

        self._defineOutputs(outputModes=nmSet)
        self._defineSourceRelation(self.inputStructure, nmSet)

PDB = 0
DCD = 1
ENS_NPZ = 2

NO_SUP = 0
YES_SUP = 1
ITERPOSE = 2

POINTER_TYPES = 'AtomStruct,SetOfAtomStructs,ProDyNpzEnsemble'

try:
    from pwchem.objects import MDSystem
    imported_chem = True
    POINTER_TYPES += ',MDSystem'
except ImportError:
    imported_chem = False

class ProDyImportEnsemble(ProtImportFiles):
    """
    This protocol will import and optionally trim a ProDy Ensemble
    """
    _label = 'Import Ensemble'
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
                      pointerClass=POINTER_TYPES,
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

        form.addParam('inputStructure', params.PointerParam, label="Input structure", allowsNull=True,
                      pointerClass='AtomStruct', condition="importType==%d or importFrom!=0" % DCD,
                      help='The input structure can be an atomic model '
                           '(true PDB) or a pseudoatomic model '
                           '(an EM volume converted into pseudoatoms). '
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

        form.addParam('writeDCDFile', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      condition='selstr!=all or importType!=%d' % DCD,
                      label="Whether to write DCD trajectory file",
                      help='This will be registered as output too')

        form.addParam('inputPsf', params.PathParam, label="Input PSF topology", allowsNull=True,
                      condition="writeDCDFile==True",
                      help='An input psf topology can also be provided. '
                           'The psf should have the same number of atoms '
                           'as the original ensemble.')

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
        fixVerbositySecondary(self)
        
        self.weights = None
        self.atoms = None

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

                if self.inputStructure.get() is not None:
                    self.atoms = prody.parsePDB(self.inputStructure.get().getFileName())

            elif self.importType == ENS_NPZ:
                if not self.pattern1.endswith('.ens.npz'):
                    self.pattern1 += '.ens.npz'
                self.outEns = prody.loadEnsemble(os.path.join(folderPath, self.pattern1))
                if not isinstance(self.outEns, prody.PDBEnsemble):
                    self.outEns = prody.PDBEnsemble(self.outEns)
                if 'size' in self.outEns.getDataLabels():
                    self.weights = self.outEns.getData('size')

                self.atoms = self.outEns.getAtoms()
        else:
            point = self.importPointer.get()

            if isinstance(point, ProDyNpzEnsemble):
                self.weights = [item.getAttributeValue(ENSEMBLE_WEIGHTS) for item in point]
                self.outEns = point.loadEnsemble()

                if self.inputStructure.get() is not None:
                    self.atoms = prody.parsePDB(self.inputStructure.get().getFileName())
                else:
                    self.atoms = self.outEns.getAtoms()

            elif isinstance(point, SetOfAtomStructs):
                self.weights = [item.getAttributeValue(ENSEMBLE_WEIGHTS) for item in point]
                self.ags = prody.parsePDB([struct.getFileName() for struct in point])
                self.outEns = prody.PDBEnsemble()
                self.outEns.setCoords(self.ags[0])
                for ag in self.ags:
                    self.outEns.addCoordset(ag)

                if self.inputStructure.get() is not None:
                    self.atoms = prody.parsePDB(self.inputStructure.get().getFileName())

            elif isinstance(point, AtomStruct):
                self.outEns = prody.PDBEnsemble(prody.parsePDB(point.getFileName()))
                if self.inputStructure.get() is not None:
                    self.atoms = prody.parsePDB(self.inputStructure.get().getFileName())

            else:
                # only activated if we can use MDSystem
                self.outEns = prody.PDBEnsemble(prody.parseDCD(point.getTrajectoryFile()))
                self.atoms = prody.parsePDB(point.getSystemFile())

        if self.weights is None or np.array_equal(self.weights, np.zeros(self.weights.shape)):
            self.weights = np.ones(self.outEns.numConfs())

        selstr = self.selstr.get()

        if self.atoms is not None:
            try:
                # setAtoms then select and trim
                self.outEns.setAtoms(self.atoms)
                self.outEns.setAtoms(self.atoms.select(selstr))
                self.outEns = prody.trimPDBEnsemble(self.outEns) # hard
            except ValueError:
                try:
                    # setAtoms with select directly then trim
                    self.outEns.setAtoms(self.atoms.select(selstr))
                    self.outEns = prody.trimPDBEnsemble(self.outEns) # hard
                except ValueError:
                    try:
                        # select then trim then setAtoms
                        self.outEns.select(selstr)
                        self.outEns = prody.trimPDBEnsemble(self.outEns) # hard
                        self.outEns.setAtoms(self.atoms.select(selstr))
                    except ValueError:
                        if hasattr(self, 'ags'):
                            try:
                                # setAtoms with first structure if available and trim
                                # then setAtoms with ref structure and trim again
                                self.outEns.setAtoms(self.ags[0])
                                self.outEns.setAtoms(self.ags[0].select(selstr))
                                self.outEns = prody.trimPDBEnsemble(self.outEns) # hard
                                self.outEns.setAtoms(self.atoms)
                                self.outEns.setAtoms(self.atoms.select(selstr))
                                self.outEns = prody.trimPDBEnsemble(self.outEns) # hard
                            except ValueError:
                                raise ValueError("Reference structure should have same number of atoms as ensemble")
                        else:
                            raise ValueError("Reference structure should have same number of atoms as ensemble")
        else:
            logger.warning("Reference structure not provided so selection will be ignored.")

        if self.superpose == YES_SUP:
            self.outEns.superpose()
        elif self.superpose == ITERPOSE:
            self.outEns.iterpose()
            
        if self.savePDBs.get():
            labels = self.outEns.getLabels()
            refLabel = labels[0]
            self.pdbs = SetOfAtomStructs().create(self._getExtraPath())
            for i, coordset in enumerate(self.outEns.getCoordsets()):
                atoms = self.atoms.select(selstr).copy()
                atoms.setCoords(coordset)
                filename = self._getExtraPath('{:s}_{:06d}_{:s}.pdb'.format(refLabel, i+1, labels[i]))
                prody.writePDB(filename, atoms)
                pdb = AtomStruct(filename)
                setattr(pdb, ENSEMBLE_WEIGHTS, pwobj.Float(self.weights[i]))
                self.pdbs.append(pdb)

        if self.writeDCDFile.get():
            prody.writeDCD(self._getPath(ENS_FILENAME), self.outEns)
            prody.writePDB(self._getPath(PDB_FILENAME), self.outEns.getAtoms())
            if self.inputPsf.get() is not None:
                psfAtoms = prody.parsePSF(self.inputPsf.get())
                if selstr=='all':
                    os.symlink(self.inputPsf.get(), self._getPath(PSF_FILENAME))
                else:
                    prody.writePSF(self._getPath(PSF_FILENAME), psfAtoms.select(selstr))

        self.filename = prody.saveEnsemble(self.outEns, self._getExtraPath('ensemble.ens.npz'))

        self.npz = ProDyNpzEnsemble().create(self._getExtraPath())
        for i in range(self.outEns.numConfs()):
            frame = TrajFrame((i+1, self.filename), objLabel=self.outEns.getLabels()[i],
                              weight=pwobj.Float(self.weights[i]))
            self.npz.append(frame)

        restoreVerbositySecondary(self)

    def createOutputStep(self):
        outputs = {"outputNpz": self.npz}

        if self.writeDCDFile.get():
            if imported_chem:
                outMDSystem = MDSystem(filename=self._getPath(PDB_FILENAME))
                if os.path.exists(self._getPath(PSF_FILENAME)):
                    outMDSystem.setTopologyFile(self._getPath(PSF_FILENAME))
                else:
                    outMDSystem.setTopologyFile(self._getPath(PDB_FILENAME))
                outMDSystem.setTrajectoryFile(self._getPath(ENS_FILENAME))
                
                outputs["outputTrajectory"] = outMDSystem
            else:
                outEMFile = EMFile(filename=self._getPath(ENS_FILENAME))
                outputs["outputTrajectory"] = outEMFile

        if self.savePDBs.get():
            outputs["outputStructures"] = self.pdbs

        self._defineOutputs(**outputs)

    def _summary(self):
        if not hasattr(self, 'outputNpz'):
            summ = ['Output ensemble not ready yet']
        else:
            if len(self.outputNpz) < 100:
                ens = self.outputNpz.loadEnsemble()
                summ = ['Ensemble imported with *{0}* structures of *{1}* atoms'.format(
                    ens.numConfs(), ens.numAtoms())]
            else:
                summ = ['Ensemble imported with *{0}* structures'.format(len(self.outputNpz))]
        return summ

    def _getImportChoices(self):
        return ['files', 'pointer']

