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
This module will provide ProDy atom tools including selection and superposition.
"""
import numpy as np
from os.path import basename, splitext

from pwem.objects import AtomStruct, SetOfAtomStructs, Transform, CsvList, Integer
from pwem.protocols import EMProtocol

from pyworkflow.utils import exists
from pyworkflow.protocol.params import (PointerParam, StringParam, FloatParam,
                                        BooleanParam, EnumParam, TextParam, IntParam,
                                        PathParam, MultiPointerParam, LEVEL_ADVANCED)

from prody2 import Plugin
from prody2.constants import (NOTHING, PWALIGN, CEALIGN, DEFAULT,  # residue mapping methods
                              BEST_MATCH, SAME_CHID, SAME_POS, CUSTOM, # chain matching
                              N_ATOMS, N_RESIDUES, N_CHAINS,
                              FIRST_RESNUM, LAST_RESNUM, MAX_RESNUM, MIN_RESNUM)

def notFoundException(inputFn):
    return Exception("Atomic structure not found at *%s*" % inputFn)

UNITE_CHAINS_LABEL = "Unite chains in mmCIF segments"
UNITE_CHAINS_HELP = ('Elect whether to unite chains in mmCIF segments for each structure like ChimeraX. '
                     'Default is **False**, which means the smaller unit IDs are used for chains like PyMOL.')

IMPORT_FROM_ID_CONDITION = 'inputPdbData == IMPORT_FROM_ID'
SUMMARY_NO_OUTPUT = 'Output structure not ready yet'

class ProDyAtomicBase(EMProtocol):
    """
    This protocol will perform atom selection
    """
    _label = 'Select'
    IMPORT_FROM_ID = 0
    IMPORT_FROM_FILES = 1
    USE_POINTER = 2

    _possibleOutputs = {'outputStructure': AtomStruct}

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form, includeSelection=True):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params
        """
        # You need a params to belong to a section:
        form.addSection(label='ProDy Select')

        form.addParam('inputPdbData', EnumParam, choices=['id', 'file', 'pointer'],
                      label="Import atomic structure from",
                      default=self.USE_POINTER,
                      display=EnumParam.DISPLAY_HLIST,
                      help='Import PDB or mmCIF data from online server or local file')
        form.addParam('pdbId', StringParam,
                      condition=IMPORT_FROM_ID_CONDITION,
                      label="Atomic structure ID ", allowsNull=True,
                      help='Type a PDB ID (four alphanumeric characters).')
        form.addParam('pdbFile', PathParam, label="File path",
                      condition='inputPdbData == IMPORT_FROM_FILES',
                      allowsNull=True,
                      help='Specify a path to desired atomic structure.')
        form.addParam('inputStructure', PointerParam, label="Input structure",
                      condition='inputPdbData == USE_POINTER',
                      pointerClass='AtomStruct',
                      help='The input structure can be an atomic model '
                           '(true PDB) or a pseudoatomic model\n'
                           '(an EM volume converted into pseudoatoms)')

        form.addParam('selection', StringParam, default="protein and name CA or nucleic and name P C4' C2",
                      label="selection string", condition=includeSelection,
                      help='This determines which atoms are selected. '
                           'There is a rich selection engine with similarities to VMD. '
                           'See http://http://www.bahargroup.org/prody/tutorials/prody_tutorial/selection.html')

        form.addParam('uniteChains', BooleanParam, default=False,
                      label=UNITE_CHAINS_LABEL,
                      help=UNITE_CHAINS_HELP)

    def getPdbFileName(self, inputFn):
        return self._getPath(splitext(basename(inputFn))[0] + '_atoms.pdb')

    def getInputFn(self):
        if self.inputPdbData == self.IMPORT_FROM_ID:
            args = '--pdb {0} --folder {1}'.format(self.pdbId.get(), self._getPath())
            self.runJob(Plugin.getProgram('parse.py', script=True), args)
            with open(self._getPath('inputFn.txt'), 'r') as fi:
                inputFn = fi.readlines()[0]
        elif self.inputPdbData == self.IMPORT_FROM_FILES:
            inputFn = self.pdbFile.get()
            if not exists(inputFn):
                raise notFoundException(inputFn)
        else:
            inputFn = self.inputStructure.get().getFileName()

        return inputFn

class ProDySelect(ProDyAtomicBase):
    """
    This protocol will perform atom selection
    """     

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.selectionStep)
        self._insertFunctionStep(self.createOutputStep)

    def selectionStep(self):
        # First handle inputs
        self.inputFn = self.getInputFn()

        # Then actually perform the selection
        self.pdbFileName = self.getPdbFileName(self.inputFn)
        args = '"{0}" {1} -o {2}'.format(str(self.selection), self.inputFn,
                                         self.pdbFileName)
        if self.uniteChains.get():
            args += '--unite-chains'
        self.runJob(Plugin.getProgram('select'), args)

    def createOutputStep(self):
        if not self.hasAttribute('inputFn'):
            self.inputFn = self.getInputFn()
        if not self.hasAttribute('pdbFileName'):
            self.pdbFileName = self.getPdbFileName(self.inputFn)
        if exists(self.pdbFileName):
            outputPdb = AtomStruct()
            outputPdb.setFileName(self.pdbFileName)
            self._defineOutputs(outputStructure=outputPdb)

    def _summary(self):
        if not hasattr(self, 'outputStructure'):
            if self.getStatus() != 'finished':
                summ = [SUMMARY_NO_OUTPUT]
            else:
                summ = ['No atoms match selection so no output structure']
        else:
            summ = ['Atom selection successful']
        return summ


class ProDyAlign(EMProtocol):
    """
    This protocol will perform atomic structure mapping and superposition
    """
    _label = 'Align'
    _possibleOutputs = {'outputStructureMob': AtomStruct,
                        'outputStructureTar': AtomStruct,
                        'outputTransform': Transform}

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params
        """
        # You need a params to belong to a section:
        form.addSection(label='ProDy Align')

        form.addParam('mobStructure', PointerParam, label="Mobile structure",
                      important=True,
                      pointerClass='AtomStruct',
                      help='The structure to be moved can be an atomic model '
                           '(true PDB) or a pseudoatomic model\n'
                           '(an EM volume converted into pseudoatoms).')

        form.addParam('tarStructure', PointerParam, label="Target structure",
                      important=True, allowsNull=True,
                      pointerClass='AtomStruct',
                      help='The target structure can be an atomic model '
                           '(true PDB) or a pseudoatomic model\n'
                           '(an EM volume converted into pseudoatoms). '
                           'If no target is provided but a transformation is, '
                           'then the whole mobile structure is moved.')

        form.addParam('uniteChains', BooleanParam, default=False,
                      label=UNITE_CHAINS_LABEL,
                      help=UNITE_CHAINS_HELP)

        form.addParam('seqid', FloatParam, default=100.,
                      expertLevel=LEVEL_ADVANCED,
                      label="Sequence Identity Cut-off (%)",
                      help='Alignment mapping with lower sequence identity will not be accepted.\n'
                           'This should be a number between 0 and 100')

        form.addParam('overlap', FloatParam, default=100.,
                      expertLevel=LEVEL_ADVANCED,
                      label="Overlap Cut-off (%)",
                      help='Alignment mapping with lower sequence coverage will not be accepted.\n'
                           'This should be a number between 0 and 100') 

        form.addParam('matchFunc', EnumParam, choices=['bestMatch', 'sameChid', 'sameChainPos', 'custom'], 
                      default=BEST_MATCH,
                      label="Chain matching function",
                      help='Chains can be matched by either trying all combinations and taking the best one '
                           'based on a number of criteria including final RMSD or by taking chains with the same ID '
                           'or position in the list of chains.\n'
                           'See http://http://www.bahargroup.org/prody/manual/release/v1.11_series.html for more details.')
        
        matchFuncCheck = 'matchFunc == %d'
        group = form.addGroup('Custom chain orders', condition=matchFuncCheck % CUSTOM)
        
        group.addParam('chainOrders', TextParam, width=30, readOnly=True,
                       condition=matchFuncCheck % CUSTOM,
                       label='Custom chain match list',
                       help='Defined order of chains from custom matching')
        
        group.addParam('insertOrder', EnumParam, choices=['1. mobile', '2. target'], default=0,
                       condition=matchFuncCheck % CUSTOM,
                       label='Insert custom match order number',
                       help='Insert the chain order with the specified index into the match list.\n'
                            'The default (when empty) is the last position')
        
        group.addParam('customOrder', StringParam, default='',
                       condition=matchFuncCheck % CUSTOM,
                       label='Custom match order to insert at the specified number',
                       help='Enter the desired chain order here.\n'
                            'The default (when empty) is the chain order in the structure file')
        
        group.addParam('label', StringParam, default='', readOnly=True,
                       condition=matchFuncCheck % CUSTOM,
                       label='Label for item with the specified number for custom match',
                       help='This cannot be changed by the user and is for display only.')

        group.addParam('recoverOrder', EnumParam, choices=['1. mobile', '2. target'], default=0,
                       condition=matchFuncCheck % CUSTOM,
                       label='Recover custom match order number',
                       help='Enter the desired chain order here.\n'
                            'Recover the chain order with the specified index from the match list.')

        form.addParam('mapping', EnumParam, choices=['Nothing',
                                                     'Biopython pwalign local sequence alignment',
                                                     'Combinatorial extension (CE) structural alignment',
                                                     'Auto (try pwalign then ce)'],
                      default=PWALIGN,
                      expertLevel=LEVEL_ADVANCED,
                      label="Residue mapping function",
                      help='This method will be used for matching residues if the residue numbers and types aren\'t identical. \n'
                           'See http://http://www.bahargroup.org/prody/manual/reference/proteins/compare.html?highlight=mapchainontochain#prody.proteins.compare.mapChainOntoChain '
                           'for more details.')

        form.addParam('rmsd_reject', FloatParam, default=15.,
                      expertLevel=LEVEL_ADVANCED,
                      label="Rejection RMSD (A)",
                      help='Alignments with worse RMSDs than this will be rejected.')

        form.addParam('use_trans', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Use existing transformation?",
                      help='Set to True to select a previously-calculated transformation.')
        form.addParam('transformation', PointerParam,
                      pointerClass='Transform',
                      expertLevel=LEVEL_ADVANCED,
                      condition="use_trans==True",
                      label="Existing transformation",
                      help='Previously-calculated transformations can be applied instead.')

        form.addParam('keepMismatching', BooleanParam, default=False,
                      label="Keep mismatching atoms?",
                      help="If this is set to True then only the ")

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('alignStep')
        self._insertFunctionStep('createOutputStep')

    def alignStep(self):
        """This step includes alignment mapping and superposition"""

        mobFn = self.mobStructure.get().getFileName()
        self.pdbFileNameMob = self.getPdbFileNameMob()
        args = '--mobFn {0} --uniteChains {1} --folder {2}'.format(
            mobFn, self.uniteChains.get(), self._getPath())
        
        if self.tarStructure.hasValue():
            tarFn = self.tarStructure.get().getFileName()
            self.pdbFileNameTar = self.getPdbFileNameTar()

            args += ' --tarFn {0} --keepMismatching {1} --seqid {2}' \
                ' --overlap {3} --rmsdReject {4}'.format(
                tarFn, self.keepMismatching.get(), self.seqid.get(),
                self.overlap.get(), self.rmsd_reject.get())

            matchFuncId = self.matchFunc.get()
            args += ' --matchFunc {0}'.format(matchFuncId)

            if matchFuncId == CUSTOM:
                chmapFn = self._getPath('chmap.txt')
                fo = open(chmapFn, 'w')
                fo.write(eval(self.chainOrders.get()))
                fo.close()
                args += ' --chmapFn {0}'.format(chmapFn)

            mappingId = self.mapping.get()
            args += ' --mapping {0}'.format(mappingId)

            if self.transformation.get() is not None:
                transformation = self.transformation.get().getMatrix()
                transFn = self._getPath('transform.txt')
                np.savetxt(transFn, transformation)
                args += ' --transformationFn {0}'.format(transFn)                
        else:
            transformation = self.transformation.get().getMatrix()
            transFn = self._getPath('transform.txt')
            np.savetxt(transFn, transformation)
            args += ' --transformationFn {0}'.format(transFn)

        self.runJob(Plugin.getProgram('align.py', script=True), args)

    def createOutputStep(self):
        if hasattr(self, "pdbFileNameMob"):
            outputPdbMob = AtomStruct()
            outputPdbMob.setFileName(self.getPdbFileNameMob())
            self._defineOutputs(outputStructureMob=outputPdbMob)
            
        if hasattr(self, "pdbFileNameTar"):
            outputPdbTar = AtomStruct()
            outputPdbTar.setFileName(self.getPdbFileNameTar())

            outputTrans = Transform()
            transMatrix = np.loadtxt(self.getTransFileName())
            outputTrans.setMatrix(transMatrix)

            self._defineOutputs(outputStructureMob=outputPdbMob,
                                outputStructureTar=outputPdbTar,
                                outputTransform=outputTrans)

    def countMatches(self):
        matchesStr = self.chainOrders.get() if self.chainOrders.get() is not None else ''
        matches = matchesStr.split('\n')
        return len(matches) - 1
    
    def getMaxMatches(self):
        return self._maxMatches
    
    def createMatchDic(self, index):

        args += ' --mobFn {0} --tarFn {1} --uniteChains {2}' \
                ' --chainOrders {3} --customOrder {4} --index {5} --folder {6}'.format(
                self.mobStructure.get().getFileName(), 
                self.tarStructure.get().getFileName(), self.uniteChains.get(),
                self.chainOrders.get(), self.customOrder.get(), int(index), self._getPath())
        
        self.runJob(Plugin.getProgram('alignMatchDic.py', script=True), args)

        with open(self._getPath('matchDic.txt')) as fo:
            self.matchDic = fo.readlines()

        return self.matchDic
    
    def getInitialMobileChainOrder(self):
        return ''.join([ch.getChid() for ch in self.mob.iterChains()])

    def getInitialTargetChainOrder(self):
        return ''.join([ch.getChid() for ch in self.tar.iterChains()])

    def getPdbFileNameMob(self):
        return self._getPath('mobile.pdb')

    def getPdbFileNameTar(self):
        return self._getPath('target.pdb')

    def getTransFileName(self):
        return self._getPath('transformation.txt')

    def _validate(self):
        errors = []
        if not (self.tarStructure.hasValue() or (self.use_trans.get() 
                                                 and self.transformation.hasValue())):
            errors.append('A target structure or transformation matrix must be provided')

        return errors


class ProDyBiomol(ProDyAtomicBase):
    """
    This protocol will extract biologically relevant molecular assemblies,
    including orientations of proteins in membranes from OPM.
    """
    _label = 'Biomol'
    IMPORT_FROM_ID = 0
    IMPORT_FROM_FILES = 1
    USE_POINTER = 2

    _possibleOutputs = {'outputStructure': AtomStruct}

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form, includeSelection=False):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params
        """
        # You need a params to belong to a section:
        ProDyAtomicBase._defineParams(self, form, includeSelection)

        form.addParam('membrane', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      condition=IMPORT_FROM_ID_CONDITION,
                      label="Download membrane placement model?",
                      help='Use the OPM database to to model placement in the membrane.')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('extractionStep')
        self._insertFunctionStep('createOutputStep')

    def extractionStep(self):
        self.inputFn = self.getInputFn()

        args = '--inputFn {0} --uniteChains {1} --folder {2}'.format(
            self.inputFn, self.uniteChains.get(), self._getPath())
        self.runJob(Plugin.getProgram('biomol.py', script=True), args)
        with open(self._getPath('pdb_data.txt'), 'r') as fi:
            lines = fi.readlines()

        self.pdbs = SetOfAtomStructs().create(self._getExtraPath())
        for line in lines:
            filename, numAtoms, numResidues, numChains = line.split('\t')
            pdb = AtomStruct(filename)
            setattr(pdb, N_ATOMS, Integer(numAtoms))
            setattr(pdb, N_RESIDUES, Integer(numResidues))
            setattr(pdb, N_CHAINS, Integer(numChains))
            self.pdbs.append(pdb)

    def createOutputStep(self):
        self._defineOutputs(outputStructures=self.pdbs)

    def _summary(self):
        if not hasattr(self, '_summ'):
            self._summ = CsvList()

        if not hasattr(self, 'outputStructures'):
            self._summ = CsvList()
            self._summ.append(SUMMARY_NO_OUTPUT)
        else:
            if len(self._summ) == 0 or not self._summ[0].startswith('Extracted'):
                self._summ = CsvList()
                numStructs = len(self.outputStructures)
                self._summ.append('Extracted *{0}* biomolecular assemblies'.format(numStructs))
        return self._summ


class ProDyAddPDBs(EMProtocol):
    """
    This protocol will add pdb/mmcif files together into a single pdb file
    """
    _label = 'Add PDBs'
    IMPORT_FROM_ID = 0
    IMPORT_FROM_FILES = 1
    USE_POINTER = 2

    _possibleOutputs = {'outputStructure': AtomStruct}

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params
        """
        # You need a params to belong to a section:
        form.addSection(label='ProDy Add PDBs')

        form.addParam('inputStructure', MultiPointerParam, label="Input structures",
                      important=True,
                      pointerClass='AtomStruct',
                      help='Each input structures should be an atomic model '
                           '(true PDB) or a pseudoatomic model\n'
                           '(an EM volume converted into pseudoatoms)')

        form.addParam('uniteChains', BooleanParam, default=False,
                      label=UNITE_CHAINS_LABEL,
                      help=UNITE_CHAINS_HELP)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):

        self._insertFunctionStep('additionStep')
        self._insertFunctionStep('createOutputStep')

    def additionStep(self):

        pdbs = ' '.join([struct.get().getFileName() for struct in self.inputStructure])
        args = '--inputFns "{0}" --uniteChains {1} --folder {2}'.format(
            pdbs, self.uniteChains.get(), self._getPath())
        self.runJob(Plugin.getProgram('add_pdbs.py', script=True), args)

    def createOutputStep(self):
        with open(self._getPath('pdb_data.txt'), 'r') as fi:
            line = fi.readlines()[0]

        self.pdbFileName, numAtoms, numResidues, numChains = line.split('\t')
        if exists(self.pdbFileName):
            outputPdb = AtomStruct()
            setattr(outputPdb, N_ATOMS, Integer(numAtoms))
            setattr(outputPdb, N_RESIDUES, Integer(numResidues))
            setattr(outputPdb, N_CHAINS, Integer(numChains))
            outputPdb.setFileName(self.pdbFileName)
            self._defineOutputs(outputStructure=outputPdb)

    def _summary(self):
        if not hasattr(self, 'outputStructure'):
            summ = [SUMMARY_NO_OUTPUT]
        else:
            summ = ['The new structure has *{0}* residues '
                     'and *{1}* atoms in *{2}* chains'.format(
                     self.outputStructure.getAttributeValue(N_RESIDUES),
                     self.outputStructure.getAttributeValue(N_ATOMS),
                     self.outputStructure.getAttributeValue(N_CHAINS))]
        return summ

class ProDyRenumber(ProDyAtomicBase):
    """
    This protocol will perform atom renumbering
    """
    _label = 'Renumber'
    IMPORT_FROM_ID = 0
    IMPORT_FROM_FILES = 1
    USE_POINTER = 2

    _possibleOutputs = {'outputStructure': AtomStruct}

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form, includeSelection=True):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params
        """
        ProDyAtomicBase._defineParams(self, form, includeSelection)
        
        form.addParam('offset', IntParam, default=0,
                      label="Renumbering offset",
                      help='This number is added to all the residue numbers of the selection')

        form.addParam('chain', StringParam, default='',
                      label="New chain ID",
                      help='This will replace the chain ID of all atoms in the selection')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):

        if self.inputPdbData == self.USE_POINTER:
            inputFn = self.inputStructure.get().getFileName()
        else:
            inputFn = self.pdbId.get()

        self._insertFunctionStep('renumStep', inputFn)
        self._insertFunctionStep('createOutputStep')

    def renumStep(self, inputFn):

        args = '--inputFn {0} --uniteChains {1} --folder {2}'.format(
            inputFn, self.uniteChains.get(), self._getPath())

        chain = self.chain.get()
        if chain == '':
            chain = ' '
        args += ' --selection "{0}" --offset {1} --chain "{2}"'.format(
            self.selection.get(), self.offset.get(), chain)

        self.runJob(Plugin.getProgram('renumber_pdbs.py', script=True), args)

    def createOutputStep(self):
        with open(self._getPath('pdb_data.txt'), 'r') as fi:
            line = fi.readlines()[0]

        (self.pdbFileName, numAtoms, numResidues, numChains,
         firstResnum, lastResnum, maxResnum, minResnum) = line.split('\t')
        if exists(self.pdbFileName):
            outputPdb = AtomStruct(self.pdbFileName)
            setattr(outputPdb, N_ATOMS, Integer(numAtoms))
            setattr(outputPdb, N_RESIDUES, Integer(numResidues))
            setattr(outputPdb, N_CHAINS, Integer(numChains))
            setattr(outputPdb, FIRST_RESNUM, Integer(firstResnum))
            setattr(outputPdb, LAST_RESNUM, Integer(lastResnum))
            setattr(outputPdb, MAX_RESNUM, Integer(maxResnum))
            setattr(outputPdb, MIN_RESNUM, Integer(minResnum))
            self._defineOutputs(outputStructure=outputPdb)

    def _summary(self):
        if not hasattr(self, 'outputStructure'):
            if self.getStatus() != 'finished':
                summ = [SUMMARY_NO_OUTPUT]
            else:
                summ = ['No atoms match selection so no output structure']
        else:
            summ = ['Selected *{0}* atoms'.format(
                self.outputStructure.getAttributeValue(N_ATOMS))]
            if self.hasAttribute('outputStructure'):
                summ.append('The new structure has *{0}* residues'.format(
                            self.outputStructure.getAttributeValue(N_RESIDUES)))
            else:
                summ.append('The new structure has *0* residues')
        return summ
