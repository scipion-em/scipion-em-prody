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
from collections import OrderedDict
from os.path import basename, splitext, abspath

from pwem.objects import AtomStruct, SetOfAtomStructs, Transform, CsvList
from pwem.protocols import EMProtocol

from pyworkflow.utils import exists
from pyworkflow.protocol.params import (PointerParam, StringParam, FloatParam,
                                        BooleanParam, EnumParam, TextParam, IntParam,
                                        PathParam, MultiPointerParam, LEVEL_ADVANCED)

import prody
from pyworkflow.utils import logger

from prody2 import Plugin
from prody2.objects import Atom, SetOfAtoms
from prody2.constants import (NOTHING, PWALIGN, CEALIGN, DEFAULT,  # residue mapping methods
                              BEST_MATCH, SAME_CHID, SAME_POS, CUSTOM) # chain matching

def notFoundException(inputFn):
    return Exception("Atomic structure not found at *%s*" % inputFn)

UNITE_CHAINS_LABEL = "Unite chains in mmCIF segments"
UNITE_CHAINS_HELP = ('Elect whether to unite chains in mmCIF segments for each structure like ChimeraX. '
                     'Default is **False**, which means the smaller unit IDs are used for chains like PyMOL.')

IMPORT_FROM_ID_CONDITION = 'inputPdbData == IMPORT_FROM_ID'
SUMMARY_NO_OUTPUT = 'Output structure not ready yet'
NOT_DUMMY_SELSTR = "not dummy"


class ProDyAtomicProtocols(EMProtocol):
    """
    Collection of ProDy-based protocols for atomic structure selection,
    alignment, biomolecular assembly extraction, structure merging,
    metadata conversion, and residue renumbering.

    AI Generated:

    ProDy Atomic Protocols — User Manual
        Overview

        This collection of protocols provides several structure-level
        operations for atomic models such as PDB or mmCIF files. These
        protocols are intended for structural biology workflows where
        users need to manipulate, compare, transform, or extract
        biologically relevant subsets of atomic coordinates.

        The main operations covered are:

            - Atom selection
            - Structure alignment and superposition
            - Biomolecular assembly extraction
            - Merging multiple structures
            - Conversion of structures into atom-level metadata
            - Residue renumbering and chain reassignment

        All protocols are designed to work with standard atomic models
        as well as pseudoatomic models derived from EM volumes.

        General Input Strategy

        Most protocols accept atomic structures through one of three
        mechanisms:

            1. Import by PDB identifier
            2. Import from a local PDB/mmCIF file
            3. Import from an existing AtomStruct object

        This flexibility allows seamless integration with both external
        structural databases and existing Scipion workflows.

        ------------------------------------------------------------------
        ProDyAtomicBase — Shared Atomic Input Definition
        ------------------------------------------------------------------

        This base protocol defines the common atomic input parameters
        shared by several derived protocols.

        It provides:

            - Atomic structure import options
            - Optional atom selection strings
            - Optional chain unification

        The selection syntax follows the ProDy selection language,
        allowing users to define subsets of atoms based on residue type,
        atom names, chains, and many other structural attributes.

        In practical biological workflows, this makes it possible to
        isolate:

            - Protein alpha carbons
            - Nucleic acid backbone atoms
            - Specific chains
            - Domain-specific regions

        ------------------------------------------------------------------
        ProDySelect — Atom Selection
        ------------------------------------------------------------------

        The ProDySelect protocol extracts a subset of atoms from an
        input atomic structure according to a user-defined selection
        string.

        Workflow

            1. Load the input structure
            2. Apply the selection expression
            3. Write the selected atoms into a new PDB file

        Biological Interpretation

        This protocol is especially useful when users want to isolate
        biologically meaningful subsets such as:

            - Protein backbone atoms
            - Ligand-free protein regions
            - Nucleic acid atoms only
            - Domain-specific structural fragments

        Outputs

            - A new atomic structure containing only the selected atoms

        Summary Information

        The protocol reports:

            - Number of selected atoms
            - Number of original atoms
            - Number of protein residues retained

        Practical Note

        If no atoms satisfy the selection criteria, no output structure
        is generated.

        ------------------------------------------------------------------
        ProDyAlign — Atomic Structure Mapping and Superposition
        ------------------------------------------------------------------

        The ProDyAlign protocol performs structural correspondence
        mapping and rigid-body superposition between a mobile structure
        and a target structure.

        Biological Purpose

        This protocol is useful for:

            - Comparing homologous proteins
            - Mapping conformational states
            - Evaluating structural similarity
            - Preparing structures for direct comparison

        Main Workflow

            1. Parse mobile and target structures
            2. Match chains using one of several strategies
            3. Map residues between structures
            4. Compute optimal rigid transformation
            5. Apply transformation to the mobile structure
            6. Save aligned structures and transformation matrix

        Chain Matching Strategies

            - bestMatch:
              explores chain combinations and chooses the best match

            - sameChid:
              matches chains with identical chain identifiers

            - sameChainPos:
              matches chains according to chain order

            - custom:
              user-defined chain correspondence

        Residue Mapping Options

            - Sequence alignment
            - Structural alignment
            - Automatic hybrid mapping
            - No mapping

        Outputs

            - Aligned mobile structure
            - Target structure
            - Transformation matrix

        RMSD Evaluation

        The protocol computes RMSD after alignment, providing an
        immediate structural similarity measure.

        Biological Considerations

        High RMSD values may indicate:

            - Low structural similarity
            - Incorrect chain mapping
            - Flexible domain differences
            - Large conformational rearrangements

        Advanced Option

        Users may optionally keep mismatching atoms, which can be useful
        when analyzing insertions, flexible regions, or incomplete
        structural correspondence.

        ------------------------------------------------------------------
        ProDyBiomol — Biomolecular Assembly Extraction
        ------------------------------------------------------------------

        This protocol extracts biologically relevant biomolecular
        assemblies from atomic structure files.

        Biological Relevance

        Structures deposited in the PDB often contain only the
        asymmetric unit. This protocol reconstructs biologically
        meaningful assemblies such as:

            - Homodimers
            - Multimeric complexes
            - Biological oligomers

        Optional OPM Support

        When importing from a PDB identifier, users may also request
        membrane placement models from the OPM database.

        Workflow

            1. Load the input structure
            2. Expand biomolecular assemblies
            3. Write each assembly as a separate structure

        Outputs

            - A set of atomic structures, one per biological assembly

        Summary Information

        For each extracted assembly, the protocol reports:

            - Number of residues
            - Number of chains

        Practical Importance

        This protocol is especially useful for membrane proteins,
        oligomeric enzymes, and complexes where the deposited asymmetric
        unit does not correspond to the functional biological state.

        ------------------------------------------------------------------
        ProDyAddPDBs — Merge Atomic Structures
        ------------------------------------------------------------------

        This protocol merges multiple atomic structures into a single
        output structure.

        Workflow

            1. Read all input structures
            2. Concatenate atom groups
            3. Write the merged structure

        Biological Applications

        Useful for:

            - Building composite assemblies
            - Combining chains from separate files
            - Preparing multicomponent models

        Output

            - A single merged atomic structure

        Summary Information

        The protocol reports:

            - Number of protein residues
            - Total number of atoms
            - Total number of chains

        Important Consideration

        The protocol concatenates structures directly and does not
        perform collision detection or spatial optimization.

        ------------------------------------------------------------------
        ProDyToBiopythonMetadata — Convert Structure to Atom Metadata
        ------------------------------------------------------------------

        This protocol converts an atomic structure into a metadata-like
        collection of atom objects.

        Workflow

            1. Parse the input structure
            2. Create one metadata atom entry per atomic coordinate
            3. Store the resulting SetOfAtoms

        Biological Use

        This protocol is useful when downstream operations require
        atom-level indexing or metadata representation rather than a
        coordinate-only atomic model.

        Output

            - SetOfAtoms metadata object

        Summary Information

        The protocol reports the total number of atoms converted.

        ------------------------------------------------------------------
        ProDyRenumber — Residue Renumbering and Chain Reassignment
        ------------------------------------------------------------------

        The ProDyRenumber protocol modifies residue numbering and
        optionally replaces chain identifiers for a selected subset of
        atoms.

        Workflow

            1. Load input structure
            2. Apply atom selection
            3. Add an integer offset to residue numbers
            4. Optionally assign a new chain identifier
            5. Write modified structure

        Biological Applications

        This protocol is particularly useful when:

            - Matching numbering between homologous structures
            - Standardizing residue indices before modeling
            - Preparing structures for comparative analysis
            - Resolving chain naming inconsistencies

        Important Practical Note

        Only atoms satisfying the selection string are modified.

        Outputs

            - A renumbered atomic structure

        Summary Information

        The protocol reports:

            - Number of selected atoms
            - Number of original atoms
            - Number of retained protein residues

        ------------------------------------------------------------------
        Final Perspective
        ------------------------------------------------------------------

        These ProDy atomic protocols provide a compact but highly useful
        toolbox for atomic-level structural manipulation.

        In practical structural biology workflows, they allow users to:

            - isolate relevant structural regions,
            - compare homologous models,
            - reconstruct biological assemblies,
            - merge independent structural components,
            - transform coordinate models into metadata representations,
            - standardize residue numbering.

        Although computationally straightforward, these operations often
        have strong biological consequences. Careful selection of chains,
        residues, mapping strategies, and biological assemblies is
        essential for obtaining meaningful downstream structural
        interpretation.
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

class ProDySelect(ProDyAtomicBase):
    """
    This protocol will perform atom selection
    """     

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):

        if self.inputPdbData == self.IMPORT_FROM_ID:
            prody.pathPDBFolder(self.getPath(""))
            inputFn = prody.fetchPDB(self.pdbId.get(), compressed=False)
            
            if inputFn == None:
                inputFn = prody.fetchPDB(self.pdbId.get(), format="cif",
                                         compressed=False)

            prody.pathPDBFolder("")

        elif self.inputPdbData == self.IMPORT_FROM_FILES:
            inputFn = self.pdbFile.get()
            if not exists(inputFn):
                raise notFoundException(inputFn)

        else:
            inputFn = self.inputStructure.get().getFileName()

        self.inputStruct = AtomStruct()
        self.inputStruct.setFileName(inputFn)

        self._insertFunctionStep('selectionStep', inputFn)
        self._insertFunctionStep('createOutputStep')

    def selectionStep(self, inputFn):
        self.pdbFileName = self._getPath(splitext(basename(inputFn))[0] + '_atoms.pdb')
        args = '"{0}" {1} -o {2}'.format(str(self.selection), inputFn,
                                         self.pdbFileName)
        if self.uniteChains.get():
            args += '--unite-chains'
        self.runJob(Plugin.getProgram('select'), args)

    def createOutputStep(self):
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
            inputAg = prody.parsePDB(self.inputStruct.getFileName(),
                                     unite_chains=self.uniteChains.get())
            outputAg = prody.parsePDB(self.outputStructure.getFileName(),
                                      unite_chains=self.uniteChains.get())

            summ = ['Selected *{0}* atoms from original *{1}* atoms'.format(
                outputAg.numAtoms(), inputAg.numAtoms())]
            if outputAg.ca is not None:
                summ.append('The new structure has *{0}* protein residues '
                            'from original *{1}* protein residues'.format(
                            outputAg.ca.numAtoms(), inputAg.ca.numAtoms()))
            else:
                summ.append('The new structure has *0* protein residues')
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
        mob = prody.parsePDB(mobFn, alt='all',
                             unite_chains=self.uniteChains.get())
        
        if self.tarStructure.hasValue():
            tarFn = self.tarStructure.get().getFileName()
            tar = prody.parsePDB(tarFn, alt='all',
                                unite_chains=self.uniteChains.get())

            if self.matchFunc.get() == BEST_MATCH:
                matchFunc = prody.bestMatch
                logger.info('\nUsing bestMatch\n')
            elif self.matchFunc.get() == SAME_CHID:
                matchFunc = prody.sameChid
                logger.info('\nUsing sameChid\n')
            elif self.matchFunc.get() == SAME_POS:
                matchFunc = prody.sameChainPos
                logger.info('\nUsing sameChainPos\n')
            else:
                chmap = eval(self.chainOrders.get())
                logger.info('\nUsing user-defined match function based on \n{0}\n'.format(self.chainOrders.get()))
                matchFunc = lambda chain1, chain2: prody.userDefined(chain1, chain2, chmap)

            if self.mapping.get() == DEFAULT:
                mapping = 'auto'
            elif self.mapping.get() == PWALIGN:
                mapping = 'pwalign'
            elif self.mapping.get() == CEALIGN:
                mapping = 'ce'
            else:
                mapping = False

            mobAmapList = prody.alignChains(mob.protein, tar.protein,
                                            seqid=self.seqid.get(),
                                            overlap=self.overlap.get(),
                                            match_func=matchFunc,
                                            mapping=mapping,
                                            rmsd_reject=self.rmsd_reject.get())
            if len(mobAmapList):
                mobAmap = mobAmapList[0]
                mobSel = mobAmap.select(NOT_DUMMY_SELSTR).copy()
                mobSel.setTitle(mob.getTitle())

                tarAmapList = prody.alignChains(tar.protein, mobSel,
                                                seqid=self.seqid.get(),
                                                overlap=self.overlap.get(),
                                                match_func=matchFunc,
                                                mapping=mapping,
                                                rmsd_reject=self.rmsd_reject.get())
                if len(tarAmapList):
                    tarAmap = tarAmapList[0]
                    tarSel = tarAmap.select(NOT_DUMMY_SELSTR).copy()
                    tarSel.setTitle(tar.getTitle())

                    if mobSel.numAtoms != tarSel.numAtoms():
                        mobAmapList = prody.alignChains(mobSel, tarSel,
                                                        seqid=self.seqid.get(),
                                                        overlap=self.overlap.get(),
                                                        match_func=matchFunc,
                                                        mapping=mapping,
                                                        rmsd_reject=self.rmsd_reject.get())
                        if len(mobAmapList):
                            mobAmap = mobAmapList[0]
                            mobSel = mobAmap.select(NOT_DUMMY_SELSTR).copy()
                            mobSel.setTitle(mob.getTitle())

                        tarAmapList = prody.alignChains(tarSel, mobSel,
                                                        seqid=self.seqid.get(),
                                                        overlap=self.overlap.get(),
                                                        match_func=matchFunc,
                                                        mapping=mapping,
                                                        rmsd_reject=self.rmsd_reject.get())
                        if len(tarAmapList):
                            tarAmap = tarAmapList[0]
                            tarSel = tarAmap.select(NOT_DUMMY_SELSTR).copy()
                            tarSel.setTitle(tar.getTitle())

                    if self.transformation.get() is None:
                        self.T = prody.calcTransformation(mobSel, tarSel)
                    else:
                        self.T = prody.Transformation(self.transformation.get().getMatrix())

                    alg = prody.applyTransformation(self.T, mobSel)

                    self.rmsd = prody.calcRMSD(mobSel, tarSel)
                    logger.info("\nRMSD = {:6.2f}\n".format(self.rmsd))

                    if self.keepMismatching.get():
                        alg = prody.applyTransformation(self.T, mob)
                        tarSel = tar

                    self.pdbFileNameMob = self._getPath('mobile.pdb')
                    prody.writePDB(self.pdbFileNameMob, alg)

                    self.pdbFileNameTar = self._getPath('target.pdb')
                    prody.writePDB(self.pdbFileNameTar, tarSel)

                    self.matrixFileName = self._getPath('transformation.txt')
                    prody.writeArray(self.matrixFileName, self.T.getMatrix())
        else:
            if self.transformation.get() is not None:
                self.T = prody.Transformation(self.transformation.get().getMatrix())
                alg = prody.applyTransformation(self.T, mob)

                self.pdbFileNameMob = self._getPath('mobile.pdb')
                prody.writePDB(self.pdbFileNameMob, alg)

    def createOutputStep(self):
        if hasattr(self, "pdbFileNameMob"):
            outputPdbMob = AtomStruct()
            outputPdbMob.setFileName(self.pdbFileNameMob)

            self._defineOutputs(outputStructureMob=outputPdbMob)
            
        if hasattr(self, "pdbFileNameTar"):
            outputPdbTar = AtomStruct()
            outputPdbTar.setFileName(self.pdbFileNameTar)

            outputTrans = Transform()
            outputTrans.setMatrix(self.T.getMatrix())

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

        index = int(index)

        self.mob = prody.parsePDB(self.mobStructure.get().getFileName(), alt='all',
                                  unite_chains=self.uniteChains.get())
        self.tar = prody.parsePDB(self.tarStructure.get().getFileName(), alt='all',
                                  unite_chains=self.uniteChains.get())
        
        try:
            self.matchDic = eval(self.chainOrders.get())
            _ = self.matchDic.keys()
        except (AttributeError, TypeError):
            self.matchDic = OrderedDict()
            self.matchDic[self.mob.getTitle()] = self.getInitialMobileChainOrder()
            self.matchDic[self.tar.getTitle()] = self.getInitialTargetChainOrder()
            
        if index == 0:
            label = self.mob.getTitle()
            if self.customOrder.get() == '':
                self.matchDic[label] = self.getInitialMobileChainOrder()
            else:
                self.matchDic[label] = self.customOrder.get()
                
        else:
            label = self.tar.getTitle()
            if self.customOrder.get() == '':
                self.matchDic[label] = self.getInitialTargetChainOrder()
            else:
                self.matchDic[label] = self.customOrder.get()
                            
        return self.matchDic
    
    def getInitialMobileChainOrder(self):
        return ''.join([ch.getChid() for ch in self.mob.iterChains()])

    def getInitialTargetChainOrder(self):
        return ''.join([ch.getChid() for ch in self.tar.iterChains()])



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

        if self.inputPdbData == self.IMPORT_FROM_ID:
            prody.pathPDBFolder(self.getPath(""))

            if self.membrane.get():
                inputFn = prody.fetchPDBfromOPM(self.pdbId.get(), filename=self.getPath(self.pdbId.get()+'-opm.pdb'))
            else:
                inputFn = prody.fetchPDB(self.pdbId.get(), compressed=False)
            
            if inputFn == None:
                inputFn = prody.fetchPDB(self.pdbId.get(), format="cif",
                                         compressed=False)

            prody.pathPDBFolder("")

        elif self.inputPdbData == self.IMPORT_FROM_FILES:
            inputFn = self.pdbFile.get()
            if not exists(inputFn):
                raise notFoundException(inputFn)

        else:
            inputFn = self.inputStructure.get().getFileName()

        self.inputStruct = AtomStruct()
        self.inputStruct.setFileName(inputFn)

        self._insertFunctionStep('extractionStep', inputFn)
        self._insertFunctionStep('createOutputStep')

    def extractionStep(self, inputFn):
        ags = prody.parsePDB(inputFn, alt='all', compressed=False,
                             biomol=True, extend_biomol=True,
                             unite_chains=self.uniteChains.get())
        if isinstance(ags, prody.AtomGroup):
            ags = [ags] 

        self.pdbs = SetOfAtomStructs().create(self._getExtraPath())
        for i, ag in enumerate(ags):
            filename = self._getPath(splitext(basename(inputFn))[0] + '_atoms_{0}.pdb'.format(i))
            prody.writePDB(filename, ag)
            pdb = AtomStruct(filename)
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

                ags = prody.parsePDB([struct.getFileName() for struct in self.outputStructures],
                                     unite_chains=self.uniteChains.get())
                if numStructs == 1:
                    ags = [ags]
                     
                for i, ag in enumerate(ags):
                    self._summ.append('New structure {0} has *{1}* residues '
                                     'across *{2}* chains'.format(i+1, ag.numResidues(), 
                                                                 ag.numChains()))
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
        pdbs = [struct.get().getFileName() for struct in self.inputStructure]

        ags = prody.parsePDB(pdbs, unite_chains=self.uniteChains.get())

        outAg = ags[0]
        for ag in ags[1:]:
            outAg += ag

        self.pdbFileName = self._getPath('joined_atoms.pdb')
        prody.writePDB(self.pdbFileName, outAg)

    def createOutputStep(self):
        if exists(self.pdbFileName):
            outputPdb = AtomStruct()
            outputPdb.setFileName(self.pdbFileName)
            self._defineOutputs(outputStructure=outputPdb)

    def _summary(self):
        if not hasattr(self, 'outputStructure'):
            summ = [SUMMARY_NO_OUTPUT]
        else:
            outputAg = prody.parsePDB(self.outputStructure.getFileName(), 
                                      unite_chains=self.uniteChains.get())

            summ = ['The new structure has *{0}* protein residues '
                     'and *{1}* atoms in *{2}* chains'.format(
                     outputAg.ca.numAtoms(), outputAg.numAtoms(),
                     outputAg.numChains())]
        return summ


class ProDyToBiopythonMetadata(EMProtocol):
    """
    This protocol will add pdb/mmcif files together into a single pdb file
    """
    _label = 'convert to metadata'
    IMPORT_FROM_ID = 0
    IMPORT_FROM_FILES = 1
    USE_POINTER = 2

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params
        """
        # You need a params to belong to a section:
        form.addSection(label='ProDy PDB to metadata')

        form.addParam('inputStructure', PointerParam, label="Input structures",
                      important=True,
                      pointerClass='AtomStruct',
                      help='Each input structures should be an atomic model '
                           '(true PDB) or a pseudoatomic model\n'
                           '(an EM volume converted into pseudoatoms)')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep')

    def createOutputStep(self):
        filename = abspath(self.inputStructure.get().getFileName())
        struct = prody.parsePDB(filename)
        outputPdb = SetOfAtoms().create(self._getExtraPath())
        for i in range(struct.numAtoms()):
            atom = Atom((i, filename))
            outputPdb.append(atom)
        self._defineOutputs(outputStructure=outputPdb)

    def _summary(self):
        if not hasattr(self, 'outputStructure'):
            summ = [SUMMARY_NO_OUTPUT]
        else:
            summ = ['The new structure has *{0}* atoms'.format(len(self.outputStructure))]
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

        if self.inputPdbData == self.IMPORT_FROM_ID:
            prody.pathPDBFolder(self.getPath(""))
            inputFn = prody.fetchPDB(self.pdbId.get(), compressed=False)
            
            if inputFn == None:
                inputFn = prody.fetchPDB(self.pdbId.get(), format="cif",
                                         compressed=False)

            prody.pathPDBFolder("")

        elif self.inputPdbData == self.IMPORT_FROM_FILES:
            inputFn = self.pdbFile.get()
            if not exists(inputFn):
                raise notFoundException(inputFn)

        else:
            inputFn = self.inputStructure.get().getFileName()

        self.inputStruct = AtomStruct()
        self.inputStruct.setFileName(inputFn)

        self._insertFunctionStep('renumStep', inputFn)
        self._insertFunctionStep('createOutputStep')

    def renumStep(self, inputFn):
        self.pdbFileName = self._getPath(splitext(basename(inputFn))[0] + '_atoms.pdb')
        ag = prody.parsePDB(inputFn)

        sel = ag.select(self.selection.get())
        sel.setResnums(sel.getResnums() + self.offset.get())
        prody.writePDB(self.pdbFileName, ag)

        chain = self.chain.get()
        if chain != '':
            sel.setChids(chain)

    def createOutputStep(self):
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
            inputAg = prody.parsePDB(self.inputStruct.getFileName(),
                                     unite_chains=self.uniteChains.get())
            outputAg = prody.parsePDB(self.outputStructure.getFileName(),
                                      unite_chains=self.uniteChains.get())

            summ = ['Selected *{0}* atoms from original *{1}* atoms'.format(
                outputAg.numAtoms(), inputAg.numAtoms())]
            if outputAg.ca is not None:
                summ.append('The new structure has *{0}* protein residues '
                            'from original *{1}* protein residues'.format(
                            outputAg.ca.numAtoms(), inputAg.ca.numAtoms()))
            else:
                summ.append('The new structure has *0* protein residues')
        return summ