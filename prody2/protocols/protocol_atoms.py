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
from os.path import basename, splitext

from pyworkflow.protocol import params

from pwem import *
from pwem.objects import AtomStruct, Transform, SetOfAtomStructs, String, CsvList
from pwem.protocols import EMProtocol

from pyworkflow.utils import *
from pyworkflow.protocol.params import (PointerParam, StringParam, FloatParam,
                                        BooleanParam, EnumParam, TextParam,
                                        PathParam, LEVEL_ADVANCED)

import prody
import logging
logger = logging.getLogger(__name__)

# chain matching methods
BEST_MATCH = 0
SAME_CHID = 1
SAME_POS = 2
CUSTOM = 3

# residue mapping methods
NOTHING = 0 # stop trivial mapping if trivial mapping fails
PWALIGN = 1 # biopython pwalign local pairwise sequence alignment after trivial mapping
CEALIGN = 2 # combinatorial extension (CE) as in PyMOL
DEFAULT = 3 # try pwalign then CE

class ProDySelect(EMProtocol):
    """
    This protocol will perform atom selection
    """
    _label = 'Select'
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
        form.addSection(label='ProDy Select')

        form.addParam('inputPdbData', EnumParam, choices=['id', 'file', 'pointer'],
                      label="Import atomic structure from",
                      default=self.USE_POINTER,
                      display=EnumParam.DISPLAY_HLIST,
                      help='Import PDB or mmCIF data from online server or local file')
        form.addParam('pdbId', StringParam,
                      condition='inputPdbData == IMPORT_FROM_ID',
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
                      label="selection string",
                      help='This determines which atoms are selected. '
                           'There is a rich selection engine with similarities to VMD. '
                           'See http://prody.csb.pitt.edu/tutorials/prody_tutorial/selection.html')

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
                raise Exception("Atomic structure not found at *%s*" % inputFn)

        else:
            inputFn = self.inputStructure.get().getFileName()

        self.inputStruct = AtomStruct()
        self.inputStruct.setFileName(inputFn)

        self._insertFunctionStep('selectionStep', inputFn)
        self._insertFunctionStep('createOutputStep')

    def selectionStep(self, inputFn):
        # configure ProDy to automatically handle secondary structure information and verbosity
        old_secondary = prody.confProDy("auto_secondary")
        old_verbosity = prody.confProDy("verbosity")
        
        from pyworkflow import Config
        prodyVerbosity =  'none' if not Config.debugOn() else 'debug'
        prody.confProDy(auto_secondary=True, verbosity='{0}'.format(prodyVerbosity))

        ag = prody.parsePDB(inputFn, alt='all', compressed=False)
        selection = ag.select(str(self.selection))

        logger.info("%d atoms selected from %d" % (selection.numAtoms(),
                                                          ag.numAtoms()))

        self.pdbFileName = self._getPath(splitext(basename(inputFn))[0] + '_atoms.pdb')
        prody.writePDB(self.pdbFileName, selection)

        # configure ProDy to restore secondary structure information and verbosity
        prody.confProDy(auto_secondary=old_secondary, verbosity='{0}'.format(old_verbosity))

    def createOutputStep(self):
        outputPdb = AtomStruct()
        outputPdb.setFileName(self.pdbFileName)
        self._defineOutputs(outputStructure=outputPdb)   

    def _summary(self):
        if not hasattr(self, 'outputStructure'):
            sum = ['Output structure not ready yet']
        else:
            input_ag = prody.parsePDB(self.inputStruct.getFileName())
            output_ag = prody.parsePDB(self.outputStructure.getFileName())
            sum = ['Selected *{0}* atoms from original *{1}* atoms'.format(
                   output_ag.numAtoms(), input_ag.numAtoms())]
            sum.append('The new structure has *{0}* protein residues '
                        'from original *{1}* protein residues'.format(
                        output_ag.ca.numAtoms(), input_ag.ca.numAtoms()))
        return sum


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
                      important=True,
                      pointerClass='AtomStruct',
                      help='The target structure can be an atomic model '
                           '(true PDB) or a pseudoatomic model\n'
                           '(an EM volume converted into pseudoatoms).')

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
                           'See http://prody.csb.pitt.edu/manual/release/v1.11_series.html for more details.')
        
        
        group = form.addGroup('custom chain orders', condition='matchFunc == %d' % CUSTOM)
        
        group.addParam('chainOrders', TextParam, width=30, readOnly=True,
                       condition='matchFunc == %d' % CUSTOM,
                       label='Custom chain match list',
                       help='Defined order of chains from custom matching. \nManual modification will have no '
                            'effect, use the wizards to add / delete the entries')
        
        group.addParam('insertOrder', EnumParam, choices=['1. mobile', '2. target'], default=0,
                       condition='matchFunc == %d' % CUSTOM,
                       label='Insert custom match order number',
                       help='Insert the chain order with the specified index into the match list.\n'
                            'The default (when empty) is the last position')
        
        group.addParam('customOrder', StringParam, default='',
                       condition='matchFunc == %d' % CUSTOM,
                       label='Custom match order to insert at the specified number',
                       help='Enter the desired chain order here.\n'
                            'The default (when empty) is the chain order in the structure file')
        
        group.addParam('label', StringParam, default='', readOnly=True,
                       condition='matchFunc == %d' % CUSTOM,
                       label='Label for item with the specified number for custom match',
                       help='This cannot be changed by the user and is for display only.')

        group.addParam('recoverOrder', EnumParam, choices=['1. mobile', '2. target'], default=0,
                       condition='matchFunc == %d' % CUSTOM,
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
                           'See http://prody.csb.pitt.edu/manual/reference/proteins/compare.html?highlight=mapchainontochain#prody.proteins.compare.mapChainOntoChain '
                           'for more details.')

        form.addParam('rmsd_reject', FloatParam, default=15.,
                      expertLevel=LEVEL_ADVANCED,
                      label="Rejection RMSD (A)",
                      help='Alignments with worse RMSDs than this will be rejected.')

        form.addParam('use_trans', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Use existing transformation?",
                      help='Select to True to select a previously-calculated transformation.')
        form.addParam('transformation', PointerParam,
                      pointerClass='Transform',
                      expertLevel=LEVEL_ADVANCED,
                      condition="use_trans==True",
                      label="Existing transformation",
                      help='Previously-calculated transformations can be applied instead.')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('alignStep')
        self._insertFunctionStep('createOutputStep')

    def alignStep(self):
        """This step includes alignment mapping and superposition"""
        # configure ProDy to automatically handle secondary structure information and verbosity
        old_secondary = prody.confProDy("auto_secondary")
        old_verbosity = prody.confProDy("verbosity")
        
        from pyworkflow import Config
        prodyVerbosity =  'none' if not Config.debugOn() else 'debug'
        prody.confProDy(auto_secondary=True, verbosity='{0}'.format(prodyVerbosity))

        mobFn = self.mobStructure.get().getFileName()
        tarFn = self.tarStructure.get().getFileName()

        mob = prody.parsePDB(mobFn, alt='all')
        tar = prody.parsePDB(tarFn, alt='all')

        if self.matchFunc.get() == BEST_MATCH:
            match_func = prody.bestMatch
            logger.info('\nUsing bestMatch\n')
        elif self.matchFunc.get() == SAME_CHID:
            match_func = prody.sameChid
            logger.info('\nUsing sameChid\n')
        elif self.matchFunc.get() == SAME_POS:
            match_func = prody.sameChainPos
            logger.info('\nUsing sameChainPos\n')
        else:
            chmap = eval(self.chainOrders.get())
            logger.info('\nUsing user-defined match function based on \n{0}\n'.format(self.chainOrders.get()))
            match_func = lambda chain1, chain2: prody.userDefined(chain1, chain2, chmap)

        if self.mapping.get() == DEFAULT:
            mapping = 'auto'
        elif self.mapping.get() == PWALIGN:
            mapping = 'pwalign'
        elif self.mapping.get() == CEALIGN:
            mapping = 'ce'
        else:
            mapping = False

        mob_amap_list = prody.alignChains(mob.protein, tar.protein,
                                          seqid=self.seqid.get(),
                                          overlap=self.overlap.get(),
                                          match_func=match_func,
                                          mapping=mapping,
                                          rmsd_reject=self.rmsd_reject.get())
        if len(mob_amap_list):
            mob_amap = mob_amap_list[0]
            mob_sel = mob_amap.select("not dummy").copy()
            mob_sel.setTitle(mob.getTitle())

            tar_amap_list = prody.alignChains(tar.protein, mob_sel,
                                              seqid=self.seqid.get(),
                                              overlap=self.overlap.get(),
                                              match_func=match_func,
                                              mapping=mapping,
                                              rmsd_reject=self.rmsd_reject.get())
            if len(tar_amap_list):
                tar_amap = tar_amap_list[0]
                tar_sel = tar_amap.select("not dummy").copy()
                tar_sel.setTitle(tar.getTitle())

                if mob_sel.numAtoms != tar_sel.numAtoms():
                    mob_amap_list = prody.alignChains(mob_sel, tar_sel,
                                                      seqid=self.seqid.get(),
                                                      overlap=self.overlap.get(),
                                                      match_func=match_func,
                                                      mapping=mapping,
                                                      rmsd_reject=self.rmsd_reject.get())
                    if len(mob_amap_list):
                        mob_amap = mob_amap_list[0]
                        mob_sel = mob_amap.select("not dummy").copy()
                        mob_sel.setTitle(mob.getTitle())

                    tar_amap_list = prody.alignChains(tar_sel, mob_sel,
                                                      seqid=self.seqid.get(),
                                                      overlap=self.overlap.get(),
                                                      match_func=match_func,
                                                      mapping=mapping,
                                                      rmsd_reject=self.rmsd_reject.get())
                    if len(tar_amap_list):
                        tar_amap = tar_amap_list[0]
                        tar_sel = tar_amap.select("not dummy").copy()
                        tar_sel.setTitle(tar.getTitle())

                if self.transformation.get() is None:
                    self.T = prody.calcTransformation(mob_sel, tar_sel)
                else:
                    self.T = prody.Transformation(self.transformation.get().getMatrix())

                alg = prody.applyTransformation(self.T, mob_sel)

                self.rmsd = prody.calcRMSD(mob_sel, tar_sel)
                logger.info("\nRMSD = {:6.2f}\n".format(self.rmsd))

                self.pdbFileNameMob = self._getPath('mobile.pdb')
                prody.writePDB(self.pdbFileNameMob, alg)

                self.pdbFileNameTar = self._getPath('target.pdb')
                prody.writePDB(self.pdbFileNameTar, tar_sel)

                self.matrixFileName = self._getPath('transformation.txt')
                prody.writeArray(self.matrixFileName, self.T.getMatrix())

        # configure ProDy to restore secondary structure information and verbosity
        prody.confProDy(auto_secondary=old_secondary, verbosity='{0}'.format(old_verbosity))

    def createOutputStep(self):
        if hasattr(self, "pdbFileNameMob"):
            outputPdbMob = AtomStruct()
            outputPdbMob.setFileName(self.pdbFileNameMob)

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

        self.mob = prody.parsePDB(self.mobStructure.get().getFileName(), alt='all')
        self.tar = prody.parsePDB(self.tarStructure.get().getFileName(), alt='all')
        
        try:
            self.matchDic = eval(self.chainOrders.get())
            keys = self.matchDic.keys()
        except:
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



class ProDyBiomol(EMProtocol):
    """
    This protocol will extract biomolecular assemblies
    """
    _label = 'extract biomol'
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
        form.addSection(label='ProDy Select')

        form.addParam('inputPdbData', EnumParam, choices=['id', 'file', 'pointer'],
                      label="Import atomic structure from",
                      default=self.USE_POINTER,
                      display=EnumParam.DISPLAY_HLIST,
                      help='Import PDB or mmCIF data from online server or local file')
        form.addParam('pdbId', StringParam,
                      condition='inputPdbData == IMPORT_FROM_ID',
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
                raise Exception("Atomic structure not found at *%s*" % inputFn)

        else:
            inputFn = self.inputStructure.get().getFileName()

        self.inputStruct = AtomStruct()
        self.inputStruct.setFileName(inputFn)

        self._insertFunctionStep('extractionStep', inputFn)
        self._insertFunctionStep('createOutputStep')

    def extractionStep(self, inputFn):
        # configure ProDy to automatically handle secondary structure information and verbosity
        old_secondary = prody.confProDy("auto_secondary")
        old_verbosity = prody.confProDy("verbosity")
        
        from pyworkflow import Config
        prodyVerbosity =  'none' if not Config.debugOn() else 'debug'
        prody.confProDy(auto_secondary=True, verbosity='{0}'.format(prodyVerbosity))

        ags = prody.parsePDB(inputFn, alt='all', compressed=False,
                             biomol=True, extend_biomol=True)
        if isinstance(ags, prody.AtomGroup):
            ags =[ags] 

        self.pdbs = SetOfAtomStructs().create(self._getExtraPath())
        for i, ag in enumerate(ags):
            filename = self._getPath(splitext(basename(inputFn))[0] + '_atoms_{0}.pdb'.format(i))
            prody.writePDB(filename, ag)
            pdb = AtomStruct(filename)
            self.pdbs.append(pdb)

        # configure ProDy to restore secondary structure information and verbosity
        prody.confProDy(auto_secondary=old_secondary, verbosity='{0}'.format(old_verbosity))

    def createOutputStep(self):
        self._defineOutputs(outputStructures=self.pdbs)

    def _summary(self):
        if not hasattr(self, '_sum'):
            self._sum = CsvList()

        if not hasattr(self, 'outputStructures'):
            self._sum = CsvList()
            self._sum.append('Output structure not ready yet')
        else:
            if len(self._sum) == 0 or not self._sum[0].startswith('Extracted'): 
                self._sum = CsvList()
                num_structs = len(self.outputStructures)
                self._sum.append('Extracted *{0}* biomolecular assemblies'.format(num_structs))

                ags = prody.parsePDB([struct.getFileName() for struct in self.outputStructures])
                if num_structs == 1:
                    ags = [ags]
                     
                for i, ag in enumerate(ags):
                    self._sum.append('New structure {0} has *{1}* residues '
                                     'across *{2}* chains'.format(i+1, ag.numResidues(), 
                                                                 ag.numChains()))
        return self._sum

