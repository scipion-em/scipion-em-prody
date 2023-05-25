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
This module will provide ProDy ensemble tools.
"""
from collections import OrderedDict
import numpy as np

from pwem import *
from pwem.objects import AtomStruct, SetOfAtomStructs, SetOfSequences
from pwem.protocols import EMProtocol

from pyworkflow.utils import *
from pyworkflow.protocol.params import (PointerParam, MultiPointerParam,
                                        StringParam, IntParam, FloatParam,
                                        EnumParam, TextParam, NumericRangeParam,
                                        BooleanParam, LEVEL_ADVANCED)

import prody
import time

from prody2.objects import ProDyNpzEnsemble, TrajFrame
from prody2.protocols.protocol_atoms import (NOTHING, PWALIGN, CEALIGN,
                                             DEFAULT)  # residue mapping methods

STRUCTURE = 0
INDEX = 1

BEST_MATCH = 0
SAME_CHID = 1
SAME_POS = 2
CUSTOM = 3

class ProDyBuildPDBEnsemble(EMProtocol):
    """
    This protocol will use ProDy's buildPDBEnsemble method to align atomic structures
    """
    _label = 'buildPDBEnsemble'
    _possibleOutputs = {'outputAtomStructs': SetOfAtomStructs,
                        'outputNpz': ProDyNpzEnsemble,
                        'outAlignment': SetOfSequences}

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label='ProDy buildPDBEnsemble')

        form.addParam('inputType', EnumParam, choices=['structures', 'id for DALI search'],
                      default=STRUCTURE, important=True,
                      label="Type of input for building the ensemble",
                      help='The input can be a SetOfAtomStructs or an indexed to search DALI')

        form.addParam('structures', MultiPointerParam, label="Set of structures",
                      condition="inputType == %d" % STRUCTURE,
                      pointerClass='AtomStruct,SetOfAtomStructs', allowsNull=True,
                      help='The structures to be aligned must be atomic models.')

        form.addParam('id', StringParam, label="PDB ID and chain ID for DALI search",
                      condition="inputType == %d" % INDEX,
                      help='This ID should be a 5-character combination of a PDB ID and chain ID e.g., 3h5vA.')

        form.addParam('cutoff_len', StringParam, label="cutoff_len for filtering DALI results",
                      condition="inputType == %d" % INDEX, default='-1',
                      expertLevel=LEVEL_ADVANCED,
                      help='Select results with length of aligned residues < cutoff_len '
                      '(must be an integer up to the number of residues or a float between 0 and 1)')
        form.addParam('cutoff_rmsd', StringParam, label="cutoff_rmsd for filtering DALI results",
                      condition="inputType == %d" % INDEX, default='-1',
                      expertLevel=LEVEL_ADVANCED,
                      help='Select results with RMSD < cutoff_rmsd (must be a positive number)')
        form.addParam('cutoff_Z', StringParam, label="cutoff_Z for filtering DALI results",
                      condition="inputType == %d" % INDEX, default='-1',
                      expertLevel=LEVEL_ADVANCED,
                      help='Select results with Z score < cutoff_Z (must be a positive number)')
        form.addParam('cutoff_identity', StringParam, label="cutoff_identity for filtering DALI results",
                      condition="inputType == %d" % INDEX, default='-1',
                      expertLevel=LEVEL_ADVANCED,
                      help='Select results sequence identity > cutoff_identity '
                      '(must be an integer up to 100 or a float between 0 and 1).')

        form.addParam('refType', EnumParam, choices=['structure', 'index'], 
                      default=STRUCTURE, condition="inputType == %d" % STRUCTURE,
                      label="Reference structure selection type",
                      help='The reference structure can be a separate structure or indexed from the set')

        form.addParam('refStructure', PointerParam, label="Reference structure",
                      condition="refType == %d" % STRUCTURE,
                      pointerClass='AtomStruct', allowsNull=True,
                      help='Select an atomic model as the reference structure. '
                      'When using Dali, this is optional and is used for selecting atoms at the end.')

        form.addParam('refIndex', IntParam, label="Reference structure index", default=1,
                      condition="refType == %d and inputType != %d" % (INDEX, INDEX),
                      help='Select the index of the reference structure in the set, starting from 1. '
                      'When using Dali, this is optional and is used for selecting atoms at the end.')
        
        form.addParam('matchFunc', EnumParam, choices=['bestMatch', 'sameChid', 'sameChainPos', 'custom'], 
                      default=SAME_CHID, condition="inputType == %d" % STRUCTURE,
                      label="Chain matching function",
                      help='See http://prody.csb.pitt.edu/manual/release/v1.11_series.html for more details.\n')

        form.addParam('seqid', FloatParam, default=0.,
                      expertLevel=LEVEL_ADVANCED,
                      label="Sequence identity cutoff",
                      help='Alignment mapping with lower sequence identity will not be accepted.\n'
                           'This can be a number between 0 and 100 or a decimal between 0 and 1')

        form.addParam('overlap', FloatParam, default=0.,
                      expertLevel=LEVEL_ADVANCED,
                      label="Overlap cutoff",
                      help='Alignment mapping with lower sequence coverage will not be accepted.\n'
                           'This can be a number between 0 and 100 or a decimal between 0 and 1')

        form.addParam('selstr', StringParam, default="name CA",
                      expertLevel=LEVEL_ADVANCED,
                      label="Selection string",
                      help='Selection string for atoms to include in the ensemble.\n'
                           'It is recommended to use "protein" or "name CA" (default)')

        form.addParam('trim', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Whether to trim away dummy atoms",
                      help='If any structure lacks some atom in the reference '
                           'then it will be replaced by a dummy atom at the average '
                           'position in the ensemble. This option allows these to be trimmed away.')
        form.addParam('trimFraction', FloatParam, default=1,
                      expertLevel=LEVEL_ADVANCED, condition="trim == True",
                      label="Occupancy fraction to trim away dummy atoms",
                      help='This option controls how many dummy atoms are trimmed away '
                           'and should take a value between 0 and 1.\n'
                           'The resulting ensemble will contain atoms whose occupancies are greater '
                           'than or equal to this value.')
        
        group = form.addGroup('custom chain orders', condition='matchFunc == %d' % CUSTOM)
        
        group.addParam('chainOrders', TextParam, width=40, readOnly=True,
                       condition='matchFunc == %d' % CUSTOM,
                       label='Custom chain match list',
                       help='Defined order of chains from custom matching. \nManual modification will have no '
                            'effect, use the wizards to add / delete the entries')
        
        group.addParam('insertOrder', NumericRangeParam, default='1',
                       condition='matchFunc == %d' % CUSTOM,
                       label='Insert custom match order number',
                       help='Insert the chain order with the specified index into the match list.\n'
                            'The default (when empty) is the last position')
        
        group.addParam('customOrder', StringParam, default='',
                       condition='matchFunc == %d' % CUSTOM,
                       label='Custom match order to insert at the specified number',
                       help='Enter the desired chain order here.\n'
                            'The default (when empty) is the current chain order in the form. '
                            'The initial value is the order in the structure file.')
        
        group.addParam('label', StringParam, default='', readOnly=True,
                       condition='matchFunc == %d' % CUSTOM,
                       label='Label for item with the specified number for recovering custom match',
                       help='This cannot be changed by the user and is for display only.')

        group.addParam('recoverOrder', StringParam, default='1',
                       condition='matchFunc == %d' % CUSTOM,
                       label='Recover custom match order number',
                       help='Enter the desired chain order here.\n'
                            'Recover the chain order with the specified index from the match list.')

        form.addParam('mapping', EnumParam, choices=['Nothing',
                                                     'Biopython pwalign local sequence alignment',
                                                     'Combinatorial extension (CE) structural alignment',
                                                     'Auto (try pwalign then ce)'],
                      default=PWALIGN, condition="inputType == %d" % STRUCTURE,
                      expertLevel=LEVEL_ADVANCED,
                      label="Residue mapping function",
                      help='This method will be used for matching residues if the residue numbers and types aren\'t identical. \n'
                           'See http://prody.csb.pitt.edu/manual/reference/proteins/compare.html?highlight=mapchainontochain#prody.proteins.compare.mapChainOntoChain '
                           'for more details.')

        form.addParam('rmsd_reject', FloatParam, default=15.,
                      expertLevel=LEVEL_ADVANCED,
                      label="Rejection RMSD (A)",
                      help='Alignments with worse RMSDs than this will be rejected.')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # handle inputs
        if self.refType.get() == STRUCTURE:
            ref = prody.parsePDB(self.refStructure.get().getFileName(), alt='all')
        else:
            ref = self.refIndex.get() - 1 # convert from Scipion (sqlite) to ProDy (python) nomenclature

        if self.inputType.get() == STRUCTURE:
            pdbs = []
            for i, obj in enumerate(self.structures):
                if isinstance(obj.get(), AtomStruct):
                    pdbs.append(obj.get().getFileName())
                else:
                    pdbs.extend([tarStructure.getFileName() for tarStructure in obj.get()])

            if self.mapping.get() == DEFAULT:
                mappings = 'auto'
            elif self.mapping.get() == PWALIGN:
                mappings = 'pwalign'
            elif self.mapping.get() == CEALIGN:
                mappings = 'ce'
            else:
                mappings = False

        else:
            idstr = self.id.get()
            dali_rec = prody.searchDali(idstr[:4], idstr[4], timeout=10000)
            while dali_rec.isSuccess != True:
                dali_rec.fetch(timeout=1000)
                time.sleep(10)

            cutoff_len = eval(str(self.cutoff_len.get()))
            if cutoff_len == -1:
                cutoff_len = None

            cutoff_rmsd = eval(str(self.cutoff_rmsd.get()))
            if cutoff_rmsd == -1:
                cutoff_rmsd = None
                
            cutoff_Z = eval(str(self.cutoff_Z.get()))
            if cutoff_Z == -1:
                cutoff_Z = None

            cutoff_identity = eval(str(self.cutoff_identity.get()))
            if cutoff_identity == -1:
                cutoff_identity = None

            pdbs = dali_rec.filter(cutoff_len=cutoff_len, cutoff_rmsd=cutoff_rmsd,
                                   cutoff_Z=cutoff_Z, cutoff_identity=cutoff_identity)
            mappings = dali_rec.getMappings()

        self.tars = prody.parsePDB(pdbs, alt='all')

        if isinstance(self.tars, prody.Atomic):
            n_models = self.tars.numCoordsets()
            self.tars = []
            for i in range(n_models):
                self.tars.append(prody.parsePDB(pdbs, alt='all',
                                                model=i+1))

        # actual steps
        self._insertFunctionStep('alignStep', ref, mappings)
        self._insertFunctionStep('createOutputStep')

    def alignStep(self, ref, mappings):
        """This step includes alignment mapping and superposition"""

        # configure ProDy to automatically handle secondary structure information and verbosity
        old_secondary = prody.confProDy("auto_secondary")
        old_verbosity = prody.confProDy("verbosity")
        
        from pyworkflow import Config
        prodyVerbosity =  'none' if not Config.debugOn() else 'debug'
        prody.confProDy(auto_secondary=True, verbosity='{0}'.format(prodyVerbosity))

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
        
        atommaps = [] # output argument for collecting atommaps

        if mappings != 'auto':
            # from dali so don't use match func or ref
            ens = prody.buildPDBEnsemble([tar.select(self.selstr.get()) for tar in self.tars],
                                          seqid=self.seqid.get(),
                                          overlap=self.overlap.get(),
                                          mapping=mappings,
                                          atommaps=atommaps,
                                          rmsd_reject=self.rmsd_reject.get())

            # instead use them later for selection
            ens_ref = ens.getAtoms()
            amap = prody.alignChains(ens_ref, ref.select(self.selstr.get()),
                                     seqid=self.seqid.get(),
                                     overlap=self.overlap.get(),
                                     match_func=match_func,
                                     rmsd_reject=self.rmsd_reject.get())[0]
            ens.setAtoms(amap)

        else:
            if self.refType.get() == STRUCTURE:
                self.tars = [ref] + self.tars
                ref=0
                
            tars = [tar.select(self.selstr.get()).copy() for tar in self.tars]

            if self.matchFunc <= SAME_POS:
                self.labels = [tar.getTitle() for tar in tars]
            else:
                self.matchDic = self.createMatchDic("1")
                self.labels = list(self.matchDic.keys())

            for i, label in enumerate(self.labels):
                if label.endswith(" Selection 'name CA'"):
                    label = label.replace(" Selection 'name CA'", "")
                    if not label.endswith("_ca"):
                        label += "_ca"

                label = label.replace("_atoms", "")

                tars[i].setTitle(label)
                self.labels[i] = label

            self.labels = list(np.array(self.labels, dtype='<U20'))

            ens = prody.buildPDBEnsemble(tars,
                                         ref=ref,
                                         seqid=self.seqid.get(),
                                         overlap=self.overlap.get(),
                                         match_func=match_func,
                                         atommaps=atommaps,
                                         rmsd_reject=self.rmsd_reject.get())

        if self.trim.get():
            ens = prody.trimPDBEnsemble(ens, self.trimFraction.get())

        indices = ens.getIndices()

        msa = ens.getMSA()
        prody.writeMSA(self._getExtraPath('ensemble.fasta'), msa)
        
        pos_aligned = prody.alignByEnsemble(self.tars, ens)

        self.pdbs = SetOfAtomStructs().create(self._getExtraPath())
        for i, ag in enumerate(pos_aligned):
            amap = atommaps[i]
            if indices is not None:
                amap = amap[indices]
            
            amap.setTitle(amap.getTitle().split('[')[0])
            filename = self._getExtraPath('{:06d}_{:s}_amap.pdb'.format(i+1, ag.getTitle()))
            prody.writePDB(filename, amap)
            pdb = AtomStruct(filename)
            self.pdbs.append(pdb)

        prody.writePDB(self._getPath('ensemble.pdb'), ens)

        self.npzFileName = self._getPath('ensemble.ens.npz')
        prody.saveEnsemble(ens, self.npzFileName)

        self.npz = ProDyNpzEnsemble().create(self._getExtraPath())
        for j in range(ens.numConfs()):
            frame = TrajFrame((j+1, self.npzFileName), objLabel=ens.getLabels()[j])
            self.npz.append(frame)

        # configure ProDy to restore secondary structure information and verbosity
        prody.confProDy(auto_secondary=old_secondary, verbosity='{0}'.format(old_verbosity))

    def createOutputStep(self):
        
        outputSeqs = SetOfSequences().create(self._getExtraPath())
        outputSeqs.importFromFile(self._getExtraPath('ensemble.fasta'))

        self._defineOutputs(outputAtomStructs=self.pdbs,
                            outputNpz=self.npz,
                            outAlignment=outputSeqs)

    def createMatchDic(self, index):

        # configure ProDy to automatically handle secondary structure information and verbosity
        old_secondary = prody.confProDy("auto_secondary")
        old_verbosity = prody.confProDy("verbosity")
        
        from pyworkflow import Config
        prodyVerbosity =  'none' if not Config.debugOn() else 'debug'
        prody.confProDy(auto_secondary=False, verbosity='{0}'.format(prodyVerbosity))
        
        pdbs = []
        for i, obj in enumerate(self.structures):
            if isinstance(obj.get(), AtomStruct):
                pdbs.append(obj.get().getFileName())
            else:
                pdbs.extend([tarStructure.getFileName() for tarStructure in obj.get()])
        
        try:
            self.matchDic = eval(self.chainOrders.get())
            self.labels = list(self.matchDic.keys())
            self.orders = list(self.matchDic.values())
            
            # reinitialise to update with new keys
            # that are still ordered correctly
            self.matchDic = OrderedDict()
        except:
            self.matchDic = OrderedDict()
            self.labels = []
            self.orders = []

            tars = prody.parsePDB(pdbs, alt='all')

            for ag in tars:
                title = ag.getTitle()
                self.labels.append(title)
                self.orders.append(self.getInitialChainOrder(ag))
                
        self.labels = np.array(self.labels, dtype='<U20')
        self.orders = np.array(self.orders)

        inds = list(np.array(getListFromRangeString(index)) - 1)
        
        for idx in inds:
            if self.customOrder.get() != '':
                self.orders[idx] = self.customOrder.get()
        
        # configure ProDy to restore secondary structure information and verbosity
        prody.confProDy(auto_secondary=old_secondary, verbosity='{0}'.format(old_verbosity))

        self.matchDic.update(zip(self.labels, self.orders))
        return self.matchDic
    
    def getInitialChainOrder(self, ag):
        return ''.join([ch.getChid() for ch in ag.protein.getHierView().iterChains()])

    def _summary(self):
        if not hasattr(self, 'outputNpz'):
            sum = ['Output ensemble not ready yet']
        else:
            ens = self.outputNpz.loadEnsemble()
            sum = ['Ensemble built with *{0}* structures of *{1}* atoms'.format(
                   ens.numConfs(), ens.numAtoms())]
        return sum
    
    