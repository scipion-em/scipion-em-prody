# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     James Krieger (jamesmkrieger@gmail.com)
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

from pwem.objects import (AtomStruct, SetOfAtomStructs, SetOfSequences,
                          EMFile)
from pwem.protocols import EMProtocol

from pyworkflow.utils import logger, getListFromRangeString, redStr
from pyworkflow.protocol.params import (PointerParam, MultiPointerParam,
                                        StringParam, IntParam, FloatParam,
                                        EnumParam, TextParam, NumericRangeParam,
                                        BooleanParam, LEVEL_ADVANCED)
from pyworkflow.object import Float

import prody
from prody2.objects import ProDyNpzEnsemble, TrajFrame
from prody2.constants import (NOTHING, PWALIGN, CEALIGN, DEFAULT,  # residue mapping methods
                              BEST_MATCH, SAME_CHID, SAME_POS, CUSTOM, # chain matching
                              ENSEMBLE_WEIGHTS)
from prody2 import parseMatchDict

import time

STRUCTURE = 0
INDEX = 1

BLAST = 0
DALI = 1

ENS_FILENAME = 'ensemble.dcd'

from prody2.objects import HAVE_CHEM
if HAVE_CHEM:
    from prody2.objects import DcdMDSystem

class ProDyBuildPDBEnsemble(EMProtocol):
    """
    Builds a structural ensemble by aligning multiple atomic models
    using ProDy's buildPDBEnsemble framework.

    AI Generated:

    Build PDB Ensemble (ProDyBuildPDBEnsemble) — User Manual
        Overview

        The Build PDB Ensemble protocol aligns a collection of atomic
        structures into a common structural frame and generates a
        consistent ensemble representation. It is intended for the
        comparative analysis of related macromolecular conformations,
        structural variability studies, or preparation of structural
        ensembles for downstream flexibility and dynamics analysis.

        From a biological perspective, this protocol is useful when
        several structures correspond to homologous proteins, distinct
        functional states, alternative experimental conditions, or
        multiple structures identified through structural similarity
        searches such as DALI.

        Inputs and General Workflow

        The protocol accepts two input modes:

        1. A set of user-provided atomic structures.
        2. A PDB identifier and chain identifier used for automatic
           DALI structural similarity search.

        In both cases, the protocol collects the target structures,
        parses atomic coordinates with ProDy, and aligns them into
        a common ensemble.

        When a set of structures is used, one structure can be chosen
        as the reference either explicitly or by index within the set.

        When DALI search is used, homologous structures are retrieved,
        filtered according to structural similarity criteria, and
        then assembled into the ensemble automatically.

        Reference Structure Selection

        The reference structure defines the coordinate system used for
        the final ensemble.

        Two reference modes are available:

        - Reference by structure:
          a separate atomic model is used as the alignment target.

        - Reference by index:
          one of the input structures is selected as reference.

        The reference can optionally be removed after alignment. This
        is useful when the reference is only needed to define the
        coordinate frame but should not be part of the final ensemble.

        Structural Matching and Alignment

        When input structures are supplied directly, the protocol
        performs structural matching before superposition.

        Several chain matching strategies are available:

        - bestMatch:
          chooses the best chain correspondence automatically.

        - sameChid:
          matches chains using identical chain identifiers.

        - sameChainPos:
          matches chains by positional order.

        - custom:
          allows user-defined chain matching order.

        For custom matching, the protocol builds a chain matching
        dictionary that specifies how chains in different structures
        correspond to each other.

        This is particularly useful for:

        - multimeric assemblies,
        - structures with inconsistent chain naming,
        - manually curated biological comparisons.

        Residue Mapping

        If residue numbering or residue identity differs between
        structures, the protocol can apply residue mapping strategies.

        Available mapping options include:

        - pairwise sequence alignment,
        - combinatorial extension structural alignment,
        - automatic mapping,
        - no mapping.

        These mapping procedures help preserve biologically meaningful
        residue correspondence across structures.

        DALI-Based Structural Search

        When a PDB identifier is used instead of explicit structures,
        the protocol queries the DALI server for structurally similar
        entries.

        Retrieved structures may be filtered using:

        - alignment length cutoff,
        - RMSD cutoff,
        - Z-score cutoff,
        - sequence identity cutoff.

        This mode is especially useful when the goal is to construct
        an ensemble of structural homologs automatically from known
        structural databases.

        Ensemble Construction

        Once target structures and mappings are defined, ProDy builds
        the structural ensemble.

        During this step the protocol:

        - aligns all selected structures,
        - stores atom mappings,
        - tracks unmapped structures,
        - preserves structure weights,
        - supports coordinate degeneracy if multiple coordinate sets
          are present.

        If a structure contains multiple coordinate sets and
        degeneracy is disabled, all coordinate sets may contribute
        independently to the final ensemble.

        Label Handling and Ensemble Consistency

        After alignment, the protocol normalizes structure labels.

        It removes unnecessary suffixes and resolves duplicated names
        by assigning unique indices.

        This guarantees that each ensemble member can be identified
        uniquely during downstream analysis.

        Trimming Dummy Atoms

        Missing atoms in some structures may generate dummy atoms in
        the ensemble.

        The protocol optionally trims these positions according to
        an occupancy threshold.

        Biologically, trimming is important because it removes poorly
        conserved or absent regions that may otherwise introduce noise
        into the final structural comparison.

        Reordering by Custom Chain Dictionary

        If custom chain matching is used, the final ensemble may also
        be reordered according to the custom chain definition.

        This preserves biologically meaningful ordering of ensemble
        members when the input order is not the desired final order.

        Optional Output Files

        The protocol can generate several optional outputs in addition
        to the main ensemble.

        Multiple aligned PDB files

            Each aligned structure can be written as an individual PDB
            file after applying the ensemble transformations.

        DCD trajectory

            The aligned ensemble can be exported as a trajectory file.
            This is useful for visualization, conformational analysis,
            or molecular dynamics–style downstream workflows.

        Sequence alignment

            A multiple sequence alignment is extracted from the final
            ensemble and written as a FASTA file.

        Main Outputs

        The protocol produces:

        outputNpz

            A ProDy ensemble stored in NPZ format. Each frame contains
            the aligned coordinates and associated structure weight.

        outAlignment

            A multiple sequence alignment derived from the final
            structural ensemble.

        outputStructures (optional)

            Individual aligned PDB files for each structure.

        outputTrajectory (optional)

            A DCD trajectory representation of the ensemble.

        Custom Chain Matching Utilities

        The protocol contains internal utilities that support custom
        matching workflows.

        createMatchDic()

            Builds the chain matching dictionary used for custom
            alignment strategies.

        getInitialChainOrder()

            Extracts the original chain order from an atomic model.

        These utilities are especially important when aligning
        assemblies whose biological chain organization must be
        preserved explicitly.

        Summary of Internal Steps

        _insertAllSteps()

            Defines protocol execution.

        alignStep()

            Main computational stage.

            - Parses structures
            - Handles references
            - Performs DALI search if needed
            - Builds mappings
            - Aligns structures
            - Creates ensemble
            - Writes ensemble files

        createOutputStep()

            Registers all generated outputs into Scipion.

        _summary()

            Reports the number of aligned structures and atoms.

        _setWeights()

            Assigns stored ensemble weights to output objects.

        Practical Biological Interpretation

        This protocol is particularly useful when the scientific goal
        is to compare multiple related structures in a unified
        coordinate system.

        Common applications include:

        - comparing conformational states,
        - studying structural heterogeneity,
        - preparing homologous ensembles,
        - generating input for normal mode analysis,
        - producing aligned datasets for visualization.

        For best biological results:

        - choose a biologically representative reference,
        - use conservative residue mapping when sequence divergence
          is high,
        - trim poorly conserved positions,
        - verify chain matching carefully for multimeric assemblies.

        Final Perspective

        Build PDB Ensemble is not simply a structural superposition
        tool.

        It provides a biologically meaningful way to transform a
        collection of related atomic models into a coherent ensemble
        that can be analyzed as a unified structural population.

        Proper reference selection, careful chain matching, and
        appropriate trimming are the key factors that determine
        whether the final ensemble accurately reflects biologically
        relevant structural variability.
    """
    _label = 'buildPDBEnsemble'
    _possibleOutputs = {'outputStructures': SetOfAtomStructs,
                        'outputNpz': ProDyNpzEnsemble,
                        'outAlignment': SetOfSequences}

    weights = []

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label='ProDy buildPDBEnsemble')

        form.addParam('inputType', EnumParam, choices=['structures', 'id for search'],
                      default=STRUCTURE, important=True,
                      label="Type of input for building the ensemble",
                      help='The input can be a SetOfAtomStructs or an ID to search the PDB')

        inputTypeCheck = "inputType == %d"
        form.addParam('structures', MultiPointerParam, label="Set of structures",
                      condition=inputTypeCheck % STRUCTURE,
                      pointerClass='AtomStruct,SetOfAtomStructs', allowsNull=True,
                      help='The structures to be aligned must be atomic models.')

        form.addParam('uniteChains', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Unite chains in mmCIF segments",
                      help='Elect whether to unite chains in mmCIF segments for each structure like ChimeraX. '
                            'Default is **False**, which means the smaller unit IDs are used for chains like PyMOL.')

        form.addParam('incrementTers', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Increment serial numbers for TER lines when writing PDB files",
                      help='Elect whether to increment serial numbers when writing PDB files. '
                            'Default is **True** as this is the normal behaviour but the **False** option is needed for Plumed.')

        form.addParam('id', StringParam, label="PDB ID and chain ID for DALI search",
                      condition=inputTypeCheck % INDEX,
                      help='This ID should be a 5-character combination of a PDB ID and chain ID e.g., 3h5vA.')

        form.addParam('degeneracy', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Take only first conformation from each structure/set",
                      help='Elect whether only the active coordinate set (degeneracy=**True**) or all the coordinate sets '
                           '(degeneracy=**False**) of each structure should be added to the ensemble. Default is **False**.')

        form.addParam('lenCutoff', StringParam, label="length cutoff for filtering DALI results",
                      condition=inputTypeCheck % INDEX, default='-1',
                      expertLevel=LEVEL_ADVANCED,
                      help='Filter out results with length of aligned residues < length cutoff'
                      '(must be an integer up to the number of residues or a float between 0 and 1)')
        form.addParam('rmsdCutoff', StringParam, label="rmsdCutoff for filtering DALI results",
                      condition=inputTypeCheck % INDEX, default='-1',
                      expertLevel=LEVEL_ADVANCED,
                      help='Filter out results with RMSD > RMSD Cutoff (must be a positive number)')
        form.addParam('zCutoff', StringParam, label="Z score cutoff for filtering DALI results",
                      condition=inputTypeCheck % INDEX, default='-1',
                      expertLevel=LEVEL_ADVANCED,
                      help='Select results with Z score < Z score cutoff (must be a positive number)')
        form.addParam('idCutoff', StringParam, label="Sequence ID cutoff for filtering DALI results",
                      condition=inputTypeCheck % INDEX, default='-1',
                      expertLevel=LEVEL_ADVANCED,
                      help='Filter out results with sequence identity < sequence ID cutoff '
                      '(must be an integer up to 100 or a float between 0 and 1).')

        form.addParam('refType', EnumParam, choices=['structure', 'index'], 
                      default=INDEX, condition=inputTypeCheck % STRUCTURE,
                      label="Reference structure type",
                      help='The reference structure can be a separate structure or indexed from the set')

        form.addParam('refStructure', PointerParam, label="Reference structure",
                      condition="refType == %d" % STRUCTURE,
                      pointerClass='AtomStruct', allowsNull=True,
                      help='Select an atomic model as the reference structure. '
                      'When using Dali, this is optional and is used for selecting atoms at the end.')

        form.addParam('delReference', BooleanParam, default=False,
                      label="Whether to delete the reference from the ensemble",
                      help='This could be useful if you just want to use the reference for alignment.')

        form.addParam('refIndex', IntParam, label="Reference structure index", default=1,
                      condition="refType == %d and inputType != %d" % (INDEX, INDEX),
                      help='Select the index of the reference structure in the set, starting from 1. '
                      'When using Dali, this is optional and is used for selecting atoms at the end.')
        
        form.addParam('matchFunc', EnumParam, choices=['bestMatch', 'sameChid', 'sameChainPos', 'custom'], 
                      default=SAME_CHID, condition=inputTypeCheck % STRUCTURE,
                      label="Chain matching function",
                      help='See http://http://www.bahargroup.org/prody/manual/release/v1.11_series.html for more details.\n')

        form.addParam('seqid', FloatParam, default=0.,
                      expertLevel=LEVEL_ADVANCED,
                      label="Sequence identity cutoff",
                      help='Alignment mapping with lower percent sequence identity will not be accepted.\n'
                           'This can be a number between 0 and 100')

        form.addParam('overlap', FloatParam, default=0.,
                      expertLevel=LEVEL_ADVANCED,
                      label="Overlap cutoff",
                      help='Alignment mapping with lower percent sequence coverage will not be accepted.\n'
                           'This can be a number between 0 and 100')

        form.addParam('rmsdReject', FloatParam, default=15.,
                      expertLevel=LEVEL_ADVANCED,
                      label="Rejection RMSD (A)",
                      help='Alignments with worse RMSDs than this will be rejected.')

        form.addParam('selstr', StringParam, default="name CA",
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
        
        matchFuncCheck = 'matchFunc == %d'
        group = form.addGroup('Custom chain orders', condition=matchFuncCheck % CUSTOM)
        
        group.addParam('chainOrders', TextParam, width=50,
                       condition=matchFuncCheck % CUSTOM, default="",
                       label='Custom chain match dictionary',
                       help='Defined order of chains from custom matching.')
        
        group.addParam('insertOrder', NumericRangeParam, default='1',
                       condition=matchFuncCheck % CUSTOM,
                       label='Insert custom match order number',
                       help='Insert the chain order with the specified index into the match list.\n'
                            'The default (when empty) is the last position.\n'
                            'If the reference is a structure then that is structure 1.')
        
        group.addParam('customOrder', StringParam, default='',
                       condition=matchFuncCheck % CUSTOM,
                       label='Custom match order to insert at the specified number',
                       help='Enter the desired chain order here.\n'
                            'The default (when empty) is the current chain order in the form. '
                            'The initial value is the order in the structure file.')
        
        group.addParam('label', StringParam, default='',
                       condition=matchFuncCheck % CUSTOM,
                       label='Label for item with the specified number for recovering custom match',
                       help='This cannot be changed by the user and is for display only.')

        group.addParam('recoverOrder', StringParam, default='1',
                       condition=matchFuncCheck % CUSTOM,
                       label='Recover custom match order number',
                       help='Enter the desired chain order index here.\n'
                            'Recover the chain order with the specified index from the match list.')

        form.addParam('mapping', EnumParam, choices=['Nothing',
                                                     'Biopython pwalign local sequence alignment',
                                                     'Combinatorial extension (CE) structural alignment',
                                                     'Auto (try pwalign then ce)'],
                      default=PWALIGN, condition=inputTypeCheck % STRUCTURE,
                      expertLevel=LEVEL_ADVANCED,
                      label="Residue mapping function",
                      help='This method will be used for matching residues if the residue numbers and types aren\'t identical. \n'
                           'See http://http://www.bahargroup.org/prody/manual/reference/proteins/compare.html?highlight=mapchainontochain#prody.proteins.compare.mapChainOntoChain '
                           'for more details.')

        form.addParam('writeDCDFile', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Whether to write DCD trajectory file",
                      help='This will be registered as output too')
        
        form.addParam('writePDBFiles', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      condition=(HAVE_CHEM==True and 'degeneracy'),
                      label="Whether to write many PDB files",
                      help='These will be registered as output too')
        
        form.addParam('doReorder', BooleanParam, default=False,
                      condition=matchFuncCheck % CUSTOM,
                      label="Whether to reorder ensemble by custom match dict",
                      help='Otherwise the order matches the input')

        form.addParam('keepAlignment', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Keep alignment",
                      help="The alternative is to realign the structures over the whole structure")

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # actual steps
        self._insertFunctionStep('alignStep')
        self._insertFunctionStep('createOutputStep')

    def alignStep(self):
        """This step includes alignment mapping and superposition"""
        degeneracy = self.degeneracy.get()

        if self.keepAlignment.get():
            superpose = False
        else:
            superpose = 'iter'

        # handle reference
        self.weights = []
        if self.refType.get() == STRUCTURE:
            ref = prody.parsePDB(self.refStructure.get().getFileName(), alt='all',
                                 unite_chains=self.uniteChains.get())
            if ref.numCoordsets() > 1 and not degeneracy:
                self.weights.extend([self.refStructure.get().getAttributeValue(ENSEMBLE_WEIGHTS,
                                                                               defaultValue=1)] * ref.numCoordsets())
            else:
                self.weights.append(self.refStructure.get().getAttributeValue(ENSEMBLE_WEIGHTS,
                                                                              defaultValue=1))
        else:
            ref = self.refIndex.get() - 1 # convert from Scipion (sqlite) to ProDy (python) nomenclature

        # handle other inputs
        if self.inputType.get() == STRUCTURE:
            self.pdbs = []
            structureObjects = []
            for i, obj in enumerate(self.structures):
                if isinstance(obj.get(), AtomStruct):
                    structureObjects.append(obj.get())
                    self.pdbs.append(obj.get().getFileName())
                else:
                    structureObjects.extend([tarStructure for tarStructure in obj.get()])
                    self.pdbs.extend([tarStructure.getFileName() for tarStructure in obj.get()])

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
            daliRec = prody.searchDali(idstr[:4], idstr[4], timeout=10000)
            while daliRec.isSuccess != True:
                daliRec.fetch(timeout=1000)
                time.sleep(10)

            lenCutoff = eval(str(self.lenCutoff.get()))
            if lenCutoff == -1:
                lenCutoff = None

            rmsdCutoff = eval(str(self.rmsdCutoff.get()))
            if rmsdCutoff == -1:
                rmsdCutoff = None

            zCutoff = eval(str(self.zCutoff.get()))
            if zCutoff == -1:
                zCutoff = None

            idCutoff = eval(str(self.idCutoff.get()))
            if idCutoff == -1:
                idCutoff = None

            self.pdbs = daliRec.filter(cutoff_len=lenCutoff, cutoff_rmsd=rmsdCutoff,
                                       cutoff_Z=zCutoff, cutoff_identity=idCutoff,
                                       stringency=True)
            mappings = daliRec.getMappings()

            if idstr not in self.pdbs:
                self.pdbs.insert(0, idstr)

        if not hasattr(self, "tars"):
            self.tars = prody.parsePDB(self.pdbs, alt='all',
                                       unite_chains=self.uniteChains.get())
            if isinstance(self.tars, prody.Atomic):
                self.tars = [self.tars]

            for i, tar in enumerate(self.tars):
                self.weights.extend([structureObjects[i].getAttributeValue(ENSEMBLE_WEIGHTS,
                                                                           defaultValue=1)] * tar.numCoordsets())

        atommaps = [] # output argument for collecting atommaps
        unmapped = []

        if self.inputType.get() != STRUCTURE:
            ens = prody.buildPDBEnsemble([tar.select(self.selstr.get()) for tar in self.tars],
                                          seqid=self.seqid.get(),
                                          overlap=self.overlap.get(),
                                          mapping=mappings,
                                          atommaps=atommaps,
                                          unmapped=unmapped,
                                          rmsd_reject=self.rmsdReject.get(),
                                          superpose=superpose)
            self.weights = list(np.ones(ens.numConfs()))
        else:
            if self.matchFunc.get() == BEST_MATCH:
                matchFunc = prody.bestMatch
                logger.info('\nUsing bestMatch\n')
            elif self.matchFunc.get() == SAME_CHID:
                matchFunc = prody.sameChid
                logger.info('\nUsing sameChid\n')
            elif self.matchFunc.get() == SAME_POS:
                matchFunc = prody.sameChainPos
                logger.info('\nUsing sameChainPos\n')

            if self.refType.get() == STRUCTURE:
                if self.matchFunc.get() < CUSTOM:
                    self.tars = [ref] + self.tars
                ref=0

            if self.matchFunc <= SAME_POS:
                tars = [tar.select(self.selstr.get()).copy() for tar in self.tars]
                self.labels = [tar.getTitle() for tar in tars]
            else:
                self.matchDic = self.createMatchDic(self.insertOrder.get())
                self.labels = list(self.matchDic.keys())
                self.orders = list(self.matchDic.values())

                if isinstance(self.labels[0], tuple):
                    self.labels = [label[1] for label in self.labels]

                self.matchDic = OrderedDict()
                self.matchDic.update(zip(self.labels, self.orders))

                logger.info('\nUsing user-defined match function based on \n{0}\n'.format(self.matchDic))
                matchFunc = lambda chain1, chain2: prody.userDefined(chain1, chain2, self.matchDic)

                tars = [tar.select(self.selstr.get()).copy() for tar in self.tars]

            if len(tars) != len(self.labels):
                logger.warn(redStr('labels e.g. from matchDic ({0}) do not match '
                            'target structures ({1})'.format(len(self.labels), len(tars))))
                
            matchDictLabels = self.labels
                
            titles = [tar.getTitle() for tar in tars]
            for i, title in enumerate(titles):
                title = title.replace(" Selection 'name CA'", "")
                title = title.replace("_atoms", "")
                tars[i].setTitle(title)

            ens = prody.buildPDBEnsemble(tars,
                                         ref=ref,
                                         seqid=self.seqid.get(),
                                         overlap=self.overlap.get(),
                                         match_func=matchFunc,
                                         atommaps=atommaps,
                                         unmapped=unmapped,
                                         rmsd_reject=self.rmsdReject.get(),
                                         degeneracy=self.degeneracy.get(),
                                         mapping=mappings,
                                         superpose=superpose)
            
            if self.delReference.get():
                ens.delCoordset(ref)
                self.tars.pop(ref)

        logger.info('\nUnmapped structures: {0}\n'.format(unmapped))

        self.labels = ens.getLabels()
        _, idx, inv, c = np.unique(self.labels, return_index=True,
                                   return_inverse=True, return_counts=True)

        for i, label in enumerate(self.labels):
            if label.endswith('_ca'):
                self.labels[i] = label[:-3]

            if i in idx:
                j = 0
            else:
                j += 1

            if c[inv][i] > 1:
                self.labels[i] = self.labels[i] + '_' + str(j)

        ens._labels = self.labels

        if self.trim.get():
            ens = prody.trimPDBEnsemble(ens, self.trimFraction.get())

        if self.doReorder.get():
            newIndices = [self.labels.index(label) for label in matchDictLabels 
                          if label in self.labels]
            ens = ens[newIndices]

        if self.writePDBFiles.get():
            indices = ens.getIndices()
            amapTitles = [amap.getAtomGroup().getTitle() for amap in atommaps]

            tars = [tar for tar in tars if tar.getTitle() in amapTitles]

            aligned = prody.alignByEnsemble(tars, ens)
            self.pdbs = SetOfAtomStructs().create(self._getExtraPath())
            for i, ag in enumerate(aligned):
                amap = atommaps[i]
                if indices is not None:
                    amap = amap[indices]

                ag_title = ag.getTitle().replace(' ', '_').replace("'", "")
                
                amap.setTitle(amap.getTitle().split('[')[0])
                filename = self._getExtraPath(
                    '{:06d}_{:s}_amap.pdb'.format(i+1, ag_title))
                prody.writePDB(filename, amap,
                               increment_ter=self.incrementTers.get())
                pdb = AtomStruct(filename)
                setattr(pdb, ENSEMBLE_WEIGHTS, Float(self.weights[i]))
                self.pdbs.append(pdb)

        prody.writePDB(self._getPath('ensemble.pdb'), ens,
                       increment_ter=self.incrementTers.get())

        self.npzFileName = self._getPath('ensemble.ens.npz')
        prody.saveEnsemble(ens, self.npzFileName)

        self.npz = ProDyNpzEnsemble().create(self._getExtraPath())
        for j in range(ens.numConfs()):
            frame = TrajFrame((j+1, self.npzFileName), objLabel=ens.getLabels()[j])
            setattr(frame, ENSEMBLE_WEIGHTS, Float(self.weights[j]))
            self.npz.append(frame)

        if self.writeDCDFile.get():
            self.pdbFilename = self._getPath('refStructure.pdb')
            prody.writeDCD(self._getPath(ENS_FILENAME), ens)
            prody.writePDB(self.pdbFilename, ens.getAtoms(),
                           increment_ter=self.incrementTers.get())

    def createOutputStep(self):
        outputs = {"outputNpz": self.npz}
        
        if self.writeDCDFile.get():
            if HAVE_CHEM:
                outMDSystem = DcdMDSystem(filename=self.pdbFilename)
                outMDSystem.setTopologyFile(self.pdbFilename)
                outMDSystem.setTrajectoryFile(self._getPath(ENS_FILENAME))
                outputs["outputTrajectory"] = outMDSystem
            else:
                outEMFile = EMFile(filename=self._getPath(ENS_FILENAME))
                outputs["outputTrajectory"] = outEMFile

        if self.writePDBFiles.get():
            outputs["outputStructures"] = self.pdbs

        self._defineOutputs(**outputs)

    def createMatchDic(self, index, label=""):     
        parseMatchDict(self)
        self.orders = list(self.matchDic.values())

        # reinitialise to update with new keys
        # that are still ordered correctly
        self.matchDic = OrderedDict()

        inds = [item-1 for item in getListFromRangeString(index)]

        if len(self.labels) == 0:
            if not hasattr(self, 'tars'):
                if self.refType.get() == STRUCTURE:
                    structures = [self.refStructure] + self.structures
                else:
                    structures = self.structures

                pdbs = []
                structureObjects = []
                for _, obj in enumerate(structures):
                    if isinstance(obj.get(), AtomStruct):
                        pdbs.append(obj.get().getFileName())
                        structureObjects.append(obj.get())
                    else:
                        pdbs.extend([tarStructure.getFileName() for tarStructure in obj.get()])
                        structureObjects.extend([tarStructure for tarStructure in obj.get()])

                self.tars = prody.parsePDB(pdbs, alt='all',
                                           unite_chains=self.uniteChains.get())

                if isinstance(self.tars, prody.Atomic):
                    self.tars = [self.tars]

                for i, tar in enumerate(self.tars):
                    self.weights.extend([structureObjects[i].getAttributeValue(ENSEMBLE_WEIGHTS, defaultValue=1)] * tar.numCoordsets())
                    if tar.numCoordsets() > 1:
                        self.haveMultipleCoordsets = True
            
            titles = [ag.getTitle() for ag in self.tars]
            _, counts = np.unique(np.array(titles), return_counts=True)

            for idx, ag in enumerate(self.tars):
                if (idx in inds or counts[idx] > 1) and label!="":
                    if len(inds) == 1:
                        ag.setTitle(label)
                    else:
                        ag.setTitle(label + str(inds.index(idx)))

                title = ag.getTitle()
                self.labels.append(title)
                self.orders.append(self.getInitialChainOrder(ag))

        self.orders = np.array(self.orders)

        if not isinstance(self.labels[0], tuple):
            self.labels = [(i+1, label) for i, label in enumerate(self.labels)]
        
        for idx in inds:
            if self.customOrder.get() != '':
                self.orders[idx] = self.customOrder.get()

        self.matchDic.update(zip(list(self.labels), list(self.orders)))
        return self.matchDic
    
    def getInitialChainOrder(self, ag):
        return ''.join([ch.getChid() for ch in ag.protein.getHierView().iterChains()])

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
    
    def _setWeights(self, item, row=None):
        weight = Float(self.weights[item.getObjId()-1])
        setattr(item, ENSEMBLE_WEIGHTS, weight)
