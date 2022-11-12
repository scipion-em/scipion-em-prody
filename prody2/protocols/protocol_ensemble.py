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
from pyworkflow.protocol import params

from pwem import *
from pwem.objects import AtomStruct, Transform, String, EMFile
from pwem.objects import SetOfAtomStructs, SetOfSequences
from pwem.protocols import EMProtocol

from pyworkflow.utils import *
from pyworkflow.protocol.params import (PointerParam, StringParam, IntParam, FloatParam,
                                        EnumParam, LEVEL_ADVANCED)

import prody
import time

STRUCTURE = 0
INDEX = 1

BEST_MATCH = 0
SAME_CHID = 1
SAME_POS = 2

class ProDyBuildPDBEnsemble(EMProtocol):
    """
    This protocol will build a PDBEnsemble
    """
    _label = 'buildPDBEnsemble'

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

        form.addParam('structures', PointerParam, label="Set of structures",
                      condition="inputType == %d" % STRUCTURE,
                      pointerClass='SetOfAtomStructs', allowsNull=True,
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

        form.addParam('refIndex', IntParam, label="Reference structure index",
                      condition="refType == %d and inputType != %d" % (INDEX, INDEX),
                      help='Select the index of the reference structure in the set, starting from 1. '
                      'When using Dali, this is optional and is used for selecting atoms at the end.')

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

        form.addParam('matchFunc', EnumParam, choices=['bestMatch', 'sameChid', 'sameChainPos'], 
                      default=SAME_CHID, condition="inputType == %d" % STRUCTURE,
                      expertLevel=LEVEL_ADVANCED,
                      label="Chain matching function",
                      help='See http://prody.csb.pitt.edu/manual/release/v1.11_series.html for more details.\n'
                           'Custom match functions will be added soon.')

        form.addParam('selstr', StringParam, default="name CA",
                      expertLevel=LEVEL_ADVANCED,
                      label="Selection string",
                      help='Selection string for atoms to include in the ensemble.\n'
                           'It is recommended to use "protein" or "name CA" (default)')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # handle inputs
        if self.refType.get() == STRUCTURE:
            ref = prody.parsePDB(self.refStructure.get().getFileName(), alt='all')
        else:
            ref = self.refIndex.get() - 1 # convert from Scipion (sqlite) to Prody (python) nomenclature

        if self.inputType.get() == STRUCTURE:
            pdbs = [tarStructure.getFileName() for tarStructure in self.structures.get()]
            mappings = 'auto'
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

        match_func = eval('prody.'+['bestMatch', 'sameChid', 'sameChainPos'][self.matchFunc.get()])
        
        atommaps = [] # output argument for collecting atommaps

        if mappings != 'auto':
            # from dali so don't use match func or ref
            ens = prody.buildPDBEnsemble([tar.select(self.selstr.get()) for tar in self.tars],
                                        seqid=self.seqid.get(),
                                        overlap=self.overlap.get(),
                                        mapping=mappings,
                                        atommaps=atommaps)

            # instead use them later for selection
            ens_ref = ens.getAtoms()
            amap = prody.alignChains(ens_ref, ref.select(self.selstr.get()),
                                     seqid=self.seqid.get(),
                                     overlap=self.overlap.get(),
                                     match_func=match_func)[0]
            ens.setAtoms(amap)

        else:
            if self.refType.get() == STRUCTURE:
                self.tars = [ref] + self.tars
                ref=0
                
            ens = prody.buildPDBEnsemble([tar.select(self.selstr.get()) for tar in self.tars],
                                          ref=ref,
                                          seqid=self.seqid.get(),
                                          overlap=self.overlap.get(),
                                          match_func=match_func,
                                          atommaps=atommaps)

        ens = prody.trimPDBEnsemble(ens, 1.)
        indices = ens.getIndices()

        msa = ens.getMSA()
        prody.writeMSA(self._getExtraPath('ensemble.fasta'), msa)

        self.pdbs = SetOfAtomStructs().create(self._getExtraPath())
        for i, ag in enumerate(self.tars):
            amap = atommaps[i]
            if indices is not None:
                amap = amap[indices]
            
            amap.setTitle(amap.getTitle().split('[')[0])
            filename = self._getExtraPath('{:06d}_{:s}_amap.pdb'.format(i+1, ag.getTitle()))
            prody.writePDB(filename, amap)
            pdb = AtomStruct(filename)
            self.pdbs.append(pdb)

        self.pdbFileName = self._getPath('ensemble.pdb')
        prody.writePDB(self.pdbFileName, ens)

        self.npzFileName = self._getPath('ensemble.ens.npz')
        prody.saveEnsemble(ens, self.npzFileName)

        self.T = [T.getMatrix() for T in ens.getTransformations()]
        self.matrixFileName = self._getPath('transformation_%05d.txt')
        [prody.writeArray(self.matrixFileName % i, T) for i, T in enumerate(self.T)]

        # configure ProDy to restore secondary structure information and verbosity
        prody.confProDy(auto_secondary=old_secondary, verbosity='{0}'.format(old_verbosity))

    def createOutputStep(self):
        
        outputSeqs = SetOfSequences().create(self._getExtraPath())
        outputSeqs.importFromFile(self._getExtraPath('ensemble.fasta'))

        #outputNpz = EMFile(filename=self.npzFileName)
        #outputTrans = [Transform(matrix=T) for T in self.T]

        self._defineOutputs(outputStructures=self.pdbs,
                            #outputNpz=outputNpz,
                            outAlignment=outputSeqs#,
                            #outputTransformations=outputTrans
                            )
