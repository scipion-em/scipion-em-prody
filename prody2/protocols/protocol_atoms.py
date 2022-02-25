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
from pyworkflow.protocol import params

from pwem import *
from pwem.objects import AtomStruct, Transform, String
from pwem.protocols import EMProtocol

from pyworkflow.utils import *
from pyworkflow.protocol.params import (PointerParam, StringParam, FloatParam,
                                        EnumParam, LEVEL_ADVANCED)

import prody

BEST_MATCH = 0
SAME_CHID = 1

class ProDySelect(EMProtocol):
    """
    This protocol will perform atom selection
    """
    _label = 'Atom Selection'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label='ProDy Select')

        form.addParam('inputStructure', PointerParam, label="Input structure",
                      important=True,
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
        inputFn = self.inputStructure.get().getFileName()
        self._insertFunctionStep('selectionStep', inputFn)
        self._insertFunctionStep('createOutputStep')

    def selectionStep(self, inputFn):
        ag = prody.parsePDB(inputFn, alt='all')
        selection = ag.select(str(self.selection))

        self.pdbFileName = self._getPath('atoms.pdb')
        prody.writePDB(self.pdbFileName, selection)

    def createOutputStep(self):
        outputPdb = AtomStruct()
        outputPdb.setFileName(self.pdbFileName)
        self._defineOutputs(outputStructure=outputPdb)   


class ProDyAlign(EMProtocol):
    """
    This protocol will perform atomic structure mapping and superposition
    """
    _label = 'Atom Alignment'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
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

        form.addParam('seqid', FloatParam, default=1.,
                      expertLevel=LEVEL_ADVANCED,
                      label="Sequence Identity Cut-off (%)",
                      help='Alignment mapping with lower sequence identity will not be accepted.')

        form.addParam('overlap', FloatParam, default=1.,
                      expertLevel=LEVEL_ADVANCED,
                      label="Overlap Cut-off (%)",
                      help='Alignment mapping with lower sequence coverage will not be accepted.\n'
                           'This can be a number between 0 and 100 or a decimal between 0 and 1') 

        form.addParam('matchFunc', EnumParam, choices=['bestMatch', 'sameChid'], 
                      default=BEST_MATCH,
                      expertLevel=LEVEL_ADVANCED,
                      label="Chain matching function",
                      help='See http://prody.csb.pitt.edu/manual/release/v1.11_series.html for more details.\n'
                           'Custom match functions will be added soon.')    

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        mobFn = self.mobStructure.get().getFileName()
        tarFn = self.tarStructure.get().getFileName()
        self._insertFunctionStep('alignStep', mobFn, tarFn)
        self._insertFunctionStep('createOutputStep')

    def alignStep(self, mobFn, tarFn):
        """This step includes alignment mapping and superposition"""

        mob = prody.parsePDB(mobFn, alt='all')
        tar = prody.parsePDB(tarFn, alt='all')

        if self.matchFunc.get() == BEST_MATCH:
            match_func = prody.bestMatch
        else:
            match_func = prody.sameChid

        mob_amap = prody.alignChains(mob, tar,
                                     seqid=self.seqid.get(),
                                     overlap=self.overlap.get(),
                                     match_func=match_func)[0]
        mob_sel = mob.select(mob_amap.getSelstr())

        tar_amap = prody.alignChains(tar, mob_sel,
                                     seqid=self.seqid.get(),
                                     overlap=self.overlap.get(),
                                     match_func=match_func)[0]
        tar_sel = tar.select(tar_amap.getSelstr())

        alg, self.T = prody.superpose(mob_sel, tar_sel)

        self.pdbFileNameMob = self._getPath('mobile.pdb')
        prody.writePDB(self.pdbFileNameMob, alg)

        self.pdbFileNameTar = self._getPath('target.pdb')
        prody.writePDB(self.pdbFileNameTar, tar_sel)

        self.matrixFileName = self._getPath('transformation.txt')
        prody.writeArray(self.matrixFileName, self.T.getMatrix())

    def createOutputStep(self):
        outputPdbMob = AtomStruct()
        outputPdbMob.setFileName(self.pdbFileNameMob)

        outputPdbTar = AtomStruct()
        outputPdbTar.setFileName(self.pdbFileNameTar)

        outputTrans = Transform()
        outputTrans.setMatrix(self.T.getMatrix())

        self._defineOutputs(outputStructureMob=outputPdbMob,
                            outputStructureTar=outputPdbTar,
                            outputTransformation=outputTrans)

