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

import numpy as np
from os.path import split

from pwem.objects import SetOfNormalModes, SetOfPrincipalComponents
from pwem.protocols import ProtImportPdb, ProtImportSetOfAtomStructs, exists
from pwem.tests.workflows import TestWorkflow
from pyworkflow.tests import setupTestProject

from prody2.protocols import (ProDySelect, ProDyBuildPDBEnsemble,
                              ProDyImportEnsemble, ProDyPCA, ProDyCompare,
                              ProDyProject, ProDyANM, ProDyImportModes, ProDyEdit,
                              ProDyMeasure)

from prody2.protocols.protocol_ensemble import INDEX, STRUCTURE
from prody2.protocols.protocol_edit import NMA_SLICE
from prody2.protocols.protocol_project import ONE, TWO, THREE
from prody2.protocols.protocol_import import SCIPION

import prody

def importSelect(cls):
    cls.protSel1 = cls.newProtocol(ProDySelect, selection="chain B D and name CA")
    cls.protSel1.inputPdbData.set(cls.protSel1.IMPORT_FROM_ID)
    cls.protSel1.pdbId.set('3p3w')
    cls.protSel1.setObjLabel('Sel 3p3w_BD_ca') # intermediate
    cls.launchProtocol(cls.protSel1)

    cls.protSel2 = cls.newProtocol(ProDySelect, selection="chain A C and name CA")
    cls.protSel2.inputPdbData.set(cls.protSel1.IMPORT_FROM_ID)
    cls.protSel2.pdbId.set('6fpj')
    cls.protSel2.setObjLabel('Sel 6fpj_AC_ca') # displaced
    cls.launchProtocol(cls.protSel2)

    cls.protSel3 = cls.newProtocol(ProDySelect, selection="name CA")
    cls.protSel3.inputPdbData.set(cls.protSel1.IMPORT_FROM_ID)
    cls.protSel3.pdbId.set('6flr')
    cls.protSel3.setObjLabel('Sel 6flr_ca') # open
    cls.launchProtocol(cls.protSel3)

    cls.protSel4 = cls.newProtocol(ProDySelect, selection="chain C D and name CA")
    cls.protSel4.inputPdbData.set(cls.protSel1.IMPORT_FROM_ID)
    cls.protSel4.pdbId.set('3o21')
    cls.protSel4.setObjLabel('Sel 3o21_CD_ca') # parallel
    cls.launchProtocol(cls.protSel4)

def importSetAS(cls):
    cls.protSetAS = cls.newProtocol(ProtImportSetOfAtomStructs, inputPdbData=1)
    cls.protSetAS.filesPath.set(split(cls.protSel2._getPath())[0])
    cls.protSetAS.filesPattern.set('*Select/*atoms.pdb')
    cls.protSetAS.setObjLabel('pwem import SetOfAS')
    cls.launchProtocol(cls.protSetAS)


class TestProDyEnsemble_1_SetAS_ASref(TestWorkflow):
    """ Test ProDy Ensemble protocol in way 1: SetOfAtomStructures with 
    atom struct ref and delete reference default (False)"""

    @classmethod
    def setUpClass(cls):
        # Create a new project
        setupTestProject(cls)
        importSelect(cls)
        importSetAS(cls)

    def testProDyEnsemble(cls):
        protEns1 = cls.newProtocol(ProDyBuildPDBEnsemble, refType=0,
                                    matchFunc=0)
        protEns1.structures.set([cls.protSetAS.outputAtomStructs])
        protEns1.refStructure.set(cls.protSel4.outputStructure)
        protEns1.setObjLabel('buildPDBEns_ref_3o21_CD')
        cls.launchProtocol(protEns1)

        numAS = len(protEns1.outputNpz)
        cls.assertSetSize(protEns1.outputNpz, 5,
                          "wrong size Npz Set object ({0} not 5) with "
                          "SetOfAtomStructs input and added atom struct ref"
                          "and delete reference default (False)".format(numAS))
        
class TestProDyEnsemble_2_SetAS_ASref_delRef(TestWorkflow):
    """ Test ProDy Ensemble protocol in way 2: SetOfAtomStructures with 
    atom struct ref and delete reference True"""

    @classmethod
    def setUpClass(cls):
        # Create a new project
        setupTestProject(cls)
        importSelect(cls)
        importSetAS(cls)

    def testProDyEnsemble(cls):
        protEns1 = cls.newProtocol(ProDyBuildPDBEnsemble, refType=0,
                                    matchFunc=0)
        protEns1.structures.set([cls.protSetAS.outputAtomStructs])
        protEns1.refStructure.set(cls.protSel4.outputStructure)
        protEns1.delReference.set(True)
        protEns1.setObjLabel('buildPDBEns_ref_3o21_CD_del')
        cls.launchProtocol(protEns1)

        numAS = len(protEns1.outputNpz)
        cls.assertSetSize(protEns1.outputNpz, 4,
                          "wrong size Npz Set object ({0} not 4) with "
                          "SetOfAtomStructs input and added atom struct ref"
                          "and delete reference True".format(numAS))


class TestProDyEnsemble_3_SetAS_IDref(TestWorkflow):
    """ Test ProDy Ensemble protocol in way 3: SetOfAtomStructures with 
    index ref default (1)"""

    @classmethod
    def setUpClass(cls):
        # Create a new project
        setupTestProject(cls)
        importSelect(cls)
        importSetAS(cls)

    def testProDyEnsemble(cls):
        protEns1 = cls.newProtocol(ProDyBuildPDBEnsemble, refType=0,
                                    matchFunc=0)
        protEns1.structures.set([cls.protSetAS.outputAtomStructs])
        protEns1.refType.set(INDEX)
        protEns1.setObjLabel('buildPDBEns_ref_id_1')
        cls.launchProtocol(protEns1)

        numAS = len(protEns1.outputNpz)
        cls.assertSetSize(protEns1.outputNpz, 4,
                          "wrong size Npz Set object ({0} not 4) with "
                          "SetOfAtomStructs input and index ref default (1)".format(numAS))

# class TestProDyImportEns_1(TestWorkflow):
#     """ Test ProDy Import Ensemble protocol with """

#     @classmethod
#     def setUpClass(cls):
#         # Create a new project
#         setupTestProject(cls)

#     def testImportEnsemble(cls):
#         protImportEns1 = cls.newProtocol(ProDyImportEnsemble)
#         #protImportEns1.inputEnsemble.set(protEns1.outputNpz)
#         protImportEns1.setObjLabel('Import ensemble')
#         cls.launchProtocol(protImportEns1)


# class TestProDyPCA_2(TestWorkflow):
#     """ Test protocol for ProDy Ensemble and Principal Component Analysis"""

#     @classmethod
#     def setUpClass(cls):
#         # Create a new project
#         setupTestProject(cls)
#         importSelect(cls)
#         importSetAS(cls)

#     def testProDyPCA(cls):
#         """buildPDBEns from SetOfAtomStructs with index ref -> PCA"""

#         ens1 = prody.loadEnsemble(protEns1._getPath("ensemble.ens.npz"))
#         idx = ens1.getLabels().index("3o21") + 1
        
#         protEns2 = cls.newProtocol(ProDyBuildPDBEnsemble, refType=1,
#                                     matchFunc=0)
#         protEns2.structures.set([protSetAS.outputAtomStructs])
#         protEns2.refIndex.set(idx)
#         protEns2.setObjLabel('buildPDBEns_set_ref_idx_{0}'.format(idx))
#         cls.launchProtocol(protEns2)

#         numAS = len(protEns2.outputNpz)
#         cls.assertSetSize(protEns2.outputNpz, 3,
#                            "wrong size SetOfAtomStructs ({0} not 3) with "
#                            "SetOfAtomStructs input and index ref".format(numAS))

#         protPca2 = cls.newProtocol(ProDyPCA, numberOfModes=2)
#         protPca2.inputEnsemble.set(protEns2.outputNpz)
#         protPca2.setObjLabel('PCA_from_set_ref_idx')
#         cls.launchProtocol(protPca2)
