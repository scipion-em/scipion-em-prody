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

from pwem.tests.workflows import TestWorkflow
from pyworkflow.tests import setupTestProject

from prody2.protocols import (ProDySelect, ProDyBlastPDB)
from pwem.protocols.protocol_import import ProtImportSequence

class TestProDyBlast(TestWorkflow):
    @classmethod
    def setUpClass(cls):
        # Create a new project
        setupTestProject(cls)
        importSelectSeq(cls)

    def testProDyBlastStruct(cls):
        protBlast = cls.newProtocol(ProDyBlastPDB)
        protBlast.inputSeqData.set(protBlast.IMPORT_FROM_STRUCT)
        protBlast.inputStructure.set(cls.protSel.outputStructure)
        protBlast.setObjLabel('blast_3hsyB')
        cls.launchProtocol(protBlast)

    def testProDyBlastSeq(cls):
        protBlast = cls.newProtocol(ProDyBlastPDB)
        protBlast.inputSeqData.set(protBlast.IMPORT_FROM_SEQ)
        protBlast.inputSequence.set(cls.protSeq.outputSequence)
        protBlast.setObjLabel('blast_3hsyB_seq')
        cls.launchProtocol(protBlast)

    def testProDyBlastText(cls):
        protBlast = cls.newProtocol(ProDyBlastPDB)
        protBlast.inputSeqData.set(protBlast.USE_TEXT)
        protBlast.inputSeqText.set("""MQKIMHISVLLSPVLWGLIFGVSSNSIQIGGLFPRGADQEYSAFRVGMVQFSTSEFRLTP
HIDNLEVANSFAVTNAFCSQFSRGVYAIFGFYDKKSVNTITSFCGTLHVSFITPSFPTDG
THPFVIQMRPDLKGALLSLIEYYQWDKFAYLYDSDRGLSTLQAVLDSAAEKKWQVTAINV
GNINNDKKDETYRSLFQDLELKKERRVILDCERDKVNDIVDQVITIGKHVKGYHYIIANL
GFTDGDLLKIQFGGANVSGFQIVDYDDSLVSKFIERWSTLEEKEYPGAHTATIKYTSALT
YDAVQVMTEAFRNLRKQRIEISRRGNAGDCLANPAVPWGQGVEIERALKQVQVEGLSGNI
KFDQNGKRINYTINIMELKTNGPRKIGYWSEVDKMVVTLTELPSGNDTSGLENKTVVVTT
ILESPYVMMKKNHEMLEGNERYEGYCVDLAAEIAKHCGFKYKLTIVGDGKYGARDADTKI
WNGMVGELVYGKADIAIAPLTITLVREEVIDFSKPFMSLGISIMIKKPQKSKPGVFSFLD
PLAYEIWMCIVFAYIGVSVVLFLVSRFSPYEWHTEEFEDGRETQSSESTNEFGIFNSLWF
SLGAFMQQGCDISPRSLSGRIVGGVWWFFTLIIISSYTANLAAFLTVERMVSPIESAEDL
SKQTEIAYGTLDSGSTKEFFRRSKIAVFDKMWTYMRSAEPSVFVRTTAEGVARVRKSKGK
YAYLLESTMNEYIEQRKPCDTMKVGGNLDSKGYGIATPKGSSLGNAVNLAVLKLNEQGLL
DKLKNKWWYDKGECGSGGGDSKEKTSALSLSNVAGVFYILVGGLGLAMLVALIEFCYKSR
AEAKRMKVAKNPQNINPSSSQNSQNFATYKEGYNVYGIESVKI""")
        protBlast.setObjLabel('blast_P19491_text')
        cls.launchProtocol(protBlast)

def importSelectSeq(cls):
    cls.protSel = cls.newProtocol(ProDySelect, selection="protein and chain B",
                                  inputPdbData=0)
    cls.protSel.pdbId.set("3hsy")
    cls.protSel.setObjLabel('sel_3hsyB')
    cls.launchProtocol(cls.protSel)

    cls.protSeq = cls.newProtocol(ProtImportSequence)
    cls.protSeq.inputProteinSequence.set(cls.protSeq.IMPORT_FROM_STRUCTURE)
    cls.protSeq.inputStructureSequence.set(cls.protSeq.IMPORT_STRUCTURE_FROM_AS)
    cls.protSeq.inputAS.set(cls.protSel.outputStructure)
    cls.protSeq.setObjLabel('seq_3hsyB')
    cls.launchProtocol(cls.protSeq)
