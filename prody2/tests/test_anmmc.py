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

from pwem.tests.workflows import TestWorkflow
from pyworkflow.tests import setupTestProject

from prody2.protocols import ProDySelect, ProDyANMMC

class TestProDyANMMC(TestWorkflow):
    @classmethod
    def setUpClass(cls):
        # Create a new project
        setupTestProject(cls)
        importSelect(cls)

    def testProDyANMMC_single(cls):
        protMC = cls.newProtocol(ProDyANMMC)
        protMC.startingStructure.set(cls.protSel.outputStructure)
        cls.launchProtocol(protMC)

    def testProDyANMMC_targeted(cls):
        protMC = cls.newProtocol(ProDyANMMC, useTarget=True)
        protMC.startingStructure.set(cls.protSel.outputStructure)
        protMC.targetStructure.set(cls.protSel2.outputStructure)
        cls.launchProtocol(protMC)


def importSelect(cls):
    cls.protSel = cls.newProtocol(ProDySelect, 
                                  selection="name CA and chain A",
                                  inputPdbData=0)
    cls.protSel.pdbId.set("4ake")
    cls.protSel.setObjLabel('sel_4akeA_ca')
    cls.launchProtocol(cls.protSel)

    cls.protSel2 = cls.newProtocol(ProDySelect, 
                                  selection="name CA and chain A",
                                  inputPdbData=0)
    cls.protSel2.pdbId.set("4ake")
    cls.protSel2.setObjLabel('sel_4akeA_ca')
    cls.launchProtocol(cls.protSel2)
