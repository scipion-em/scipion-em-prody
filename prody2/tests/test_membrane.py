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

from prody2.protocols import (ProDyBiomol, ProDySelect, ProDyANM, ProDyGNM)

class TestProDyMembraneANM(TestWorkflow):
    @classmethod
    def setUpClass(cls):
        # Create a new project
        setupTestProject(cls)
        importBiomols(cls)
        selectCalpha(cls)

    def testProDyMembraneANM(cls):
        protANM = cls.newProtocol(ProDyANM, membrane=True)
        protANM.inputStructure.set(cls.protSel.outputStructure)
        protANM.setObjLabel('exANM_2NWL')
        cls.launchProtocol(protANM)

class TestProDyMembraneGNM(TestWorkflow):
    @classmethod
    def setUpClass(cls):
        # Create a new project
        setupTestProject(cls)
        importBiomols(cls)
        selectCalpha(cls)

    def testMemGNM(cls):
        protGNM = cls.newProtocol(ProDyGNM, membrane=True)
        protGNM.inputStructure.set(cls.protSel.outputStructure)
        protGNM.setObjLabel('exGNM_2NWL')
        cls.launchProtocol(protGNM)


def importBiomols(cls):
    cls.protBM = cls.newProtocol(ProDyBiomol, inputPdbData=0, membrane=True)
    cls.protBM.pdbId.set("2NWL")
    cls.protBM.setObjLabel('Biomol_2NWL-opm')
    cls.launchProtocol(cls.protBM)

def selectCalpha(cls):
    cls.protSel = cls.newProtocol(ProDySelect, selection="name CA",
                                inputPdbData=2) # import from pointer
    cls.protSel.inputStructure.set(cls.protBM.outputStructures)
    cls.protSel.inputStructure.setExtended(1)
    cls.protSel.setObjLabel('Sel_2NWL-opm_CA')
    cls.launchProtocol(cls.protSel)
