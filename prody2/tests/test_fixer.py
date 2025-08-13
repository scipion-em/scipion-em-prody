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

from prody2.protocols import (ProDySelect, ProDyPDBFixer)
from prody2.constants import N_ATOMS

class TestProDyFixer(TestWorkflow):
    @classmethod
    def setUpClass(cls):
        # Create a new project
        setupTestProject(cls)
        importSelect(cls)

    def testProDyFixer(cls):
        protFix = cls.newProtocol(ProDyPDBFixer)
        protFix.inputStructure.set(cls.protSel.outputStructure)
        protFix.setObjLabel('fix_3hsyB')
        cls.launchProtocol(protFix)

        struct1 = protFix.outputStructure
        numAtoms = struct1.getAttributeValue(N_ATOMS)
        cls.assertTrue(numAtoms == 5956,
                       "After fixing, 3hsy B should have 5956 atoms, not {0}".format(numAtoms))

def importSelect(cls):
    cls.protSel = cls.newProtocol(ProDySelect, selection="protein and chain B",
                                  inputPdbData=0)
    cls.protSel.pdbId.set("3hsy")
    cls.protSel.setObjLabel('sel_3hsyB')
    cls.launchProtocol(cls.protSel)
