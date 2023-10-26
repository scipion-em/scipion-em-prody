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

class TestProDyFixer(TestWorkflow):
    """ Test protocol for ProDy Normal Mode Analysis and Deformation Analysis. """

    @classmethod
    def setUpClass(cls):
        # Create a new project
        setupTestProject(cls)

    def testProDyFixer(cls):
        """ Run membrane ANM and GNM for 3kg2 downloaded from OPM"""

        # --------------------------------------------------------------
        # Step 1. Import 3hsy and select chain A
        # --------------------------------------------------------------

        # Select Chain A
        protSel = cls.newProtocol(ProDySelect, selection="protein and chain A",
                                    inputPdbData=0)
        protSel.pdbId.set("3hsy")
        protSel.setObjLabel('sel_3hsyA')
        cls.launchProtocol(protSel)

        # --------------------------------------------------------------
        # Step 3. Membrane ANM
        # --------------------------------------------------------------
        protFix = cls.newProtocol(ProDyPDBFixer)
        protFix.inputStructure.set(protSel.outputStructure)
        protFix.setObjLabel('fix_3hsyA')
        cls.launchProtocol(protFix)

