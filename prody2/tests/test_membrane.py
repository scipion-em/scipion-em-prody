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

class TestProDyMembrane(TestWorkflow):
    """ Test protocol for ProDy Normal Mode Analysis and Deformation Analysis. """

    @classmethod
    def setUpClass(cls):
        # Create a new project
        setupTestProject(cls)

    def testProDyMembrane(cls):
        """ Run membrane ANM and GNM for 3kg2 downloaded from OPM"""

        # --------------------------------------------------------------
        # Step 1. Import 3kg2 biologically relevant assembly in membrane
        # --------------------------------------------------------------

        # Select Chain A
        protBM = cls.newProtocol(ProDyBiomol, inputPdbData=0, membrane=True)
        protBM.pdbId.set("2NWL")
        protBM.setObjLabel('Biomol_2NWL-opm')
        cls.launchProtocol(protBM)

        # ------------------------------------------------
        # Step 2. Select CA
        # ------------------------------------------------
        # Select Calpha atoms
        protSel = cls.newProtocol(ProDySelect, selection="name CA",
                                   inputPdbData=2) # import from pointer
        protSel.inputStructure.set(protBM.outputStructures)
        protSel.inputStructure.setExtended(1)
        protSel.setObjLabel('Sel_2NWL-opm_CA')
        cls.launchProtocol(protSel)

        # --------------------------------------------------------------
        # Step 3. Membrane ANM
        # --------------------------------------------------------------
        protANM = cls.newProtocol(ProDyANM, membrane=True)
        protANM.inputStructure.set(protSel.outputStructure)
        protANM.setObjLabel('exANM_2NWL')
        cls.launchProtocol(protANM)

        # --------------------------------------------------------------
        # Step 4. Membrane GNM
        # --------------------------------------------------------------
        protGNM = cls.newProtocol(ProDyGNM, membrane=True)
        protGNM.inputStructure.set(protSel.outputStructure)
        protGNM.setObjLabel('exGNM_2NWL')
        cls.launchProtocol(protGNM)
