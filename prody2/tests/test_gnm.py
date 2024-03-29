# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     James Krieger (jmkrieger@cnb.csic.es)
# *              Ricardo Serrano Gutiérrez (rserranogut@hotmail.com)   
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

from pwem.protocols import ProtImportPdb, exists
from pwem.tests.workflows import TestWorkflow
from pyworkflow.tests import setupTestProject

from prody2.protocols import (ProDySelect, ProDyGNM, ProDyEdit, ProDyCompare,
                              ProDyDomainDecomp)

from prody2.protocols.protocol_edit import NMA_SLICE, NMA_REDUCE, NMA_EXTEND
from prody2.objects import SetOfGnmModes

import prody

gnmModesTypeWarning = "GNM modes should be parsed as a SetOfGnmModes, not {0}"

class TestProDyGNM(TestWorkflow):
    """ Test protocol for ProDy Gaussian Normal Model Analysis. """

    @classmethod
    def setUpClass(cls):
        # Create a new project
        setupTestProject(cls)

    def testProDyGNM(self):
        """ Run GNM simple workflow for two Atomic structures. """
        
        oldVerbosity = prody.confProDy("verbosity")
        oldSecondary = prody.confProDy("auto_secondary")

        # ------------------------------------------------
        # Step 1. Import a Pdb -> Select chain A -> GNM
        # ------------------------------------------------
        # Import a PDB
        protImportPdb1 = self.newProtocol(ProtImportPdb, inputPdbData=0,
                                          pdbId="4ake")
        protImportPdb1.setObjLabel('pwem import 4ake')
        self.launchProtocol(protImportPdb1)

        # Select Chain A
        protSel1 = self.newProtocol(ProDySelect, selection="protein and chain A")
        protSel1.inputStructure.set(protImportPdb1.outputPdb)
        protSel1.setObjLabel('Sel_4akeA_all')
        self.launchProtocol(protSel1)

        # Launch GNM NMA for chain A (all atoms)
        protGNM1 = self.newProtocol(ProDyGNM, cutoff=4)
        protGNM1.inputStructure.set(protSel1.outputStructure)
        protGNM1.setObjLabel('GNM_all')
        self.launchProtocol(protGNM1)

        self.assertTrue(isinstance(protGNM1.outputModes, SetOfGnmModes),
                        gnmModesTypeWarning.format(type(protGNM1.outputModes)))

        # ------------------------------------------------
        # Step 2. Select CA -> GNM NMA
        # ------------------------------------------------
        # Select Calpha atoms
        protSel2 = self.newProtocol(ProDySelect, selection="name CA")
        protSel2.inputStructure.set(protSel1.outputStructure)
        protSel2.setObjLabel('Sel_4akeA_CA')
        self.launchProtocol(protSel2)

        # Launch GNM NMA for selected atoms (CA) with zero mode (default)
        protGNM2 = self.newProtocol(ProDyGNM, cutoff=10)
        protGNM2.inputStructure.set(protSel2.outputStructure)
        protGNM2.setObjLabel('GNM_CA_z')
        self.launchProtocol(protGNM2)

        self.assertFalse(exists(protGNM2._getExtraPath("animations/animated_mode_001.pdb")))
        self.assertFalse(exists(protGNM2._getExtraPath("animations/animated_mode_002.pdb")))
        # (no animations from GNM)

        self.assertFalse(exists(protGNM2._getExtraPath("distanceProfiles/vec1.xmd")))
        self.assertTrue(exists(protGNM2._getExtraPath("distanceProfiles/vec2.xmd")))

        # Launch ANM NMA for selected atoms (CA) without zeros
        protGNM2b = self.newProtocol(ProDyGNM)
        protGNM2b.inputStructure.set(protSel2.outputStructure)
        protGNM2b.zeros.set(False)
        protGNM2b.setObjLabel('GNM_CA_n-z')
        self.launchProtocol(protGNM2b)

        self.assertFalse(exists(protGNM2b._getExtraPath("animations/animated_mode_001.pdb")))
        # (no animations from GNM)

        self.assertTrue(exists(protGNM2b._getExtraPath("distanceProfiles/vec1.xmd")))

        # ------------------------------------------------
        # Step 3. Slice -> Compare
        # ------------------------------------------------
        # Slice all-atom GNM NMA to CA
        protEdit1 = self.newProtocol(ProDyEdit, edit=NMA_SLICE)
        protEdit1.modes.set(protGNM1.outputModes)
        protEdit1.newNodes.set(protSel2.outputStructure)
        protEdit1.setObjLabel('Slice_to_CA')
        self.launchProtocol(protEdit1)

        self.assertTrue(isinstance(protEdit1.outputModes, SetOfGnmModes),
                        gnmModesTypeWarning.format(type(protEdit1.outputModes)))

        # Compare sliced and original CA NMA
        protComp1 = self.newProtocol(ProDyCompare)
        protComp1.modes1.set(protGNM2.outputModes)
        protComp1.modes2.set(protEdit1.outputModes)
        protComp1.setObjLabel('Compare_slcAA_to_CA')
        self.launchProtocol(protComp1)

        # ------------------------------------------------
        # Step 4. Reduce -> Compare
        # ------------------------------------------------
        # Reduce all-atom GNM NMA to CA
        protEdit2 = self.newProtocol(ProDyEdit, edit=NMA_REDUCE)
        protEdit2.modes.set(protGNM1.outputModes)
        protEdit2.newNodes.set(protSel2.outputStructure)
        protEdit2.setObjLabel('Reduce_to_CA')
        self.launchProtocol(protEdit2)

        # Compare reduced and original CA NMA
        protComp2 = self.newProtocol(ProDyCompare)
        protComp2.modes1.set(protGNM2.outputModes)
        protComp2.modes2.set(protEdit2.outputModes)
        protComp2.setObjLabel('Compare_redAA_to_CA')
        self.launchProtocol(protComp2)

        # ------------------------------------------------
        # Step 5. Extend -> Compare
        # ------------------------------------------------
        # Extend CA NMA to all-atoms
        protEdit3 = self.newProtocol(ProDyEdit, edit=NMA_EXTEND)
        protEdit3.modes.set(protGNM2.outputModes)
        protEdit3.newNodes.set(protSel1.outputStructure)
        protEdit3.setObjLabel('Extend_to_AA')
        self.launchProtocol(protEdit3)

        # Compare original AA GNM NMA and extended CA GNM NMA
        # Test matching too
        protComp3 = self.newProtocol(ProDyCompare, match=True)
        protComp3.modes1.set(protGNM1.outputModes)
        protComp3.modes2.set(protEdit3.outputModes)
        protComp3.setObjLabel('Compare_AA_to_extCA')
        self.launchProtocol(protComp3)

        self.assertTrue(isinstance(protComp3.outputModes, SetOfGnmModes),
                        gnmModesTypeWarning.format(type(protComp3.outputModes)))

        # ------------------------------------------------
        # Step 6. CA -> Domain Decomposition
        # ------------------------------------------------
        # Domain Decomposition on CA GNM
        protDomDec1 =  self.newProtocol(ProDyDomainDecomp)
        protDomDec1.modesGNM.set(protGNM2.outputModes)
        protDomDec1.setObjLabel('DomainDecomp_CA')
        self.launchProtocol(protDomDec1)

        self.assertTrue(prody.confProDy("verbosity") == oldVerbosity, 
                        "prody verbosity changed")

        self.assertTrue(prody.confProDy("auto_secondary") == oldSecondary, 
                        "prody auto_secondary changed")
