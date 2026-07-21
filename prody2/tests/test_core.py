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

import os

from pwem.protocols import ProtImportPdb, exists
from pwem.tests.workflows import TestWorkflow
from pyworkflow.tests import setupTestProject

from prody2.protocols import (ProDySelect, ProDyAlign,
                              ProDyBiomol, ProDyRenumber, ProDyAddPDBs,
                              ProDyANM, ProDyRTB, ProDyDefvec,
                              ProDyEdit, ProDyCompare,
                              ProDyImportModes, ProDyAlgebra)

from prody2.protocols.protocol_edit import NMA_SLICE, NMA_REDUCE, NMA_EXTEND, NMA_INTERP
from prody2.protocols.protocol_rtb import BLOCKS_FROM_RES, BLOCKS_FROM_SECSTR
from prody2.protocols.protocol_import import MODES_NPZ, SCIPION, NMD

from prody2.constants import (PRODY_TEST_PDB_FILE, N_RESIDUES, N_CHAINS,
                              FIRST_RESNUM, LAST_RESNUM, MAX_RESNUM,
                              PRODY_TEST_ALG_PDB_FILE,
                              PRODY_TEST_TAR_PDB_FILE,
                              TESTDIR)

import numpy as np

animationsFile7 = "animations/animated_mode_007.pdb"
animationsFile1 = "animations/animated_mode_001.pdb"
distProfile1 = "distanceProfiles/vec1.xmd"
distProfile7 = "distanceProfiles/vec7.xmd"

renumFilename = "renum_atoms.pdb"

class TestProDyDefvec(TestWorkflow):
    """ Test protocol for ProDy Deformation Vector Analysis. """

    @classmethod
    def setUpClass(cls):
        # Create a new project
        setupTestProject(cls)
        importOnly4akeA(cls)
        importAligned1akeA(cls)


    def testProDyDefvec(cls):
        """ Run deformation vector calculation and confirm if it works """

        # Defvec from same starting point as NMA
        protDefvec5 = cls.newProtocol(ProDyDefvec, rmsd=5)
        protDefvec5.mobStructure.set(cls.protImportPdb4akeA.outputPdb)
        protDefvec5.tarStructure.set(cls.protImportPdb1akeA.outputPdb)
        protDefvec5.setObjLabel('Defvec_5A_4akeA_1akeA_CA')
        cls.launchProtocol(protDefvec5)

class TestProDyCore1(TestWorkflow):
    """ Test protocol for ProDy Anisotropic Network Model (ANM) Normal Mode Analysis (NMA) and Deformation Analysis. """

    @classmethod
    def setUpClass(cls):
        # Create a new project
        setupTestProject(cls)
        importSelect4ake(cls)
        importSelect1ake(cls)

    def testProDySelect(cls):
        """ Run different selection options and confirm if it works """

        # ----------------------------------------------------------------------
        # Step 1a. Import a Pdb -> Select chain A and other options from Pointer
        # ----------------------------------------------------------------------
        # Import a PDB
        protImportPdb1 = cls.newProtocol(ProtImportPdb, inputPdbData=0,
                                         pdbId="4ake",
                                         skipChimera=True)
        protImportPdb1.setObjLabel('pwem import 4ake')
        cls.launchProtocol(protImportPdb1)

        chainAselstr = "protein and chain A"
        outputFilename = "4ake_atoms.pdb"
        chainCselstr = "chain C"

        # Select Chain A
        protSel1 = cls.newProtocol(ProDySelect, selection=chainAselstr)
        protSel1.inputStructure.set(protImportPdb1.outputPdb)
        protSel1.setObjLabel('Sel_4akeA_all_pointer')
        cls.launchProtocol(protSel1)

        cls.assertTrue(exists(protSel1._getPath(outputFilename)))
        cls.assertTrue(hasattr(protSel1, "outputStructure"))

        # Select chain C CA to show it doesn't work
        protSel1a = cls.newProtocol(ProDySelect, selection="chain C and name CA")
        protSel1a.inputStructure.set(protImportPdb1.outputPdb)
        protSel1a.setObjLabel('Sel 4ake_C_pointer unite_False')
        cls.launchProtocol(protSel1a)

        cls.assertFalse(exists(protSel1a._getPath(outputFilename)))
        cls.assertFalse(hasattr(protSel1a, "outputStructure"))

        # Select whole chain C with uniteChains False (default) to show it does work
        protSel1a = cls.newProtocol(ProDySelect, selection=chainCselstr)
        protSel1a.inputStructure.set(protImportPdb1.outputPdb)
        protSel1a.setObjLabel('Sel 4ake_C_pointer_1')
        cls.launchProtocol(protSel1a)

        cls.assertTrue(exists(protSel1a._getPath(outputFilename)))
        cls.assertTrue(hasattr(protSel1a, "outputStructure"))

        # Select whole chain C with uniteChains True to show it doesn't work
        protSel1a = cls.newProtocol(ProDySelect, selection=chainCselstr,
                                    uniteChains=True)
        protSel1a.inputStructure.set(protImportPdb1.outputPdb)
        protSel1a.setObjLabel('Sel 4ake_C_pointer unite_True')
        cls.launchProtocol(protSel1a)

        cls.assertFalse(exists(protSel1a._getPath(outputFilename)))
        cls.assertFalse(hasattr(protSel1a, "outputStructure"))

        # Select chain C whole to show it does work with uniteChains False (default)
        protSel1a2 = cls.newProtocol(ProDySelect, selection=chainCselstr)
        protSel1a2.inputStructure.set(protImportPdb1.outputPdb)
        protSel1a2.setObjLabel('Sel 4ake_C_pointer_2')
        cls.launchProtocol(protSel1a2)

        cls.assertTrue(exists(protSel1a2._getPath(outputFilename)))
        cls.assertTrue(hasattr(protSel1a2, "outputStructure"))

        # Select chain C whole to show it does work with uniteChains True
        protSel1a3 = cls.newProtocol(ProDySelect, selection=chainCselstr,
                                     uniteChains=True)
        protSel1a3.inputStructure.set(protImportPdb1.outputPdb)
        protSel1a3.setObjLabel('Sel 4ake_C_pointer_3')
        cls.launchProtocol(protSel1a3)

        cls.assertFalse(exists(protSel1a3._getPath(outputFilename)))
        cls.assertFalse(hasattr(protSel1a3, "outputStructure"))

        # ------------------------------------------------
        # Step 1b. Select chain A from Filename
        # ------------------------------------------------
        protSel1b = cls.newProtocol(ProDySelect, selection=chainAselstr,
                                     inputPdbData=1)
        protSel1b.pdbFile.set(protImportPdb1.outputPdb.getFileName())
        protSel1b.setObjLabel('Sel_4akeA_all_file')
        cls.launchProtocol(protSel1b)

        # ----------------------------------------------------
        # Step 1c. Select chain A from PDB id (difficult case)
        # ----------------------------------------------------
        protSel1c = cls.newProtocol(ProDySelect, selection=chainAselstr,
                                     inputPdbData=0)
        protSel1c.pdbId.set("6xr8")
        protSel1c.setObjLabel('Sel_6xr8_A_all_id')
        cls.launchProtocol(protSel1c)

    def testProDyCore(cls):
        """ Run NMA simple workflow for two Atomic structures. """
        # --------------------------------------------------------------
        # Step 1a. Import a Pdb with 4ake chain A -> ANM NMA (all-atoms)
        # --------------------------------------------------------------

        # Import a PDB
        protImportPdb1 = cls.newProtocol(ProtImportPdb, inputPdbData=1,
                                         pdbFile=PRODY_TEST_PDB_FILE,
                                         skipChimera=True)
        protImportPdb1.setObjLabel('pwem import 4ake')
        cls.launchProtocol(protImportPdb1)

        # Launch ANM NMA for chain A (all atoms) with zero modes (default)
        protANM1 = cls.newProtocol(ProDyANM, cutoff=8)
        protANM1.inputStructure.set(protImportPdb1.outputPdb)
        protANM1.setObjLabel('ANM_all')
        cls.launchProtocol(protANM1)

        # check that we have animations and distance profiles for non-zero modes from 7
        cls.assertFalse(exists(protANM1._getExtraPath(animationsFile1)))
        cls.assertTrue(exists(protANM1._getExtraPath(animationsFile7)))

        cls.assertFalse(exists(protANM1._getExtraPath(distProfile1)))
        cls.assertTrue(exists(protANM1._getExtraPath(distProfile7)))

        # ------------------------------------------------
        # Step 2. CA ANM NMA
        # ------------------------------------------------

        # Launch ANM NMA for selected atoms (CA) with zeros (default)
        protANM2 = cls.newProtocol(ProDyANM)
        protANM2.inputStructure.set(cls.protSel.outputStructure)
        protANM2.setObjLabel('ANM_CA_z')
        cls.launchProtocol(protANM2)

        # check that we have animations and distance profiles for non-zero modes from 7
        cls.assertFalse(exists(protANM2._getExtraPath(animationsFile1)))
        cls.assertTrue(exists(protANM2._getExtraPath(animationsFile7)))

        cls.assertFalse(exists(protANM2._getExtraPath(distProfile1)))
        cls.assertTrue(exists(protANM2._getExtraPath(distProfile7)))

        # Launch ANM NMA for selected atoms (CA) without zeros
        protANM2b = cls.newProtocol(ProDyANM)
        protANM2b.inputStructure.set(cls.protSel.outputStructure)
        protANM2b.zeros.set(False)
        protANM2b.setObjLabel('ANM_CA_n-z')
        cls.launchProtocol(protANM2b)

        # check that we have animations and distance profiles for non-zero modes from 1
        cls.assertTrue(exists(protANM2b._getExtraPath(animationsFile1)))
        cls.assertTrue(exists(protANM2b._getExtraPath(distProfile1)))

        # Launch ed-ENM NMA for selected atoms (CA) without zeros
        protANM3 = cls.newProtocol(ProDyANM)
        protANM3.inputStructure.set(cls.protSel.outputStructure)
        protANM3.zeros.set(False)
        protANM3.gamma.set("GammaED")
        protANM3.cutoff.set("2.9 * math.log(214) - 2.9")
        protANM3.setObjLabel('edENM_CA_n-z')
        cls.launchProtocol(protANM3)

        # check that we have animations and distance profiles for non-zero modes from 1
        cls.assertTrue(exists(protANM3._getExtraPath(animationsFile1)))
        cls.assertTrue(exists(protANM3._getExtraPath(distProfile1)))

        # ------------------------------------------------
        # Step 3. Slice -> Compare
        # ------------------------------------------------
        # Slice all-atom ANM NMA to CA (no animations by default)
        protEdit1 = cls.newProtocol(ProDyEdit, edit=NMA_SLICE)
        protEdit1.modes.set(protANM1.outputModes)
        protEdit1.newNodes.set(cls.protSel.outputStructure)
        protEdit1.setObjLabel('Slice_to_CA')
        cls.launchProtocol(protEdit1)
        cls.assertFalse(exists(protEdit1._getExtraPath(animationsFile7)))

        # Compare sliced and original CA NMA
        protComp1 = cls.newProtocol(ProDyCompare)
        protComp1.modes1.set(protANM2.outputModes)
        protComp1.modes2.set(protEdit1.outputModes)
        protComp1.setObjLabel('Compare_slcAA_to_CA')
        cls.launchProtocol(protComp1)

        # Slice all-atom ANM NMA to CA (with animations)
        protEdit1a = cls.newProtocol(ProDyEdit, edit=NMA_SLICE,
                                     doAnimation=True)
        protEdit1a.modes.set(protANM1.outputModes)
        protEdit1a.newNodes.set(cls.protSel.outputStructure)
        protEdit1a.setObjLabel('Slice_to_CA_anim')
        cls.launchProtocol(protEdit1a)
        cls.assertTrue(exists(protEdit1a._getExtraPath(animationsFile7)))

        # ------------------------------------------------
        # Step 4. Reduce -> Compare
        # ------------------------------------------------
        # Reduce all-atom ANM NMA to CA
        protEdit2 = cls.newProtocol(ProDyEdit, edit=NMA_REDUCE)
        protEdit2.modes.set(protANM1.outputModes)
        protEdit2.newNodes.set(cls.protSel.outputStructure)
        protEdit2.setObjLabel('Reduce_to_CA')
        cls.launchProtocol(protEdit2)

        # Compare reduced and original CA NMA
        protComp2 = cls.newProtocol(ProDyCompare)
        protComp2.modes1.set(protANM2.outputModes)
        protComp2.modes2.set(protEdit2.outputModes)
        protComp2.setObjLabel('Compare_redAA_to_CA')
        cls.launchProtocol(protComp2)

        # ------------------------------------------------
        # Step 5. Extend -> Compare
        # ------------------------------------------------
        # Extend CA NMA to all-atoms
        protEdit3 = cls.newProtocol(ProDyEdit, edit=NMA_EXTEND)
        protEdit3.modes.set(protANM2.outputModes)
        protEdit3.newNodes.set(protImportPdb1.outputPdb)
        protEdit3.setObjLabel('Extend_to_AA')
        cls.launchProtocol(protEdit3)

        # Compare original AA ANM NMA and extended CA ANM NMA
        # Test matching too
        protComp3 = cls.newProtocol(ProDyCompare, match=True)
        protComp3.modes1.set(protANM1.outputModes)
        protComp3.modes2.set(protEdit3.outputModes)
        protComp3.setObjLabel('Compare_AA_to_extCA')
        cls.launchProtocol(protComp3)

        # ------------------------------------------------
        # Step 6. Interpolate -> Compare
        # ------------------------------------------------
        # Interpolate CA NMA to all-atoms
        protEdit4 = cls.newProtocol(ProDyEdit, edit=NMA_INTERP)
        protEdit4.modes.set(protANM2.outputModes)
        protEdit4.newNodes.set(protImportPdb1.outputPdb)
        protEdit4.setObjLabel('Interp_to_AA')
        cls.launchProtocol(protEdit4)

        # Compare original AA ANM NMA and interpolated CA ANM NMA
        protComp4 = cls.newProtocol(ProDyCompare)
        protComp4.modes1.set(protANM1.outputModes)
        protComp4.modes2.set(protEdit4.outputModes)
        protComp4.setObjLabel('Compare_AA_to_intCA')
        cls.launchProtocol(protComp4)

        # --------------------------------------------------
        # Step 7. Import other Pdb -> Select chain A and CA
        # -> align -> defvec -> compare
        # --------------------------------------------------

        # Align new pdb from protSel3 to match old pdb with NMA
        protAlign1 = cls.newProtocol(ProDyAlign)
        protAlign1.mobStructure.set(cls.protSel3.outputStructure)
        protAlign1.tarStructure.set(cls.protSel.outputStructure)
        protAlign1.setObjLabel('Align_1akeA_4akeA_CA')
        cls.launchProtocol(protAlign1)

        # Defvec from same starting point as NMA
        protDefvec1 = cls.newProtocol(ProDyDefvec)
        protDefvec1.mobStructure.set(protAlign1.outputStructureTar)
        protDefvec1.tarStructure.set(protAlign1.outputStructureMob)
        protDefvec1.setObjLabel('Defvec_4akeA_1akeA_CA')
        cls.launchProtocol(protDefvec1)

        # Compare original CA NMA to defvec with default overlaps
        protComp5 = cls.newProtocol(ProDyCompare)
        protComp5.modes1.set(protANM2.outputModes)
        protComp5.modes2.set(protDefvec1.outputModes)
        protComp5.setObjLabel('Compare_ANM_to_Defvec')
        cls.launchProtocol(protComp5)

        compMatrix5 = np.loadtxt(protComp5._getPath('matrix.txt'))
        cls.assertTrue(max(compMatrix5) <= 1, "Default defvec comparison is normalised")

        # Compare original CA NMA to defvec with raw overlaps
        protComp6 = cls.newProtocol(ProDyCompare)
        protComp6.norm.set(False)
        protComp6.modes1.set(protANM2.outputModes)
        protComp6.modes2.set(protDefvec1.outputModes)
        protComp6.setObjLabel('Compare_ANM_to_Defvec_raw')
        cls.launchProtocol(protComp6)

        compMatrix6 = np.loadtxt(protComp6._getPath('matrix.txt'))
        cls.assertTrue(max(compMatrix6) > 1, "Raw defvec comparison generated larger numbers")

        # ------------------------------------------------
        # Step 8. Import ANM & compare scipion vs prody npz
        # -> import -> import -> compare
        # ------------------------------------------------
        # Define path
        modes = protANM2.outputModes
        modesPath = os.path.dirname(os.path.dirname(
            modes._getMapper().selectFirst().getModeFile()))

        # Import modes from prody npz
        protImportModes1 = cls.newProtocol(ProDyImportModes)
        protImportModes1.importType.set(MODES_NPZ)
        protImportModes1.filesPath.set(modesPath)
        protImportModes1.filesPattern.set("modes.anm.npz")
        protImportModes1.inputStructure.set(cls.protSel.outputStructure)
        protImportModes1.setObjLabel('import_npz_ANM_CA')
        cls.launchProtocol(protImportModes1)

        # Import scipion modes
        protImportModes2 = cls.newProtocol(ProDyImportModes)
        protImportModes2.importType.set(SCIPION)
        protImportModes2.filesPath.set(protANM2.outputModes.getFileName())
        protImportModes2.inputStructure.set(cls.protSel.outputStructure)
        protImportModes2.setObjLabel('import_scipion_ANM_CA')
        cls.launchProtocol(protImportModes2)  

        # Compare two imported ANMs
        protComp6 = cls.newProtocol(ProDyCompare)
        protComp6.modes1.set(protImportModes1.outputModes)
        protComp6.modes2.set(protImportModes2.outputModes)
        protComp6.setObjLabel('Compare_imported_ANMs')
        cls.launchProtocol(protComp6)

        # ------------------------------------------------
        # Step 9. Confirm mode algebra works for ANM
        # ------------------------------------------------
        protAlgebra = cls.newProtocol(ProDyAlgebra)
        protAlgebra.modes.set(protANM1.outputModes)
        protAlgebra.coeffString.set("1,2")
        protAlgebra.setObjLabel('mode algebra')
        cls.launchProtocol(protAlgebra)

class TestProDyRTB(TestWorkflow):
    """ Test protocol for ProDy Rotating and Translating Blocks (RTB) Normal Mode Analysis (NMA)"""

    @classmethod
    def setUpClass(cls):
        # Create a new project
        setupTestProject(cls)
        importSelect4ake(cls)

    def testProDyRTB(cls):
        # -------------------------------------------------------
        # Step 9. RTB in 2 ways -> Compare to each other and ANM
        # -------------------------------------------------------

        # Launch RTB NMA for selected atoms (CA) with 10 res per block
        #  with zeros (default)
        protRTB1 = cls.newProtocol(ProDyRTB, blockDef=BLOCKS_FROM_RES)
        protRTB1.inputStructure.set(cls.protSel.outputStructure)
        protRTB1.setObjLabel('RTB_CA_10_res')
        cls.launchProtocol(protRTB1)

        # check that we have animations and distance profiles for non-zero modes from 7
        cls.assertFalse(exists(protRTB1._getExtraPath(animationsFile1)))
        cls.assertTrue(exists(protRTB1._getExtraPath(animationsFile7)))

        cls.assertFalse(exists(protRTB1._getExtraPath(distProfile1)))
        cls.assertTrue(exists(protRTB1._getExtraPath(distProfile7)))

        # Launch RTB NMA for selected atoms (CA) with secstr block
        #  without zeros
        protRTB2 = cls.newProtocol(ProDyRTB, blockDef=BLOCKS_FROM_SECSTR)
        protRTB2.inputStructure.set(cls.protSel.outputStructure)
        protRTB2.zeros.set(False)
        protRTB2.setObjLabel('RTB_CA_secstr')
        cls.launchProtocol(protRTB2)

        # check that we have animations and distance profiles for non-zero modes from 1
        cls.assertTrue(exists(protRTB2._getExtraPath(animationsFile1)))
        cls.assertTrue(exists(protRTB2._getExtraPath(distProfile1)))

        # Compare RTB1 and RTB2
        protComp6 = cls.newProtocol(ProDyCompare)
        protComp6.modes1.set(protRTB1.outputModes)
        protComp6.modes2.set(protRTB2.outputModes)
        protComp6.setObjLabel('Compare_RTB1_to_RTB2')
        cls.launchProtocol(protComp6)

class TestProDyAtomic(TestWorkflow):
    """ Test protocol for other ProDy atomic operations"""

    @classmethod
    def setUpClass(cls):
        # Create a new project
        setupTestProject(cls)
        importSelect4ake(cls)
        importSelect1ake(cls)

    def testProDyAdd(cls):
        # extract biomol from 4ake (dimer) from id
        protAdd = cls.newProtocol(ProDyAddPDBs)
        protAdd.inputStructure.set([cls.protSel.outputStructure,
                                    cls.protSel3.outputStructure])
        protAdd.setObjLabel('Add_4akeA_1akeA')
        cls.launchProtocol(protAdd)

        struct1 = protAdd.outputStructure
        nResidues = struct1.getAttributeValue(N_RESIDUES)
        nChains = struct1.getAttributeValue(N_CHAINS)
        cls.assertTrue(nResidues == 428,
                       "AddPDBs output should have 428 residues, not {0}".format(nResidues))
        cls.assertTrue(nChains == 2,
                       "AddPDBs output should have 2 chains, not {0}".format(nChains))

    def testProDyBiomol(cls):
        # extract biomol from 4ake (dimer) from id
        protBm1 = cls.newProtocol(ProDyBiomol)
        protBm1.inputPdbData.set(0)
        protBm1.pdbId.set("4ake")
        protBm1.setObjLabel('Biomol_4ake_id')
        cls.launchProtocol(protBm1)

        numStructs = len(protBm1.outputStructures)
        cls.assertTrue(numStructs == 1, "Failed to extract 1 biomol from 4ake (dimer)")

        struct1 = protBm1.outputStructures.getFirstItem()
        nResidues = struct1.getAttributeValue(N_RESIDUES)
        nChains = struct1.getAttributeValue(N_CHAINS)
        cls.assertTrue(nResidues == 575,
                       "1ake biomol 1 should have 575 residues, not {0}".format(nResidues))
        cls.assertTrue(nChains == 2,
                       "1ake biomol 1 should have 2 chains, not {0}".format(nChains))

        # extract biomols from 1ake (2 monomers) from pointer with uniteChains False (default)
        protBm2 = cls.newProtocol(ProDyBiomol)
        protBm2.inputPdbData.set(2)
        protBm2.inputStructure.set(cls.protImportPdb2.outputPdb)
        protBm2.setObjLabel('Biomol_1ake_pointer_uniteChains_False')
        cls.launchProtocol(protBm2)

        numStructs = len(protBm2.outputStructures)
        cls.assertTrue(numStructs == 2, "Failed to extract 2 biomols from 1ake (no dimer)")

        struct1 = protBm2.outputStructures.getFirstItem()
        nResidues = struct1.getAttributeValue(N_RESIDUES)
        nChains = struct1.getAttributeValue(N_CHAINS)
        cls.assertTrue(nResidues == 456,
                       "1ake biomol 1 should have 456 residues, not {0}".format(nResidues))
        cls.assertTrue(nChains == 3,
                       "1ake biomol 1 should have 3 chains with uniteChains False, not {0}".format(nChains))

        # extract biomols from 1ake (2 monomers) from pointer with uniteChains True
        protBm2b = cls.newProtocol(ProDyBiomol, uniteChains=True)
        protBm2b.inputPdbData.set(2)
        protBm2b.inputStructure.set(cls.protImportPdb2.outputPdb)
        protBm2b.setObjLabel('Biomol_1ake_pointer_uniteChains_True')
        cls.launchProtocol(protBm2b)

        numStructs = len(protBm2b.outputStructures)
        cls.assertTrue(numStructs == 2, "Failed to extract 2 biomols from 1ake (no dimer)")

        struct1 = protBm2b.outputStructures.getFirstItem()
        nResidues = struct1.getAttributeValue(N_RESIDUES)
        nChains = struct1.getAttributeValue(N_CHAINS)
        cls.assertTrue(nResidues == 456,
                       "1ake biomol 1 should have 456 residues, not {0}".format(nResidues))
        cls.assertTrue(nChains == 1,
                       "1ake biomol 1 should have 1 chains with uniteChains True, not {0}".format(nChains))

    def testProDyRenumberAll(cls):
        """ Run different selection options and confirm if it works """

        # ----------------------------------------------------------------------
        # Step 1a. Renumber imported selected 4akeA_ca to add 100
        # ----------------------------------------------------------------------
        protRenum = cls.newProtocol(ProDyRenumber, selection='all', offset=100)
        protRenum.inputStructure.set(cls.protSel.outputStructure)
        protRenum.setObjLabel('Renum_all_4akeA_ca_100')
        cls.launchProtocol(protRenum)

        outputFilename = renumFilename
        cls.assertTrue(exists(protRenum._getPath(outputFilename)))
        cls.assertTrue(hasattr(protRenum, "outputStructure"))

        struct1 = protRenum.outputStructure
        cls.assertTrue(struct1.getAttributeValue(FIRST_RESNUM) == 101,
                        "renumbered 4ake should have first resnum 101, not {0}".format(
                            struct1.getAttributeValue(FIRST_RESNUM)))
        cls.assertTrue(struct1.getAttributeValue(LAST_RESNUM) == 314,
                        "renumbered 4ake should have last resnum 314, not {0}".format(
                            struct1.getAttributeValue(LAST_RESNUM)))

    def testProDyRenumberSome(cls):
        """ Run different selection options and confirm if it works """

        # ----------------------------------------------------------------------
        # Step 1a. Renumber imported selected 4akeA_ca to add 100
        # ----------------------------------------------------------------------
        protRenum = cls.newProtocol(ProDyRenumber, selection='resnum 100 to 150', offset=1000)
        protRenum.inputStructure.set(cls.protSel.outputStructure)
        protRenum.setObjLabel('Renum_some_4akeA_ca_1000')
        cls.launchProtocol(protRenum)

        outputFilename = renumFilename
        cls.assertTrue(exists(protRenum._getPath(outputFilename)))
        cls.assertTrue(hasattr(protRenum, "outputStructure"))

        struct1 = protRenum.outputStructure

        cls.assertTrue(struct1.getAttributeValue(N_RESIDUES) == 214,
                       "Partially renumbered 4ake should still have 214 residues, not {0}".format(
                           struct1.getAttributeValue(N_RESIDUES)))
        cls.assertTrue(struct1.getAttributeValue(N_CHAINS) == 1,
                       "Partially renumbered 4ake should still have 1 chain, not {0}".format(
                           struct1.getAttributeValue(N_CHAINS)))

        # check that most of the structure stays the same
        cls.assertTrue(struct1.getAttributeValue(FIRST_RESNUM) == 1,
                        "Partially renumbered 4ake should have first resnum 1, not {0}".format(
                            struct1.getAttributeValue(FIRST_RESNUM)))
        cls.assertTrue(struct1.getAttributeValue(LAST_RESNUM) == 214,
                        "Partially renumbered 4ake should have last resnum 214, not {0}".format(
                            struct1.getAttributeValue(LAST_RESNUM)))

        # check that the part got renumbered
        cls.assertTrue(struct1.getAttributeValue(MAX_RESNUM) == 1150,
                        "Partially renumbered 4ake should have max resnum 1150, not {0}".format(
                            struct1.getAttributeValue(MAX_RESNUM)))

    def testProDyRenumberSomeChid(cls):
        """ Run different selection options and confirm if it works """

        # ----------------------------------------------------------------------
        # Step 1a. Renumber imported selected 4akeA_ca to add 100
        # ----------------------------------------------------------------------
        protRenum = cls.newProtocol(ProDyRenumber, selection='resnum 100 to 150',
                                    offset=1000, chain='B')
        protRenum.inputStructure.set(cls.protSel.outputStructure)
        protRenum.setObjLabel('Renum_some_4akeA_ca_1000-B')
        cls.launchProtocol(protRenum)

        outputFilename = renumFilename
        cls.assertTrue(exists(protRenum._getPath(outputFilename)))
        cls.assertTrue(hasattr(protRenum, "outputStructure"))

        struct1 = protRenum.outputStructure

        cls.assertTrue(struct1.getAttributeValue(N_RESIDUES) == 214,
                       "Partially renumbered 4ake should still have 214 residues, not {0}".format(
                           struct1.getAttributeValue(N_RESIDUES)))

        # check that most of the structure stays the same
        cls.assertTrue(struct1.getAttributeValue(FIRST_RESNUM) == 1,
                        "Partially renumbered 4ake should have first resnum 1, not {0}".format(
                            struct1.getAttributeValue(FIRST_RESNUM)))
        cls.assertTrue(struct1.getAttributeValue(LAST_RESNUM) == 214,
                        "Partially renumbered 4ake should have last resnum 214, not {0}".format(
                            struct1.getAttributeValue(LAST_RESNUM)))

        # check that the part got renumbered
        cls.assertTrue(struct1.getAttributeValue(MAX_RESNUM) == 1150,
                        "Partially renumbered 4ake should have max resnum 1150, not {0}".format(
                            struct1.getAttributeValue(MAX_RESNUM)))
        cls.assertTrue(struct1.getAttributeValue(N_CHAINS) == 2,
                       "Partially renumbered and rechained 4ake should now have 2 chain, not {0}".format(
                           struct1.getAttributeValue(N_CHAINS)))

class TestProDyCompareModes(TestWorkflow):
    """ Test protocol for comparing modes. """

    @classmethod
    def setUpClass(cls):
        # Create a new project
        setupTestProject(cls)
        importSelect4ake(cls)
        importANM2(cls)
        importDefvec(cls)

    def testProDyCompareANMvsDefvec(cls):
        # Compare original CA ANM NMA to defvec with default overlaps
        protComp5 = cls.newProtocol(ProDyCompare)
        protComp5.modes1.set(cls.protANM2.outputModes)
        protComp5.modes2.set(cls.protDefvec1.outputModes)
        protComp5.setObjLabel('Compare_ANM_to_Defvec_overlap')
        cls.launchProtocol(protComp5)

        matrix = np.loadtxt(protComp5.matrixFile.getFileName())
        cls.assertTrue(matrix.shape == (14,)) # excluding 6 zero modes

    def testProDyCompareANMvsANMdefaultOverlap(cls):
        # Compare original CA ANM NMA to itself with default overlaps
        protComp5 = cls.newProtocol(ProDyCompare)
        protComp5.modes1.set(cls.protANM2.outputModes)
        protComp5.modes2.set(cls.protANM2.outputModes)
        protComp5.setObjLabel('Compare_ANM_to_ANM_overlap')
        cls.launchProtocol(protComp5)

        matrix = np.loadtxt(protComp5.matrixFile.getFileName())
        cls.assertTrue(matrix.shape == (14,14)) # excluding 6 zero modes

    def testProDyCompareANMvsANMcovOverlap(cls):
        # Compare original CA ANM NMA to itself with default overlaps
        protComp5 = cls.newProtocol(ProDyCompare)
        protComp5.modes1.set(cls.protANM2.outputModes)
        protComp5.modes2.set(cls.protANM2.outputModes)
        protComp5.setObjLabel('Compare_ANM_to_ANM_covOverlap')
        protComp5.metric.set(1)
        cls.launchProtocol(protComp5)

        matrix = np.loadtxt(protComp5.matrixFile.getFileName())
        cls.assertTrue(matrix.size == 1) # cov overlap collapses

    def testProDyCompareANMvsANMrwsip(cls):
        # Compare original CA ANM NMA to itself with default overlaps
        protComp5 = cls.newProtocol(ProDyCompare)
        protComp5.modes1.set(cls.protANM2.outputModes)
        protComp5.modes2.set(cls.protANM2.outputModes)
        protComp5.setObjLabel('Compare_ANM_to_ANM_rwsip')
        protComp5.metric.set(2)
        cls.launchProtocol(protComp5)

        matrix = np.loadtxt(protComp5.matrixFile.getFileName())
        cls.assertTrue(matrix.size == 1) # rwsip collapses

    def testProDyCompareANMvsANMdefaultOverlapMatch(cls):
        # Compare original CA ANM NMA to itself with default overlaps
        protComp5 = cls.newProtocol(ProDyCompare)
        protComp5.modes1.set(cls.protANM2.outputModes)
        protComp5.modes2.set(cls.protANM2.outputModes)
        protComp5.setObjLabel('Compare_ANM_to_ANM_overlap_match')
        cls.launchProtocol(protComp5)

        matrix = np.loadtxt(protComp5.matrixFile.getFileName())
        cls.assertTrue(matrix.shape == (14,14)) # excluding 6 zero modes

    def testProDyCompareANMvsANMcovOverlapMatch(cls):
        # Compare original CA ANM NMA to itself with default overlaps
        protComp5 = cls.newProtocol(ProDyCompare, match=True)
        protComp5.modes1.set(cls.protANM2.outputModes)
        protComp5.modes2.set(cls.protANM2.outputModes)
        protComp5.setObjLabel('Compare_ANM_to_ANM_covOverlap_match')
        protComp5.metric.set(1)
        cls.launchProtocol(protComp5)

        matrix = np.loadtxt(protComp5.matrixFile.getFileName())
        cls.assertTrue(matrix.size == 1) # cov overlap collapses

    def testProDyCompareANMvsANMrwsipMatch(cls):
        # Compare original CA ANM NMA to itself with default overlaps
        protComp5 = cls.newProtocol(ProDyCompare, match=True)
        protComp5.modes1.set(cls.protANM2.outputModes)
        protComp5.modes2.set(cls.protANM2.outputModes)
        protComp5.setObjLabel('Compare_ANM_to_ANM_rwsip_match')
        protComp5.metric.set(2)
        cls.launchProtocol(protComp5)

        matrix = np.loadtxt(protComp5.matrixFile.getFileName())
        cls.assertTrue(matrix.size == 1) # rwsip collapses

def importSelect4ake(cls):
    cls.protSel = cls.newProtocol(ProDySelect, 
                                  selection="name CA and chain A",
                                  inputPdbData=0)
    cls.protSel.pdbId.set("4ake")
    cls.protSel.setObjLabel('sel_4akeA_ca')
    cls.launchProtocol(cls.protSel)

def importSelect1ake(cls):
    # Import a PDB
    cls.protImportPdb2 = cls.newProtocol(ProtImportPdb, inputPdbData=0,
                                         pdbId="1ake",
                                         skipChimera=True)
    cls.protImportPdb2.setObjLabel('pwem import 1ake')
    cls.launchProtocol(cls.protImportPdb2)

    # Select Chain A
    cls.protSel3 = cls.newProtocol(ProDySelect,
        selection="protein and chain A and name CA")
    cls.protSel3.inputStructure.set(cls.protImportPdb2.outputPdb)
    cls.protSel3.setObjLabel('Sel_1akeA_CA')
    cls.launchProtocol(cls.protSel3)

def importAligned1akeA(cls):
    # Import the already processed PDB
    cls.protImportPdb1akeA = cls.newProtocol(ProtImportPdb, inputPdbData=1,
                                             pdbFile=PRODY_TEST_ALG_PDB_FILE,
                                             skipChimera=True)
    cls.protImportPdb1akeA.setObjLabel('pwem import 1akeA_ca')
    cls.launchProtocol(cls.protImportPdb1akeA)

def importOnly4akeA(cls):
    # Import the already processed PDB
    cls.protImportPdb4akeA = cls.newProtocol(ProtImportPdb, inputPdbData=1,
                                             pdbFile=PRODY_TEST_TAR_PDB_FILE,
                                             skipChimera=True)
    cls.protImportPdb4akeA.setObjLabel('pwem import 4akeA_ca')
    cls.launchProtocol(cls.protImportPdb4akeA)

def importANM2(cls):
    # Import modes from prody npz
    cls.protANM2 = cls.newProtocol(ProDyImportModes)
    cls.protANM2.importType.set(MODES_NPZ)
    cls.protANM2.filesPath.set(TESTDIR)
    cls.protANM2.filesPattern.set("modes.anm.npz")
    cls.protANM2.inputStructure.set(cls.protSel.outputStructure)
    cls.protANM2.setObjLabel('import_npz_ANM_CA')
    cls.launchProtocol(cls.protANM2)

def importDefvec(cls):
    # Import modes from prody npz
    cls.protDefvec1 = cls.newProtocol(ProDyImportModes)
    cls.protDefvec1.importType.set(NMD)
    cls.protDefvec1.filesPath.set(TESTDIR)
    cls.protDefvec1.filesPattern.set("defvec.nmd")
    cls.protDefvec1.setObjLabel('import_nmd_Defvec_CA')
    cls.launchProtocol(cls.protDefvec1)
