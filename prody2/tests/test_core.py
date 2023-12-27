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

import math
import os

from pwem.protocols import ProtImportPdb, exists
from pwem.tests.workflows import TestWorkflow
from pyworkflow.tests import setupTestProject

from prody2.protocols import (ProDySelect, ProDyAlign, ProDyBiomol, ProDyANM, ProDyRTB,
                              ProDyDefvec, ProDyEdit, ProDyCompare, ProDyImportModes)

from prody2.protocols.protocol_edit import NMA_SLICE, NMA_REDUCE, NMA_EXTEND, NMA_INTERP
from prody2.protocols.protocol_rtb import BLOCKS_FROM_RES, BLOCKS_FROM_SECSTR
from prody2.protocols.protocol_import import MODES_NPZ, SCIPION

import prody
from prody.tests.datafiles import pathDatafile

oldVerbosity = prody.confProDy("verbosity")
oldSecondary = prody.confProDy("auto_secondary")

class TestProDyCore(TestWorkflow):
    """ Test protocol for ProDy Normal Mode Analysis and Deformation Analysis. """

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
                                          pdbId="4ake")
        protImportPdb1.setObjLabel('pwem import 4ake')
        cls.launchProtocol(protImportPdb1)

        chainAselstr = "protein and chain A"
        outputFilename = "4ake_atoms.pdb"

        # Select Chain A
        protSel1 = cls.newProtocol(ProDySelect, selection=chainAselstr)
        protSel1.inputStructure.set(protImportPdb1.outputPdb)
        protSel1.setObjLabel('Sel_4akeA_all_pointer')
        cls.launchProtocol(protSel1)

        cls.assertTrue(exists(protSel1._getPath(outputFilename)))
        cls.assertTrue(hasattr(protSel1, "outputStructure"))

        # Select chain C and name CA to show it doesn't work
        protSel1a = cls.newProtocol(ProDySelect, selection="chain C and name CA")
        protSel1a.inputStructure.set(protImportPdb1.outputPdb)
        protSel1a.setObjLabel('Sel 4ake_C_pointer_1')
        cls.launchProtocol(protSel1a)

        cls.assertFalse(exists(protSel1a._getPath(outputFilename)))
        cls.assertFalse(hasattr(protSel1a, "outputStructure"))

        # Select chain C whole to show it does work with uniteChains False (default)
        protSel1a2 = cls.newProtocol(ProDySelect, selection="chain C")
        protSel1a2.inputStructure.set(protImportPdb1.outputPdb)
        protSel1a2.setObjLabel('Sel 4ake_C_pointer_2')
        cls.launchProtocol(protSel1a2)

        cls.assertTrue(exists(protSel1a2._getPath(outputFilename)))
        cls.assertTrue(hasattr(protSel1a2, "outputStructure"))

        # Select chain C whole to show it does work with uniteChains True
        protSel1a3 = cls.newProtocol(ProDySelect, selection="chain C",
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
                                         pdbFile=pathDatafile("pdb4ake_fixed"))
        protImportPdb1.setObjLabel('pwem import 4ake')
        cls.launchProtocol(protImportPdb1)

        # Launch ANM NMA for chain A (all atoms)
        protANM1 = cls.newProtocol(ProDyANM, cutoff=8)
        protANM1.inputStructure.set(protImportPdb1.outputPdb)
        protANM1.setObjLabel('ANM_all')
        cls.launchProtocol(protANM1)

        # ------------------------------------------------
        # Step 2. CA ANM NMA
        # ------------------------------------------------

        animationsFile7 = "animations/animated_mode_007.pdb"

        # Launch ANM NMA for selected atoms (CA) with zeros (default)
        protANM2 = cls.newProtocol(ProDyANM)
        protANM2.inputStructure.set(cls.protSel.outputStructure)
        protANM2.setObjLabel('ANM_CA_z')
        cls.launchProtocol(protANM2)

        cls.assertFalse(exists(protANM1._getExtraPath("animations/animated_mode_001.pdb")))
        cls.assertTrue(exists(protANM1._getExtraPath(animationsFile7)))

        cls.assertFalse(exists(protANM1._getExtraPath("distanceProfiles/vec1.xmd")))
        cls.assertTrue(exists(protANM1._getExtraPath("distanceProfiles/vec7.xmd")))

        # Launch ANM NMA for selected atoms (CA) without zeros
        protANM2b = cls.newProtocol(ProDyANM)
        protANM2b.inputStructure.set(cls.protSel.outputStructure)
        protANM2b.zeros.set(False)
        protANM2b.setObjLabel('ANM_CA_n-z')
        cls.launchProtocol(protANM2b)

        cls.assertTrue(exists(protANM2b._getExtraPath("animations/animated_mode_001.pdb")))
        cls.assertTrue(exists(protANM2b._getExtraPath("distanceProfiles/vec1.xmd")))

        # Launch ed-ENM NMA for selected atoms (CA) without zeros
        protANM3 = cls.newProtocol(ProDyANM)
        protANM3.inputStructure.set(cls.protSel.outputStructure)
        protANM3.zeros.set(False)
        protANM3.gamma.set("GammaED")
        protANM3.cutoff.set("2.9 * math.log(214) - 2.9")
        protANM3.setObjLabel('edENM_CA_n-z')
        cls.launchProtocol(protANM3)

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

        compMatrix5 = prody.parseArray(protComp5._getExtraPath('matrix.txt'))
        cls.assertTrue(max(compMatrix5) <= 1, "Default defvec comparison didn't normalise")

        # Compare original CA NMA to defvec with raw overlaps
        protComp6 = cls.newProtocol(ProDyCompare)
        protComp6.norm.set(False)
        protComp6.modes1.set(protANM2.outputModes)
        protComp6.modes2.set(protDefvec1.outputModes)
        protComp6.setObjLabel('Compare_ANM_to_Defvec_raw')
        cls.launchProtocol(protComp6)

        compMatrix6 = prody.parseArray(protComp6._getExtraPath('matrix.txt'))
        cls.assertTrue(max(compMatrix6) > 1, "Raw defvec comparison didn't generate large numbers")

        # ------------------------------------------------
        # Step 8. Import ANM & compare scipion vs prody npz
        # -> import -> import -> compare
        # ------------------------------------------------
        # Define path
        modes = protANM2.outputModes
        modesPath = os.path.dirname(os.path.dirname(modes._getMapper().selectFirst().getModeFile()))

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

    def testProDyRTB(cls):
        # -------------------------------------------------------
        # Step 9. RTB in 2 ways -> Compare to each other and ANM
        # -------------------------------------------------------

        # Launch RTB NMA for selected atoms (CA) with 10 res per block
        protRTB1 = cls.newProtocol(ProDyRTB, blockDef=BLOCKS_FROM_RES)
        protRTB1.inputStructure.set(cls.protSel.outputStructure)
        protRTB1.setObjLabel('RTB_CA_10_res')
        cls.launchProtocol(protRTB1)

        # Launch RTB NMA for selected atoms (CA) with secstr block
        protRTB2 = cls.newProtocol(ProDyRTB, blockDef=BLOCKS_FROM_SECSTR)
        protRTB2.inputStructure.set(cls.protSel.outputStructure)
        protRTB2.setObjLabel('RTB_CA_secstr')
        cls.launchProtocol(protRTB2)

        # Compare RTB1 and RTB2
        protComp6 = cls.newProtocol(ProDyCompare)
        protComp6.modes1.set(protRTB1.outputModes)
        protComp6.modes2.set(protRTB2.outputModes)
        protComp6.setObjLabel('Compare_RTB1_to_RTB2')
        cls.launchProtocol(protComp6)

        cls.assertTrue(prody.confProDy("verbosity") == oldVerbosity,
                        "prody verbosity changed")
        
        cls.assertTrue(prody.confProDy("auto_secondary") == oldSecondary,
                        "prody auto_secondary changed")

    def testProDyBiomol(cls):
        # extract biomol from 4ake (dimer) from id
        protBm1 = cls.newProtocol(ProDyBiomol)
        protBm1.inputPdbData.set(0)
        protBm1.pdbId.set("4ake")
        protBm1.setObjLabel('Biomol_4ake_id')
        cls.launchProtocol(protBm1)

        numStructs = len(protBm1.outputStructures)
        cls.assertTrue(numStructs == 1, "Failed to extract 1 biomol from 4ake (dimer)")

        ag = prody.parsePDB([struct.getFileName() for struct in protBm1.outputStructures])
        cls.assertTrue(ag.numResidues() == 575,
                        "4ake biomol 1 should have 575 residues, not {0}".format(ag.numResidues()))
        cls.assertTrue(ag.numChains() == 2,
                        "4ake biomol 1 should have 2 chains, not {0}".format(ag.numChains()))

        # extract biomols from 1ake (2 monomers) from pointer with uniteChains False (default)
        protBm2 = cls.newProtocol(ProDyBiomol)
        protBm2.inputPdbData.set(2)
        protBm2.inputStructure.set(cls.protImportPdb2.outputPdb)
        protBm2.setObjLabel('Biomol_1ake_pointer')
        cls.launchProtocol(protBm2)

        numStructs = len(protBm2.outputStructures)
        cls.assertTrue(numStructs == 2, "Failed to extract 2 biomols from 1ake (no dimer)")

        ag = prody.parsePDB([struct.getFileName() for struct in protBm2.outputStructures])[0]
        cls.assertTrue(ag.numResidues() == 456,
                       "1ake biomol 1 should have 456 residues, not {0}".format(ag.numResidues()))
        cls.assertTrue(ag.numChains() == 3,
                       "1ake biomol 1 should have 3 chains, not {0}".format(ag.numChains()))

        # extract biomols from 1ake (2 monomers) from pointer with uniteChains True
        protBm2b = cls.newProtocol(ProDyBiomol, uniteChains=True)
        protBm2b.inputPdbData.set(2)
        protBm2b.inputStructure.set(cls.protImportPdb2.outputPdb)
        protBm2b.setObjLabel('Biomol_1ake_pointer')
        cls.launchProtocol(protBm2b)

        numStructs = len(protBm2b.outputStructures)
        cls.assertTrue(numStructs == 2, "Failed to extract 2 biomols from 1ake (no dimer)")

        ag = prody.parsePDB([struct.getFileName() for struct in protBm2.outputStructures])[0]
        cls.assertTrue(ag.numResidues() == 456,
                       "1ake biomol 1 should have 456 residues, not {0}".format(ag.numResidues()))
        cls.assertTrue(ag.numChains() == 3,
                       "1ake biomol 1 should have 1 chains, not {0}".format(ag.numChains()))

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
                                        pdbId="1ake")
    cls.protImportPdb2.setObjLabel('pwem import 1ake')
    cls.launchProtocol(cls.protImportPdb2)

    # Select Chain A
    cls.protSel3 = cls.newProtocol(ProDySelect,
                                selection="protein and chain A and name CA")
    cls.protSel3.inputStructure.set(cls.protImportPdb2.outputPdb)
    cls.protSel3.setObjLabel('Sel_1akeA_CA')
    cls.launchProtocol(cls.protSel3)