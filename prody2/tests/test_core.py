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
from prody2.protocols.protocol_import import NMD, MODES_NPZ, SCIPION, GROMACS

import prody

class TestProDyCore(TestWorkflow):
    """ Test protocol for ProDy Normal Mode Analysis and Deformation Analysis. """

    @classmethod
    def setUpClass(cls):
        # Create a new project
        setupTestProject(cls)

    def testProDyCore(self):
        """ Run NMA simple workflow for two Atomic structures. """

        oldVerbosity = prody.confProDy("verbosity")
        oldSecondary = prody.confProDy("auto_secondary")

        # --------------------------------------------------------------
        # Step 1a. Import a Pdb -> Select chain A from Pointer -> NMA
        # --------------------------------------------------------------
        # Import a PDB
        protImportPdb1 = self.newProtocol(ProtImportPdb, inputPdbData=0,
                                          pdbId="4ake")
        protImportPdb1.setObjLabel('pwem import 4ake')
        self.launchProtocol(protImportPdb1)

        chainAselstr = "protein and chain A"

        # Select Chain A
        protSel1 = self.newProtocol(ProDySelect, selection=chainAselstr)
        protSel1.inputStructure.set(protImportPdb1.outputPdb)
        protSel1.setObjLabel('Sel_4akeA_all_pointer')
        self.launchProtocol(protSel1)

        self.assertTrue(exists(protSel1._getPath("4ake_atoms.pdb")))
        self.assertTrue(hasattr(protSel1, "outputStructure"))

        # Select chain C to show it doesn't work
        protSel1a = self.newProtocol(ProDySelect, selection="chain C")
        protSel1a.inputStructure.set(protImportPdb1.outputPdb)
        protSel1a.setObjLabel('Sel 4ake_C_pointer')
        self.launchProtocol(protSel1a)

        self.assertFalse(exists(protSel1a._getPath("4ake_atoms.pdb")))
        self.assertFalse(hasattr(protSel1a, "outputStructure"))

        # Launch ANM NMA for chain A (all atoms)
        protANM1 = self.newProtocol(ProDyANM, cutoff=8)
        protANM1.inputStructure.set(protSel1.outputStructure)
        protANM1.setObjLabel('ANM_all')
        self.launchProtocol(protANM1)

        # ------------------------------------------------
        # Step 1b. Select chain A from Filename
        # ------------------------------------------------
        protSel1b = self.newProtocol(ProDySelect, selection=chainAselstr,
                                     inputPdbData=1)
        protSel1b.pdbFile.set(protImportPdb1.outputPdb.getFileName())
        protSel1b.setObjLabel('Sel_4akeA_all_file')
        self.launchProtocol(protSel1b)

        # ----------------------------------------------------
        # Step 1c. Select chain A from PDB id (difficult case)
        # ----------------------------------------------------
        protSel1c = self.newProtocol(ProDySelect, selection=chainAselstr,
                                     inputPdbData=0)
        protSel1c.pdbId.set("6xr8")
        protSel1c.setObjLabel('Sel_6xr8_A_all_id')
        self.launchProtocol(protSel1c)

        # ------------------------------------------------
        # Step 2. Select CA -> ANM NMA
        # ------------------------------------------------
        # Select Calpha atoms
        protSel2 = self.newProtocol(ProDySelect, selection="name CA")
        protSel2.inputStructure.set(protSel1.outputStructure)
        protSel2.setObjLabel('Sel_4akeA_CA')
        self.launchProtocol(protSel2)

        # Launch ANM NMA for selected atoms (CA) with zeros (default)
        protANM2 = self.newProtocol(ProDyANM)
        protANM2.inputStructure.set(protSel2.outputStructure)
        protANM2.setObjLabel('ANM_CA_z')
        self.launchProtocol(protANM2)

        self.assertFalse(exists(protANM1._getExtraPath("animations/animated_mode_001.pdb")))
        self.assertTrue(exists(protANM1._getExtraPath("animations/animated_mode_007.pdb")))

        self.assertFalse(exists(protANM1._getExtraPath("distanceProfiles/vec1.xmd")))
        self.assertTrue(exists(protANM1._getExtraPath("distanceProfiles/vec7.xmd")))

        # Launch ANM NMA for selected atoms (CA) without zeros
        protANM2b = self.newProtocol(ProDyANM)
        protANM2b.inputStructure.set(protSel2.outputStructure)
        protANM2b.zeros.set(False)
        protANM2b.setObjLabel('ANM_CA_n-z')
        self.launchProtocol(protANM2b)

        self.assertTrue(exists(protANM2b._getExtraPath("animations/animated_mode_001.pdb")))
        self.assertTrue(exists(protANM2b._getExtraPath("distanceProfiles/vec1.xmd")))

        # Launch ed-ENM NMA for selected atoms (CA) without zeros
        protANM3 = self.newProtocol(ProDyANM)
        protANM3.inputStructure.set(protSel2.outputStructure)
        protANM3.zeros.set(False)
        protANM3.gamma.set("GammaED")
        protANM3.cutoff.set("2.9 * math.log(214) - 2.9")
        protANM3.setObjLabel('edENM_CA_n-z')
        self.launchProtocol(protANM3)

        # ------------------------------------------------
        # Step 3. Slice -> Compare
        # ------------------------------------------------
        # Slice all-atom ANM NMA to CA
        protEdit1 = self.newProtocol(ProDyEdit, edit=NMA_SLICE)
        protEdit1.modes.set(protANM1.outputModes)
        protEdit1.newNodes.set(protSel2.outputStructure)
        protEdit1.setObjLabel('Slice_to_CA')
        self.launchProtocol(protEdit1)        

        # Compare sliced and original CA NMA
        protComp1 = self.newProtocol(ProDyCompare)
        protComp1.modes1.set(protANM2.outputModes)
        protComp1.modes2.set(protEdit1.outputModes)
        protComp1.setObjLabel('Compare_slcAA_to_CA')
        self.launchProtocol(protComp1)           

        # ------------------------------------------------
        # Step 4. Reduce -> Compare
        # ------------------------------------------------
        # Reduce all-atom ANM NMA to CA
        protEdit2 = self.newProtocol(ProDyEdit, edit=NMA_REDUCE)
        protEdit2.modes.set(protANM1.outputModes)
        protEdit2.newNodes.set(protSel2.outputStructure)
        protEdit2.setObjLabel('Reduce_to_CA')
        self.launchProtocol(protEdit2)        

        # Compare reduced and original CA NMA
        protComp2 = self.newProtocol(ProDyCompare)
        protComp2.modes1.set(protANM2.outputModes)
        protComp2.modes2.set(protEdit2.outputModes)
        protComp2.setObjLabel('Compare_redAA_to_CA')
        self.launchProtocol(protComp2)           

        # ------------------------------------------------
        # Step 5. Extend -> Compare
        # ------------------------------------------------
        # Extend CA NMA to all-atoms
        protEdit3 = self.newProtocol(ProDyEdit, edit=NMA_EXTEND)
        protEdit3.modes.set(protANM2.outputModes)
        protEdit3.newNodes.set(protSel1.outputStructure)
        protEdit3.setObjLabel('Extend_to_AA')
        self.launchProtocol(protEdit3)        

        # Compare original AA ANM NMA and extended CA ANM NMA
        # Test matching too
        protComp3 = self.newProtocol(ProDyCompare, match=True)
        protComp3.modes1.set(protANM1.outputModes)
        protComp3.modes2.set(protEdit3.outputModes)
        protComp3.setObjLabel('Compare_AA_to_extCA')
        self.launchProtocol(protComp3)           

        # ------------------------------------------------
        # Step 6. Interpolate -> Compare
        # ------------------------------------------------
        # Interpolate CA NMA to all-atoms
        protEdit4 = self.newProtocol(ProDyEdit, edit=NMA_INTERP)
        protEdit4.modes.set(protANM2.outputModes)
        protEdit4.newNodes.set(protSel1.outputStructure)
        protEdit4.setObjLabel('Interp_to_AA')
        self.launchProtocol(protEdit4)

        # Compare original AA ANM NMA and interpolated CA ANM NMA
        protComp4 = self.newProtocol(ProDyCompare)
        protComp4.modes1.set(protANM1.outputModes)
        protComp4.modes2.set(protEdit4.outputModes)
        protComp4.setObjLabel('Compare_AA_to_intCA')
        self.launchProtocol(protComp4)

        # ------------------------------------------------
        # Step 7. Import other Pdb -> Select chain A and CA
        # -> align -> defvec -> compare
        # ------------------------------------------------
        # Import a PDB
        protImportPdb2 = self.newProtocol(ProtImportPdb, inputPdbData=0,
                                         pdbId="1ake")
        protImportPdb2.setObjLabel('pwem import 1ake')
        self.launchProtocol(protImportPdb2)

        # Select Chain A
        protSel3 = self.newProtocol(ProDySelect,
                                    selection="protein and chain A and name CA")
        protSel3.inputStructure.set(protImportPdb2.outputPdb)
        protSel3.setObjLabel('Sel_1akeA_CA')
        self.launchProtocol(protSel3)

        # Align new pdb from protSel3 to match old pdb with NMA
        protAlign1 = self.newProtocol(ProDyAlign)
        protAlign1.mobStructure.set(protSel3.outputStructure)
        protAlign1.tarStructure.set(protSel2.outputStructure)
        protAlign1.setObjLabel('Align_1akeA_4akeA_CA')
        self.launchProtocol(protAlign1) 

        # Defvec from same starting point as NMA
        protDefvec1 = self.newProtocol(ProDyDefvec)
        protDefvec1.mobStructure.set(protAlign1.outputStructureTar)
        protDefvec1.tarStructure.set(protAlign1.outputStructureMob)
        protDefvec1.setObjLabel('Defvec_4akeA_1akeA_CA')
        self.launchProtocol(protDefvec1) 

        # Compare original CA NMA to defvec with default overlaps
        protComp5 = self.newProtocol(ProDyCompare)
        protComp5.modes1.set(protANM2.outputModes)
        protComp5.modes2.set(protDefvec1.outputModes)
        protComp5.setObjLabel('Compare_ANM_to_Defvec')
        self.launchProtocol(protComp5)  

        compMatrix5 = prody.parseArray(protComp5._getExtraPath('matrix.txt'))
        self.assertTrue(max(compMatrix5) <= 1, "Default defvec comparison didn't normalise")

        # Compare original CA NMA to defvec with raw overlaps
        protComp6 = self.newProtocol(ProDyCompare)
        protComp6.norm.set(False)
        protComp6.modes1.set(protANM2.outputModes)
        protComp6.modes2.set(protDefvec1.outputModes)
        protComp6.setObjLabel('Compare_ANM_to_Defvec_raw')
        self.launchProtocol(protComp6)  

        compMatrix6 = prody.parseArray(protComp6._getExtraPath('matrix.txt'))
        self.assertTrue(max(compMatrix6) > 1, "Raw defvec comparison didn't generate large numbers")

        # ------------------------------------------------
        # Step 8. Import ANM & compare scipion vs prody npz
        # -> import -> import -> compare
        # ------------------------------------------------
        # Define path
        modes = protANM2.outputModes
        modesPath = os.path.dirname(os.path.dirname(modes._getMapper().selectFirst().getModeFile()))

        # Import modes from prody npz
        protImportModes1 = self.newProtocol(ProDyImportModes)
        protImportModes1.importType.set(MODES_NPZ)
        protImportModes1.filesPath.set(modesPath)
        protImportModes1.filesPattern.set("modes.anm.npz")
        protImportModes1.inputStructure.set(protSel2.outputStructure)
        protImportModes1.setObjLabel('import_npz_ANM_CA')
        self.launchProtocol(protImportModes1)   

        # Import scipion modes
        protImportModes2 = self.newProtocol(ProDyImportModes)
        protImportModes2.importType.set(SCIPION)
        protImportModes2.filesPath.set(protANM2.outputModes.getFileName())
        protImportModes2.inputStructure.set(protSel2.outputStructure)
        protImportModes2.setObjLabel('import_scipion_ANM_CA')
        self.launchProtocol(protImportModes2)  

        # Compare two imported ANMs
        protComp6 = self.newProtocol(ProDyCompare)
        protComp6.modes1.set(protImportModes1.outputModes)
        protComp6.modes2.set(protImportModes2.outputModes)
        protComp6.setObjLabel('Compare_imported_ANMs')
        self.launchProtocol(protComp6)  

        # -------------------------------------------------------
        # Step 9. RTB in 2 ways -> Compare to each other and ANM
        # -------------------------------------------------------

        # Launch RTB NMA for selected atoms (CA) with 10 res per block
        protRTB1 = self.newProtocol(ProDyRTB, blockDef=BLOCKS_FROM_RES)
        protRTB1.inputStructure.set(protSel2.outputStructure)
        protRTB1.setObjLabel('RTB_CA_10_res')
        self.launchProtocol(protRTB1)

        # Launch RTB NMA for selected atoms (CA) with secstr block
        protRTB2 = self.newProtocol(ProDyRTB, blockDef=BLOCKS_FROM_SECSTR)
        protRTB2.inputStructure.set(protSel2.outputStructure)
        protRTB2.setObjLabel('RTB_CA_secstr')
        self.launchProtocol(protRTB2)

        # Compare RTB1 and RTB2
        protComp6 = self.newProtocol(ProDyCompare)
        protComp6.modes1.set(protRTB1.outputModes)
        protComp6.modes2.set(protRTB2.outputModes)
        protComp6.setObjLabel('Compare_RTB1_to_RTB2')
        self.launchProtocol(protComp6)

        # Compare CA ANM and RTB1
        protComp7 = self.newProtocol(ProDyCompare)
        protComp7.modes1.set(protANM2.outputModes)
        protComp7.modes2.set(protRTB1.outputModes)
        protComp7.setObjLabel('Compare_ANM_to_RTB1')
        self.launchProtocol(protComp7)

        # Compare CA ANM and RTB2
        protComp8 = self.newProtocol(ProDyCompare)
        protComp8.modes1.set(protANM2.outputModes)
        protComp8.modes2.set(protRTB2.outputModes)
        protComp8.setObjLabel('Compare_ANM_to_RTB2')
        self.launchProtocol(protComp8)

        self.assertTrue(prody.confProDy("verbosity") == oldVerbosity, 
                        "prody verbosity changed")
        
        self.assertTrue(prody.confProDy("auto_secondary") == oldSecondary, 
                        "prody auto_secondary changed")


        # extract biomol from 4ake (dimer) from id
        protBm1 = self.newProtocol(ProDyBiomol)
        protBm1.inputPdbData.set(0)
        protBm1.pdbId.set("4ake")
        protBm1.setObjLabel('Biomol_4ake_id')
        self.launchProtocol(protBm1)

        numStructs = len(protBm1.outputStructures)
        self.assertTrue(numStructs == 1, "Failed to extract 1 biomol from 4ake (dimer)")

        ag = prody.parsePDB([struct.getFileName() for struct in protBm1.outputStructures])
        self.assertTrue(ag.numResidues() == 575, 
                        "4ake biomol 1 should have 575 residues, not {0}".format(ag.numResidues()))
        self.assertTrue(ag.numChains() == 2, 
                        "4ake biomol 1 should have 2 chains, not {0}".format(ag.numChains()))

        # extract biomols from 1ake (2 monomers) from pointer
        protBm2 = self.newProtocol(ProDyBiomol)
        protBm2.inputPdbData.set(2)
        protBm2.inputStructure.set(protImportPdb2.outputPdb)
        protBm2.setObjLabel('Biomol_1ake_pointer')
        self.launchProtocol(protBm2)

        numStructs = len(protBm2.outputStructures)
        self.assertTrue(numStructs == 2, "Failed to extract 2 biomols from 1ake (no dimer)")

        ag = prody.parsePDB([struct.getFileName() for struct in protBm2.outputStructures])[0]
        self.assertTrue(ag.numResidues() == 456, 
                        "1ake biomol 1 should have 456 residues, not {0}".format(ag.numResidues()))
        self.assertTrue(ag.numChains() == 1, 
                        "1ake biomol 1 should have 1 chain, not {0}".format(ag.numChains()))
