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

from pwem.protocols import *
from pwem.tests.workflows import TestWorkflow
from pyworkflow.tests import setupTestProject

from prody2.protocols import (ProDySelect, ProDyAlign, ProDyANM, ProDyRTB,
                              ProDyDefvec, ProDyEdit, ProDyCompare, ProDyImportModes)

from prody2.protocols.protocol_edit import NMA_SLICE, NMA_REDUCE, NMA_EXTEND, NMA_INTERP
from prody2.protocols.protocol_rtb import BLOCKS_FROM_RES, BLOCKS_FROM_SECSTR
from prody2.protocols.protocol_import import NMD, NPZ, SCIPION, GROMACS

import prody

class TestProDy_1(TestWorkflow):
    """ Test protocol for ProDy Normal Mode Analysis and Deformation Analysis. """

    @classmethod
    def setUpClass(cls):
        # Create a new project
        setupTestProject(cls)

    def test_ProDy_core(self):
        """ Run NMA simple workflow for two Atomic structures. """
        # ------------------------------------------------
        # Step 1. Import a Pdb -> Select chain A -> NMA
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

        # Launch ANM NMA for chain A (all atoms)
        protANM1 = self.newProtocol(ProDyANM, cutoff=8)
        protANM1.inputStructure.set(protSel1.outputStructure)
        protANM1.setObjLabel('ANM_all')
        self.launchProtocol(protANM1)

        # ------------------------------------------------
        # Step 2. Select CA -> ANM NMA
        # ------------------------------------------------
        # Select Calpha atoms
        protSel2 = self.newProtocol(ProDySelect, selection="name CA")
        protSel2.inputStructure.set(protSel1.outputStructure)
        protSel2.setObjLabel('Sel_4akeA_CA')
        self.launchProtocol(protSel2)

        # Launch ANM NMA for selected atoms (CA)
        protANM2 = self.newProtocol(ProDyANM)
        protANM2.inputStructure.set(protSel2.outputStructure)
        protANM2.setObjLabel('ANM_CA')
        self.launchProtocol(protANM2)        

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

        comp5_matrix = prody.parseArray(protComp5._getExtraPath('matrix.txt'))
        self.assertTrue(max(comp5_matrix) <= 1, "Default defvec comparison didn't normalise")

        # Compare original CA NMA to defvec with raw overlaps
        protComp6 = self.newProtocol(ProDyCompare)
        protComp6.norm.set(False)
        protComp6.modes1.set(protANM2.outputModes)
        protComp6.modes2.set(protDefvec1.outputModes)
        protComp6.setObjLabel('Compare_ANM_to_Defvec_raw')
        self.launchProtocol(protComp6)  

        comp6_matrix = prody.parseArray(protComp6._getExtraPath('matrix.txt'))
        self.assertTrue(max(comp6_matrix) > 1, "Raw defvec comparison didn't generate large numbers")

        # ------------------------------------------------
        # Step 8. Import ANM & compare scipion vs prody npz
        # -> import -> import -> compare
        # ------------------------------------------------
        # Define path
        modes = protANM2.outputModes
        modes_path = os.path.dirname(os.path.dirname(modes[1].getModeFile()))

        # Import modes from prody npz
        protImportModes1 = self.newProtocol(ProDyImportModes)
        protImportModes1.importType.set(NPZ)
        protImportModes1.filesPath.set(modes_path)
        protImportModes1.filesPattern.set("modes.anm.npz")
        protImportModes1.inputStructure.set(protSel2.outputStructure)
        protImportModes1.setObjLabel('import_npz_ANM_CA')
        self.launchProtocol(protImportModes1)   

        # Import scipion modes
        protImportModes2 = self.newProtocol(ProDyImportModes)
        protImportModes2.importType.set(SCIPION)
        protImportModes2.filesPath.set(modes_path)
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
        