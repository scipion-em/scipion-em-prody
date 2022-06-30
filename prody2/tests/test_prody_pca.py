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

import numpy as np
from os.path import basename, split

from pwem.objects import SetOfNormalModes, SetOfPrincipalComponents
from pwem.protocols import *
from pwem.tests.workflows import TestWorkflow
from pyworkflow.tests import setupTestProject

from prody2.protocols import (ProDySelect, ProDyBuildPDBEnsemble,
                              ProDyImportEnsemble, ProDyPCA, ProDyCompare,
                              ProDyProject, ProDyANM, ProDyImportModes, ProDyEdit)

from prody2.protocols.protocol_project import ONE, TWO, THREE
from prody2.protocols.protocol_import import SCIPION

import prody

class TestProDy_pca(TestWorkflow):
    """ Test protocol for ProDy Normal Mode Analysis and Deformation Analysis. """

    @classmethod
    def setUpClass(cls):
        # Create a new project
        setupTestProject(cls)

    def test_ProDy_pca(self):
        """ Run PCA simple workflow for two ways of building ensembles. """
        # ----------------------------------------------------
        # Step 1. Import some structures -> Select CA 
        # ---------------------------------------------------
        
        # Import PDB 3o21
        protImportPdb1 = self.newProtocol(ProtImportPdb, inputPdbData=0,
                                          pdbId="3o21")
        protImportPdb1.setObjLabel('pwem import 3o21')
        self.launchProtocol(protImportPdb1)

        # Import PDB 6fpj
        protImportPdb2 = self.newProtocol(ProtImportPdb, inputPdbData=0,
                                          pdbId="6fpj")
        protImportPdb2.setObjLabel('pwem import 6fpj')
        self.launchProtocol(protImportPdb2)

        # Import PDB 6flr
        protImportPdb3 = self.newProtocol(ProtImportPdb, inputPdbData=0,
                                          pdbId="6flr")
        protImportPdb3.setObjLabel('pwem import 6flr')
        self.launchProtocol(protImportPdb3)

        # Import PDB 3o21
        protImportPdb4 = self.newProtocol(ProtImportPdb, inputPdbData=0,
                                          pdbId="3p3w")
        protImportPdb4.setObjLabel('pwem import 3p3w')
        self.launchProtocol(protImportPdb4)

        # Select dimer and CA
        protSel1 = self.newProtocol(ProDySelect, selection="chain C D and name CA")
        protSel1.inputStructure.set(protImportPdb1.outputPdb)
        protSel1.setObjLabel('Sel 3o21_CD_ca') # parallel
        self.launchProtocol(protSel1)

        protSel2 = self.newProtocol(ProDySelect, selection="chain A C and name CA")
        protSel2.inputStructure.set(protImportPdb2.outputPdb)
        protSel2.setObjLabel('Sel 6fpj_AC_ca') # displaced
        self.launchProtocol(protSel2)

        protSel3 = self.newProtocol(ProDySelect, selection="name CA")
        protSel3.inputStructure.set(protImportPdb3.outputPdb)
        protSel3.setObjLabel('Sel 6flr_ca') # open
        self.launchProtocol(protSel3)

        protSel4 = self.newProtocol(ProDySelect, selection="chain B D and name CA")
        protSel4.inputStructure.set(protImportPdb4.outputPdb)
        protSel4.setObjLabel('Sel 3p3w_BD_ca') # intermediate
        self.launchProtocol(protSel4)

        # ------------------------------------------------
        # Step 2. Import set of structs
        # ------------------------------------------------
        protSetAS = self.newProtocol(ProtImportSetOfAtomStructs, inputPdbData=1)
        protSetAS.filesPath.set(split(protSel2._getPath())[0])
        protSetAS.filesPattern.set('*Select/*atoms.pdb')
        protSetAS.setObjLabel('pwem import SetOfAS')
        self.launchProtocol(protSetAS)   

        # --------------------------------------------------
        # Step 3a. buildPDBEns with atom struct ref -> PCA 1
        # --------------------------------------------------
        protEns1 = self.newProtocol(ProDyBuildPDBEnsemble, refType=0,
                                    matchFunc=0)
        protEns1.structures.set(protSetAS.outputAtomStructs)
        protEns1.refStructure.set(protSel1.outputStructure)
        protEns1.setObjLabel('buildPDBEns_ref_3o21_CD')
        self.launchProtocol(protEns1)   

        protPca1 = self.newProtocol(ProDyPCA, numberOfModes=3)
        protPca1.inputEnsemble.set(protEns1.outputNpz)
        protPca1.setObjLabel('PCA_ref_3o21_CD')
        self.launchProtocol(protPca1)  

        # --------------------------------------------------
        # Step 3a. buildPDBEns with index ref -> PCA 2
        # --------------------------------------------------

        ens1 = prody.loadEnsemble(protEns1._getPath("ensemble.ens.npz"))
        
        protEns2 = self.newProtocol(ProDyBuildPDBEnsemble, refType=1,
                                    matchFunc=0)
        protEns2.structures.set(protSetAS.outputAtomStructs)
        protEns2.refIndex.set(ens1.getLabels().index("3o21_atoms")+1)
        protEns2.setObjLabel('buildPDBEns_ref_idx_4')
        self.launchProtocol(protEns2)   

        protPca2 = self.newProtocol(ProDyPCA, numberOfModes=3)
        protPca2.inputEnsemble.set(protEns2.outputNpz)
        protPca2.setObjLabel('PCA_ref_idx_4')
        self.launchProtocol(protPca2)  

        # ------------------------------------------------
        # Step 4. Compare
        # ------------------------------------------------    
        protComp1 = self.newProtocol(ProDyCompare,
                                     diag=True, match=True)
        protComp1.modes1.set(protPca1.outputModes)
        protComp1.modes2.set(protPca2.outputModes)
        protComp1.setObjLabel('Compare_A3_PCAs')
        self.launchProtocol(protComp1)     

        comp_matrix = prody.parseArray(protComp1._getExtraPath('matrix.txt'))
        self.assertTrue(np.allclose(comp_matrix, np.ones(3)), "The modes aren't identical")      

        # ------------------------------------------------
        # Step 5. Project 1D, 2D and 3D
        # ------------------------------------------------    
        protProj1 = self.newProtocol(ProDyProject)
        protProj1.inputEnsemble.set(protEns2.outputNpz)
        protProj1.inputModes.set(protPca2.outputModes)
        protProj1.numModes.set(ONE)
        protProj1.setObjLabel('Project 1D')
        self.launchProtocol(protProj1)   

        protProj2 = self.newProtocol(ProDyProject,
                                      byFrame=True)
        protProj2.inputEnsemble.set(protEns2.outputNpz)
        protProj2.inputModes.set(protPca2.outputModes)
        protProj2.numModes.set(TWO)
        protProj2.setObjLabel('Project 2D')
        self.launchProtocol(protProj2)   

        protProj3 = self.newProtocol(ProDyProject,
                                      byFrame=True)
        protProj3.inputEnsemble.set(protEns2.outputNpz)
        protProj3.inputModes.set(protPca2.outputModes)
        protProj3.numModes.set(THREE)
        protProj3.setObjLabel('Project 3D')
        self.launchProtocol(protProj3)  

        # ------------------------------------------------
        # Step 6. Import 2k39 NMR ensemble -> select N+CA -> PCA
        # ------------------------------------------------  
        protImportPdb4 = self.newProtocol(ProtImportPdb, inputPdbData=0,
                                          pdbId="2k39")
        protImportPdb4.setObjLabel('pwem import 2k39')
        self.launchProtocol(protImportPdb4)

        protSel4 = self.newProtocol(ProDySelect, selection="name N CA")
        protSel4.inputStructure.set(protImportPdb4.outputPdb)
        protSel4.setObjLabel('Sel 2k39_n_ca')
        self.launchProtocol(protSel4)        

        protImportEns = self.newProtocol(ProDyImportEnsemble, superpose=2)
        protImportEns.filesPath.set(protSel4._getPath("2k39_atoms.pdb"))
        protImportEns.setObjLabel('prody import ens - 2k39_n_ca')
        self.launchProtocol(protImportEns)

        protPca3 = self.newProtocol(ProDyPCA)
        protPca3.inputEnsemble.set(protImportEns.outputNpz)
        protPca3.setObjLabel('PCA_2k39')
        self.launchProtocol(protPca3)  

        # ------------------------------------------------
        # Step 7. ANM -> Import ANM -> import PCA 
        # -> check types -> compare
        # ------------------------------------------------  

        # Launch ANM NMA for selected atoms (2k39 CA)
        protANM1 = self.newProtocol(ProDyANM)
        protANM1.inputStructure.set(protSel4.outputStructure)
        protANM1.setObjLabel('ANM_2k39_CA')
        self.launchProtocol(protANM1)   

        # Import scipion ANM modes
        protImportModes1 = self.newProtocol(ProDyImportModes)
        protImportModes1.importType.set(SCIPION)
        protImportModes1.filesPath.set(protANM1.outputModes.getFileName())
        protImportModes1.inputStructure.set(protSel4.outputStructure)
        protImportModes1.setObjLabel('import_scipion_ANM_CA')
        self.launchProtocol(protImportModes1)

        # Import scipion PCA modes
        protImportModes2 = self.newProtocol(ProDyImportModes)
        protImportModes2.importType.set(SCIPION)
        protImportModes2.filesPath.set(protPca3.outputModes.getFileName())
        protImportModes2.inputStructure.set(protSel4.outputStructure)
        protImportModes2.setObjLabel('import_scipion_PCA_CA')
        self.launchProtocol(protImportModes2)  

        # Check types
        self.assertTrue(isinstance(protImportModes1.outputModes, SetOfNormalModes), 
                        "ANM modes should be parsed as a SetOfNormalModes, not {0}".format(
                            type(protImportModes1.outputModes)
                        )) 

        self.assertTrue(isinstance(protImportModes2.outputModes, SetOfPrincipalComponents), 
                        "PCA modes should be parsed as a SetOfPrincipalComponents, not {0}".format(
                            type(protImportModes2.outputModes)
                        )) 

        # compare
        protComp2 = self.newProtocol(ProDyCompare,
                                     diag=True, match=True)
        protComp2.modes1.set(protImportModes1.outputModes)
        protComp2.modes2.set(protImportModes2.outputModes)
        protComp2.setObjLabel('Compare_imported_2k39')
        self.launchProtocol(protComp2)   
