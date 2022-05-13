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

from os.path import basename, split

from pwem.protocols import *
from pwem.tests.workflows import TestWorkflow
from pyworkflow.tests import setupTestProject

from prody2.protocols import ProDySelect, ProDyBuildPDBEnsemble, ProDyPCA

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
        
        # Import PDB 7kec
        protImportPdb1 = self.newProtocol(ProtImportPdb, inputPdbData=0,
                                          pdbId="7kec")
        protImportPdb1.setObjLabel('pwem import 7kec')
        self.launchProtocol(protImportPdb1)

        # Import PDB 7bnm
        protImportPdb2 = self.newProtocol(ProtImportPdb, inputPdbData=0,
                                          pdbId="7bnm")
        protImportPdb2.setObjLabel('pwem import 7bnm')
        self.launchProtocol(protImportPdb2)

        # Import PDB 7bnn
        protImportPdb3 = self.newProtocol(ProtImportPdb, inputPdbData=0,
                                          pdbId="7bnn")
        protImportPdb3.setObjLabel('pwem import 7bnn')
        self.launchProtocol(protImportPdb3)

        # Import PDB 7bno
        protImportPdb4 = self.newProtocol(ProtImportPdb, inputPdbData=0,
                                          pdbId="7bno")
        protImportPdb4.setObjLabel('pwem import 7bno')
        self.launchProtocol(protImportPdb4)

        # Select CA
        protSel1 = self.newProtocol(ProDySelect, selection="name CA")
        protSel1.inputStructure.set(protImportPdb1.outputPdb)
        protSel1.setObjLabel('Sel 7kec_ca')
        self.launchProtocol(protSel1)

        protSel2 = self.newProtocol(ProDySelect, selection="name CA")
        protSel2.inputStructure.set(protImportPdb2.outputPdb)
        protSel2.setObjLabel('Sel 7bnm_ca')
        self.launchProtocol(protSel2)

        protSel3 = self.newProtocol(ProDySelect, selection="name CA")
        protSel3.inputStructure.set(protImportPdb3.outputPdb)
        protSel3.setObjLabel('Sel 7bnn_ca')
        self.launchProtocol(protSel3)

        protSel4 = self.newProtocol(ProDySelect, selection="name CA")
        protSel4.inputStructure.set(protImportPdb4.outputPdb)
        protSel4.setObjLabel('Sel 7bno_ca')
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
        protEns1 = self.newProtocol(ProDyBuildPDBEnsemble, refType=0)
        protEns1.structures.set(protSetAS.outputAtomStructs)
        protEns1.refStructure.set(protSel4.outputStructure)
        protEns1.setObjLabel('buildPDBEns_ref_str_4')
        self.launchProtocol(protEns1)   

        protPca1 = self.newProtocol(ProDyPCA, numberOfModes=3)
        protPca1.inputEnsemble.set(protEns1.outputNpz)
        protPca1.setObjLabel('PCA_ref_str_4')
        self.launchProtocol(protPca1)  

        # --------------------------------------------------
        # Step 3a. buildPDBEns with index ref -> PCA 2
        # --------------------------------------------------
        protEns2 = self.newProtocol(ProDyBuildPDBEnsemble, refType=1)
        protEns2.structures.set(protSetAS.outputAtomStructs)
        protEns2.refIndex.set(4)
        protEns2.setObjLabel('buildPDBEns_ref_idx_4')
        self.launchProtocol(protEns2)   

        protPca2 = self.newProtocol(ProDyPCA, numberOfModes=3)
        protPca2.inputEnsemble.set(protEns2.outputNpz)
        protPca2.setObjLabel('PCA_ref_idx_4')
        self.launchProtocol(protPca2)  

        # ------------------------------------------------
        # Step 4. Compare
        # ------------------------------------------------    
        protComp1 = self.newProtocol(ProDyCompare)
        protComp1.modes1.set(protPca1.outputModes)
        protComp1.modes2.set(protPca2.outputModes)
        protComp1.setObjLabel('Compare')
        self.launchProtocol(protComp1)           
