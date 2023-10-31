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
from os.path import split

from pwem.objects import SetOfNormalModes, SetOfPrincipalComponents
from pwem.protocols import ProtImportPdb, ProtImportSetOfAtomStructs, exists
from pwem.tests.workflows import TestWorkflow
from pyworkflow.tests import setupTestProject

from prody2.protocols import (ProDySelect, ProDyBuildPDBEnsemble,
                              ProDyImportEnsemble, ProDyPCA, ProDyCompare,
                              ProDyProject, ProDyANM, ProDyImportModes, ProDyEdit,
                              ProDyMeasure)

from prody2.protocols.protocol_edit import NMA_SLICE
from prody2.protocols.protocol_project import ONE, TWO, THREE
from prody2.protocols.protocol_import import SCIPION

import prody

class TestProDyPCA(TestWorkflow):
    """ Test protocol for ProDy Ensemble and Principal Component Analysis"""

    @classmethod
    def setUpClass(cls):
        # Create a new project
        setupTestProject(cls)

    def testProDyPCA(cls):
        """ Run PCA simple workflow for two ways of building ensembles. """
        
        oldVerbosity = prody.confProDy("verbosity")
        oldSecondary = prody.confProDy("auto_secondary")

        # ---------------------------------------------------------------
        # Step 1. Import some structures -> Select CA from all but one
        # --------------------------------------------------------------
        
        # Import PDB 3o21
        protImportPdb1 = cls.newProtocol(ProtImportPdb, inputPdbData=0,
                                          pdbId="3o21")
        protImportPdb1.setObjLabel('pwem import 3o21')
        cls.launchProtocol(protImportPdb1)

        # Import PDB 6fpj
        protImportPdb2 = cls.newProtocol(ProtImportPdb, inputPdbData=0,
                                          pdbId="6fpj")
        protImportPdb2.setObjLabel('pwem import 6fpj')
        cls.launchProtocol(protImportPdb2)

        # Import PDB 6flr
        protImportPdb3 = cls.newProtocol(ProtImportPdb, inputPdbData=0,
                                          pdbId="6flr")
        protImportPdb3.setObjLabel('pwem import 6flr')
        cls.launchProtocol(protImportPdb3)

        # Import PDB 3o21
        protImportPdb4 = cls.newProtocol(ProtImportPdb, inputPdbData=0,
                                          pdbId="3p3w")
        protImportPdb4.setObjLabel('pwem import 3p3w')
        cls.launchProtocol(protImportPdb4)

        # Select dimer and CA
        protSel1 = cls.newProtocol(ProDySelect, selection="chain B D and name CA")
        protSel1.inputStructure.set(protImportPdb4.outputPdb)
        protSel1.setObjLabel('Sel 3p3w_BD_ca') # parallel
        cls.launchProtocol(protSel1)

        protSel2 = cls.newProtocol(ProDySelect, selection="chain A C and name CA")
        protSel2.inputStructure.set(protImportPdb2.outputPdb)
        protSel2.setObjLabel('Sel 6fpj_AC_ca') # displaced
        cls.launchProtocol(protSel2)

        protSel3 = cls.newProtocol(ProDySelect, selection="name CA")
        protSel3.inputStructure.set(protImportPdb3.outputPdb)
        protSel3.setObjLabel('Sel 6flr_ca') # open
        cls.launchProtocol(protSel3)

        # -----------------------------------------------------------
        # Step 2. Import set of atom structs from existing selections
        # -----------------------------------------------------------
        protSetAS = cls.newProtocol(ProtImportSetOfAtomStructs, inputPdbData=1)
        protSetAS.filesPath.set(split(protSel2._getPath())[0])
        protSetAS.filesPattern.set('*Select/*atoms.pdb')
        protSetAS.setObjLabel('pwem import SetOfAS')
        cls.launchProtocol(protSetAS)

        # -------------------------------------------------------------------------
        # Step 3a. last CA selection (3o21_CD) -> buildPDBEns from SetOfAtomStructs
        # with atom struct ref from last CA selection -> PCA 1
        # -------------------------------------------------------------------------

        protSel4 = cls.newProtocol(ProDySelect, selection="chain C D and name CA")
        protSel4.inputStructure.set(protImportPdb1.outputPdb)
        protSel4.setObjLabel('Sel 3o21_CD_ca') # intermediate
        cls.launchProtocol(protSel4)

        protEns1 = cls.newProtocol(ProDyBuildPDBEnsemble, refType=0,
                                    matchFunc=0)
        protEns1.structures.set([protSetAS.outputAtomStructs])
        protEns1.refStructure.set(protSel4.outputStructure)
        protEns1.setObjLabel('buildPDBEns_ref_3o21_CD')
        cls.launchProtocol(protEns1)

        numAS = len(protEns1.outputNpz)
        cls.assertSetSize(protEns1.outputNpz, 4,
                           "wrong size SetOfAtomStructs ({0} not 4) with "
                           "SetOfAtomStructs input and added atom struct ref".format(numAS))

        protPca1 = cls.newProtocol(ProDyPCA, numberOfModes=3)
        protPca1.inputEnsemble.set(protEns1.outputNpz)
        protPca1.setObjLabel('PCA_from_set_ref_3o21_CD')
        cls.launchProtocol(protPca1)


        # -------------------------------------------------------------------
        # Step 3b. buildPDBEns from SetOfAtomStructs with index ref -> PCA 2
        # -------------------------------------------------------------------

        ens1 = prody.loadEnsemble(protEns1._getPath("ensemble.ens.npz"))
        idx = ens1.getLabels().index("3o21") + 1
        
        protEns2 = cls.newProtocol(ProDyBuildPDBEnsemble, refType=1,
                                    matchFunc=0)
        protEns2.structures.set([protSetAS.outputAtomStructs])
        protEns2.refIndex.set(idx)
        protEns2.setObjLabel('buildPDBEns_set_ref_idx_{0}'.format(idx))
        cls.launchProtocol(protEns2)

        numAS = len(protEns2.outputNpz)
        cls.assertSetSize(protEns2.outputNpz, 3,
                           "wrong size SetOfAtomStructs ({0} not 3) with "
                           "SetOfAtomStructs input and index ref".format(numAS))

        protPca2 = cls.newProtocol(ProDyPCA, numberOfModes=2)
        protPca2.inputEnsemble.set(protEns2.outputNpz)
        protPca2.setObjLabel('PCA_from_set_ref_idx')
        cls.launchProtocol(protPca2)

        # -------------------------------------------------------------------
        # Step 3c. buildPDBEns from SetOfAtomStructs and last CA selection
        # with index ref --> PCA with 3 and 2 components
        # -------------------------------------------------------------------

        protEns3 = cls.newProtocol(ProDyBuildPDBEnsemble, refType=1,
                                    matchFunc=0)
        protEns3.structures.set([protSel4.outputStructure,
                                 protSetAS.outputAtomStructs])
        protEns3.refIndex.set(1)
        protEns3.setObjLabel('buildPDBEns_set_plus_sel_ref_idx_{0}'.format(idx))
        cls.launchProtocol(protEns3)

        numAS = len(protEns3.outputNpz)
        cls.assertSetSize(protEns3.outputNpz, 4,
                           "wrong size SetOfAtomStructs ({0} not 4) with "
                           "SetOfAtomStructs and AtomStruct input and index ref".format(numAS))

        protPca3 = cls.newProtocol(ProDyPCA, numberOfModes=3)
        protPca3.inputEnsemble.set(protEns3.outputNpz)
        protPca3.setObjLabel('PCA_3_from_set_plus_sel_ref_idx')
        cls.launchProtocol(protPca3)

        cls.assertSetSize(protPca3.outputModes, 3,
                           "wrong size SetOfPrincipalComponents ({0} not 3)".format(
                               len(protPca3.outputModes)))

        protPca4 = cls.newProtocol(ProDyPCA, numberOfModes=2)
        protPca4.inputEnsemble.set(protEns3.outputNpz)
        protPca4.setObjLabel('PCA_2_from_set_plus_sel_ref_idx')
        cls.launchProtocol(protPca4)

        cls.assertSetSize(protPca4.outputModes, 2,
                           "wrong size SetOfPrincipalComponents ({0} not 2)".format(
                               len(protPca4.outputModes)))

        # --------------------------------------------------
        # Step 3d. buildPDBEns from list of AtomStruct
        # objects with index ref
        # --------------------------------------------------

        protEns4 = cls.newProtocol(ProDyBuildPDBEnsemble, refType=1,
                                    matchFunc=0)
        protEns4.structures.set([protSel4.outputStructure,
                                 protSel2.outputStructure])
        protEns4.setObjLabel('buildPDBEns_2_structs_ref_idx_{0}'.format(idx))
        cls.launchProtocol(protEns4)

        numAS = len(protEns4.outputNpz)
        cls.assertSetSize(protEns4.outputNpz, 2,
                           "wrong size SetOfAtomStructs ({0} not 2) with "
                           "multiple AtomStructs input and index ref".format(numAS))

        # ------------------------------------------------
        # Step 4. Compare
        # ------------------------------------------------
        # compare equivalent modes
        protComp1 = cls.newProtocol(ProDyCompare,
                                     diag=True, match=True)
        protComp1.modes1.set(protPca1.outputModes)
        protComp1.modes2.set(protPca3.outputModes)
        protComp1.setObjLabel('Compare_A3_PCAs_same')
        cls.launchProtocol(protComp1)

        compMatrix1 = prody.parseArray(protComp1._getExtraPath('matrix.txt'))
        cls.assertTrue(np.allclose(compMatrix1, np.ones(3)), "The modes aren't identical")

        # compare modes from different size sets
        protComp2 = cls.newProtocol(ProDyCompare)
        protComp2.modes1.set(protPca1.outputModes)
        protComp2.modes2.set(protPca4.outputModes)
        protComp2.setObjLabel('Compare_A3_PCAs_2_vs_3')
        cls.launchProtocol(protComp2)

        compMatrix2 = prody.parseArray(protComp2._getExtraPath('matrix.txt'))
        cls.assertEqual(compMatrix2.shape, (3,2),
                         "Comparing 3 and 2 modes doesn't give the 3x2 matrix")

        # ------------------------------------------------
        # Step 5. Project 1D, 2D and 3D
        # ------------------------------------------------
        protProj1 = cls.newProtocol(ProDyProject)
        protProj1.inputEnsemble.set([protEns2.outputNpz])
        protProj1.inputModes.set(protPca2.outputModes)
        protProj1.numModes.set(ONE)
        protProj1.setObjLabel('Project 1D')
        cls.launchProtocol(protProj1)

        cls.assertTrue(hasattr(protProj1.outputEns1.getFirstItem(), "_prodyProjCoefficients"),
                        "1D Project protocol didn't add coefficients to SetOfAtomStructs")

        cls.assertEqual(len(protProj1.outputEns1.getFirstItem()._prodyProjCoefficients), 1,
                        "1D Project protocol didn't add 1 coefficient to first item")

        protProj2 = cls.newProtocol(ProDyProject,
                                     byFrame=True)
        protProj2.inputEnsemble.set([protEns2.outputNpz])
        protProj2.inputModes.set(protPca2.outputModes)
        protProj2.numModes.set(TWO)
        protProj2.setObjLabel('Project 2D')
        cls.launchProtocol(protProj2)

        cls.assertTrue(hasattr(protProj2.outputEns1.getFirstItem(), "_prodyProjCoefficients"),
                        "2D Project protocol didn't add coefficients to SetOfAtomStructs")

        cls.assertEqual(len(protProj2.outputEns1.getFirstItem()._prodyProjCoefficients), 2,
                        "2D Project protocol didn't add 2 coefficient to first item")

        protProj3 = cls.newProtocol(ProDyProject,
                                     byFrame=True)
        protProj3.inputEnsemble.set([protEns2.outputNpz])
        protProj3.inputModes.set(protPca2.outputModes)
        protProj3.numModes.set(THREE)
        protProj3.setObjLabel('Project 3D')
        cls.launchProtocol(protProj3)

        cls.assertTrue(hasattr(protProj3.outputEns1.getFirstItem(), "_prodyProjCoefficients"),
                        "3D Project protocol didn't add coefficients to SetOfAtomStructs")

        cls.assertEqual(len(protProj3.outputEns1.getFirstItem()._prodyProjCoefficients), 2,
                        "3D Project protocol didn't add 2 coefficients to first item when given 2 components")

        protProj4 = cls.newProtocol(ProDyProject,
                                     byFrame=True)
        protProj4.inputEnsemble.set([protEns3.outputNpz])
        protProj4.inputModes.set(protPca3.outputModes)
        protProj4.numModes.set(THREE)
        protProj4.setObjLabel('Project 3D')
        cls.launchProtocol(protProj4)

        cls.assertTrue(hasattr(protProj4.outputEns1.getFirstItem(), "_prodyProjCoefficients"),
                        "3D Project protocol didn't add coefficients to SetOfAtomStructs")

        cls.assertEqual(len(protProj4.outputEns1.getFirstItem()._prodyProjCoefficients), 3,
                        "3D Project protocol didn't add 3 coefficients to first item when given 3 components")

        protProj5 = cls.newProtocol(ProDyProject,
                                     byFrame=True)
        protProj5.inputEnsemble.set([protEns3.outputNpz])
        protProj5.inputModes.set(protPca3.outputModes)
        protProj5.modeList.set("2")
        protProj5.numModes.set(TWO)
        protProj5.setObjLabel('Project 2D sel')
        cls.launchProtocol(protProj5)

        cls.assertTrue(hasattr(protProj5.outputEns1.getFirstItem(), "_prodyProjCoefficients"),
                        "2D Project protocol didn't add coefficients to SetOfAtomStructs")

        cls.assertEqual(len(protProj5.outputEns1.getFirstItem()._prodyProjCoefficients), 2,
                        "2D Project protocol didn't add 2 coefficient to first item")


        # -------------------------------------------------------
        # Step 6. Import 2k39 NMR ensemble -> select N+CA -> PCA
        # -------------------------------------------------------
        protImportPdb4 = cls.newProtocol(ProtImportPdb, inputPdbData=0,
                                          pdbId="2k39")
        protImportPdb4.setObjLabel('pwem import 2k39')
        cls.launchProtocol(protImportPdb4)

        protSel5 = cls.newProtocol(ProDySelect, selection="name N CA")
        protSel5.inputStructure.set(protImportPdb4.outputPdb)
        protSel5.setObjLabel('Sel 2k39_n_ca')
        cls.launchProtocol(protSel5)

        protImportEns = cls.newProtocol(ProDyImportEnsemble, superpose=2) # iterpose
        protImportEns.filesPath.set(protSel5._getPath("2k39_atoms.pdb"))
        protImportEns.setObjLabel('prody import ens - 2k39_n_ca')
        cls.launchProtocol(protImportEns)

        protPca5 = cls.newProtocol(ProDyPCA)
        protPca5.inputEnsemble.set(protImportEns.outputNpz)
        protPca5.setObjLabel('PCA_2k39')
        cls.launchProtocol(protPca5)

        # ------------------------------------------------
        # Step 7. ANM -> Import ANM -> import PCA 
        # -> check types -> compare
        # ------------------------------------------------

        # Launch ANM NMA for selected atoms (2k39 CA)
        protANM1 = cls.newProtocol(ProDyANM)
        protANM1.inputStructure.set(protSel5.outputStructure)
        protANM1.setObjLabel('ANM_2k39_CA')
        cls.launchProtocol(protANM1)

        # Import scipion ANM modes
        protImportModes1 = cls.newProtocol(ProDyImportModes)
        protImportModes1.importType.set(SCIPION)
        protImportModes1.filesPath.set(protANM1.outputModes.getFileName())
        protImportModes1.inputStructure.set(protSel5.outputStructure)
        protImportModes1.setObjLabel('import_scipion_ANM_CA')
        cls.launchProtocol(protImportModes1)

        # Import scipion PCA modes
        protImportModes2 = cls.newProtocol(ProDyImportModes)
        protImportModes2.importType.set(SCIPION)
        protImportModes2.filesPath.set(protPca5.outputModes.getFileName())
        protImportModes2.inputStructure.set(protSel5.outputStructure)
        protImportModes2.setObjLabel('import_scipion_PCA_CA')
        cls.launchProtocol(protImportModes2)

        # Check types
        cls.assertTrue(isinstance(protImportModes1.outputModes, SetOfNormalModes),
                        "ANM modes should be parsed as a SetOfNormalModes, not {0}".format(
                            type(protImportModes1.outputModes)
                        )) 

        cls.assertTrue(isinstance(protImportModes2.outputModes, SetOfPrincipalComponents),
                        "PCA modes should be parsed as a SetOfPrincipalComponents, not {0}".format(
                            type(protImportModes2.outputModes)
                        )) 

        # compare
        protComp3 = cls.newProtocol(ProDyCompare)
        protComp3.modes1.set(protImportModes1.outputModes)
        protComp3.modes2.set(protImportModes2.outputModes)
        protComp3.setObjLabel('Compare_imported_2k39')
        cls.launchProtocol(protComp3)

        cls.assertTrue(prody.confProDy("verbosity") == oldVerbosity,
                        "prody verbosity changed")

        cls.assertTrue(prody.confProDy("auto_secondary") == oldSecondary,
                        "prody auto_secondary changed")

        # ------------------------------------------------
        # Step 8. Select chain B from PCA -> slice
        # ------------------------------------------------
        protSel6 = cls.newProtocol(ProDySelect, selection="chain C")
        protSel6.inputStructure.set(protPca4.refPdb)
        protSel6.setObjLabel('Sel ref_C')
        cls.launchProtocol(protSel6)

        cls.assertTrue(exists(protSel6._getPath("atoms_atoms.pdb")))

        protEdit1 = cls.newProtocol(ProDyEdit, edit=NMA_SLICE)
        protEdit1.modes.set(protPca4.outputModes)
        protEdit1.newNodes.set(protSel6.outputStructure)
        protEdit1.setObjLabel('Slice_to_C')
        cls.launchProtocol(protEdit1)

        cls.assertTrue(exists(protEdit1._getExtraPath("animations/animated_mode_001.pdb")),
                        'slicing a SetOfPrincipalComponents should create animations')
        cls.assertTrue(exists(protEdit1._getExtraPath("distanceProfiles/vec1.xmd")),
                        'slicing a SetOfPrincipalComponents should create distance profiles')
        
        # ------------------------------------------------
        # Step 9. Test measures
        # ------------------------------------------------
        measuresFilename = "measures_1.csv"

        protMeasure1 = cls.newProtocol(ProDyMeasure)
        protMeasure1.inputEnsemble.set([protEns1.outputNpz]) # default: distance
        protMeasure1.selection1.set('chain C and resnum 1 to 114 249 to 350')
        protMeasure1.selection2.set('chain C and resnum 117 to 243 355 to 400')
        protMeasure1.setObjLabel('measure cleft dist')
        cls.launchProtocol(protMeasure1)

        cls.assertTrue(exists(protMeasure1._getPath(measuresFilename)),
                        'measuring distances should create measures_1.csv')

        dists = prody.parseArray(protMeasure1._getPath(measuresFilename), delimiter=',')
        cls.assertEqual(np.round(dists[0], 1), 30.5,
                        'cleft distances should measure 30.5 A')
        
        
        protMeasure2 = cls.newProtocol(ProDyMeasure, measureType=1) # angle
        protMeasure2.inputEnsemble.set([protEns1.outputNpz])
        protMeasure2.selection1.set('chain C and resnum 117 to 243 355 to 400')
        protMeasure2.selection2.set('resnum 1 to 114 249 to 350')
        protMeasure2.selection3.set('chain D and resnum 117 to 243 355 to 400')
        protMeasure2.setObjLabel('measure LL angle')
        cls.launchProtocol(protMeasure2)

        angles = prody.parseArray(protMeasure2._getPath(measuresFilename), delimiter=',')
        cls.assertEqual(np.round(angles[0], 0), 58,
                        '1st LL angle should measure 58 degrees')


        protMeasure3 = cls.newProtocol(ProDyMeasure, measureType=2) # dihedral
        protMeasure3.inputEnsemble.set([protEns1.outputNpz])
        protMeasure3.selection1.set('chain C and resnum 117 to 243 355 to 400')
        protMeasure3.selection2.set('chain C and resnum 1 to 114 249 to 350')
        protMeasure3.selection3.set('chain D and resnum 1 to 114 249 to 350')
        protMeasure3.selection4.set('chain D and resnum 117 to 243 355 to 400')
        protMeasure3.setObjLabel('measure displacement dihedral')
        cls.launchProtocol(protMeasure3)

        dihedrals = prody.parseArray(protMeasure3._getPath(measuresFilename), delimiter=',')
        cls.assertEqual(np.round(dihedrals[0], 0), -18,
                        '1st displacement dihedral should measure -18 degrees')
