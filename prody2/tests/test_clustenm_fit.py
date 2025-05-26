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

from pwem.protocols import ProtImportPdb, ProtImportVolumes
from pwem.tests.workflows import TestWorkflow
from pyworkflow.tests import setupTestProject

from prody2.protocols import ProDyClustENM
from prody2.constants import (PRODY_TEST_PDB_FILE,
                              PRODY_TEST_MRC_FILE,
                              ENSEMBLE_CCS)

import prody

class TestProDyClustenmFit(TestWorkflow):
    """ Test protocol for ProDy Normal Mode Analysis and Deformation Analysis. """

    @classmethod
    def setUpClass(cls):
        # Create a new project
        setupTestProject(cls)

        # Import starting structure
        cls.protPdb4ake = cls.newProtocol(ProtImportPdb, inputPdbData=1,
                                          pdbFile=PRODY_TEST_PDB_FILE)
        cls.protPdb4ake.setObjLabel('Input PDB')
        cls.launchProtocol(cls.protPdb4ake)

        # Import target EM map
        cls.protImportVol = cls.newProtocol(ProtImportVolumes,
                                            importFrom=ProtImportVolumes.IMPORT_FROM_FILES,
                                            filesPath=PRODY_TEST_MRC_FILE,
                                            samplingRate=2.0)
        cls.protImportVol.setObjLabel('EM map')
        cls.launchProtocol(cls.protImportVol)

    def testProDyClustENMFitting(cls):

        # Run ClustENM fitting in with replace filtered False (default)
        protClustenm3 = cls.newProtocol(ProDyClustENM, n_gens=1, numberOfModes=3,
                                        clusterMode=0, maxclust=2, rmsd=5,
                                        n_confs=10, sim=False, doFitting=True)
        protClustenm3.inputStructures.set([cls.protPdb4ake.outputPdb])
        protClustenm3.inputVolumes.set([cls.protImportVol.outputVolume])
        protClustenm3.setObjLabel('ClustENM_fitting_4akeA')
        cls.launchProtocol(protClustenm3)

        cc = [struct.getAttributeValue(ENSEMBLE_CCS)
              for struct in protClustenm3.outputStructures1]
        cls.assertTrue(cc[-1] > cc[0],
                       "Last CC should be more than starting CC when filtering and clustering")
        cls.assertTrue(len(cc) == 3,
                       "Number of structures should be 3 (1+2) when filtering and clustering to max 2")

    def testProDyClustENMFittingReplace(cls):

        # Run ClustENM fitting in with replace filtered True
        protClustenm3 = cls.newProtocol(ProDyClustENM, n_gens=1, numberOfModes=3,
                                        clusterMode=0, maxclust=10, rmsd=5,
                                        n_confs=10, sim=False, doFitting=True,
                                        replaceFiltered=True)
        protClustenm3.inputStructures.set([cls.protPdb4ake.outputPdb])
        protClustenm3.inputVolumes.set([cls.protImportVol.outputVolume])
        protClustenm3.setObjLabel('ClustENM_fitting_4akeA_replace')
        cls.launchProtocol(protClustenm3)

        cc = [struct.getAttributeValue(ENSEMBLE_CCS)
              for struct in protClustenm3.outputStructures1]
        cls.assertTrue(cc[-1] > cc[0],
                       "Last CC should be more than starting CC")
        cls.assertTrue(len(cc) == 11,
                       "Number of structures should be 11 (1+10) when filtering and replacing to 10")
