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

from prody2.protocols import (ProDySelect, ProDyClustENM)

from prody.tests.datafiles import pathDatafile
import prody

class TestProDyClustenmFit(TestWorkflow):
    """ Test protocol for ProDy Normal Mode Analysis and Deformation Analysis. """

    @classmethod
    def setUpClass(cls):
        # Create a new project
        setupTestProject(cls)

    def testProDyClustENMFitting(cls):

        # Import starting structure
        cls.protPdb4ake = cls.newProtocol(ProtImportPdb, inputPdbData=1,
                                          pdbFile=pathDatafile('pdb4ake_fixed'))
        cls.protPdb4ake.setObjLabel('Input PDB')
        cls.launchProtocol(cls.protPdb4ake)

        # Import target EM map
        cls.protImportVol = cls.newProtocol(ProtImportVolumes, importFrom=ProtImportVolumes.IMPORT_FROM_FILES,
                                            filesPath=pathDatafile('mrc1ake'),  samplingRate=2.0)
        cls.protImportVol.setObjLabel('EM map')
        cls.launchProtocol(cls.protImportVol)

        # Run ClustENM fitting
        protClustenm3 = cls.newProtocol(ProDyClustENM, n_gens=3, numberOfModes=32,
                                        clusterMode=1, threshold='1.5',
                                        n_confs=20, sim=False, doFitting=True)
        protClustenm3.inputStructures.set([cls.protPdb4ake.outputPdb])
        protClustenm3.inputVolumes.set([cls.protImportVol.outputVolume])
        protClustenm3.setObjLabel('ClustENM_fitting_4akeA')
        cls.launchProtocol(protClustenm3)

        stderr = open(protClustenm3._getLogsPath('run.stderr'), 'r')
        lines = stderr.readlines()
        stderr.close()

        cc = []
        for line in lines:
            if line.find('CC') != -1:
                cc.append(float(line.split()[4]))

        cls.assertTrue(cc[-1] > cc[0],
                       "Best last CC should be more than starting CC")
