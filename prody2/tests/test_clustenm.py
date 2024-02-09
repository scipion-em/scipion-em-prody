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

from pwem.tests.workflows import TestWorkflow
from pyworkflow.tests import setupTestProject

from pwem.protocols import ProtImportPdb, ProtImportVolumes
from prody2.protocols import (ProDySelect, ProDyClustENM)

from prody.tests.datafiles import pathDatafile

class TestProDyClustENMsingle(TestWorkflow):
    @classmethod
    def setUpClass(cls):
        # Create a new project
        setupTestProject(cls)
        importSelect4akeA(cls)

    def testProDyClustENMsingle(cls):
        """Run ClustENM for chain A from 4ake to test single structure option"""
        protClustenm1 = cls.newProtocol(ProDyClustENM, n_gens=2,
                                        clusterMode=0, maxclust='(2, 3)',
                                        n_confs=5, sim=False, outlier=True)
        protClustenm1.inputStructures.set([cls.protSelA.outputStructure])
        protClustenm1.setObjLabel('ClustENM_4akeA')
        cls.launchProtocol(protClustenm1)


class TestProDyClustENMmulti(TestWorkflow):
    @classmethod
    def setUpClass(cls):
        # Create a new project
        setupTestProject(cls)
        importSelect4akeA(cls)
        importSelect1akeA(cls)

    def testProDyClustENMmulti(cls):
        """Run ClustENM for chain A from 4ake and 1ake to test multi input option"""
        protClustenm2 = cls.newProtocol(ProDyClustENM, n_gens=1,
                                        clusterMode=1, threshold='1.',
                                        n_confs=2, sim=False, outlier=True)
        protClustenm2.inputStructures.set([cls.protSelA.outputStructure,
                                           cls.protSelB.outputStructure])
        protClustenm2.setObjLabel('ClustENM_2_structs')
        cls.launchProtocol(protClustenm2)


class TestProDyClustenmFit(TestWorkflow):
    @classmethod
    def setUpClass(cls):
        # Create a new project
        setupTestProject(cls)
        importPdbVol(cls)

    def testProDyClustENMFitting(cls):
        protClustenm3 = cls.newProtocol(ProDyClustENM, n_gens=3, numberOfModes=32,
                                        clusterMode=1, threshold='1.5',
                                        n_confs=20, sim=False, doFitting=True)
        protClustenm3.inputStructures.set([cls.protPdb4ake.outputPdb])
        protClustenm3.inputVolumes.set([cls.protImportVol.outputVolume])
        protClustenm3.setObjLabel('ClustENM_fitting_4akeA')
        cls.launchProtocol(protClustenm3)


def importSelect4akeA(cls):
    cls.protSelA = cls.newProtocol(ProDySelect, selection="protein and chain A",
                               inputPdbData=0)
    cls.protSelA.pdbId.set("4ake")
    cls.protSelA.setObjLabel('Sel_4akeA')
    cls.launchProtocol(cls.protSelA)

def importSelect1akeA(cls):
    cls.protSelB = cls.newProtocol(ProDySelect, selection="protein and chain A",
                               inputPdbData=0)
    cls.protSelB.pdbId.set("1ake")
    cls.protSelB.setObjLabel('Sel_1akeA')
    cls.launchProtocol(cls.protSelB)

def importPdbVol(cls):
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
