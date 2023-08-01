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

from prody2.protocols import (ProDySelect, ProDyBuildPDBEnsemble, ProDyClustENM)

class TestProDyClustENM(TestWorkflow):
    """ Test protocol for ProDy Normal Mode Analysis and Deformation Analysis. """

    @classmethod
    def setUpClass(cls):
        # Create a new project
        setupTestProject(cls)

    def testProDyClustENM(self):
        """ Run ClustENM for one and both chains of 4ake, test single and multi input options """

        # --------------------------------------------------------------
        # Step 1. Import 4ake and 1ake and select chain A from each
        # --------------------------------------------------------------

        # Select Chain A
        protSelA = self.newProtocol(ProDySelect, selection="protein and chain A", 
                                    inputPdbData=0)
        protSelA.pdbId.set("4ake")
        protSelA.setObjLabel('Sel_4akeA')
        self.launchProtocol(protSelA)

        # Select Chain A
        protSelB = self.newProtocol(ProDySelect, selection="protein and chain A", 
                                    inputPdbData=0)
        protSelB.pdbId.set("1ake")
        protSelB.setObjLabel('Sel_1akeA')
        self.launchProtocol(protSelB)

        # --------------------------------------------------------------
        # Step 2. Build PDB ensemble to combine 
        # --------------------------------------------------------------
        protEns = self.newProtocol(ProDyBuildPDBEnsemble, refType=1,
                                    matchFunc=0, selstr='all', degeneracy=True)
        protEns.structures.set([protSelA.outputStructure,
                                protSelB.outputStructure])
        protEns.setObjLabel('buildPDBEns_2_structs')
        self.launchProtocol(protEns)

        # --------------------------------------------------------------
        # Step 3. ClustENM with 4ake A
        # --------------------------------------------------------------
        protClustenm1 = self.newProtocol(ProDyClustENM, n_gens=1, 
                                         clusterMode=1, threshold=1,
                                         n_confs=5)
        protClustenm1.inputStructure.set(protSelA.outputStructure)
        protClustenm1.setObjLabel('ClustENM_4akeA')
        self.launchProtocol(protClustenm1)

        # --------------------------------------------------------------
        # Step 4. ClustENM with ensemble
        # --------------------------------------------------------------
        protClustenm2 = self.newProtocol(ProDyClustENM, n_gens=1, 
                                         clusterMode=1, threshold=1,
                                         n_confs=5)
        protClustenm2.inputStructure.set(protEns.outputStructures)
        protClustenm2.setObjLabel('ClustENM_2_structs')
        self.launchProtocol(protClustenm2)