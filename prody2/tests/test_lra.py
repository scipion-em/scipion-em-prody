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

from pwem.protocols import ProtImportSetOfAtomStructs
from pwem.tests.workflows import TestWorkflow
from pyworkflow.tests import setupTestProject

from prody2.protocols import (ProDyBuildPDBEnsemble,
                              ProDyLRA, ProDyProject)

from prody2.protocols.protocol_project import ONE, TWO
from prody2.constants import CUSTOM

class TestProDyLRA(TestWorkflow):
    """ Test protocol for ProDy Ensemble and Logistic Regression Analysis"""

    @classmethod
    def setUpClass(cls):
        # Create a new project
        setupTestProject(cls)
        importStructs(cls)

    def testProDyLRA_1(cls):
        """ Run LRA simple workflows for two ways of building ensembles with A3 NTDs."""

        protEns1 = cls.newProtocol(ProDyBuildPDBEnsemble, refType=1,
                                   matchFunc=CUSTOM, chainOrders="""OrderedDict([((1, '7bnm'), 'BCA'),
                      ((2, '7bnn'), 'BCA'),
                      ((4, '7krq'), 'ABC'),
                      ((5, '7krr'), 'ABC')])""")
        protEns1.structures.set([cls.protImportPdbs.outputAtomStructs])
        protEns1.setObjLabel('buildPDBEns')
        cls.launchProtocol(protEns1)

        numAS = len(protEns1.outputNpz)
        cls.assertSetSize(protEns1.outputNpz, 4,
                           "wrong size SetOfAtomStructs ({0} not 4) with "
                           "SetOfAtomStructs input and added atom struct ref".format(numAS))

        protLra1 = cls.newProtocol(ProDyLRA, chainOrders="""OrderedDict([((1, '7bnm'), '1'),
                      ((2, '7krq'), '1'),
                      ((3, '7bnn'), '2'),
                      ((4, '7krr'), '2')])""")
        protLra1.inputEnsemble.set([protEns1.outputNpz])
        protLra1.setObjLabel('LRA_3d_vs_1u')
        cls.launchProtocol(protLra1)

        # ------------------------------------------------
        # Step 5. Project 1D and 2D (check it's still 1D)
        # ------------------------------------------------
        protProj1 = cls.newProtocol(ProDyProject)
        protProj1.inputEnsemble.set([protEns1.outputNpz])
        protProj1.inputModes.set(protLra1.outputModes)
        protProj1.numModes.set(ONE)
        protProj1.setObjLabel('Project 1D')
        cls.launchProtocol(protProj1)

        cls.assertTrue(hasattr(protProj1.outputEns1.getFirstItem(), "_prodyProjCoefficients"),
                        "1D Project protocol didn't add coefficients to SetOfAtomStructs")

        cls.assertEqual(len(protProj1.outputEns1.getFirstItem()._prodyProjCoefficients), 1,
                        "1D Project protocol didn't add 1 coefficient to first item")

        protProj2 = cls.newProtocol(ProDyProject,
                                     byFrame=True)
        protProj2.inputEnsemble.set([protEns1.outputNpz])
        protProj2.inputModes.set(protLra1.outputModes)
        protProj2.numModes.set(TWO)
        protProj2.setObjLabel('Project 2D')
        cls.launchProtocol(protProj2)

        cls.assertTrue(hasattr(protProj2.outputEns1.getFirstItem(), "_prodyProjCoefficients"),
                        "2D Project protocol didn't add coefficients to SetOfAtomStructs")

        cls.assertEqual(len(protProj2.outputEns1.getFirstItem()._prodyProjCoefficients), 1,
                        "2D Project protocol didn't add only 1 coefficient to first item for LRA")

def importStructs(cls):
    # ---------------------------------------------------------------
    # Step 1. Import some structures -> Select CA from all but one
    # --------------------------------------------------------------
    
    # Import PDB 3-down 7bnm and 7krq and 1-up 7bnn and 7krr
    cls.protImportPdbs = cls.newProtocol(ProtImportSetOfAtomStructs, inputPdbData=0,
                                         pdbIds="7bnm, 7krq, 7bnn, 7krr")
    cls.protImportPdbs.setObjLabel('pwem import Set AS')
    cls.launchProtocol(cls.protImportPdbs)
