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
# *  e-mail address 'you@yourinstitution.email'
# *
# **************************************************************************


"""
This module will provide ProDy normal mode analysis using the anisotropic network model (ANM).
"""
from pyworkflow.protocol import Protocol, params, Integer
from pwem.protocols import EMProtocol

import os
from os.path import basename, exists, join
import math
import numpy as np

from pwem import *
from pwem.convert.atom_struct import cifToPdb
from pwem.emlib import MetaData, MDL_NMA_ATOMSHIFT, MDL_NMA_MODEFILE
from pyworkflow.utils import *
from pyworkflow.utils.path import copyFile, createLink, makePath, cleanPath, moveFile
from pyworkflow.protocol.params import (PointerParam, IntParam, FloatParam, StringParam,
                                        BooleanParam, LEVEL_ADVANCED)
from pwem.objects import NormalMode, SetOfNormalModes

from xmipp3.base import XmippMdRow
import prody

class prodyANM(EMProtocol):
    """
    This protocol will perform normal mode analysis using the anisotropic network model (ANM)
    """
    _label = 'ANM analysis'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label='ProDy ANM NMA')

        form.addParam('inputStructure', PointerParam, label="Input structure",
                      important=True,
                      pointerClass='AtomStruct',
                      help='The input structure can be an atomic model '
                           '(true PDB) or a pseudoatomic model\n'
                           '(an EM volume converted into pseudoatoms)')

        form.addParam('numberOfModes', IntParam, default=20,
                      label='Number of modes',
                      help='The maximum number of modes allowed by the method for '
                           'atomic normal mode analysis is 6 times the number of '
                           'RTB blocks and for pseudoatomic normal mode analysis 3 '
                           'times the number of pseudoatoms. However, the protocol '
                           'allows only up to 200 modes as 20-100 modes are usually '
                           'enough. The number of modes given here should be below '
                           'the minimum between these two numbers.')

        form.addParam('cutoff', FloatParam, default=15.,
                      expertLevel=LEVEL_ADVANCED,
                      label="Cut-off distance (A)",
                      help='Atoms or pseudoatoms beyond this distance will not interact. \n'
                           'For Calpha atoms, the distance of 15 Angstroms works well in the majority of cases. \n'
                           'For pseudoatoms, please set this according to the level of coarse-graining '
                           'following Doruker et al., J Comput Chem 2002.')

        form.addParam('gamma', StringParam, default=1.,
                      expertLevel=LEVEL_ADVANCED,
                      label="Spring constant",
                      help='This determines the strength of the strings and can be a function. \n'
                           'This allows the use of ENMs based on pseudoatom, atom or residue properties including radii, '
                           'secondary structure, distance-dependent potentials, and other custom functions. \n'
                           'See http://prody.csb.pitt.edu/tutorials/enm_analysis/gamma.html')

        form.addParam('sparse', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Use sparse matrices",
                      help='Elect whether to use Scipy sparse matrices. \n'
                           'This allows efficient usage of memory at the cost of computation speed')

        form.addParam('kdtree', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Use KDTree for building Hessian matrix",
                      help='Elect whether to use KDTree for building Hessian matrix. This is usually slower.')

        form.addParam('zeros', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="include zero eigvals",
                      help='Elect whether modes with zero eigenvalues will be kept.')

        form.addParam('turbo', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Use turbo mode",
                      help='Elect whether to use a memory intensive, but faster way to calculate modes.')


    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps

        # Link the input
        inputFn = self.inputStructure.get().getFileName()

        # Compute modes
        self._insertFunctionStep('computeModesStep', inputFn)
        self._insertFunctionStep('createOutputStep')

    def computeModesStep(self, inputFn):
        self.structureEM = False #self.inputStructure.get().getPseudoAtoms()

        ag = prody.parsePDB(inputFn)
        
        if self.structureEM:
            self.atoms = ag
        else:
            self.atoms = ag.ca

        self.anm = prody.ANM(ag.getTitle())
        self.anm.buildHessian(self.atoms, cutoff=self.cutoff.get(),
                              gamma=eval(self.gamma.get()), 
                              sparse=self.sparse.get(), 
                              kdtree=self.kdtree.get())
        self.anm.calcModes(n_modes=self.numberOfModes.get(), 
                           zeros=self.zeros.get(),
                           turbo=self.turbo.get())

        prody.saveModel(self.anm, self._getExtraPath('{0}.anm.npz'.format(basename(inputFn))))
        prody.writeNMD(self._getExtraPath('{0}.anm.nmd'.format(basename(inputFn))), self.anm, self.atoms)

        prody.writeCFlexModes(self._getPath(), self.anm)

    def createOutputStep(self):
        """Copied from continuous_flex protocol_nma"""
        fnSqlite = self._getPath('modes.sqlite')
        nmSet = SetOfNormalModes(filename=fnSqlite)

        md = MetaData(self._getPath('modes.xmd'))
        row = XmippMdRow()
        
        for objId in md:
            row.readFromMd(md, objId)
            nmSet.append(rowToMode(row))
        inputPdb = self.inputStructure.get()
        nmSet.setPdb(inputPdb)
        self._defineOutputs(outputModes=nmSet)
        self._defineSourceRelation(self.inputStructure, nmSet)

