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


"""
This module will provide ProDy deformation vector analysis.
"""
from pwem import *
from pwem.emlib import MetaData, MDL_NMA_MODEFILE, MDL_NMA_ATOMSHIFT
from pwem.objects import AtomStruct, SetOfNormalModes, String
from pwem.protocols import EMProtocol

from pyworkflow.utils import *
from pyworkflow.protocol.params import PointerParam, StringParam

import prody

import math

class ProDyDefvec(EMProtocol):
    """
    This protocol will perform deformation vector analysis
    """
    _label = 'Deformation'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label='ProDy Defvec')

        form.addParam('mobStructure', PointerParam, label="Mobile structure",
                      important=True,
                      pointerClass='AtomStruct',
                      help='The structure to be moved can be an atomic model '
                           '(true PDB) or a pseudoatomic model\n'
                           '(an EM volume converted into pseudoatoms).'
                           'The two structures should have the same number of nodes.')

        form.addParam('tarStructure', PointerParam, label="Target structure",
                      important=True,
                      pointerClass='AtomStruct',
                      help='The target structure can be an atomic model '
                           '(true PDB) or a pseudoatomic model\n'
                           '(an EM volume converted into pseudoatoms)'
                           'The two structures should have the same number of nodes.')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self.mobFn = self.mobStructure.get().getFileName()
        self.tarFn = self.tarStructure.get().getFileName()
        self._insertFunctionStep('defvecStep', self.mobFn, self.tarFn)
        self._insertFunctionStep('computeAtomShiftsStep')
        self._insertFunctionStep('createOutputStep')

    def defvecStep(self, mobFn, tarFn):
        mob = prody.parsePDB(mobFn, alt='all')
        tar = prody.parsePDB(tarFn, alt='all')

        rmsd = prody.calcRMSD(mob, tar)
        defvec = prody.calcDeformVector(mob, tar)

        self.outAtoms = prody.traverseMode(defvec, mob, rmsd=rmsd, neg=False) # morph
        
        self.outModes = prody.NMA('defvec')
        self.outModes.setEigens(defvec.getArray().reshape(-1, 1))

        prody.writePDB(self._getPath('atoms.pdb'), self.outAtoms)
        prody.writeScipionModes(self._getPath(), self.outModes, write_star=True)
        prody.writeNMD(self._getPath('modes.nmd'), self.outModes, mob)

    def computeAtomShiftsStep(self):
        fnOutDir = self._getExtraPath("distanceProfiles")
        makePath(fnOutDir)
        maxShift=[]
        maxShiftMode=[]
        
        n = 1
        fnVec = self._getPath("modes", "vec.%d" % n)
        fhIn = open(fnVec)
        md = MetaData()
        for line in fhIn:
            x, y, z = map(float, line.split())
            d = math.sqrt(x*x+y*y+z*z)
            maxShift.append(d)
            maxShiftMode.append(1)
            md.setValue(MDL_NMA_ATOMSHIFT,d,md.addObject())
        md.write(join(fnOutDir,"vec%d.xmd" % n))
        fhIn.close()

        md = MetaData()
        for i, _ in enumerate(maxShift):
            objId = md.addObject()
            md.setValue(MDL_NMA_ATOMSHIFT, maxShift[i],objId)
            md.setValue(MDL_NMA_MODEFILE, fnVec, objId)
        md.write(self._getExtraPath('maxAtomShifts.xmd'))

    def createOutputStep(self):
        fnSqlite = self._getPath('modes.sqlite')
        nmSet = SetOfNormalModes(filename=fnSqlite)
        nmSet._nmdFileName = String(self._getPath('modes.nmd'))

        outputPdb = AtomStruct()
        outputPdb.setFileName(self._getPath('atoms.pdb'))

        inputPdb = AtomStruct()
        inputPdb.setFileName(self.mobFn)
        nmSet.setPdb(inputPdb.get())

        self._defineOutputs(outputModes=nmSet, outputMorph=outputPdb)
        self._defineSourceRelation(outputPdb, nmSet)

