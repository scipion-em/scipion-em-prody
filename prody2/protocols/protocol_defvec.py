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
from pyworkflow.protocol.params import (PointerParam, StringParam,
                                        FloatParam, IntParam, 
                                        BooleanParam, LEVEL_ADVANCED)

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

        form.addSection(label='Animation')        
        form.addParam('rmsd', FloatParam, default=0,
                      label='RMSD Amplitude (A)',
                      help='Used only for animations of computed normal modes. '
                      'This is the maximal amplitude with which atoms or pseudoatoms are moved '
                      'along the deformation vector in the animations. \n'
                      'The default value of 0 means use the actual RMSD between the structures.')
        form.addParam('n_steps', IntParam, default=10,
                      expertLevel=LEVEL_ADVANCED,
                      label='Number of frames',
                      help='Number of frames used in each direction of animations.')
        form.addParam('pos', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Include positive direction",
                      help='Elect whether to animate in the positive direction '
                           'from mobile to target.')
        form.addParam('neg', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Include negative direction",
                      help='Elect whether to animate in the negative direction, '
                           'extrapolating from mobile further away from target.')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self.mobFn = self.mobStructure.get().getFileName()
        self.tarFn = self.tarStructure.get().getFileName()
        self._insertFunctionStep('defvecStep', self.mobFn, self.tarFn)
        self._insertFunctionStep('animateModesStep', self.n_steps.get(),
                                 self.neg.get(), self.pos.get())
        self._insertFunctionStep('computeAtomShiftsStep')
        self._insertFunctionStep('createOutputStep')

    def defvecStep(self, mobFn, tarFn):
        self.mob = prody.parsePDB(mobFn, alt='all')
        self.tar = prody.parsePDB(tarFn, alt='all')

        if self.rmsd.get() == 0:
            self.rmsd = prody.calcRMSD(self.mob, self.tar)
        else:
            self.rmsd = self.rmsd.get()

        self.defvec = prody.calcDeformVector(self.mob, self.tar)

        self.outModes = prody.NMA('defvec')
        self.outModes.setEigens(self.defvec.getArray().reshape(-1, 1))
        prody.writeScipionModes(self._getPath(), self.outModes, write_star=True)
        prody.writeNMD(self._getPath('modes.nmd'), self.outModes, self.mob)

    def animateModesStep(self, n_steps, pos, neg):
        animations_dir = self._getExtraPath('animations')
        makePath(animations_dir)

        fnAnimation = join(animations_dir, "animated_mode_001")

        self.outAtoms = prody.traverseMode(self.defvec, self.mob, rmsd=self.rmsd,
                                           n_steps=self.n_steps.get(),
                                           pos=self.pos.get(), neg=self.neg.get())
        prody.writePDB(fnAnimation+".pdb", self.outAtoms)

        fhCmd=open(fnAnimation+".vmd",'w')
        fhCmd.write("mol new %s.pdb\n" % fnAnimation)
        fhCmd.write("animate style Rock\n")
        fhCmd.write("display projection Orthographic\n")
        fhCmd.write("mol modcolor 0 0 Index\n")
        if self.mob.ca.numAtoms() == self.mob.numAtoms():
            fhCmd.write("mol modstyle 0 0 Beads 1.000000 8.000000\n")
            # fhCmd.write("mol modstyle 0 0 Beads 1.800000 6.000000 "
            #         "2.600000 0\n")
        else:
            fhCmd.write("mol modstyle 0 0 NewRibbons 1.800000 6.000000 "
                    "2.600000 0\n")
        fhCmd.write("animate speed 0.5\n")
        fhCmd.write("animate forward\n")
        fhCmd.close()   

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
        nmSet.setPdb(self.mobStructure.get())

        self._defineOutputs(outputModes=nmSet)
        self._defineSourceRelation(self.mobStructure, nmSet)

