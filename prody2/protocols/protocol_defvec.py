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
from pwem.emlib import MetaData, MDL_NMA_MODEFILE, MDL_NMA_ATOMSHIFT
from pwem.objects import AtomStruct, SetOfNormalModes, String
from pwem.protocols import EMProtocol

from pyworkflow.utils import join, makePath
from pyworkflow.protocol.params import (PointerParam, FloatParam, IntParam, 
                                        BooleanParam, LEVEL_ADVANCED)

import prody
import math

class ProDyDefvec(EMProtocol):
    """
    Deformation Vector Analysis (ProDyDefvec) — User Manual

    Overview

    The ProDyDefvec protocol computes a deformation vector between two
    structures. It describes the coordinate displacement required to move
    a mobile structure toward a target structure.

    In structural biology, this protocol is useful when comparing two
    conformational states of the same macromolecule. Rather than describing
    motion through multiple normal modes, it directly captures the observed
    structural transition as a single deformation vector.

    This makes the protocol particularly valuable for studying domain
    rearrangements, hinge motions, ligand-induced changes, or conformational
    differences between experimentally determined structures.

    Inputs and General Workflow

    The protocol requires two structures:

    - A mobile structure, which serves as the starting conformation.
    - A target structure, which defines the destination conformation.

    Both structures may be true atomic models (PDB files) or pseudoatomic
    models derived from EM maps.

    A critical requirement is that both structures contain the same number
    of nodes in corresponding order. The protocol assumes a one-to-one
    positional correspondence between atoms or pseudoatoms.

    During execution, the protocol:

    1. Loads both structures.
    2. Computes the RMSD between them (unless the user explicitly provides
       a different amplitude for animation).
    3. Calculates the deformation vector connecting the mobile and target
       conformations.
    4. Stores the result as a single-mode NMA object.
    5. Generates animations and atom-shift profiles.

    Biological Interpretation

    Unlike classical normal mode analysis, where motion is predicted from
    intrinsic flexibility, this protocol measures an actual structural
    displacement observed between two conformations.

    Biologically, this allows direct interpretation of:

    - Which regions move the most
    - Whether motion is localized or distributed
    - Which domains behave as rigid bodies
    - Whether structural change is collective or highly local

    Large coherent displacements often indicate biologically meaningful
    transitions such as domain closure, subunit rotation, or large-scale
    functional rearrangements.

    RMSD Amplitude and Animation

    The RMSD parameter controls the amplitude used for animation.

    If RMSD is set to zero, the protocol automatically uses the true RMSD
    measured between the mobile and target structures.

    This is usually the most biologically meaningful setting because the
    animation then reflects the experimentally observed displacement.

    Alternatively, the user may provide a custom RMSD amplitude.

    This is useful when:

    - Exaggerating subtle motions for visualization
    - Normalizing amplitudes across multiple comparisons
    - Producing presentation-quality animations

    Animation Parameters

    The protocol generates an animation along the deformation vector.

    Important parameters include:

    - Number of frames:
      Controls smoothness of the animation. Higher values produce smoother
      transitions but require more storage.

    - Positive direction:
      Animates the motion from the mobile structure toward the target.

    - Negative direction:
      Extrapolates motion beyond the mobile structure, effectively extending
      the deformation past the starting conformation.

    For most biological applications, the positive direction is the most
    relevant because it represents the actual observed conformational change.

    The negative direction can be useful for exploratory analysis, especially
    when investigating whether the detected motion belongs to a larger
    structural trajectory.

    Structural Representation in Visualization

    The protocol automatically writes VMD visualization commands.

    It adapts the display style depending on the type of atoms detected:

    - If the structure mainly contains CA atoms or phosphate atoms,
      a bead representation is used.

    - Otherwise, ribbon-like representations are generated.

    This provides immediate qualitative visualization of the conformational
    transition without additional manual setup.

    Atom Shift Profiles

    A particularly informative output is the per-atom displacement profile.

    For every atom (or pseudoatom), the protocol computes the magnitude of
    displacement induced by the deformation vector.

    Biologically, this helps identify:

    - Flexible loops
    - Hinge regions
    - Stable structural cores
    - Domains undergoing the largest displacement

    These profiles are often more informative than global RMSD values,
    because they reveal where structural change is concentrated.

    Outputs and Their Interpretation

    The protocol produces:

    - A SetOfNormalModes object containing a single deformation mode
    - An NMD file for visualization
    - Animated PDB trajectories
    - Per-atom shift metadata
    - Maximum displacement profiles

    The resulting deformation mode can be interpreted exactly like a
    normal mode in downstream ProDy workflows, although biologically
    it represents an observed transition rather than a predicted intrinsic
    fluctuation.

    Practical Recommendations

    For best results:

    - Ensure both structures are already aligned before running the protocol.
      Misalignment will strongly distort the deformation vector.

    - Use structures with identical atom correspondence.
      Missing residues or inconsistent atom ordering can invalidate results.

    - When comparing different conformational states, inspect atom-shift
      profiles rather than relying only on RMSD.

    - Use the automatically estimated RMSD for biologically realistic
      animations.

    In practice, this protocol is especially useful when studying
    experimentally observed structural transitions where direct geometric
    interpretation is more important than harmonic flexibility analysis.

    Final Perspective

    For structural biologists, ProDyDefvec provides a direct way to
    characterize conformational transitions.

    Rather than asking how a molecule could move, this protocol asks how
    it did move between two experimentally observed states.

    That distinction makes it particularly powerful for mechanistic
    interpretation, visualization of functional transitions, and the
    identification of structurally important moving regions.
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
        self._insertFunctionStep('defvecStep')
        self._insertFunctionStep('animateModesStep', self.rmsd.get(), self.n_steps.get(),
                                 self.neg.get(), self.pos.get())
        self._insertFunctionStep('computeAtomShiftsStep')
        self._insertFunctionStep('createOutputStep')

    def defvecStep(self):
        mobStruct = self.mobStructure.get()
        self.mobFn = mobStruct.getFileName()

        tarStruct = self.tarStructure.get()
        self.tarFn = tarStruct.getFileName()

        self.mob = prody.parsePDB(self.mobFn, alt='all')
        self.tar = prody.parsePDB(self.tarFn, alt='all')

        if self.rmsd.get() == 0:
            self.rmsd = prody.calcRMSD(self.mob, self.tar)
        else:
            self.rmsd = self.rmsd.get()

        self.defvec = prody.calcDeformVector(self.mob, self.tar)

        self.outModes = prody.NMA('defvec')
        self.outModes.setEigens(self.defvec.getArray().reshape(-1, 1))
        prody.writeScipionModes(self._getPath(), self.outModes, write_star=True)
        prody.writeNMD(self._getPath('modes.nmd'), self.outModes, self.mob)

    def animateModesStep(self, rmsd, nSteps, pos, neg):
        animationsDir = self._getExtraPath('animations')
        makePath(animationsDir)

        fnAnimation = join(animationsDir, "animated_mode_001")

        self.outAtoms = prody.traverseMode(self.defvec, self.mob, rmsd=rmsd,
                                           n_steps=nSteps,
                                           pos=pos, neg=neg)
        prody.writePDB(fnAnimation+".pdb", self.outAtoms)

        fhCmd=open(fnAnimation+".vmd",'w')
        fhCmd.write("mol new %s.pdb\n" % fnAnimation)
        fhCmd.write("animate style Rock\n")
        fhCmd.write("display projection Orthographic\n")
        fhCmd.write("mol modcolor 0 0 Index\n")

        if self.mob.select('name P') is not None:
            numAtomsP = self.mob.select('name P').numAtoms()
        else:
            numAtomsP = 0

        if self.mob.ca is not None:
            numAtomsCA = self.mob.ca.numAtoms()
        else:
            numAtomsCA = 0

        numAtomsRep = numAtomsCA + numAtomsP
        if numAtomsRep == self.mob.numAtoms():
            fhCmd.write("mol modstyle 0 0 Beads 2.000000 8.000000\n")
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

