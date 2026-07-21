# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     James Krieger (jamesmkrieger@gmail.com)
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
This module will provide ProDy mode algebra tools (linear combinations).
"""
import os
import numpy as np

from pwem.objects import SetOfNormalModes, String, CsvList

from pyworkflow.protocol.params import (PointerParam, EnumParam, IntParam,
                                        StringParam)

from prody2.protocols.protocol_modes_base import ProDyModesBase
from prody2 import Plugin

COEFF_POINTER = 0
COEFF_STRING = 1

class ProDyAlgebra(ProDyModesBase):
    """
    Combines several normal modes into a single collective mode by applying
    user-defined coefficients.

    AI Generated:

    ProDy Algebra (ProDyAlgebra) — User Manual
        Overview

        The ProDy Algebra protocol creates a new motion mode by linearly
        combining modes from an existing set of normal modes or principal
        components. Instead of analyzing each mode independently, this
        protocol allows the user to build a custom collective deformation
        by assigning coefficients to selected components.

        In structural biology, this is useful when several modes contribute
        jointly to a biologically relevant conformational change. For example,
        a transition between open and closed states may not be explained by
        a single low-frequency mode, but rather by a weighted combination of
        several motions.

        Inputs and General Workflow

        The protocol requires a SetOfNormalModes as input. These modes may
        originate from:

        - A normal mode analysis of an atomic structure (PDB-based model)
        - A pseudoatomic model derived from an EM volume
        - A set of principal components obtained from previous analyses

        The protocol first reads both the normal modes and the associated
        atomic coordinates. These coordinates are necessary because the final
        output must remain associated with a structural model.

        Coefficient Definition

        The key operation is the assignment of coefficients that determine
        how strongly each mode contributes to the final combined vector.

        Two coefficient input methods are available:

        1. Pointer input
           Coefficients are read from an external file. This is useful when
           coefficients have been generated automatically, for example from
           overlap analysis or external numerical processing.

        2. String input
           Coefficients are typed manually by the user as a list of numeric
           values separated by commas or spaces.

        The order of coefficients follows the order of the modes in the input
        set. If more modes are available than coefficients, only the modes
        with assigned coefficients are used.

        Number of Components

        The parameter "Number of components" controls how many modes are
        included in the linear combination.

        - If set to -1, all provided coefficients are used.
        - If a positive value is given, only the first N modes are included.

        This is useful when the user wants to restrict the combination to
        only the most relevant low-frequency modes, which are often the most
        biologically meaningful.

        Mathematical Operation

        The protocol performs a simple weighted sum:

            combined_mode = c1*m1 + c2*m2 + ... + cn*mn

        where each coefficient c multiplies its corresponding mode m.

        The resulting vector is stored as a new normal mode containing a
        single eigendirection.

        Biological Interpretation

        The generated output should not be interpreted as a new independent
        normal mode in the strict physical sense. Instead, it represents a
        synthetic collective motion built from existing components.

        This makes the protocol especially valuable for:

        - Reconstructing experimentally observed conformational transitions
        - Building custom trajectories for visualization
        - Exploring hypotheses about coupled domain motions
        - Combining principal components into interpretable structural changes

        In practice, this protocol is often used to generate motions that
        better match experimental observations than any individual mode alone.

        Output Files

        After execution, the protocol produces:

        - A new SetOfNormalModes object containing the combined mode
        - A list of the coefficients that were actually used
        - A ProDy NMD file for visualization of the resulting motion

        The original atomic model is linked to the output so the combined
        motion can be visualized directly in molecular viewers.

        Practical Recommendations

        In most biological applications, low-frequency modes should be
        prioritized because they usually describe large-scale collective
        motions such as hinge bending, domain rearrangements, or breathing
        motions.

        Large coefficients can amplify unrealistic distortions, especially
        when many modes are combined. It is therefore good practice to begin
        with small coefficients and visually inspect the resulting motion.

        If the goal is to reproduce a known conformational change, overlap-
        derived coefficients are often more reliable than arbitrary manual
        values.

        Final Perspective

        For structural interpretation, ProDy Algebra is best understood as
        a flexible mode-combination tool rather than a strict physical
        normal mode calculation.

        Its strength lies in allowing users to construct biologically
        meaningful collective motions from existing dynamical components,
        making it especially useful for hypothesis-driven conformational
        analysis and visualization.
    """
    _label = 'Vector alegebra'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form, besidesAnimation=False):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label='ProDy linear algebra')
        form.addParam('modes', PointerParam, label='Input set of modes',
                      pointerClass='SetOfNormalModes',
                      help='The input modes can be a SetOfNormalModes '
                           'from an atomic model (true PDB) or a pseudoatomic model '
                           '(an EM volume compared into pseudoatoms)'
                           'or a SetOfPrincipalComponents.')

        form.addParam('coeffSource', EnumParam, choices=['Pointer', 'String'],
                      default=COEFF_STRING,
                      label='Type of edit',
                      help='Modes will be added together with these coefficients. If there are more modes '
                           'than coefficients then the remaining modes will be ignored.')

        form.addParam('coeffPointer', PointerParam,
                      label='Coefficients',
                      condition='coeffSource==%d' % COEFF_POINTER,
                      pointerClass='EMFile',
                      help='Atoms or pseudoatoms to use as new nodes.')   

        form.addParam('coeffString', StringParam, default='', 
                      condition='coeffSource==%d' % COEFF_STRING,
                      label='Coefficients',
                      help='Modes will be added together with these coefficients. If there are more modes '
                           'than coefficients then the remaining modes will be ignored.')

        form.addParam('numCoeffs', IntParam, default='-1',
                      label='Number of components',
                      help='This number of modes will be added together with coefficients. '
                           'The remaining modes will be ignored.')

        ProDyModesBase._defineParams(self, form, besidesAnimation=besidesAnimation)

    # --------------------------- STEPS functions ------------------------------
    # This is inherited from modes base protocol
    def _insertAllSteps(self):
        super(ProDyAlgebra, self)._insertAllSteps()

    def computeModesStep(self):
        if self.coeffSource == COEFF_STRING:
            sep = ''
            coeffsString = self.coeffString.get()
            if ',' in coeffsString:
                sep += ','
            if ' ' in coeffsString:
                sep += ' '

            coeffs = np.array(coeffsString.split(sep), dtype=float)
            coeffsFn = self._getExtraPath('coeffs.txt')
            np.savetxt(coeffsFn, coeffs)
        else:
            coeffsFn = self.coeffPointer.get().getFileName()
            coeffs = np.loadtxt(coeffsFn)

        self.modesFn = self.modes.get().getFileName()
        self.atomsFn = self.modes.get().getPdb().getFileName()

        numCoeffs = min(self.numCoeffs.get(), len(coeffs))
        if numCoeffs == -1:
            numCoeffs = len(coeffs)

        args = '--modesFn {0} --atomsFn {1} --coeffsFn {2} --numCoeffs {3} ' \
            '--folder {4} --nmdFileName {5} --npzFileName {6}'.format(
                self.modesFn, self.atomsFn, coeffsFn, numCoeffs,
                self._getPath(), self.getNmdFileName(), self.getNpzFileName())

        self.runJob(Plugin.getProgram('algebra.py', script=True), args)

        self.outModesFn = self._getPath("modes.sqlite")
        self.coeffs = CsvList()
        self.coeffs._convertValue(["{:18.15f}".format(x) for x in coeffs[:numCoeffs]])

    def createOutputStep(self):
        fnSqlite = self._getPath('modes.sqlite')

        nmSet = SetOfNormalModes(filename=fnSqlite)
        nmSet._nmdFileName = String(self.getNmdFileName())

        pdb = self.modes.get().getPdb()
        nmSet.setPdb(pdb)
        os.symlink(os.path.abspath(pdb.getFileName()), 
                   self._getPath('atoms.pdb'))

        self._defineOutputs(outputModes=nmSet, coeffs=self.coeffs)
        self._defineSourceRelation(pdb, nmSet)

    def _summary(self):
        if not hasattr(self, 'outputModes'):
            summ = ['Output modes not ready yet']
        elif len(self.coeffs) == 1:
            summ = ['*1* mode scaled with coefficient *{0}*'.format(self.coeffs)]
        else:
            summ = ['*{0}* modes added with coefficients *{1}*'.format(
                    len(self.coeffs), self.coeffs)]
        return summ

    def getNmdFileName(self):
        return self._getPath('modes.nmd')
    
    def getNpzFileName(self):
        return self._getPath('modes.nma.npz')
