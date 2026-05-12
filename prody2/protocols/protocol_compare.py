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
This module will provide ProDy normal mode analysis using the anisotropic network model (ANM).
"""
import os
import numpy as np

from pwem.objects import AtomStruct, EMFile, String
from pwem.protocols import EMProtocol

from pyworkflow.utils import glob
from pyworkflow.protocol.params import PointerParam, EnumParam, BooleanParam

import prody

NMA_METRIC_OVERLAP = 0
NMA_METRIC_COV_OVERLAP = 1
NMA_METRIC_RWSIP = 2

class ProDyCompare(EMProtocol):
    """
    Compares two sets of normal modes using overlap-based or
    ensemble-based similarity metrics.

    AI Generated:

    ProDy Compare (ProDyCompare) — User Manual
        Overview

        The ProDy Compare protocol evaluates the similarity between two
        sets of normal modes or principal components. Its main goal is to
        determine whether two dynamical descriptions capture comparable
        structural motions.

        In structural biology, this is especially useful when comparing
        dynamical behavior derived from different experimental conditions,
        different structural models, different force fields, or different
        computational methods. Rather than comparing static structures,
        this protocol compares the intrinsic motions encoded in mode sets.

        Inputs and General Workflow

        The protocol requires two input sets of normal modes.

        Each input may come from:

        - A normal mode analysis of an atomic structure
        - A pseudoatomic model derived from an EM map
        - A principal component analysis of structural ensembles

        In general, both mode sets should describe systems with the same
        number of nodes. The only exception is when one of the inputs
        contains exactly one mode, in which case single-mode comparison
        is allowed.

        The protocol reads both mode sets together with their associated
        structural coordinates when available.

        Comparison Metrics

        Three comparison metrics are available.

        Overlap

        Overlap is the most direct comparison. It measures the cosine
        similarity between individual modes.

        Biologically, overlap answers the question:

        “Does a motion from one model point in the same direction as a
        motion from another model?”

        Values close to 1 indicate highly similar motions.
        Values close to 0 indicate unrelated motions.
        Negative values indicate opposite directions.

        This metric is most useful when comparing individual modes or
        mode-by-mode correspondence.

        Covariance Overlap

        Covariance overlap compares whole mode subspaces rather than
        individual modes.

        It measures how similarly two mode sets describe collective
        fluctuations.

        This is often more biologically meaningful when the exact order
        of modes differs but the overall dynamical space remains similar.

        RWSIP

        The Root Weighted Square Inner Product (RWSIP) also compares
        dynamical subspaces.

        It provides a robust measure of global similarity between two
        sets of collective motions.

        Like covariance overlap, it is especially useful when comparing
        the overall flexibility landscape rather than direct mode pairing.

        Zero Modes and Physical Interpretation

        For covariance overlap and RWSIP calculations, the first six
        trivial rigid-body modes are excluded.

        These modes correspond to global translations and rotations and
        generally do not carry biologically relevant internal dynamics.

        This makes the comparison focus on genuine internal structural
        flexibility.

        Pairwise Overlap Options

        When the overlap metric is selected, additional controls are
        available.

        Diagonal Only

        If enabled, only corresponding mode pairs are compared.

        This is useful when the user already expects a one-to-one
        correspondence between modes.

        Normalization

        By default, overlaps are normalized.

        This produces cosine-like similarities independent of vector
        magnitude.

        If normalization is disabled, the protocol computes raw dot
        products instead. This can be useful in advanced analyses but is
        usually less intuitive biologically.

        Mode Matching

        For overlap and covariance overlap, the protocol can attempt to
        match modes before comparison.

        This is especially useful when the two mode sets contain similar
        motions but not in the same order.

        The protocol internally searches for the best correspondence
        between modes and produces:

        - A reordered matched mode set
        - A file listing the matched mode indices

        Biologically, this can help identify equivalent collective
        motions across different models or experimental states.

        Output Matrix

        The main output is a numerical matrix stored as a file.

        Its interpretation depends on the selected metric.

        For overlap:
        - The matrix contains pairwise mode-to-mode similarities.

        For covariance overlap or RWSIP:
        - The matrix contains similarity values computed over growing
          mode subsets.

        This means the protocol progressively evaluates how similarity
        evolves as more modes are included.

        Additional Outputs

        When mode matching is enabled, the protocol also generates:

        - A matched mode set
        - A file containing the matched mode indices

        These outputs are useful for downstream inspection, visualization,
        and interpretation of mode correspondence.

        Practical Recommendations

        For most biological applications, overlap is the best starting
        point when the user wants to compare specific individual modes.

        Covariance overlap and RWSIP are usually better when the goal is
        to compare the global dynamical behavior of two systems.

        Mode matching should generally be enabled when the two mode sets
        are expected to describe similar motions but may not preserve the
        same ordering.

        If the biological question concerns functional flexibility rather
        than exact mode identity, subspace-based metrics often provide a
        more robust interpretation.

        Final Perspective

        In biological terms, ProDy Compare does not ask whether two
        structures look similar.

        Instead, it asks whether they tend to move in similar ways.

        This makes it especially valuable for studying conserved
        flexibility, conformational transitions, and the dynamical
        consequences of mutations, ligand binding, or alternative
        structural models.
    """
    _label = 'Compare modes'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label='ProDy compare')

        form.addParam('modes1', PointerParam, label="Input modes set 1",
                      important=True,
                      pointerClass='SetOfNormalModes',
                      help='The input modes can be a SetOfNormalModes '
                           'from an atomic model (true PDB) or a pseudoatomic model '
                           '(an EM volume compared into pseudoatoms) '
                           'or a SetOfPrincipalComponents.\n'
                           'The two sets should have the same number of nodes '
                           'unless one of them has exactly 1 mode in it.')

        form.addParam('modes2', PointerParam, label="Input modes set 2",
                      important=True,
                      pointerClass='SetOfNormalModes',
                      help='The input modes can be a SetOfNormalModes '
                           'from an atomic model (true PDB) or a pseudoatomic model '
                           '(an EM volume compared into pseudoatoms) '
                           'or a SetOfPrincipalComponents.\n'
                           'The two sets should have the same number of nodes '
                           'unless one of them has exactly 1 mode in it.')

        form.addParam('metric', EnumParam, choices=['Overlap', 'Covariance Overlap', 'RWSIP'],
                      default=NMA_METRIC_OVERLAP,
                      label='Comparison metric',
                      help='Modes can be compared pairwise using correlation cosine overlaps (inner products), '
                      'or in sets using either covariance overlap (Hess, Phys Rev E 2002; aka spectral overlap) '
                      'or the root weighted square inner product (RWSIP; Carnevale et al., J Phys Condens Matter 2007). \n'
                      'Covariance overlaps and RWSIPs are calculated over growing mode sets.\n'
                      'The six zero modes are excluded from covariance overlap and RWSIP calculations.')

        form.addParam('diag', BooleanParam, default=False, 
                      condition='metric==%d' % NMA_METRIC_OVERLAP,
                      label='Calculate diagonal values only',
                      help='Elect whether to calculate diagonal values only.')      

        form.addParam('match', BooleanParam, default=False, 
                      condition='metric!=%d' % NMA_METRIC_RWSIP,
                      label='Match modes',
                      help='Elect whether to match modes.')     

        form.addParam('norm', BooleanParam, default=True, 
                      condition='metric==%d' % NMA_METRIC_OVERLAP,
                      label='Normalise overlaps',
                      help='Elect whether to normalise vectors for overlaps or calculate raw dot products.')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('compareModesStep')
        self._insertFunctionStep('createOutputStep')

    def compareModesStep(self):
        modesPath1 = os.path.dirname(os.path.dirname(
            self.modes1.get()._getMapper().selectFirst().getModeFile()))

        pdb1 = glob(modesPath1+"/*atoms.pdb")
        if len(pdb1) == 0:
            pdb1 = None

        modes1 = prody.parseScipionModes(self.modes1.get().getFileName(), pdb=pdb1)

        modesPath2 = os.path.dirname(os.path.dirname(
            self.modes2.get()._getMapper().selectFirst().getModeFile()))
            
        pdb2 = glob(modesPath2+"/*atoms.pdb")
        if len(pdb2) == 0:
            pdb2 = None

        modes2 = prody.parseScipionModes(self.modes2.get().getFileName(), pdb=pdb2)

        nModes = np.max([modes1.numModes(), modes2.numModes()])
        nModesMin = np.min([modes1.numModes(), modes2.numModes()])

        if nModesMin != 1 and self.match:
            modeEns = prody.ModeEnsemble()
            modeEns.addModeSet(modes1)
            modeEns.addModeSet(modes2)
            modeEns.match()

            matchInds = prody.matchModes(modes1, modes2, index=True)

            self.matchIndsFileName = prody.writeArray(self._getExtraPath('matchInds.txt'),
                                                      np.array(matchInds, dtype=int)[1]+1,
                                                      format='%3d')

            pdb = self.modes1.get().getPdb()
            if pdb is not None:
                atoms = prody.parsePDB(pdb.getFileName())
            else:
                atoms = prody.parsePDB(pdb1)

            typeStr = str(type(modes2)).lower().split('.')[-1].split("'")[0]
            self.nmdFileName = self._getPath('matched_modes.{0}.nmd'.format(typeStr))

            prody.writeNMD(self.nmdFileName, modeEns[1], atoms)
            prody.writeScipionModes(self._getPath(), modeEns[1], write_star=True)
        else:
            modeEns = [modes1, modes2]
        
        if self.metric == NMA_METRIC_OVERLAP:
            if self.norm:
                self.matrix = prody.calcOverlap(modeEns[0], modeEns[1], diag=self.diag)
            else:
                # Calculate direct dot product without vector normalisation found in calcOverlap
                self.matrix = modes1.getEigvecs().T @ modes2.getEigvecs()

            if self.matrix.ndim == 1:
                self.matrix.reshape(-1, 1)

        else:
            self.matrix = np.empty((nModes-6, 1))
            for i in range(6, nModes):
                if self.metric == NMA_METRIC_COV_OVERLAP:
                    self.matrix[i-6, 0] = prody.calcEnsembleSpectralOverlaps(modeEns[:, 6:i+1])[0, 1]
                else:
                    self.matrix[i-6, 0] = prody.calcRWSIP(modeEns[0, 6:i+1], modeEns[1, 6:i+1])

        prody.writeArray(self._getExtraPath('matrix.txt'), self.matrix,
                         format='%' + str(max([len(str(int(np.max(self.matrix)))),
                                               len(str(int(np.min(self.matrix))))]) + 4) + '.2f')

    def createOutputStep(self):
        outputMatrix = EMFile(filename=self._getExtraPath('matrix.txt'))

        if self.match:
            fnSqlite = self._getPath('modes.sqlite')
            inputClass = type(self.modes1.get())
            nmSet = inputClass(filename=fnSqlite)
            nmSet._nmdFileName = String(self.nmdFileName)

            outputPdb = AtomStruct()
            outputPdb.setFileName(self._getPath('atoms.pdb'))
            nmSet.setPdb(outputPdb.get())

            outputMatch = EMFile(filename=self.matchIndsFileName)
            nmSet._indsFileName = String(self.matchIndsFileName)

            self._defineOutputs(matrixFile=outputMatrix, matchFile=outputMatch, outputModes=nmSet)
        else:
            self._defineOutputs(matrixFile=outputMatrix)
