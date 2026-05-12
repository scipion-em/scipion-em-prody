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
This module will provide ProDy projection of structural ensembles on principal component or normal modes
"""
import numpy as np
import os

from pwem.objects import SetOfAtomStructs
from pwem.protocols import EMProtocol

import pyworkflow.object as pwobj
from pyworkflow.protocol.params import (PointerParam, EnumParam, BooleanParam,
                                        MultiPointerParam, NumericRangeParam)
from pyworkflow.utils import getListFromRangeString, glob

import prody
from prody2.constants import PROJ_COEFFS

ONE = 0
TWO = 1
THREE = 2

class ProDyProject(EMProtocol):
    """
    Projects structural ensembles onto normal modes or principal
    components computed previously with ProDy.

    The protocol converts high-dimensional structural variability into a
    reduced coordinate representation.

    AI Generated:

    ProDy Projection (ProDyProject) — User Manual
        Overview

        The ProDyProject protocol projects structural ensembles onto a
        selected set of normal modes or principal components.

        Its main purpose is to describe each conformation of an ensemble
        in terms of coordinates along a reduced number of collective
        structural directions.

        Instead of analyzing every atomic coordinate directly, the
        protocol asks:

            "How much does each conformation move along these selected
            collective modes?"

        This is especially useful for simplifying structural variability
        and revealing dominant conformational trends.

        Input Data

        The protocol requires two kinds of input.

        Structural Ensembles

        One or more structural ensembles can be provided as:

            - SetOfAtomStructs
            - ProDyNpzEnsemble

        All conformations within an ensemble must contain equivalent
        atoms so that projections remain structurally meaningful.

        Input Modes

        The protocol also requires a previously computed set of modes.

        These can come from:

            - ProDy normal mode analysis
            - ProDy principal component analysis
            - Continuous-Flex NMA

        Only the selected modes are used for projection.

        Mode Selection

        The parameter:

            modeList

        allows the user to choose which modes are included.

        If left empty, the protocol uses all available modes starting
        from the first one.

        The user may also provide ranges or explicit lists of mode
        numbers.

        Examples include:

            - 7,8-10
            - 8,10,12
            - 8-12

        This is biologically important because different modes may
        capture different structural phenomena.

        Number of Projection Dimensions

        The parameter:

            numModes

        determines whether projection is performed onto:

            - 1 mode
            - 2 modes
            - 3 modes

        This effectively defines the dimensionality of the reduced
        conformational space.

        In practical structural analysis:

            - 1D projection reveals a single dominant structural trend
            - 2D projection often reveals conformational landscapes
            - 3D projection allows richer exploration of structural
              heterogeneity

        Projection Scaling

        Normalize

        The parameter:

            norm

        determines whether projections are normalized.

        Normalization is useful when comparing relative positions along
        the selected collective coordinates.

        RMSD Scaling

        The parameter:

            rmsd

        determines whether projection amplitudes are scaled to RMSD-like
        units.

        This often makes the projected coordinates easier to interpret in
        structural terms.

        Computational Workflow

        The protocol performs the following steps.

        Mode Preparation

        First, the selected input modes are loaded and optionally
        filtered according to the requested mode list.

        A new reduced mode set is then written to disk.

        The selected modes are also exported in NMD format for
        visualization.

        Ensemble Loading

        For each input ensemble, the protocol loads the conformations.

        If the input consists of atomic structures, it builds a ProDy
        ensemble directly from the structures.

        If the input is already a ProDy ensemble, it is loaded directly.

        Projection Calculation

        Each conformation is projected onto the selected modes.

        The result is a low-dimensional vector describing that
        conformation in the chosen collective coordinate system.

        For each structure, projection coefficients are associated with
        the original object identifier.

        Output Files

        For each input ensemble, the protocol writes a CSV file
        containing the projection coordinates.

        Each row corresponds to one conformation.

        The protocol also exports the weights associated with the
        ensemble entries.

        These additional files can be useful for downstream numerical
        analysis, plotting, or statistical interpretation.

        Output Ensembles

        For every input ensemble, the protocol creates an output
        ensemble.

        Each structure in the output receives a new attribute containing
        its projection coefficients.

        This preserves the identity of each conformation while enriching
        it with reduced-dimensional structural descriptors.

        Output Modes

        The protocol also creates an output mode set corresponding only
        to the modes used in the projection.

        This ensures that the reduced coordinates remain directly linked
        to the structural directions that define them.

        Biological Interpretation

        Projection is one of the most useful tools for understanding
        conformational landscapes.

        Biologically, projection allows the user to see whether
        conformations cluster, separate into states, or populate
        continuous transitions.

        Typical biological applications include:

            - detecting conformational substates
            - comparing functional structural states
            - mapping molecular dynamics trajectories
            - identifying transition pathways

        A projection does not define new motions.

        Instead, it quantifies how much each conformation expresses
        already defined collective motions.

        Practical Recommendations

        In most structural biology applications, projecting onto the
        first two or three most informative modes provides the clearest
        interpretation.

        Projection becomes especially powerful when combined with:

            - PCA
            - normal mode analysis
            - clustering
            - structural visualization

        Interpreting projections together with the original modes often
        reveals whether conformational variability corresponds to:

            - domain closure
            - hinge bending
            - twisting motions
            - continuous structural transitions

        Summary Information

        Once execution is complete, the protocol reports how many
        components were used for projection.

        If output is not yet available, the summary indicates that the
        projection is still pending.

        Final Perspective

        ProDyProject is best understood as a structural dimensionality
        reduction tool.

        Rather than asking:

            "What collective motions exist?"

        it asks:

            "Where does each conformation lie within the space defined
            by those collective motions?"

        This makes it especially useful for interpreting structural
        ensembles in a compact and biologically meaningful way.
    """
    _label = 'Projection'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        form.addSection(label='ProDy Projection')
        form.addParam('inputEnsemble', MultiPointerParam, label="Input ensemble(s)",
                      important=True,
                      pointerClass='SetOfAtomStructs,ProDyNpzEnsemble',
                      help='The input ensemble should be SetOfAtomStructs or ProDyNpzEnsemble '
                      'objects where all structures have the same number of atoms.')

        form.addParam('inputModes', PointerParam, label="Input set of modes",
                      important=True,
                      pointerClass='SetOfNormalModes,SetOfPrincipalComponents',
                      help='The input modes can come from Continuous-Flex NMA, ProDy NMA, or ProDy PCA.\n'
                           'The first modes from this set will be used. To use other modes, make a subset.')

        form.addParam('modeList', NumericRangeParam,
                      label="Modes selection", allowsNull=True, default="",
                      help='Select the normal modes that will be used for analysis.\n'
                           'If you leave this field empty, all the computed modes will be selected from.\n'
                           'If you only enter one number, all the computed modes will be selected from starting with that one.\n'
                           'You have several ways to specify the modes.\n'
                           '   Examples:\n'
                           ' "7,8-10" -> [7,8,9,10]\n'
                           ' "8, 10, 12" -> [8,10,12]\n'
                           ' "8 9, 10-12" -> [8,9,10,11,12])\n')
        
        form.addParam('numModes', EnumParam, choices=['1', '2', '3'],
                      label='Number of modes', default=TWO,
                      display=EnumParam.DISPLAY_HLIST,
                      help='1, 2 or 3 modes can be used for projection')

        form.addParam('norm', BooleanParam, label="Normalize?", default=False,
                      help='Select whether to normalise projections.')

        form.addParam('rmsd', BooleanParam, label="RMSD scale?", default=True,
                      help='Select whether to scale projections to RMSDs.')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('computeStep')
        self._insertFunctionStep('createOutputStep')

    def computeStep(self):
        inputModes = self.inputModes.get()
        modesPath = inputModes.getFileName()
        modes = prody.parseScipionModes(modesPath)

        if not self.modeList.empty():
            modeSelection = list(np.array(getListFromRangeString(self.modeList.get())) - 1)
            if len(modeSelection) == 1:
                modes = modes[modeSelection[0]:]
            else:
                modes = modes[modeSelection]

        modes = modes[:self.numModes.get()+1]

        prody.writeScipionModes(self._getPath(), modes, write_star=True)
        fnSqlite = self._getPath('modes.sqlite')
        inputClass = type(inputModes)
        self.outputModes = inputClass(filename=fnSqlite)

        atoms = prody.parsePDB(glob(os.path.dirname(modesPath)+"/*atoms.pdb")[0],
                               altloc="all")
        self.nmdFileName = self._getPath('modes.nmd')
        prody.writeNMD(self.nmdFileName, modes, atoms)

        self.outputModes._nmdFileName = pwobj.String(self.nmdFileName)

        self.proj = []
        for i, inputEnsemble in enumerate(self.inputEnsemble):
            ensGot = inputEnsemble.get()
            idSet = ensGot.getIdSet()
            if isinstance(ensGot, SetOfAtomStructs):
                ags = prody.parsePDB([tarStructure.getFileName() for tarStructure in ensGot])
                ens = prody.buildPDBEnsemble(ags, match_func=prody.sameChainPos, seqid=0., overlap=0., superpose=False, mapping=None)
                # the ensemble gets built exactly as the input is setup and nothing gets rejected
            else:
                ens = ensGot.loadEnsemble()

            projection = prody.calcProjection(ens, modes[:self.numModes.get()+1], rmsd=self.rmsd.get(),
                                              norm=self.norm.get())
            projDict = dict()
            for j, idx in enumerate(idSet):
                proj = projection[j]
                if isinstance(proj, float):
                    proj = [proj]

                projDict[idx] = proj
            self.proj.append(projDict)
            prody.writeArray(self._getPath('projection_{0}.csv'.format(i+1)), projection, 
                             format='%8.5f', delimiter=',')

            weights = np.array([np.array(item._prodyWeights, dtype=float) for item in ensGot])
            prody.writeArray(self._getPath('weights_{0}.csv'.format(i+1)), weights,
                             format='%8.5f', delimiter=',')

    def createOutputStep(self):
        args = {}
        for self.ensId, inputEnsemble in enumerate(self.inputEnsemble): 
            ensGot = inputEnsemble.get()

            suffix = str(self.ensId+1)

            inputClass = type(ensGot)
            outSet = inputClass().create(self._getExtraPath(), suffix=suffix)
            outSet.copyItems(ensGot, updateItemCallback=self._setCoeffs)
            name = "outputEns" + suffix
            args[name] = outSet

        args["outputModes"] = self.outputModes

        self._defineOutputs(**args)

    # --------------------------- UTILS functions --------------------------------------------
    def _setCoeffs(self, item, row=None):
        # We provide data directly so don't need a row
        vector = pwobj.CsvList()
        vector._convertValue(["{:18.15f}".format(x) for x in (self.proj[self.ensId][item.getObjId()])])
        setattr(item, PROJ_COEFFS, vector)

    def _summary(self):
        if not hasattr(self, 'outputEns1'):
            summ = ['Projection not ready yet']
        else:
            summ = ['Projected structures onto *{0}* components'.format(self.numModes.get()+1)]
        return summ
        