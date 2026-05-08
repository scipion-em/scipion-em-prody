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
This module will provide ProDy mode editing tools.
"""
import os
import numpy as np

from pwem.objects import AtomStruct, SetOfNormalModes, SetOfPrincipalComponents, String

from pyworkflow.utils import glob, logger
from pyworkflow.protocol.params import (PointerParam, EnumParam, BooleanParam,
                                        FloatParam, IntParam, LEVEL_ADVANCED)

import prody
from prody2.protocols.protocol_modes_base import ProDyModesBase

NMA_SLICE = 0
NMA_REDUCE = 1
NMA_EXTEND = 2
NMA_INTERP = 3

class ProDyEdit(ProDyModesBase):
    """
    This protocol edits a SetOfNormalModes object so that the modes are
    represented on a different set of nodes (atoms or pseudoatoms).

    AI Generated:

    Edit Modes (ProDyEdit) — User Manual
        Overview

        The Edit Modes protocol transforms an existing normal mode model so that
        it can be expressed on a new structural representation. In practical
        terms, it allows the user to decrease, increase, or remap the number of
        nodes that define the normal modes.

        In structural biology, this is especially useful when one wants to move
        between different levels of representation. For example, a mode
        calculated on a coarse pseudoatomic model may need to be transferred to
        an atomic model for interpretation, or an atomic normal mode analysis
        may need to be simplified to focus only on a specific structural region.

        The protocol preserves the dynamical information as much as possible,
        while adapting the eigenvectors to the geometry of a new structure.

        Inputs and General Workflow

        The protocol requires two main inputs.

        The first input is a SetOfNormalModes. These modes may originate from
        standard atomic normal mode analysis, pseudoatomic models generated from
        EM maps, or principal component analysis.

        The second input is a new set of nodes, provided as an atomic structure
        or pseudoatomic structure. These nodes define the target representation
        where the edited modes will be expressed.

        During execution, the protocol first establishes a correspondence
        between the original nodes and the new nodes. This alignment is based on
        chain matching and structural correspondence. Once this mapping is
        determined, the selected editing strategy is applied.

        Choosing the Editing Method

        The biological meaning of the result depends strongly on the selected
        editing strategy.

        Slice

        Slice is the simplest option. It keeps only the subset of nodes that
        correspond to the target structure.

        This option is particularly useful when the user wants to isolate a
        domain, chain, or structural fragment from a larger normal mode model.

        Biologically, slicing is often appropriate when the removed regions are
        not central to the motion of interest. It is computationally fast and
        usually the safest first option.

        A normalization option is available for Slice. When enabled, the output
        vectors are normalized after slicing. This is generally recommended when
        the edited modes will later be compared quantitatively.

        Reduce

        Reduce performs a more physically meaningful reduction of the model
        using Hessian-based reduction (vibrational subsystem analysis).

        Rather than simply cutting away nodes, the method attempts to preserve
        the dynamical coupling between retained and removed regions.

        This option is biologically preferable when the removed part of the
        structure may still influence the internal dynamics of the retained
        region.

        Because it relies on the original ProDy Hessian representation, Reduce
        may not always be available for modes generated in other workflows. In
        those cases, the protocol automatically falls back to slicing.

        Extend

        Extend increases the number of nodes by propagating the existing mode
        values from the original representation onto a denser target model.

        This is useful when a coarse-grained model needs to be interpreted at
        higher structural resolution.

        Biologically, extension is often applied when modes computed on CA-only
        or pseudoatomic models need to be visualized on a full atomic structure.

        The resulting motion preserves the coarse dynamical character but should
        not be interpreted as a newly computed all-atom normal mode analysis.

        Interpolate

        Interpolate also increases the number of nodes, but instead of directly
        copying values, it estimates new motions through thin-plate spline
        interpolation.

        This usually produces smoother spatial deformations than Extend.

        Biologically, interpolation is often preferred when transferring
        low-resolution motions onto high-resolution structural models for
        visualization, flexible fitting, or morph generation.

        However, interpolation introduces a stronger geometric assumption, so
        users should interpret local motions cautiously.

        Selecting the New Nodes

        The new node structure should correspond biologically to the original
        system.

        Best results are obtained when the original and target structures share
        chain identity, residue correspondence, and broadly similar overall
        geometry.

        Large differences in sequence, missing domains, or substantial
        structural rearrangements can make the mapping ambiguous and may lead to
        biologically misleading edited modes.

        As a practical rule, the new node model should represent the same
        macromolecular system, ideally in a nearby conformational state.

        Animation Options

        The protocol optionally generates animations for the ContinuousFlex
        viewer.

        These animations do not change the computed modes. They simply visualize
        the motion along the edited eigenvectors.

        RMSD amplitude controls how far atoms move along the mode. Larger values
        make motions easier to visualize but can exaggerate structural changes
        beyond physically realistic amplitudes.

        Number of frames determines the smoothness of the animation.

        Positive and negative directions allow the user to visualize motion in
        either direction along the mode vector.

        From a biological perspective, animations are especially useful for
        identifying collective motions such as hinge bending, domain rotations,
        breathing motions, or interface opening and closing.

        Outputs and Their Interpretation

        The main output is a new SetOfNormalModes expressed on the target node
        structure.

        The output preserves the mode framework of the original analysis but is
        now compatible with the geometry of the new structural model.

        A new PDB file describing the target nodes is also produced, together
        with NMD files for visualization in ProDy-compatible viewers.

        These outputs can be used directly in downstream analyses such as mode
        comparison, overlap calculations, animation, or flexible fitting.

        Practical Recommendations

        In most biological applications, Slice is the best starting point when
        the goal is simply to focus on a subregion of the structure.

        Reduce should be preferred when preserving dynamical coupling is
        important and the original modes come from a full ProDy model.

        Extend is useful for straightforward transfer from coarse to finer
        resolution representations.

        Interpolate is often the most visually appealing option when generating
        smooth deformations on high-resolution structures.

        Before interpreting the results biologically, it is always advisable to
        inspect the edited modes visually and verify that the node mapping makes
        structural sense.

        Final Perspective

        Edit Modes is fundamentally a transfer-of-representation protocol.

        Its biological usefulness lies in connecting dynamical information
        computed in one structural representation with another representation
        better suited for interpretation, visualization, or downstream
        analysis.

        The reliability of the result depends less on the mathematics of the
        transformation than on whether the chosen target nodes genuinely
        represent the same underlying biological system.
    """
    _label = 'Edit modes'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form, besidesAnimation=False):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label='ProDy edit')

        form.addParam('modes', PointerParam, label='Input set of modes',
                      pointerClass='SetOfNormalModes',
                      help='The input modes can be a SetOfNormalModes '
                           'from an atomic model (true PDB) or a pseudoatomic model '
                           '(an EM volume compared into pseudoatoms)'
                           'or a SetOfPrincipalComponents.')

        form.addParam('edit', EnumParam, choices=['Slice', 'Reduce', 'Extend', 'Interpolate'],
                    default=NMA_SLICE,
                    label='Type of edit',
                    help='Modes can have the number of nodes decreased using either eigenvector slicing '
                    'or the slower but often more meaningful Hessian reduction method (aka vibrational subsystem '
                    'analysis; Hinsen et al., Chem Phys 2000; Woodcock et al., J Chem Phys 2008) for ProDy vectors. \n'
                    'The number of nodes can be increased by extending (copying) eigenvector values '
                    'from nodes of the same residue or by through-space thin plate splines interpolation')

        form.addParam('newNodes', PointerParam,
                      label='New nodes',
                      pointerClass='AtomStruct',
                      help='Atoms or pseudoatoms to use as new nodes.')   

        form.addParam('norm', BooleanParam, default=True, 
                      condition='edit==%d' % NMA_SLICE,
                      label='Normalise sliced vectors',
                      help='Elect whether to normalise vectors.')
                      
        form.addSection(label='Animation')
        form.addParam('doAnimation', BooleanParam, default=False,
                      label='Make animations for ContinuousFlex viewer')
        animCheck = 'doAnimation == True'
        form.addParam('rmsd', FloatParam, default=5,
                      condition=animCheck,
                      label='RMSD Amplitude (A)',
                      help='Used only for animations of computed normal modes. '
                      'This is the maximal amplitude with which atoms or pseudoatoms are moved '
                      'along normal modes in the animations. \n')
        form.addParam('n_steps', IntParam, default=10,
                      condition=animCheck,
                      expertLevel=LEVEL_ADVANCED,
                      label='Number of frames',
                      help='Number of frames used in each direction of animations.')
        form.addParam('pos', BooleanParam, default=True,
                      condition=animCheck,
                      expertLevel=LEVEL_ADVANCED,
                      label="Include positive direction",
                      help='Elect whether to animate in the positive mode direction.')
        form.addParam('neg', BooleanParam, default=True,
                      condition=animCheck,
                      expertLevel=LEVEL_ADVANCED,
                      label="Include negative direction",
                      help='Elect whether to animate in the negative mode direction.')

    # --------------------------- STEPS functions ------------------------------
    # This is inherited from modes base protocol
    def _insertAllSteps(self):
        modes = prody.parseScipionModes(self.modes.get().getFileName())
        self.nzero = len(np.nonzero(modes.getEigvals() < prody.utilities.ZERO)[0])

        super(ProDyEdit, self)._insertAllSteps(len(self.modes.get()), self.nzero)

    def computeModesStep(self):
        self.inputStructure = self.modes.get().getPdb()
        modes = prody.parseScipionModes(self.modes.get().getFileName(),
                                        pdb=self.inputStructure.getFileName())

        oldNodes = prody.parsePDB(self.inputStructure.getFileName(), altloc="all")
        newNodes = prody.parsePDB(self.newNodes.get().getFileName(), altloc="all")

        nodesList = [oldNodes, newNodes]
        numAtomsArr = np.array([nodes.numAtoms() for nodes in nodesList])
        smaller = nodesList[np.argmin(numAtomsArr)]
        bigger = nodesList[np.argmax(numAtomsArr)]

        amap = prody.alignChains(bigger, smaller, match_func=prody.sameChid)[0]
        
        if self.edit == NMA_SLICE:
            self.outModes, self.atoms = prody.sliceModel(modes, bigger, amap, norm=self.norm)

        elif self.edit == NMA_REDUCE:
            modesPath = os.path.dirname(os.path.dirname(
                self.modes.get()._getMapper().selectFirst().getModeFile()))

            fromPrody = len(glob(modesPath+"/*npz"))
            if fromPrody:
                modes = prody.loadModel(glob(modesPath+"/*npz")[0])
                self.outModes, self.atoms = prody.reduceModel(modes, bigger, amap)
                zeros = bool(np.any(modes.getEigvals() < prody.utilities.ZERO))
                self.outModes.calcModes(modes.numModes(), zeros=zeros)
            else:
                logger.warn('ContinuousFlex modes cannot be reduced at this time. Slicing instead')
                self.outModes, self.atoms = prody.sliceModel(modes, bigger, amap, norm=self.norm)

        elif self.edit == NMA_EXTEND:
            self.outModes, self.atoms = prody.extendModel(modes, amap, bigger, norm=True)

        else:
            self.outModes, self.atoms = prody.interpolateModel(modes, amap, bigger, norm=True)

        prody.writePDB(self._getPath('atoms.pdb'), self.atoms)
        prody.writeScipionModes(self._getPath(), self.outModes, write_star=True)

        typeStr = str(type(self.outModes)).lower().split('.')[-1].split("'")[0]
        self.nmdFileName = self._getPath('modes.{0}.nmd'.format(typeStr))
        prody.writeNMD(self.nmdFileName, self.outModes, self.atoms)

        if isinstance(self.outModes, prody.GNM):
            self.gnm = True

    def createOutputStep(self):
        fnSqlite = self._getPath('modes.sqlite')

        inputClass = type(self.modes.get())
        nmSet = inputClass(filename=fnSqlite)
        nmSet._nmdFileName = String(self.nmdFileName)
        nmSet.setPdb(self.newNodes.get())

        self._defineOutputs(outputModes=nmSet)
        self._defineSourceRelation(self.newNodes, nmSet)
