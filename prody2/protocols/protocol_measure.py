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
This module will provide ProDy distance and angle measurement for structural ensembles
"""
import numpy as np

from pwem.objects import SetOfAtomStructs
from pwem.protocols import EMProtocol

import pyworkflow.object as pwobj
from pyworkflow.protocol import params

import prody
from prody2.constants import MEASURES, ENSEMBLE_WEIGHTS

RMSD = 0
DISTANCE = 1
ANGLE = 2
DIHEDRAL = 3

selstrHelp = '''The distance, angle or dihedral will be calculated between the centers of 2, 3 or 4 selections.
There is a rich selection engine with similarities to VMD. 
See http://http://www.bahargroup.org/prody/tutorials/prody_tutorial/selection.html'''

defaultSelstr = "protein and name CA or nucleic and name P C4' C2"
measureTypeCheck = "measureType>%d"

class ProDyMeasure(EMProtocol):
    """
    Performs geometric measurements on structural ensembles using ProDy.

    The protocol computes distances, angles, or dihedral angles between
    atom selections across all conformations of one or more ensembles.

    AI Generated:

    ProDy Measure (ProDyMeasure) — User Manual
        Overview

        The ProDyMeasure protocol performs geometric measurements on
        structural ensembles.

        Its main purpose is to quantify how specific regions of a
        molecular structure move relative to one another across multiple
        conformations.

        This protocol is particularly useful when studying structural
        variability, conformational transitions, domain rearrangements,
        or flexible motions in proteins and other macromolecular systems.

        Instead of focusing on global motions, it extracts specific
        geometric descriptors from selected regions.

        Input Data

        The protocol accepts one or more input ensembles provided as:

            - SetOfAtomStructs
            - ProDyNpzEnsemble

        All structures within each ensemble must contain the same number
        of atoms so that equivalent atom selections can be compared
        consistently across conformations.

        If the input consists of atomic structures, the protocol builds
        a ProDy ensemble directly from the input files without rejecting
        any conformations.

        If the input is already a ProDy ensemble, it is loaded directly.

        Measurement Types

        The protocol supports three kinds of geometric measurements:

            - distance
            - angle
            - dihedral

        The selected measure determines how many atom selections are
        required.

        Distance

        Distance requires two atom selections.

        For each conformation, the protocol calculates the geometric
        center of both selected atom groups and measures the distance
        between those centers.

        Biologically, this is useful for monitoring:

            - domain opening and closing
            - inter-subunit separation
            - ligand-induced displacement
            - motion between flexible structural regions

        Angle

        Angle requires three atom selections.

        The protocol computes the centers of the three selected regions
        and measures the angle formed by those centers.

        This can help characterize hinge motions or bending events
        involving three structural regions.

        Dihedral

        Dihedral requires four atom selections.

        The protocol computes the centers of four selected regions and
        calculates the dihedral angle for each conformation.

        This is particularly useful when studying torsional rearrangements
        or rotational motions involving multiple domains.

        Atom Selections

        The protocol uses atom selection strings to define the regions
        involved in the measurement.

        Each selection may contain any valid ProDy atom selection syntax.

        The measurement is not performed on individual atoms directly,
        but rather on the geometric centers of the selected groups.

        This design is especially useful in biological systems because
        it reduces local atomic noise and captures collective positional
        behavior.

        Computational Workflow

        For each input ensemble, the protocol performs the following
        steps:

            1. Load or construct the structural ensemble.
            2. Preserve the original atom set.
            3. Apply each atom selection independently.
            4. Compute the geometric center for each selection.
            5. Calculate the requested geometric measurement.
            6. Store one measurement value per conformation.

        Each measurement is associated with the corresponding object ID
        of the original structure.

        Output Files

        For every input ensemble, the protocol writes a CSV file
        containing the measured values.

        Each file stores one numerical measurement per conformation.

        This makes the results easy to inspect, plot, or use in
        downstream structural analysis.

        Output Ensembles

        The protocol also creates output ensembles that preserve the
        identity of the original structures.

        Each item in the output ensemble receives an additional attribute
        containing the computed measurement value.

        This means the measured geometry remains linked to the original
        structural objects, which is especially useful for filtering,
        sorting, or correlating structural states.

        Biological Interpretation

        The biological value of this protocol lies in converting complex
        structural variability into simple interpretable descriptors.

        For example:

            - a changing distance may indicate domain separation
            - a changing angle may reveal hinge flexibility
            - a changing dihedral may uncover rotational transitions

        Because measurements are computed across all conformations, the
        protocol can reveal continuous trends, state-dependent changes,
        or structural heterogeneity within the ensemble.

        Practical Recommendations

        The biological interpretation depends strongly on the quality of
        the atom selections.

        It is generally preferable to select structurally meaningful
        groups such as:

            - protein domains
            - helices
            - loops
            - active-site regions
            - subunit interfaces

        Very small selections may become noisy, whereas larger coherent
        structural regions often produce more robust measurements.

        Summary Information

        Once execution is complete, the protocol reports that the
        measurements have been calculated.

        If outputs are not yet available, the summary indicates that the
        calculation is still pending.

        Final Perspective

        ProDyMeasure is best understood as a targeted structural analysis
        tool.

        Rather than asking:

            "What are the dominant global motions?"

        it asks:

            "How does a specific geometric relationship change across
            the ensemble?"

        This makes it especially useful when the biological question
        focuses on specific structural rearrangements rather than global
        conformational modes.
    """
    _label = 'Measure'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        form.addSection(label='ProDy Measure')
        form.addParam('inputEnsemble', params.MultiPointerParam, label="Input ensemble(s)",
                      important=True,
                      pointerClass='SetOfAtomStructs,ProDyNpzEnsemble',
                      help='The input ensembles should be SetOfAtomStructs or ProDyNpzEnsemble '
                      'objects where all structures have the same number of atoms.')
        
        form.addParam('measureType', params.EnumParam,
                      choices=['rmsd', 'distance', 'angle', 'dihedral'], default=DISTANCE,
                      label='Measure type',
                      help='Select the type of measure.')

        form.addParam('selection1', params.StringParam, default=defaultSelstr,
                      label="selection string 1",
                      help=selstrHelp)
        
        form.addParam('selection2', params.StringParam, default=defaultSelstr,
                      label="selection string 2", condition=measureTypeCheck % RMSD,
                      help=selstrHelp)
        
        form.addParam('selection3', params.StringParam, default=defaultSelstr,
                      label="selection string 3", condition=measureTypeCheck % DISTANCE,
                      help=selstrHelp)

        form.addParam('selection4', params.StringParam, default=defaultSelstr,
                      label="selection string 4", condition=measureTypeCheck % ANGLE,
                      help=selstrHelp)


    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('computeStep')
        self._insertFunctionStep('createOutputStep')

    def computeStep(self):
        selstr1 = self.selection1.get()

        measureType = self.measureType.get()
        if measureType > RMSD:
            selstr2 = self.selection2.get()
        if measureType > DISTANCE:
            selstr3 = self.selection3.get()
        if measureType > ANGLE:
            selstr4 = self.selection4.get()

        self.measures = []
        for i, inputEnsemble in enumerate(self.inputEnsemble):
            ensGot = inputEnsemble.get()
            idSet = ensGot.getIdSet()
            if isinstance(ensGot, SetOfAtomStructs):
                ags = prody.parsePDB([tarStructure.getFileName() for tarStructure in ensGot])
                ens = prody.buildPDBEnsemble(ags, match_func=prody.sameChainPos, seqid=0., overlap=0., superpose=False, mapping=None)
                # the ensemble gets built exactly as the input is setup and nothing gets rejected
            else:
                ens = ensGot.loadEnsemble()

            atomsCopy = ens.getAtoms().copy()

            try:
                ens.setAtoms(atomsCopy)
            except ValueError:
                ens = prody.trimPDBEnsemble(ens)
                ens.setAtoms(atomsCopy)

            ens.setAtoms(atomsCopy.select(selstr1))
            centers1 = prody.calcCenter(ens.getCoordsets())

            if measureType != RMSD:
                ens.setAtoms(atomsCopy)
                ens.setAtoms(atomsCopy.select(selstr2))
                centers2 = prody.calcCenter(ens.getCoordsets())

            if measureType > DISTANCE:
                ens.setAtoms(atomsCopy)
                ens.setAtoms(atomsCopy.select(selstr3))
                centers3 = prody.calcCenter(ens.getCoordsets())

            if measureType == DIHEDRAL:
                ens.setAtoms(atomsCopy)
                ens.setAtoms(atomsCopy.select(selstr4))
                centers4 = prody.calcCenter(ens.getCoordsets())

            ens.setAtoms(atomsCopy)

            if measureType == RMSD:
                measures = ens.getRMSDs() # from start
            elif measureType == DISTANCE:
                measures = prody.calcDistance(centers1, centers2)
            elif measureType == ANGLE:
                measures = prody.measure.getAngle(centers1, centers2, centers3)
            else:
                measures = np.zeros(len(centers1))
                for j in range(len(centers1)):
                    measures[j] = prody.measure.getDihedral(centers1[j], centers2[j],
                                                            centers3[j], centers4[j])

            measuresDict = dict()
            for j, idx in enumerate(idSet):
                measuresDict[idx] = measures[j]
            self.measures.append(measuresDict)
            prody.writeArray(self._getPath('measures_{0}.csv'.format(i+1)), measures, 
                             format='%8.5f', delimiter=',')

    def createOutputStep(self):
        args = {}
        for self.ensId, inputEnsemble in enumerate(self.inputEnsemble): 
            ensGot = inputEnsemble.get()

            suffix = str(self.ensId+1)

            inputClass = type(ensGot)
            outSet = inputClass().create(self._getExtraPath(), suffix=suffix)
            outSet.copyItems(ensGot, updateItemCallback=self._setMeasures)

            name = "outputEns" + suffix
            
            args[name] = outSet

        self._defineOutputs(**args)

    # --------------------------- UTILS functions --------------------------------------------
    def _setMeasures(self, item, row=None):
        # We provide data directly so don't need a row
        measure = pwobj.Float(self.measures[self.ensId][item.getObjId()])
        setattr(item, MEASURES, measure)
        if not hasattr(item, ENSEMBLE_WEIGHTS):
            setattr(item, ENSEMBLE_WEIGHTS, pwobj.Float(1.))

    def _summary(self):
        if not hasattr(self, 'outputEns1'):
            summ = ['Measures not ready yet']
        else:
            summ = ['Measures calculated']
        return summ
        