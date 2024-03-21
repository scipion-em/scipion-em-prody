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
This module will provide ProDy distance and angle measurement for structural ensembles
"""
import numpy as np

from pwem.objects import SetOfAtomStructs
from pwem.protocols import EMProtocol

import pyworkflow.object as pwobj
from pyworkflow.protocol import params

import prody
from prody2.constants import MEASURES
from prody2 import fixVerbositySecondary, restoreVerbositySecondary

DISTANCE = 0
ANGLE = 1
DIHEDRAL = 2

selstrHelp = '''The distance, angle or dihedral will be calculated between the centers of 2, 3 or 4 selections.
There is a rich selection engine with similarities to VMD. 
See http://prody.csb.pitt.edu/tutorials/prody_tutorial/selection.html'''

defaultSelstr = "protein and name CA or nucleic and name P C4' C2"

class ProDyMeasure(EMProtocol):
    """
    This module will provide ProDy distance and angle measurement for structural ensembles
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
                      choices=['distance', 'angle', 'dihedral'], default=DISTANCE,
                      label='Measure type',
                      help='Select the type of measure.')

        form.addParam('selection1', params.StringParam, default=defaultSelstr,
                      label="selection string 1",
                      help=selstrHelp)
        
        form.addParam('selection2', params.StringParam, default=defaultSelstr,
                      label="selection string 2",
                      help=selstrHelp)
        
        form.addParam('selection3', params.StringParam, default=defaultSelstr,
                      label="selection string 3", condition="measureType>%d" % DISTANCE,
                      help=selstrHelp)

        form.addParam('selection4', params.StringParam, default=defaultSelstr,
                      label="selection string 4", condition="measureType>%d" % ANGLE,
                      help=selstrHelp)


    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('computeStep')
        self._insertFunctionStep('createOutputStep')

    def computeStep(self):
        fixVerbositySecondary(self)

        selstr1 = self.selection1.get()
        selstr2 = self.selection2.get()

        measureType = self.measureType.get()
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

            if measureType == DISTANCE:
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

        restoreVerbositySecondary(self)

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

    def _summary(self):
        if not hasattr(self, 'outputEns1'):
            summ = ['Measures not ready yet']
        else:
            summ = ['Measures calculated']
        return summ
        