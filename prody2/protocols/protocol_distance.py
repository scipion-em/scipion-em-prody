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
This module will provide ProDy distance calculation for structural ensembles
"""

from pwem.objects import SetOfAtomStructs
from pwem.protocols import EMProtocol

import pyworkflow.object as pwobj
from pyworkflow.protocol.params import StringParam, MultiPointerParam

import prody
from prody2.constants import DISTANCES
from prody2.objects import ProDyNpzEnsemble

ONE = 0
TWO = 1
THREE = 2

class ProDyDistance(EMProtocol):
    """
    This module will provide ProDy distance calculation for structural ensembles
    """
    _label = 'Distance'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        form.addSection(label='ProDy Distance')
        form.addParam('inputEnsemble', MultiPointerParam, label="Input ensemble",
                      important=True,
                      pointerClass='SetOfAtomStructs,ProDyNpzEnsemble',
                      help='The input ensemble should be a SetOfAtomStructs '
                      'where all structures have the same number of atoms.')

        form.addParam('selection1', StringParam, default="protein and name CA or nucleic and name P C4' C2",
                      label="selection string 1",
                      help='The distance will be calculated between the average positions of the two selections. '
                           'There is a rich selection engine with similarities to VMD. '
                           'See http://prody.csb.pitt.edu/tutorials/prody_tutorial/selection.html')
        
        form.addParam('selection2', StringParam, default="protein and name CA or nucleic and name P C4' C2",
                      label="selection string 2",
                      help='The distance will be calculated between the average positions of the two selections. '
                           'There is a rich selection engine with similarities to VMD. '
                           'See http://prody.csb.pitt.edu/tutorials/prody_tutorial/selection.html')


    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('computeStep')
        self._insertFunctionStep('createOutputStep')

    def computeStep(self):
        # configure ProDy to automatically handle secondary structure information and verbosity
        oldSecondary = prody.confProDy("auto_secondary")
        oldVerbosity = prody.confProDy("verbosity")
        from pyworkflow import Config
        prodyVerbosity =  'none' if not Config.debugOn() else 'debug'
        prody.confProDy(auto_secondary=True, verbosity='{0}'.format(prodyVerbosity))

        selstr1 = self.selection1.get()
        selstr2 = self.selection2.get()

        self.distances = []
        for _, inputEnsemble in enumerate(self.inputEnsemble):
            ensGot = inputEnsemble.get()
            if isinstance(ensGot, SetOfAtomStructs):
                ags = prody.parsePDB([tarStructure.getFileName() for tarStructure in ensGot])
                ens = prody.buildPDBEnsemble(ags, match_func=prody.sameChainPos, seqid=0., overlap=0., superpose=False)
                # the ensemble gets built exactly as the input is setup and nothing gets rejected
            else:
                ens = ensGot.loadEnsemble()

            atoms = ens.getAtoms()

            ens.select(selstr1)
            centers1 = prody.calcCenter(ens.getCoordsets())
            ens.setAtoms(atoms)

            ens.select(selstr2)
            centers2 = prody.calcCenter(ens.getCoordsets())
            ens.setAtoms(atoms)

            self.distances.append(prody.calcDistance(centers1, centers2))

        # configure ProDy to restore secondary structure information and verbosity
        prody.confProDy(auto_secondary=oldSecondary, verbosity='{0}'.format(oldVerbosity))

    def createOutputStep(self):

        args = {}
        for self.ensId, inputEnsemble in enumerate(self.inputEnsemble): 
            ensGot = inputEnsemble.get()

            suffix = str(self.ensId+1)

            inputClass = type(ensGot)
            outSet = inputClass().create(self._getExtraPath(), suffix=suffix)
            outSet.copyItems(ensGot, updateItemCallback=self._setDistances)

            name = "outputEns" + suffix
            
            args[name] = outSet

        self._defineOutputs(**args)

    # --------------------------- UTILS functions --------------------------------------------
    def _setDistances(self, item, row=None):
        # We provide data directly so don't need a row
        vector = pwobj.CsvList()
        vector._convertValue(["{:6.3f}".format(self.distances[self.ensId][item.getObjId()-1])])
        setattr(item, DISTANCES, vector)

    def _summary(self):
        if not hasattr(self, 'outputStructures'):
            summ = ['Distances not ready yet']
        else:
            summ = ['Distances calculated for *{0}* structures'.format(len(self.outputStructures))]
        return summ
        