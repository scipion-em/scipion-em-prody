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
This module will provide ProDy linear discriminant analysis (LRA) using atomic structures
"""
from pwem.objects import Float, String

from prody2.protocols.protocol_lda import ProDyLDA
from prody2.objects import SetOfLogisticModes, loadAndWriteEnsemble
from prody2.constants import PRODY_FRACT_VARS
from prody2 import fixVerbositySecondary, restoreVerbositySecondary

import prody


class ProDyLRA(ProDyLDA):
    """
    This protocol will perform ProDy logistic regression analysis (LRA) using atomic structures
    """
    _label = 'LRA'
    _possibleOutputs = {'outputModes': SetOfLogisticModes}

    def computeModesStep(self, n=1):
        # configure ProDy to automatically handle secondary structure information and verbosity
        fixVerbositySecondary(self)

        loadAndWriteEnsemble(self)
        self.atoms = self.ens.getAtoms()

        restoreVerbositySecondary(self)

        self.outModes = prody.LRA()
        self.outModes.calcModes(self.ens, self.classes,
                                n_shuffles=self.numberOfShuffles.get())

        prody.writeScipionModes(self._getPath(), self.outModes)
        self._nmdFileName = String(self._getPath('modes.logreg.nmd'))
        prody.writeNMD(self._nmdFileName.get(), self.outModes, self.atoms)
        prody.saveModel(self.outModes, self._getPath('modes.logreg.npz'), matrices=True)

    def createOutputStep(self):
        fnSqlite = self._getPath('modes.sqlite')
        nmSet = SetOfLogisticModes(filename=fnSqlite)
        nmSet._nmdFileName = self._nmdFileName

        self.fractVarsDict = {}
        for _, item in enumerate(nmSet):
            self.fractVarsDict[item.getObjId()] = 1

        outSet = SetOfLogisticModes().create(self._getPath())
        outSet.copyItems(nmSet, updateItemCallback=self._setFractVars)
        outSet._nmdFileName = self._nmdFileName

        inputPdb = self.averageStructure
        self._defineOutputs(refPdb=inputPdb)
        outSet.setPdb(inputPdb)

        self._defineOutputs(outputModes=outSet, outputEnsemble=self.npz)
        self._defineSourceRelation(inputPdb, outSet)

    def _summary(self):
        if not hasattr(self, 'outputModes'):
            summ = ['Output modes not ready yet']
        else:
            modes = prody.parseScipionModes(self.outputModes.getFileName())
            ens = self.outputEnsemble.loadEnsemble()

            summ = ['*{0}* LRA components calculated from *{1}* structures of *{2}* atoms'.format(
                    modes.numModes(), ens.numConfs(), ens.numAtoms())]
        return summ

    def _setFractVars(self, item, row=None):
        # We provide data directly so don't need a row
        fractVar = Float(self.fractVarsDict[item.getObjId()])
        setattr(item, PRODY_FRACT_VARS, fractVar)
