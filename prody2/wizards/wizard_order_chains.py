# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:  James M. Krieger (jamesmkrieger@cnb.csic.es)
# *
# * Biocomputing Unit, CNB-CSIC
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
These wizards will extract the chains from a atomic structure (pdb) file in
order to select it in the protocol.
Then, it will load the structure and will take all chain related
information such as name and number of residues.
"""

# Imports
from collections import OrderedDict
import re

from ..protocols.protocol_atoms import ProDyAlign
from ..protocols.protocol_ensemble import ProDyBuildPDBEnsemble
from ..protocols.protocol_lda import ProDyLDA

from pwem.wizards import VariableWizard

class ProDyAddChainOrderWizard(VariableWizard):
    """Add a step of the workflow in the defined position"""
    _targets, _inputs, _outputs = [], {}, {}

    def show(self, form, *params):
        inputParam, outputParam = self.getInputOutput(form)
        protocol = form.protocol
        index = getattr(protocol, inputParam[0]).get()

        noLabel = True
        try:
            label = getattr(protocol, inputParam[1]).get()
            matchDic = protocol.createMatchDic(index, label)
        except IndexError:
            raise IndexError("no input param 1")
        else:
            noLabel = False

        if noLabel:
            matchDic = protocol.createMatchDic(index)

        output = str(matchDic)
        cnt = 0
        for i in re.finditer(" ", output):
            cnt=cnt+1
            if cnt%3==0:
                output = output[:i.start()] + "\n" + output[i.start() + 1:]
            
        if output.find('\n' + " "*22) == -1:
            output = output.replace('\n', '\n' + " "*22)

        form.setVar(outputParam[0], output)


ProDyAddChainOrderWizard().addTarget(protocol=ProDyAlign,
                                     targets=['insertOrder'],
                                     inputs=['insertOrder'],
                                     outputs=['chainOrders'])

ProDyAddChainOrderWizard().addTarget(protocol=ProDyBuildPDBEnsemble,
                                     targets=['insertOrder'],
                                     inputs=['insertOrder', 'label'],
                                     outputs=['chainOrders'])

ProDyAddChainOrderWizard().addTarget(protocol=ProDyLDA,
                                         targets=['insertOrder'],
                                         inputs=['insertOrder', 'label'],
                                         outputs=['chainOrders'])

class ProDyRecoverChainOrderWizard(VariableWizard):
    """Watch the parameters of the step of the workflow defined by the index"""
    _targets, _inputs, _outputs = [], {}, {}

    def show(self, form, *params):
        inputParam, outputParam = self.getInputOutput(form)
        protocol = form.protocol
        
        index = int(getattr(protocol, inputParam[0]).get()) - 1
        matchDic = eval(protocol.chainOrders.get())
        form.setVar(outputParam[0], list(matchDic.keys())[index])
        form.setVar(outputParam[1], list(matchDic.values())[index])


ProDyRecoverChainOrderWizard().addTarget(protocol=ProDyAlign,
                                         targets=['recoverOrder'],
                                         inputs=['recoverOrder'],
                                         outputs=['label', 'customOrder'])

ProDyRecoverChainOrderWizard().addTarget(protocol=ProDyBuildPDBEnsemble,
                                         targets=['recoverOrder'],
                                         inputs=['recoverOrder'],
                                         outputs=['label', 'customOrder'])

ProDyRecoverChainOrderWizard().addTarget(protocol=ProDyLDA,
                                         targets=['recoverOrder'],
                                         inputs=['recoverOrder'],
                                         outputs=['label', 'customOrder'])
