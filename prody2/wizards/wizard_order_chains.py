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
from ..protocols.protocol_atoms import *
from ..protocols.protocol_ensemble import *

import pyworkflow.wizard as pwizard
from pwchem.wizards import VariableWizard

class ProDyAddChainOrderWizard(VariableWizard):
    """Add a step of the workflow in the defined position"""
    _targets, _inputs, _outputs = [], {}, {}

    def show(self, form, *params):
        inputParam, outputParam = self.getInputOutput(form)
        protocol = form.protocol
        index = getattr(protocol, inputParam[0]).get()
        matchDic = protocol.createMatchDic(index)
        form.setVar(outputParam[0], str(matchDic).replace('),', '),\n' + ' '*20))


ProDyAddChainOrderWizard().addTarget(protocol=ProDyAlign,
                                     targets=['insertOrder'],
                                     inputs=['insertOrder'],
                                     outputs=['chainOrders'])

ProDyAddChainOrderWizard().addTarget(protocol=ProDyBuildPDBEnsemble,
                             targets=['insertOrder'],
                             inputs=['insertOrder'],
                             outputs=['chainOrders'])

class ProDyRecoverChainOrderWizard(VariableWizard):
    """Watch the parameters of the step of the workflow defined by the index"""
    #_targets = [(ProDyAlign, ['recoverOrder'])]
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
