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
This module implements wrappers around the ProDy and Xmipp NMA protocol
visualization program and the normal mode wizard NMWiz.
"""

from pyworkflow.viewer import Viewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.utils import *

from pwem.objects import SetOfNormalModes
from pwem.viewers import VmdView

from prody2.objects import ProDyNpzEnsemble
from prody2.protocols import (ProDyBuildPDBEnsemble, ProDyImportEnsemble,
                              ProDyLDA, ProDyPCA)

import os
import prody

class ProDyEnsembleViewer(Viewer):
    """ Visualization of an ensemble from the ProDy protocol or elsewhere
    """    
    _label = 'ProDy ensemble viewer'
    _targets = [ProDyNpzEnsemble, ProDyBuildPDBEnsemble, ProDyImportEnsemble,
                ProDyLDA, ProDyPCA]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    def _visualize(self, obj, **kwargs):
        """visualisation for ensembles"""

        if isinstance(obj, ProDyNpzEnsemble):
            ensemble = obj
        else:
            try:
                ensemble = obj.outputNpz
            except AttributeError:
                ensemble = obj.outputEnsemble

        ensPath = os.path.dirname(os.path.dirname(list(ensemble.getFiles())[0])) + "/"
        ensFn = ensPath + "ensemble.dcd"
        atomsFn = ensPath + "atoms.pdb"

        if not os.path.isfile(ensPath + "ensemble.dcd"):
            ens = ensemble.loadEnsemble()
            atoms = ens.getAtoms()
            prody.writePDB(atomsFn, atoms)
            prody.writeDCD(ensFn, ens)
            
        cmdFn = ensPath + "ensemble.vmd"
        fhCmd = open(cmdFn, "w")
        fhCmd.write("mol new %s\n" % atomsFn)
        fhCmd.write("mol addfile %s\n" % ensFn)
        fhCmd.write("mol modcolor 0 0 Chain\n")
        fhCmd.write("mol modstyle 0 0 Tube\n")
        fhCmd.close()
        
        return [VmdView('-e "%s"' % cmdFn)]

