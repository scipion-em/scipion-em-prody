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

from prody2.protocols import (ProDyANM, ProDyDefvec, ProDyEdit,
                              ProDyImportModes, ProDyRTB, ProDyPCA)

import os
import prody

class ProDyModeViewer(Viewer):
    """ Visualization of a SetOfNormalModes from the ProDy ANM NMA protocol or elsewhere.    
        Normally, modes with high collectivity and low NMA score are preferred.
    """    
    _label = 'ProDy mode viewer'
    _targets = [SetOfNormalModes, ProDyANM, ProDyRTB, ProDyPCA,
                ProDyDefvec, ProDyEdit, ProDyImportModes]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    def _visualize(self, obj, **kwargs):
        """visualisation for mode sets"""
        # configure ProDy to automatically handle secondary structure information and verbosity
        old_secondary = prody.confProDy("auto_secondary")
        old_verbosity = prody.confProDy("verbosity")
        from pyworkflow import Config
        prodyVerbosity =  'none' if not Config.debugOn() else 'debug'
        prody.confProDy(auto_secondary=True, verbosity='{0}'.format(prodyVerbosity))

        type_ = type(obj)

        if isinstance(obj, SetOfNormalModes):
            modes = obj
        else:
            modes = obj.outputModes

        if not os.path.isfile(self.protocol._getPath("modes.nmd")):
            prody_modes = prody.parseScipionModes(modes.getFileName())
            modes_path = os.path.dirname(os.path.dirname(modes._getMapper().selectFirst().getModeFile()))

            atoms = prody.parsePDB(glob(modes_path+"/*atoms.pdb"), altloc="all")
            if isinstance(atoms, list):
                for atoms_i in atoms:
                    if atoms_i.numAtoms() == prody_modes.numAtoms():
                        prody_atoms = atoms_i
                        break
            else:
                prody_atoms = atoms

            prody.writeNMD(modes_path+"/modes.nmd", prody_modes, prody_atoms)

        # configure ProDy to restore secondary structure information and verbosity
        prody.confProDy(auto_secondary=old_secondary, verbosity='{0}'.format(old_verbosity))
        
        return [VmdView('-e "%s"' % self.protocol._getPath("modes.nmd"))]

