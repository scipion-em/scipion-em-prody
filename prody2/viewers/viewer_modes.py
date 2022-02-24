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
This module implements wrappers around the Xmipp NMA protocol
visualization program and the normal mode wizard NMWiz.
"""

from pyworkflow.viewer import Viewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.utils import *

from pwem.objects import SetOfNormalModes
from pwem.viewers import VmdView

from prody2.protocols import ProDyANM, ProDyDefvec, ProDyEdit

import os
import prody

class ProDyModeViewer(Viewer):
    """ Visualization of a SetOfNormalModes from the ProDy ANM NMA protocol or elsewhere.    
        Normally, modes with high collectivity and low NMA score are preferred.
    """    
    _label = 'ProDy mode viewer'
    _targets = [SetOfNormalModes, ProDyANM, ProDyDefvec, ProDyEdit]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    def _visualize(self, obj, **kwargs):
        """visualisation for mode sets"""

        type_ = type(obj)

        if isinstance(obj, SetOfNormalModes):
            modes = obj
        else:
            modes = obj.outputModes

        modes_path = os.path.dirname(os.path.dirname(modes[1].getModeFile()))

        if not os.path.isfile(self.protocol._getPath("modes.nmd")):
            prody_modes = prody.parseScipionModes(modes_path)
            atoms = prody.parsePDB(glob(modes_path+"/*atoms.pdb"), altloc="all")
            prody.writeNMD(modes_path+"/modes.nmd", prody_modes, atoms)

        return [VmdView('-e "%s"' % self.protocol._getPath("modes.nmd"))]

