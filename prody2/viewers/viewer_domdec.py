# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     James Krieger (jmkrieger@cnb.csic.es)
# *              Ricardo Serrano Gutiérrez (rserranogut@hotmail.com)                 
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
This module implements the dynamical domains decomposition protocol
visualization program using VMD.
"""
from pyworkflow.viewer import Viewer, DESKTOP_TKINTER, WEB_DJANGO
from pwem.viewers import VmdView
from prody2.protocols import ProDyDomainDecomp

class ProDyDomainViewer(Viewer):
    """ Visualization of domains from GNM domain decomposition
    """    
    _label = 'Dynamical domain viewer'
    _targets = [ProDyDomainDecomp]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    def _visualize(self, obj, **kwargs):
        """visualisation for mode dynamical domains"""

        return [VmdView('-e "%s"' % self.protocol._getPath("domains.vmd"))]

