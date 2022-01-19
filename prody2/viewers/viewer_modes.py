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

from pyworkflow.gui.project import ProjectWindow
from pyworkflow.protocol.params import LabelParam, IntParam
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pwem.viewers import ObjectView, VmdView
from pwem.emlib import MDL_NMA_ATOMSHIFT

from continuousflex.viewers.nma_plotter import FlexNmaPlotter
from continuousflex.viewers.viewer_nma import createShiftPlot

from prody2.protocols import ProDyANM, ProDyWriteNMD

import os

class ProDyModeViewer(ProtocolViewer):
    """ Visualization of results from the ProDy ANM NMA protocol.    
        Normally, modes with high collectivity and low NMA score are preferred.
    """    
    _label = 'ProDy mode viewer'
    _targets = [ProDyANM, ProDyWriteNMD]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    def _defineParams(self, form):
        form.addSection(label='Visualization')
  
        form.addParam('displayModes', LabelParam,
                      label="Display output Normal Modes?", important=True)
        form.addParam('displayMaxDistanceProfile', LabelParam,
                      label="Plot max distance profile?",
                      help="Maximum unitary shift of each atom or pseudoatom over all computed modes.")
    
        # VMD display is now a property for the whole set of modes not single modes
        form.addParam('displayVmd', LabelParam,
                       label='Display mode animations with VMD NMWiz?') 
        
        group = form.addGroup('Single mode')  
        group.addParam('modeNumber', IntParam, default=7,
              label='Mode number')
        group.addParam('displayDistanceProfile', LabelParam, default=False,
                      label="Plot mode distance profile?",
                      help="Unitary shift of each atom or pseudoatom along the mode that is requested to be animated.")

    def _getVisualizeDict(self):
        return {'displayModes': self._viewParam,
                'displayMaxDistanceProfile': self._viewParam,
                'displayVmd': self._viewParam,
                'displayDistanceProfile': self._viewSingleMode,
                } 

    def _viewParam(self, paramName):
        """ visualisation for mode sets"""
        modes = self.protocol.outputModes
        modes_path = os.path.dirname(os.path.dirname(modes[1].getModeFile()))

        if paramName == 'displayModes':    
            return [ObjectView(self._project, modes.strId(), self.protocol.outputModes.getFileName())]
        elif paramName == 'displayMaxDistanceProfile':
            fn = modes_path + "/extra/maxAtomShifts.xmd"
            return [createShiftPlot(fn, "Maximum atom shifts", "maximum shift")]
        elif paramName == 'displayVmd':
            return [VmdView('-e "%s"' % self.protocol._getPath("modes.nmd"))]

    def _viewSingleMode(self, paramName):
        """ visualization for a selected mode. """
        modes =  self.protocol.outputModes
        modes_path = os.path.dirname(os.path.dirname(modes[1].getModeFile()))

        modeNumber = self.modeNumber.get()
        mode = modes[modeNumber]
        
        if mode is None:
            return [self.errorMessage("Invalid mode number *%d*\n"
                                      "Display the output Normal Modes to see "
                                      "the availables ones." % modeNumber,
                                      title="Invalid input")]
        elif paramName == 'displayDistanceProfile':
            return [createDistanceProfilePlot(modes_path, modeNumber)]

def createDistanceProfilePlot(modes_path, modeNumber):
    vectorMdFn = modes_path + "/extra/distanceProfiles/vec%d.xmd" % modeNumber
    plotter = createShiftPlot(vectorMdFn, "Atom shifts for mode %d"
                              % modeNumber, "shift")
    return plotter

