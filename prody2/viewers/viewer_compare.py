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
This module implements wrappers around the ProDy tools 
for plotting mode comparison matrices and graphs.
"""

from pyworkflow.protocol.params import LabelParam, IntParam, BooleanParam
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO

from pwem.viewers.plotter import EmPlotter

from prody2.protocols import ProDyCompare
from prody2.protocols.protocol_compare import NMA_METRIC_OVERLAP

import numpy as np
import os
import prody

class ProDyComparisonsViewer(ProtocolViewer):
    """ Visualization of results from the ProDy mode comparison protocol.    
    """    
    _label = 'ProDy mode comparisons viewer'
    _targets = [ProDyCompare]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    def _defineParams(self, form):
        metric = self.protocol.metric.get()
        diag = self.protocol.diag.get()

        form.addSection(label='Visualization')
  
        form.addParam('displayWhole', LabelParam,
                      label='Display whole output matrix?', important=True,
                      help='Truly 2D matrices are shown as heatmaps. '
                           'Diagonal or cumulative values are shown as bars.')
        
        group = form.addGroup('Single mode', 
                              condition=(metric==NMA_METRIC_OVERLAP and diag==False))      
        group.addParam('modeNumber', IntParam, default=7,
                       label='Mode number')
        group.addParam('displayBars', LabelParam, default=False,
                       label='Display matrix row?',
                       help='Matrix rows are shown as bars.')
        group.addParam('displayCumulOverlap', BooleanParam, default=False, 
                       condition=metric==NMA_METRIC_OVERLAP,
                       label='Display cumulative overlap?',
                       help='Matrix rows are shown as bars.')

        form.addParam('abs', BooleanParam, default=True,
                       label='Show absolute values?',
                       help='Otherwise signed values are shown.', 
                       condition=metric==NMA_METRIC_OVERLAP)

    def _getVisualizeDict(self):
        matrix = prody.parseArray(self.protocol.matrixFile.getFileName())

        if matrix.ndim == 2:
            return {'displayWhole': self._viewMatrix,
                    'displayBars': self._viewSingleMode} 
        else:
            return {'displayWhole': self._viewSingleMode}             

    def _viewMatrix(self, paramName):
        """ visualisation for 2D mode comparisons""" 
        matrix = prody.parseArray(self.protocol.matrixFile.getFileName())

        if self.abs:
            matrix = abs(matrix)
            vmin = 0
            vmax = 1
        else:
            vmin = -1
            vmax = 1

        plotter = EmPlotter()
        plotter.createSubPlot('', 'mode index', 'mode index')
        ax = plotter.figure.gca()
        plotter.plotMatrix(ax, matrix, vmin, vmax, origin='lower')
        return [plotter]

    def _viewSingleMode(self, paramName):
        """ visualization for a selected mode. """
        matrix = prody.parseArray(self.protocol.matrixFile.getFileName())

        if matrix.ndim == 2:
            modes =  self.protocol.modes1.get()
            modes_path = os.path.dirname(os.path.dirname(modes[1].getModeFile()))

            modeNumber = self.modeNumber.get()
            mode = modes[modeNumber]
            
            if modeNumber < 7 or mode is None:
                return [self.errorMessage("Invalid mode number *%d*\n"
                                        "Display the output Normal Modes to see "
                                        "the availables ones." % modeNumber,
                                        title="Invalid input")]

            row = matrix[modeNumber]
        else:
            row = matrix

        inds = np.arange(len(row))

        if self.abs:
            row = abs(row)

        plotter = EmPlotter()
        plotter.createSubPlot('', 'mode index', 'overlap') # needed for self.bar
        plotter.plotDataBar(inds, row, width=0.8)

        if matrix.ndim == 2:
            ca = plotter.figure.gca()
            if self.abs:
                ca.set_ylim([0, 1])
            else:
                ca.set_ylim([-1, 1])

        if self.displayCumulOverlap:
            cum_overlaps = np.sqrt(np.power(row, 2).cumsum(axis=row.ndim-1))
            plotter.plotData(inds, cum_overlaps, 'r')
            
        return [plotter]

