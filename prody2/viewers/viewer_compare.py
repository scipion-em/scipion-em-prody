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
from prody.utilities.drawtools import IndexFormatter

class ProDyComparisonsViewer(ProtocolViewer):
    """ Visualization of results from the ProDy mode comparison protocol.    
    """    
    _label = 'ProDy mode comparisons viewer'
    _targets = [ProDyCompare]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    def _defineParams(self, form):
        metric = self.protocol.metric.get()
        diag = self.protocol.diag.get()

        self.modes1 = self.protocol.modes1.get()
        self.modes2 = self.protocol.modes2.get()
        many_modes = len(self.modes1) > 1 and len(self.modes2) > 1

        have_matrix = (many_modes and metric==NMA_METRIC_OVERLAP and diag==False)

        form.addSection(label='Visualization')
  
        form.addParam('displayMatrix', LabelParam,
                      condition=have_matrix == True,
                      label='Display matrix?',
                      help='Overlap matrices are shown as heatmaps.')

        form.addParam('displayBarsMain', LabelParam,
                      condition=have_matrix == False,
                      label='Display bar graphs?',
                      help='Values are shown as bars.')

        group = form.addGroup('Single mode from set 1',
                              condition=have_matrix == True)
        group.addParam('modeNumSet1', IntParam, default=7,
                       label='Mode number')
        group.addParam('displayBarsSet1', LabelParam, default=False,
                       label='Display bar graphs?',
                       help='Matrix rows are shown as bars.')
        group.addParam('cumulOverlapSet1', BooleanParam, default=False,
                       condition=metric==NMA_METRIC_OVERLAP,
                       label='Display cumulative overlap?',
                       help='Matrix rows are shown as bars.')

        group = form.addGroup('Single mode from set 2',
                              condition=have_matrix == True)
        group.addParam('modeNumSet2', IntParam, default=7,
                       label='Mode number')
        group.addParam('displayBarsSet2', LabelParam, default=False,
                       label='Display bar graphs?',
                       help='Matrix rows are shown as bars.')
        group.addParam('cumulOverlapSet2', BooleanParam, default=False,
                       condition=metric==NMA_METRIC_OVERLAP,
                       label='Display cumulative overlap?',
                       help='Matrix rows are shown as bars.')

        form.addParam('cumulOverlapMain', BooleanParam, default=False,
                      condition=(metric==NMA_METRIC_OVERLAP and have_matrix==False),
                      label='Display cumulative overlap?',
                      help='Matrix rows are shown as bars.')

        form.addParam('abs', BooleanParam, default=True,
                      label='Show absolute values?',
                      help='Otherwise signed values are shown.',
                      condition=metric == NMA_METRIC_OVERLAP)

    def _getVisualizeDict(self):
        return {'displayMatrix': self._viewMatrix,
                'displayBarsMain': self._viewSingleMode,
                'displayBarsSet1': self._viewSingleMode,
                'displayBarsSet2': self._viewSingleMode}            

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
        plotter.createSubPlot('', 'mode index from set 1', 'mode index from set 2')
        ax = plotter.figure.gca()
        plot = plotter.plotMatrix(ax, matrix, vmin, vmax, origin='lower')

        x_ticklabels = np.arange(1, len(self.modes1)+1)
        ax.xaxis.set_major_formatter(IndexFormatter(x_ticklabels))

        y_ticklabels = np.arange(1, len(self.modes2)+1) 
        ax.yaxis.set_major_formatter(IndexFormatter(y_ticklabels))

        plotter.getColorBar(plot)
        
        return [plotter]

    def _viewSingleMode(self, paramName):
        """ visualization for a selected mode. """
        matrix = prody.parseArray(self.protocol.matrixFile.getFileName())

        if paramName == 'displayBarsSet1':
            modes =  self.modes1
            modes_path = os.path.dirname(os.path.dirname(modes[1].getModeFile()))

            modeNumber = self.modeNumSet1.get()
            mode = modes[modeNumber]
            
            if modeNumber < 7 or mode is None:
                return [self.errorMessage("Invalid mode number *%d*\n"
                                        "Display the output Normal Modes to see "
                                        "the availables ones." % modeNumber,
                                        title="Invalid input")]

            row = matrix[modeNumber-1]
            cumulOverlap = self.cumulOverlapSet1.get()

            x_ticklabels = np.arange(1, len(modes)+1)

        elif paramName == 'displayBarsSet2':
            modes =  self.modes2
            modes_path = os.path.dirname(os.path.dirname(modes[1].getModeFile()))

            modeNumber = self.modeNumSet2.get()
            mode = modes[modeNumber]
            
            if modeNumber < 7 or mode is None:
                return [self.errorMessage("Invalid mode number *%d*\n"
                                        "Display the output Normal Modes to see "
                                        "the availables ones." % modeNumber,
                                        title="Invalid input")]

            row = matrix[:, modeNumber-1]
            cumulOverlap = self.cumulOverlapSet2.get()

            x_ticklabels = np.arange(1, len(self.modes2)+1) 

        else:
            modes_sets = [self.modes1, self.modes2]
            modes = modes_sets[np.argmax([len(nmset) for nmset in modes_sets])]

            row = matrix
            cumulOverlap = self.cumulOverlapMain.get()

            x_ticklabels = np.arange(1, len(modes)+1) 

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

        if cumulOverlap:
            cum_overlaps = np.sqrt(np.power(row, 2).cumsum(axis=row.ndim-1))
            plotter.plotData(inds, cum_overlaps, 'r')

        ax = plotter.figure.gca()
        ax.xaxis.set_major_formatter(IndexFormatter(x_ticklabels))
            
        return [plotter]

