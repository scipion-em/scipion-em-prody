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

from matplotlib import ticker
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

        matched = self.protocol.match.get()

        self.modes1 = self.protocol.modes1.get()

        if not matched:
            self.modes2 = self.protocol.modes2.get()
        else:
            self.modes2 = self.protocol.outputModes

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

        group = form.addGroup('Single mode from set 1 (row)',
                              condition=have_matrix == True)
        group.addParam('displayBarsSet1', LabelParam, default=False,
                       label='Display bar graphs?',
                       help='Matrix rows are shown as bars.')
        group.addParam('modeNumSet1', IntParam, default=7,
                       label='Mode number')
        group.addParam('cumulOverlapSet1', BooleanParam, default=False,
                       condition=metric==NMA_METRIC_OVERLAP,
                       label='Display cumulative overlap?',
                       help='Cumulative overlaps from matrix rows are shown as lines.')

        group = form.addGroup('Single mode from set 2 (column)',
                              condition=have_matrix == True)
        group.addParam('displayBarsSet2', LabelParam, default=False,
                       label='Display bar graphs?',
                       help='Matrix columns are shown as bars.')
        group.addParam('modeNumSet2', IntParam, default=7,
                       label='Mode number')
        group.addParam('cumulOverlapSet2', BooleanParam, default=False,
                       condition=metric==NMA_METRIC_OVERLAP,
                       label='Display cumulative overlap?',
                       help='Cumulative overlaps from matrix rows are shown as lines.')

        form.addParam('cumulOverlapMain', BooleanParam, default=False,
                      condition=(metric==NMA_METRIC_OVERLAP and have_matrix==False),
                      label='Display cumulative overlap?',
                      help='Matrix rows are shown as bars.')

        form.addParam('matchedModeNum', BooleanParam, default=False,
                      condition=matched==True,
                      label='Use reordered mode numbers',
                      help='If true, modes are indexed by their numbers '
                            'in the original set before matching. '
                            'Otherwise, they are indexed directly using '
                            'the new order')

        form.addParam('allticks', BooleanParam, default=True,
                      label='Show all axis tick labels',
                      help='This could get too crowded with many modes, '
                            'but is helpful for seeing matches.')

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

        # configure ProDy to automatically handle secondary structure information and verbosity
        self.old_secondary = prody.confProDy("auto_secondary")
        self.old_verbosity = prody.confProDy("verbosity")
        from pyworkflow import Config
        prodyVerbosity =  'none' if not Config.debugOn() else 'debug'
        prody.confProDy(auto_secondary=True, verbosity='{0}'.format(prodyVerbosity))

        matrix = prody.parseArray(self.protocol.matrixFile.getFileName())

        if self.abs:
            matrix = abs(matrix)
            vmin = 0
            vmax = 1
        else:
            vmin = -1
            vmax = 1

        plotter = EmPlotter()
        plotter.createSubPlot('', 'mode index from set 2', 'mode index from set 1')
        ax = plotter.figure.gca()
        plot = prody.showMatrix(matrix, vmin=vmin, vmax=vmax, origin='lower')

        y_ticklabels = np.arange(1, len(self.modes1)+1)

        if os.path.isfile(self.protocol._getExtraPath('match_inds.txt')) and self.matchedModeNum.get():
            x_ticklabels = prody.parseArray(self.protocol._getExtraPath('match_inds.txt'), dtype=int)
        else:
            x_ticklabels = np.arange(1, len(self.modes2)+1) 

        ax.xaxis.set_major_formatter(IndexFormatter(x_ticklabels))
        ax.yaxis.set_major_formatter(IndexFormatter(y_ticklabels))

        if self.allticks.get():
            ax.xaxis.set_major_locator(ticker.IndexLocator(offset=0.5, base=1.))
            ax.yaxis.set_major_locator(ticker.IndexLocator(offset=0.5, base=1.))
        else:
            locator = ticker.AutoLocator()
            locator.set_params(integer=True)
            minor_locator = ticker.AutoMinorLocator()

            ax.xaxis.set_major_locator(locator)
            ax.xaxis.set_minor_locator(minor_locator)
            
            ax.yaxis.set_major_locator(locator)
            ax.yaxis.set_minor_locator(minor_locator)

        # configure ProDy to restore secondary structure information and verbosity
        prody.confProDy(auto_secondary=self.old_secondary, 
                        verbosity='{0}'.format(self.old_verbosity))
        
        return [plotter]

    def _viewSingleMode(self, paramName):
        """ visualization for a selected mode. """

        # configure ProDy to automatically handle secondary structure information and verbosity
        self.old_secondary = prody.confProDy("auto_secondary")
        self.old_verbosity = prody.confProDy("verbosity")
        from pyworkflow import Config
        prodyVerbosity =  'none' if not Config.debugOn() else 'debug'
        prody.confProDy(auto_secondary=True, verbosity='{0}'.format(prodyVerbosity))

        matrix = prody.parseArray(self.protocol.matrixFile.getFileName())

        if paramName == 'displayBarsSet1':
            modes =  self.modes1

            modeNumber = self.modeNumSet1.get()
            mode = modes[modeNumber]
            
            if mode is None:
                return [self.errorMessage("Invalid mode number *%d*\n"
                                        "Display the output Normal Modes to see "
                                        "the availables ones." % modeNumber,
                                        title="Invalid input")]

            # row
            row = matrix[modeNumber-1]
            cumulOverlap = self.cumulOverlapSet1.get()

            if os.path.isfile(self.protocol._getExtraPath('match_inds.txt')) and self.matchedModeNum.get():
                x_ticklabels = prody.parseArray(self.protocol._getExtraPath('match_inds.txt'), dtype=int)
            else:
                x_ticklabels = np.arange(1, len(self.modes2)+1) 

        elif paramName == 'displayBarsSet2':
            modes =  self.modes2
            
            modeNumber = self.modeNumSet2.get()
            mode = modes[modeNumber]
            
            if mode is None:
                return [self.errorMessage("Invalid mode number *%d*\n"
                                        "Display the output Normal Modes to see "
                                        "the availables ones." % modeNumber,
                                        title="Invalid input")]

            if os.path.isfile(self.protocol._getExtraPath('match_inds.txt')) and self.matchedModeNum.get():
                match_inds = list(prody.parseArray(self.protocol._getExtraPath('match_inds.txt')))
                modeNumber = match_inds.index(modeNumber)+1

            # column
            row = matrix[:, modeNumber-1]
            cumulOverlap = self.cumulOverlapSet2.get()

            x_ticklabels = np.arange(1, len(self.modes1)+1) 

        else:
            modes_sets = [self.modes1, self.modes2]
            modes = modes_sets[np.argmax([len(nmset) for nmset in modes_sets])]

            row = matrix
            cumulOverlap = self.cumulOverlapMain.get()

            if os.path.isfile(self.protocol._getExtraPath('match_inds.txt')) and modes == self.modes2 and self.matchedModeNum.get():
                x_ticklabels = prody.parseArray(self.protocol._getExtraPath('match_inds.txt'))
            else:
                x_ticklabels = np.arange(1, len(modes)+1) 

        inds = np.arange(len(row))

        if self.abs:
            row = abs(row)

        plotter = EmPlotter()
        plotter.createSubPlot('', 'mode index', 'overlap') # needed for self.bar
        plotter.plotDataBar(inds, row, width=0.8)

        ax = plotter.figure.gca()

        if np.max(matrix) <= 1:
            if self.abs:
                ax.set_ylim([0, 1])
            else:
                ax.set_ylim([-1, 1])

        if cumulOverlap:
            cum_overlaps = np.sqrt(np.power(row, 2).cumsum(axis=row.ndim-1))
            plotter.plotData(inds, cum_overlaps, 'r')
        
        ax.xaxis.set_major_formatter(IndexFormatter(x_ticklabels))

        if self.allticks.get():
            ax.xaxis.set_major_locator(ticker.IndexLocator(offset=0.5, base=1.))
        else:
            locator = ticker.AutoLocator()
            locator.set_params(integer=True)
            minor_locator = ticker.AutoMinorLocator()

            ax.xaxis.set_major_locator(locator)
            ax.xaxis.set_minor_locator(minor_locator)

        # configure ProDy to restore secondary structure information and verbosity
        prody.confProDy(auto_secondary=self.old_secondary, 
                        verbosity='{0}'.format(self.old_verbosity))

        return [plotter]

