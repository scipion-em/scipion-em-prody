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
for plotting projections of ensembles onto modes.
"""

from pyworkflow.protocol.params import LabelParam, IntParam, BooleanParam
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO

from pwem.viewers.plotter import EmPlotter

from prody2.protocols.protocol_project import ProDyProject, ONE, TWO, THREE

from matplotlib import ticker
import numpy as np
import os

import prody
from prody.utilities.drawtools import IndexFormatter

class ProDyProjectionsViewer(ProtocolViewer):
    """Visualization of results from the ProDy mode projection protocol.    
    """    
    _label = 'ProDy mode comparisons viewer'
    _targets = [ProDyProject]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    def _defineParams(self, form):
                
        # configure ProDy to automatically handle secondary structure information and verbosity
        from pyworkflow import Config
        global old_secondary; old_secondary = prody.confProDy("auto_secondary")
        global old_verbosity; old_verbosity = prody.confProDy("verbosity")
        prodyVerbosity =  'none' if not Config.debugOn() else 'debug'
        prody.confProDy(auto_secondary=True, verbosity='{0}'.format(prodyVerbosity))

        self.numModes = self.protocol.numModes.get()

        form.addSection(label='Visualization')

        form.addParam('showProjection', LabelParam,
                      label='Show projection?',
                      help='Projections are shown in various ways depending on the options selected')

        form.addParam('byFrame', BooleanParam, label="Plot by frame?",
                      condition=self.numModes==ONE, default=False,
                      help='Select whether to plot a line graph with the frame number on the x-axis '
                           'and the mode projection on the y-axis.\n'
                           'The alternative is to show a histogram for 1D.')   

        form.addParam('norm', BooleanParam, label="Normalize?", default=False,
                      help='Select whether to normalise projections.')   

        form.addParam('rmsd', BooleanParam, label="RMSD scale?", default=True,
                      help='Select whether to scale projections to RMSDs.')
        
        form.addParam('label', BooleanParam, label="Label points?", default=True,
                      help='Select whether to label points.',
                      condition=self.numModes!=ONE)   

        # form.addParam('kde', BooleanParam, label="Use kernel density estimation?",
        #               condition=(numModes != THREE and not byFrame),
        #               help='Select whether to use kernel density estimation from seaborn.\n'
        #                    'The alternative is to show points for 2D and a regular histogram for 1D.')

    def _getVisualizeDict(self):
        return {'showProjection': self._viewProjection}            

    def _viewProjection(self, paramName):
        """visualisation for all projections""" 
        #proj = prody.parseArray(self.protocol.outputProjection.getFileName())

        ags = prody.parsePDB([tarStructure.getFileName() for tarStructure in self.protocol.inputEnsemble.get()])
        ensemble = prody.buildPDBEnsemble(ags, match_func=prody.sameChainPos, seqid=0., overlap=0.)
        # the ensemble gets built exactly as the input is setup and nothing gets rejected
        
        if ensemble.getLabels()[0].endswith('_atoms_amap'):
            ensemble._labels = [label[:-11] for label in ensemble.getLabels()]
            
        if ensemble.getLabels()[0][:6].isnumeric():
            ensemble._labels = [str(int(label[:6])) for label in ensemble.getLabels()]

        modes_path = self.protocol.inputModes.get().getFileName()
        modes = prody.parseScipionModes(modes_path)

        plotter = EmPlotter()
        plotter.createSubPlot('', '', '')
        ax = plotter.figure.gca()

        if self.numModes == ONE:
            plot = prody.showProjection(ensemble, modes[:self.protocol.numModes.get()+1], 
                                        by_time=self.byFrame.get(),
                                        rmsd=self.rmsd.get(), norm=self.norm.get())#, kde=self.kde.get())
        else:
            if self.label.get():
                plot = prody.showProjection(ensemble, modes[:self.protocol.numModes.get()+1], 
                                            text=ensemble.getLabels(), 
                                            rmsd=self.rmsd.get(), norm=self.norm.get())#, kde=self.kde.get())
            else:
                plot = prody.showProjection(ensemble, modes[:self.protocol.numModes.get()+1], 
                                            rmsd=self.rmsd.get(), norm=self.norm.get())#, kde=self.kde.get())
                
        # configure ProDy to restore secondary structure information and verbosity
        prody.confProDy(auto_secondary=old_secondary, verbosity='{0}'.format(old_verbosity))

        return [plotter]
