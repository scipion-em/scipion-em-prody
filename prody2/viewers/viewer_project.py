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
for plotting projections of ensembles onto modes or distance distributions.
"""
import matplotlib.pyplot as plt
import numpy as np

import os

from pyworkflow.protocol.params import LabelParam, BooleanParam, FloatParam, IntParam
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO

from pwem.viewers.plotter import EmPlotter
from pwem.objects import SetOfAtomStructs, Set

from prody2.protocols.protocol_project import ProDyProject, ONE, TWO, THREE
from prody2.protocols.protocol_measure import ProDyMeasure

import prody

class ProDyProjectionsViewer(ProtocolViewer):
    """Visualization of results from the ProDy mode projection protocol.    
    """    
    _label = 'Projection viewer'
    _targets = [ProDyProject, ProDyMeasure]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    def _defineParams(self, form):
                
        # configure ProDy to automatically handle secondary structure information and verbosity
        from pyworkflow import Config
        global oldSecondary; oldSecondary = prody.confProDy("auto_secondary")
        global oldVerbosity; oldVerbosity = prody.confProDy("verbosity")
        prodyVerbosity =  'none' if not Config.debugOn() else 'debug'
        prody.confProDy(auto_secondary=True, verbosity='{0}'.format(prodyVerbosity))


        self.isProjection = isinstance(self.protocol, ProDyProject)
        if self.isProjection:
            self.numModes = self.protocol.numModes.get()
        else:
            self.numModes = ONE

        form.addSection(label='Visualization')

        form.addParam('showPlot', LabelParam,
                      label='Show plot?',
                      help='Projections or measures are shown in various ways depending on the options selected')

        form.addParam('norm', BooleanParam, label="Normalize?", default=False,
                      help='Select whether to normalise projections.',
                      condition=self.isProjection)

        form.addParam('rmsd', BooleanParam, label="RMSD scale?", default=True,
                      help='Select whether to scale projections to RMSDs.',
                      condition=self.isProjection)
        
        form.addParam('label', BooleanParam, label="Label points?", default=False,
                      help='Select whether to label points.',
                      condition=(self.numModes!=ONE and self.isProjection))
        
        form.addParam('adjustText', BooleanParam, label="Adjust labels?", default=False,
                      help='Select whether to adjust labels on points to not overlap.',
                      condition="label==True")

        form.addParam('density', BooleanParam, label="Show density?",
                      default=False, condition=(self.numModes != THREE and self.isProjection),
                      help='Select whether to use a 1D histogram or 2D kernel density estimation from seaborn.\n'
                           'The alternative is to show points for 2D and an ordered series in 1D.')
        
        form.addParam('points', BooleanParam, label="Show points?",
                      default=False, condition='density == True and numModes == %d and isProjection == True' % TWO,
                      help='Select whether to show points too.')
        
        form.addParam('useWeights', BooleanParam, label="Use cluster weights?",
                      default=False, condition=self.numModes != THREE,
                      help='Select whether to use cluster weights to rescale 1D histogram or 2D kernel density estimate'
                           ' or point size.\n')

        form.addParam('separatePlots', BooleanParam, label="Separate plots?",
                      default=False, condition=self.numModes != THREE,
                      help='Select whether to use a 1D histogram or 2D kernel density estimation from seaborn.\n'
                           'The alternative is to show points for 2D and an ordered series in 1D.')
        
        form.addParam('bins', IntParam, label="number of bins", default=-1,
                        help='Enter a number of bins here, which should be positive.\n'
                             '-1 is the dummy value and needs changing to have an effect.\n'
                             'This option will not do anything for line plots.',
                        condition="numModes==%d" % ONE)
        
        groupX = form.addGroup('xrange')
        groupX.addParam('xrange1', FloatParam, label="x-axis range limit 1 for bins", default=-1,
                        help='Enter values here and below to specify x-axis limits for bins.\n'
                             '-1 is the dummy value and needs changing to e.g. -1.1 '
                             'in both places to have an effect',
                        condition=self.numModes==ONE)
        groupX.addParam('xrange2', FloatParam, label="x-axis range limit 2 for bins", default=-1,
                        help='Enter values here and above to specify x-axis limits for bins.\n'
                             '-1 is the dummy value and needs changing to e.g. -1.1 '
                             'in both places to have an effect',
                        condition=self.numModes==ONE)

        groupX = form.addGroup('xlim')
        groupX.addParam('xlim1', FloatParam, label="x-axis limit 1", default=-1,
                        help='Enter values here and below to specify x-axis limits.\n'
                             '-1 is the dummy value and needs changing to e.g. -1.1 '
                             'to have an effect',
                        condition=self.numModes!=ONE)
        groupX.addParam('xlim2', FloatParam, label="x-axis limit 2", default=-1,
                        help='Enter values here and above to specify x-axis limits.\n'
                             '-1 is the dummy value and needs changing to e.g. -1.1 '
                             'to have an effect',
                        condition=self.numModes!=ONE)

        groupY = form.addGroup('ylim')
        groupY.addParam('ylim1', FloatParam, label="y-axis limit 1", default=-1,
                        help='Enter values here and below to specify y-axis limits.\n'
                             '-1 is the dummy value and needs changing to e.g. -1.1 '
                             'to have an effect',
                        condition=self.numModes!=ONE)
        groupY.addParam('ylim2', FloatParam, label="y-axis limit 2", default=-1,
                        help='Enter values here and above to specify y-axis limits.\n'
                             '-1 is the dummy value and needs changing to e.g. -1.1 '
                             'to have an effect',
                        condition=self.numModes!=ONE)

        groupZ = form.addGroup('zlim')
        groupZ.addParam('zlim1', FloatParam, label="z-axis limit 1", default=-1,
                        help='Enter values here and below to specify z-axis limits.\n'
                             '-1 is the dummy value and needs changing to e.g. -1.1 '
                             'to have an effect',
                        condition=self.numModes==THREE)
        groupZ.addParam('zlim2', FloatParam, label="z-axis limit 2", default=-1,
                        help='Enter values here and above to specify z-axis limits.\n'
                             '-1 is the dummy value and needs changing to e.g. -1.1 '
                             'to have an effect',
                        condition=self.numModes==THREE)

    def _getVisualizeDict(self):
        return {'showPlot': self._viewProjection}            

    def _viewProjection(self, paramName):
        """visualisation for all projections or measures"""

        inputEnsemble = self.protocol.inputEnsemble

        if isinstance(inputEnsemble.get(), Set):
            inputEnsemble = [inputEnsemble]

        if self.isProjection:
            modesPath = self.protocol.inputModes.get().getFileName()
            modes = prody.parseScipionModes(modesPath)
        else:
            extraPath = self.protocol._getExtraPath()
            outFiles = os.listdir(extraPath)

        for i, ensPointer in enumerate(inputEnsemble):
            ens = ensPointer.get()
            
            if not self.isProjection:
                inputClass = type(ens)
                ens = inputClass(filename=extraPath+'/'+outFiles[i])
                measures = np.array([float(item._prodyMeasures) for item in ens], dtype=float)

            if isinstance(ens, SetOfAtomStructs):
                # the ensemble gets built exactly as the input is setup and nothing gets rejected
                ags = prody.parsePDB([tarStructure.getFileName() for tarStructure in ens])
                ensemble = prody.buildPDBEnsemble(ags, match_func=prody.sameChainPos, seqid=0., overlap=0., superpose=False, mapping=None)
            else:
                ensemble = ens.loadEnsemble()

            ensemble = self._cleanLabels(ensemble)

            if i == 0 or self.separatePlots.get():
                plotter = EmPlotter()

            c = plt.rcParams['axes.prop_cycle'].by_key()['color'][i]

            sizes = ensemble.getData('size')
            if sizes is None:
                weights = np.ones(ensemble.numCoordsets())
            else:
                if self.density.get() or sizes.max() == 1 or not self.isProjection:
                    weights = sizes
                else:
                    weights = sizes * 100/sizes.max()

            if self.numModes == ONE:
                bins = self.bins.get()
                if bins < 1:
                    bins = None

                if self.xrange1.get() != -1 and self.xrange2.get() != -1:
                    xrange = (self.xrange1.get(), self.xrange2.get())
                else:
                    xrange = None

                if self.isProjection:
                    density = self.density.get()
                    if density:
                        prody.showProjection(ensemble, modes[:self.protocol.numModes.get()+1],
                                            rmsd=self.rmsd.get(), norm=self.norm.get(),
                                            show_density=True, c=c,
                                            use_weights=self.useWeights.get(), weights=weights,
                                            bins=bins, range=xrange)
                    else:
                        prody.showProjection(ensemble, modes[:self.protocol.numModes.get()+1],
                                            rmsd=self.rmsd.get(), norm=self.norm.get(),
                                            show_density=False, c=c,
                                            use_weights=self.useWeights.get(), weights=weights)
                else:
                    if not self.useWeights.get():
                        weights = None
                        
                    plt.hist(measures, weights=weights, bins=bins, range=xrange)
            else:
                if self.label.get():
                    prody.showProjection(ensemble, modes[:self.protocol.numModes.get()+1],
                                         text=ensemble.getLabels(),
                                         rmsd=self.rmsd.get(), norm=self.norm.get(),
                                         show_density=self.density.get(), 
                                         adjust=self.adjustText.get(), c=c,
                                         use_weights=self.useWeights.get(), weights=weights)
                else:
                    prody.showProjection(ensemble, modes[:self.protocol.numModes.get()+1],
                                         rmsd=self.rmsd.get(), norm=self.norm.get(),
                                         show_density=self.density.get(), 
                                         adjust=self.adjustText.get(), c=c,
                                         use_weights=self.useWeights.get(), weights=weights)
                    
                if self.points.get():
                    prody.showProjection(ensemble, modes[:self.protocol.numModes.get()+1],
                                         rmsd=self.rmsd.get(), norm=self.norm.get(),
                                         show_density=False, 
                                         adjust=self.adjustText.get(), c=c,
                                         use_weights=self.useWeights.get(), weights=weights)                    

            ax = plotter.figure.gca()
            
            if self.xlim1.get() != -1 or self.xlim2.get() != -1:
                xlims = ax.get_xlim()
                if self.xlim1.get() != -1 and self.xlim2.get() == -1:
                    ax.set_xlim([self.xlim1.get(), xlims[1]])
                elif self.xlim1.get() == -1 and self.xlim2.get() != -1:
                    ax.set_xlim([xlims[0], self.xlim2.get()])
                else:
                    ax.set_xlim([self.xlim1.get(), self.xlim2.get()])

            if self.ylim1.get() != -1 or self.ylim2.get() != -1:
                ylims = ax.get_ylim()
                if self.ylim1.get() != -1 and self.ylim2.get() == -1:
                    ax.set_ylim([self.ylim1.get(), ylims[1]])
                elif self.ylim1.get() == -1 and self.ylim2.get() != -1:
                    ax.set_ylim([ylims[0], self.ylim2.get()])
                else:
                    ax.set_ylim([self.ylim1.get(), self.ylim2.get()])

            if self.zlim1.get() != -1 or self.zlim2.get() != -1:
                zlims = ax.get_zlim()
                if self.zlim1.get() != -1 and self.zlim2.get() == -1:
                    ax.set_zlim([self.zlim1.get(), zlims[1]])
                elif self.zlim1.get() == -1 and self.zlim2.get() != -1:
                    ax.set_zlim([zlims[0], self.zlim2.get()])
                else:
                    ax.set_zlim([self.zlim1.get(), self.zlim2.get()])
                
        # configure ProDy to restore secondary structure information and verbosity
        prody.confProDy(auto_secondary=oldSecondary, verbosity='{0}'.format(oldVerbosity))

        return [plotter]

    def _cleanLabels(self, ensemble):

        if ensemble.getLabels()[0].find('Selection') != -1:
            ensemble._labels = [label.split('Selection')[0] for label in ensemble.getLabels()]

        if ensemble.getLabels()[0].endswith('_atoms_amap'):
            ensemble._labels = [label[:-11] for label in ensemble.getLabels()]

        if ensemble.getLabels()[0].endswith('_ca'):
            ensemble._labels = [label[:-3] for label in ensemble.getLabels()]

        if ensemble.getLabels()[0][:6].isnumeric():
            ensemble._labels = [str(int(label[:6])) for label in ensemble.getLabels()]

        if ensemble.getLabels()[0].startswith('Unknown_m'):
            ensemble._labels = [label.split('Unknown_m')[-1] for label in ensemble.getLabels()]

        if ensemble.getLabels()[0][5:12] == 'atoms_m' and ensemble.getLabels()[1][5:12] == 'atoms_m':
            ensemble._labels = [label.split('atoms_m')[-1] for label in ensemble.getLabels()]

        return ensemble