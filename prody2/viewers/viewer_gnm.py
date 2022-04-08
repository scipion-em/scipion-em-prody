# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Slavica Jonic  (slavica.jonic@upmc.fr)
# *              James Krieger (jmkrieger@cnb.csic.es)
# *              Ricardo Serrano Gutiérrez (rserranogut@hotmail.com)  
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
This module implement the wrappers around ProDy GNM 
visualization programs.
"""

from pyworkflow.gui.project import ProjectWindow
from pyworkflow.protocol.params import LabelParam, IntParam, BooleanParam
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO

from pwem.viewers.plotter import EmPlotter

from pwem.viewers import ObjectView, VmdView, DataView
from pwem.emlib import MDL_NMA_ATOMSHIFT
from pwem.objects import SetOfNormalModes

from continuousflex.protocols import FlexProtNMA
from continuousflex.viewers.nma_plotter import FlexNmaPlotter

from prody2.protocols import ProDyGNM

import os

import prody
from prody.utilities.drawtools import IndexFormatter
from matplotlib.pyplot import *

OBJCMD_NMA_PLOTDIST = "Plot distance profile"
OBJCMD_NMA_VMD = "Display VMD animation"


class ProDyGNMViewer(ProtocolViewer):
    """ Visualization of results from the GNM-NMA protocol.    
        Normally, modes with high collectivity and low 
        score are preferred.
    """
    _label = 'viewer gnm'
    _targets = [ProDyGNM]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

#     def setProtocol(self, protocol):
#         ProtocolViewer.setProtocol(self, protocol)
#         inputPdb = protocol.inputStructure.get()
#         self.isEm.set(inputPdb.getPseudoAtoms())
        
    def _defineParams(self, form):

        if isinstance(self.protocol, SetOfNormalModes):
            protocol_path = os.path.dirname(os.path.dirname(self.protocol[1].getModeFile()))
            vmdFiles = protocol_path + "/extra/animations/"
            nmdFile = protocol_path + "/modes.nmd"
        else:
            vmdFiles = self.protocol._getExtraPath("animations")
            nmdFile = self.protocol._getPath("modes.nmd")

        form.addSection(label='Visualization')
  
        form.addParam('displayModes', LabelParam,
                      label="Display output Normal Modes?", important=True)
        form.addParam('displayMaxDistanceProfile', LabelParam,
                      label="Plot max distance profile?",
                      help="Maximum unitary shift of each atom or pseudoatom over all computed modes.") 

        form.addParam('displaySqFlucts', LabelParam,
                      label="Plot MSF from all the computed modes?",
                      help="The MSF residues are found on the diagonal of "
                           "the covariance matrix and describe the residue mobility or flexibility.")

        group = form.addGroup('Single mode')  
        group.addParam('modeNumber', IntParam, default=7,
                      label='Mode number')
        group.addParam('displayVmd', LabelParam,
                      condition=os.path.isdir(vmdFiles),
                      label='Display mode animation with VMD?') 
        group.addParam('displayDistanceProfile', LabelParam, default=False,
                      label="Plot mode distance profile?",
                      help="Unitary shift of each atom or pseudoatom along the mode that is requested to be animated.")
        group.addParam('displaySingleSqFluct', LabelParam, default=False,
                      label="Plot mode square fluctuation?",
                      help="Shows the single mode Square Fluctuation")
        group.addParam('cumulSqFluct', BooleanParam, default=False,
                      #condition=metric==NMA_METRIC_OVERLAP,
                      label='Display cumulative Square Fluctuation?',
                      help='Show the cumulative MSF of the modes selected or show the single MSF of the mode selected')
        
        form.addParam('displayVmd2', LabelParam,
                      condition=os.path.isfile(nmdFile),
                      label="Display mode animations with VMD NMWiz?",
                      help="Use ProDy Normal Mode Wizard to view all modes in a more interactive way. "
                           "See http://prody.csb.pitt.edu/tutorials/nmwiz_tutorial/nmwiz.html") 
        
        group = form.addGroup('GNM covariance and Cross-correlation')
        group.addParam ('displayCovMatrix', LabelParam,
                      #condition=have_matrix == True,
                      label='Display Covariance matrix?',
                      help='Covariance matrices are shown as heatmaps.')
        group.addParam ('displayCrossCorrMatrix', LabelParam,
                      #condition=have_matrix == True,
                      label='Display Cross Correlation matrix?',
                      help='Cross Correlation matrices are shown as heatmaps.')
        

    def _getVisualizeDict(self):
        return {'displayModes': self._viewParam,
                'displayMaxDistanceProfile': self._viewParam,
                'displaySqFlucts': self._viewSQF,
                'displayVmd': self._viewSingleMode,
                'displayVmd2': self._viewSingleMode,
                'displayDistanceProfile': self._viewSingleMode,
                'displaySingleSqFluct': self._viewSQF,
                'displayCovMatrix': self._viewMatrix,
                'displayCrossCorrMatrix': self._viewMatrix
                } 

    def _viewMatrix(self, paramName):
        """ visualisation for 2D covariance and Cross correlation matrices""" 
        if paramName=='displayCovMatrix':
            matrix = prody.parseArray(self.protocol.matrixFileCV.getFileName())
            title = 'Covariance matrix'
        else:
            matrix = prody.parseArray(self.protocol.matrixFileCC.getFileName())
            title = 'Cross Correlation matrix'

        plotter = EmPlotter(mainTitle=title)
        plot = prody.showMatrix(matrix, vmax=1, vmin=-1, origin='lower')
        
        return [plotter] 

    def _viewParam(self, paramName):
        if paramName == 'displayModes':
            # The following two lines display modes.sqlite file
            # modes =  self.protocol.outputModes
            # return [ObjectView(self._project, modes.strId(), modes.getFileName())]
            # The following two lines display modes.xmd file
            if isinstance(self.protocol, SetOfNormalModes):
                modes = os.path.dirname(self.protocol[1].getModeFile()) + ".xmd"
            else:
                modes =  self.protocol._getPath("modes.xmd")
            return [DataView(modes)]

        elif paramName == 'displayMaxDistanceProfile':
            if isinstance(self.protocol, SetOfNormalModes):
                fn = os.path.dirname(os.path.dirname(self.protocol[1].getModeFile())) + "/extra/maxAtomShifts.xmd"
            else:
                fn = self.protocol._getExtraPath("maxAtomShifts.xmd")
            return [createShiftPlot(fn, "Maximum atom shifts", "maximum shift")]

    def _viewSQF(self, paramName):
        """ visualization of square fluctuation for all/single mode/s. """   
        modes_path = self.protocol._getPath("modes.gnm.npz")
        modes = prody.loadModel(modes_path)
        
        cumulSqFluct = self.cumulSqFluct.get()
        modeNumber = self.modeNumber.get()

        if paramName == 'displaySqFlucts':
            plotter = EmPlotter()
            plot = prody.showSqFlucts(modes)
        
        elif paramName == 'displaySingleSqFluct':
            if cumulSqFluct:
                plotter = EmPlotter(mainTitle="Qumulative")
                plot = prody.showSqFlucts(modes[:modeNumber])
            else:
                plotter = EmPlotter()
                plot = prody.showSqFlucts(modes[modeNumber-1])

        return [plotter]

    
    def _viewSingleMode(self, paramName):
        """ visualization for a selected mode. """
        if isinstance(self.protocol, SetOfNormalModes):
            modes = self.protocol
        else:
            modes =  self.protocol.outputModes

        modeNumber = self.modeNumber.get()
        mode = modes[modeNumber]
        
        if mode is None:
            return [self.errorMessage("Invalid mode number *%d*\n"
                                      "Display the output Normal Modes to see "
                                      "the availables ones." % modeNumber,
                                      title="Invalid input")]
                                      
        elif paramName == 'displayVmd':
            return [createVmdView(self.protocol, modeNumber)]
        elif paramName == 'displayVmd2':
            return [createVmdNmwizView(self.protocol, modeNumber)]
        elif paramName == 'displayDistanceProfile':
            return [createDistanceProfilePlot(self.protocol, modeNumber)]


def createShiftPlot(mdFn, title, ylabel):
    plotter = FlexNmaPlotter()
    plotter.createSubPlot(title, 'atom index', ylabel)
    plotter.plotMdFile(mdFn, None, MDL_NMA_ATOMSHIFT)
    return plotter


def createDistanceProfilePlot(protocol, modeNumber):
    if isinstance(protocol, SetOfNormalModes):
        vectorMdFn = os.path.dirname(os.path.dirname(protocol[1].getModeFile(
        ))) + "/extra/distanceProfiles/vec%d.xmd" % modeNumber
    else:
        vectorMdFn = protocol._getExtraPath("distanceProfiles","vec%d.xmd"
                                            % modeNumber)
    plotter = createShiftPlot(vectorMdFn, "Atom shifts for mode %d"
                              % modeNumber, "shift")
    return plotter


def createVmdView(protocol, modeNumber):
    if isinstance(protocol, SetOfNormalModes):
        vmdFile = os.path.dirname(os.path.dirname(protocol[1].getModeFile())) + "/extra/animations/animated_mode_%03d.vmd" % modeNumber
    else:
        vmdFile = protocol._getExtraPath("animations", "animated_mode_%03d.vmd"
                                        % modeNumber)
    return VmdView('-e "%s"' % vmdFile)

def createVmdNmwizView(protocol, modeNumber):
    if isinstance(protocol, SetOfNormalModes):
        nmdFile = os.path.dirname(os.path.dirname(protocol[1].getModeFile())) + "/modes.nmd"
    else:
        nmdFile = protocol._getPath("modes.nmd")
    return VmdView('-e %s' % nmdFile)

def showDistanceProfilePlot(protocol, modeNumber):
    createDistanceProfilePlot(protocol, modeNumber).show()


def showVmdView(protocol, modeNumber):
    createVmdView(protocol, modeNumber).show()


ProjectWindow.registerObjectCommand(OBJCMD_NMA_PLOTDIST,
                                    showDistanceProfilePlot)
ProjectWindow.registerObjectCommand(OBJCMD_NMA_VMD,
                                    showVmdView)

