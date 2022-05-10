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

from pyworkflow.utils import *
from pyworkflow.gui.project import ProjectWindow
from pyworkflow.protocol.params import LabelParam, IntParam, BooleanParam
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO

from pwem import *
from pwem.viewers.plotter import EmPlotter
from pwem.viewers import ObjectView, VmdView, DataView
from pwem.objects import SetOfNormalModes
from pwem.emlib import MetaData, MDL_NMA_ATOMSHIFT

from prody2.protocols import ProDyGNM

import os

import prody
from prody.utilities.drawtools import IndexFormatter
from matplotlib.pyplot import *

OBJCMD_NMA_PLOTDIST = "Plot distance profile"


class ProDyGNMViewer(ProtocolViewer):
    """ Visualization of results from the GNM-NMA protocol.    
        Normally, modes with high collectivity and low 
        score are preferred.
    """
    _label = 'viewer gnm'
    _targets = [ProDyGNM]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
        
    def _defineParams(self, form):

        if isinstance(self.protocol, SetOfNormalModes):
            protocol_path = os.path.dirname(os.path.dirname(self.protocol[1].getModeFile()))
            nmdFile = protocol_path + "/modes.nmd"
        else:
            nmdFile = self.protocol._getPath("modes.nmd")

        form.addSection(label='Visualization')

        group = form.addGroup('All non zero modes')
        group.addParam('displayModes', LabelParam,
                      label="Display output Normal Modes?", important=True)
        group.addParam('displayMaxDistanceProfile', LabelParam,
                      label="Plot max distance profile?",
                      help="Maximum unitary shift of each atom or pseudoatom over all computed modes.") 
        group.addParam('displaySqFlucts', LabelParam,
                      label="Plot MSF from all the computed modes?",
                      help="The mean square fluctuations (MSF) for each residue are found on the diagonal of "
                           "the covariance matrix and describe the residue mobility or flexibility.")
        group.addParam('displayRMSFlucts', LabelParam,
                      label="Plot RMSF from all the computed modes?",
                      help="Plot the root mean square fluctuations (RMSF) for all the computed modes.")
        group.addParam ('displayCovMatrix', LabelParam,
                      label='Display Covariance matrix?',
                      help='Raw covariance matrices are shown as heatmaps.')
        group.addParam ('displayCrossCorrMatrix', LabelParam,
                      label='Display Cross-Correlation matrix?',
                      help='Normalised, orientational cross-correlation matrices are shown as heatmaps.')
        
        group = form.addGroup('Single mode')  
        group.addParam('modeNumber', IntParam, default=2,
                      label='Mode number')
        group.addParam('displaySingleMode', LabelParam, default=False,
                      label="Plot mode shape?",
                      help="Shows the mode shape for a single GNM mode. Residues with positive and negative "
                            "values move in opposite directions.")
        group.addParam('overlaychains', BooleanParam, default=True,
                      label="Show overlay chains",
                      help="Choose whether to show chains as overlaid curves of different colours or one after the other with bars underneath to "
                        "indicate them. Different options may be better for different data.")
        

        group = form.addGroup('Modes range')  
        group.addParam('modeNumber1', IntParam, default=2,
                      label='Initial mode number')
        group.addParam('modeNumber2', IntParam, default=2,
                      label='Final mode number')
        group.addParam('displayRangeSqFluct', LabelParam, default=False,
                      label="Plot range mean square fluctuation?",
                      help="Shows the cumulative Mean Square Fluctuations of the range of modes selected, to calculate the Square fluctuations "
                      "of a single mode, just put the same number in the range boxes.")
        group.addParam('displayRangeRMSFluct', LabelParam, default=False,
                      label="Plot range root mean square fluctuation?",
                      help="Shows the cumulative Root Mean Square Fluctuations of the range of modes selected. To calculate the Root Square fluctuations "
                      "of a single mode, just put the same number in the range boxes.")

        form.addParam('displayVmd2', LabelParam,
                      condition=os.path.isfile(nmdFile),
                      label="Display mode color structures with VMD NMWiz?",
                      help="Use ProDy Normal Mode Wizard to view all modes in a more interactive way. "
                           "See http://prody.csb.pitt.edu/tutorials/nmwiz_tutorial/nmwiz.html") 
        
    def _getVisualizeDict(self):
        modes =  self.protocol.outputModes
        modes_path = os.path.dirname(os.path.dirname(modes[1].getModeFile()))
        self.modes = prody.parseScipionModes(modes_path, pdb=glob(modes_path+"/*atoms.pdb"))
        self.atoms = prody.parsePDB(glob(modes_path+"/*atoms.pdb"))
        return {'displayModes': self._viewParam,
                'displayMaxDistanceProfile': self._viewParam,
                'displaySqFlucts': self._viewSQF,
                'displayRMSFlucts': self._viewSQF,
                'displayVmd2': self._viewAllModes,
                'displayRangeSqFluct': self._viewSQF,
                'displayRangeRMSFluct': self._viewSQF,
                'displayCovMatrix': self._viewAllModes,
                'displayCrossCorrMatrix': self._viewAllModes,
                'displaySingleMode': self._viewSingleMode
                } 

    def _viewAllModes(self, paramName):
        """ visualisation for 2D covariance and Cross correlation matrices""" 
        if paramName == 'displayVmd2':
            return [createVmdNmwizView(self.protocol)]
        
        elif paramName=='displayCovMatrix':
            matrix = prody.parseArray(self.protocol.matrixFileCV.getFileName())
            title = 'Covariance matrix'

        else:
            matrix = prody.parseArray(self.protocol.matrixFileCC.getFileName())
            title = 'Cross Correlation matrix'

        plotter = EmPlotter(mainTitle=title)
        plot = prody.showAtomicMatrix(matrix, origin='lower', atoms=self.atoms)
        
        return [plotter] 

    def _viewParam(self, paramName):
        if paramName == 'displayModes':
            modes =  self.protocol._getPath("modes.xmd")
            return [DataView(modes)]

        elif paramName == 'displayMaxDistanceProfile':
            fn = self.protocol._getExtraPath("maxAtomShifts.xmd")
            return [createShiftPlot(fn, "Maximum atom shifts", "maximum shift")]

    def _viewSQF(self, paramName):
        """ visualization of square fluctuation for all the modes or the range of modes selected. """   
        modeNumber1 = self.modeNumber1.get()-1
        modeNumber2 = self.modeNumber2.get() # no subtraction as end of range

        if modeNumber1+1 > modeNumber2:
            return [self.errorMessage("Invalid mode range\n"
                                      "Initial mode number can not be " 
                                      "bigger than the final one.", title="Invalid input")]
        elif modeNumber1 < 1:
            return [self.errorMessage("Invalid mode range\n"
                                      "Initial mode number can not be " 
                                      "smaller than 2.", title="Invalid input")]

        elif modeNumber2 < 2:
            return [self.errorMessage("Invalid mode range\n"
                                      "Final mode number can not be " 
                                      "smaller than 2.", title="Invalid input")]
        
        try:
            mode1 = self.modes[modeNumber1] 
        except IndexError:
            return [self.errorMessage("Invalid initial mode number *%d*\n"
                                    "Display the output Normal Modes to see "
                                    "the availables ones." % modeNumber1+1,
                                    title="Invalid input")] 
        try:
            mode2 = self.modes[modeNumber2-1] 
        except IndexError:
            return [self.errorMessage("Invalid final mode number *%d*\n"
                                    "Display the output Normal Modes to see "
                                    "the availables ones." % modeNumber2,
                                    title="Invalid input")]  
                
        plotter = EmPlotter()

        if paramName == 'displaySqFlucts':
            plot = prody.showSqFlucts(self.modes[1:], atoms=self.atoms)
        elif paramName == 'displayRMSFlucts':
            plot = prody.showRMSFlucts(self.modes[1:], atoms=self.atoms)
        else:            
            if modeNumber1 == modeNumber2:
                if paramName == 'displayRangeSqFluct':
                    plot = prody.showSqFlucts(self.modes[modeNumber1], atoms=self.atoms)
                elif paramName == 'displayRangeRMSFluct':
                    plot = prody.showRMSFlucts(self.modes[modeNumber1], atoms=self.atoms)
            else:
                if paramName == 'displayRangeSqFluct':
                    plot = prody.showSqFlucts(self.modes[modeNumber1:modeNumber2], atoms=self.atoms)
                elif paramName == 'displayRangeRMSFluct':
                    plot = prody.showRMSFlucts(self.modes[modeNumber1:modeNumber2], atoms=self.atoms)

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
        
        if paramName == 'displaySingleMode':
            mode = self.modes[modeNumber-1]

            plotter = EmPlotter()
            plot = prody.showMode(mode, atoms=self.atoms, overlay_chains=self.overlaychains)        
            return [plotter]

def createShiftPlot(mdFn, title, ylabel):
    plotter = EmPlotter()
    plotter.createSubPlot(title, 'atom index', ylabel)
    
    md = MetaData(mdFn)
    yy = []
    for objId in md:
        yy.append(md.getValue(MDL_NMA_ATOMSHIFT, objId))

    xx = list(range(len(yy)))
        
    plotter.plotData(xx, yy)

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


def createVmdNmwizView(protocol):
    if isinstance(protocol, SetOfNormalModes):
        nmdFile = os.path.dirname(os.path.dirname(protocol[1].getModeFile())) + "/modes.nmd"
    else:
        nmdFile = protocol._getPath("modes.nmd")
    return VmdView('-e %s' % nmdFile)

def showDistanceProfilePlot(protocol, modeNumber):
    createDistanceProfilePlot(protocol, modeNumber).show()


ProjectWindow.registerObjectCommand(OBJCMD_NMA_PLOTDIST,
                                    showDistanceProfilePlot)

