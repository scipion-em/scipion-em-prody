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

from pyworkflow.protocol import params
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO

from pwem.objects import SetOfNormalModes
from pwem.viewers import VmdView
from pwem.viewers.plotter import EmPlotter
from prody2.protocols.protocol_lda import ProDyLDA

import prody

class ProDyLDAViewer(ProtocolViewer):
    """Visualization of results from the ProDy mode projection protocol.    
    """    
    _label = 'Projection viewer'
    _targets = [ProDyLDA]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('showPlot', params.LabelParam,
                      label='Show root mean square fluctuations (RMSFs)?',
                      help='RMSFs are plotted with key residues highlighted')
        form.addParam('percentile', params.FloatParam, default=99,
                      label='Percentile cutoff for best residues',
                      help='Values above this percentile (float from 0 to 100) of the shuffled '
                           'LDAs will be used to select most mobile residues.')

        if isinstance(self.protocol, SetOfNormalModes):
            modes = self.protocol
        else:
            modes =  self.protocol.outputModes

        form.addParam('displayVmd', params.LabelParam,
                      condition=os.path.isfile(modes._nmdFileName.get()),
                      label="Display mode color structures with VMD NMWiz?",
                      help="Use ProDy Normal Mode Wizard to view all modes in a more interactive way. "
                           "See http://prody.csb.pitt.edu/tutorials/nmwiz_tutorial/nmwiz.html")

    def _getVisualizeDict(self):
        return {'showPlot': self._viewFlucts,
                'displayVmd': self._viewAllModes}

    def _viewFlucts(self, paramName):
        """visualisation for LDA fluctuations and residues"""

        atoms = self.protocol.outputEnsemble.loadEnsemble().getAtoms()
        lda = prody.loadModel(self.protocol._getPath('modes.lda.npz'))

        plotter = EmPlotter()

        shuffled = lda.getShuffledEigvecs()
        rmsf = prody.calcRMSFlucts(lda)

        percentile = np.percentile(shuffled, self.percentile.get())
        inds = prody.calcMostMobileNodes(lda, cutoff=percentile)

        prody.showRMSFlucts(lda, atoms=atoms)
        plt.plot(inds, rmsf[inds], 'r*')
        resnums = atoms.getResnums()
        for i, ind in enumerate(inds):
            plt.text(ind - np.mod(i,2)*2, 
                     rmsf[ind], resnums[ind])
            
        plt.title('LDA RMSF with pecentile {0}'.format(self.percentile.get()))
        plt.savefig(self.protocol._getPath('lda_rmsf_best_{0}.png'.format(self.percentile.get())))

        fo = open(self.protocol._getExtraPath('lda_best_residues_{0}.txt'.format(self.percentile.get())), 'w')
        for ind in inds:
            fo.write('\t'.join(['{:5d}'.format(ind), 
                                '{:5.3f}'.format(rmsf[ind]),
                                str(atoms[ind].getResname()),
                                str(atoms[ind].getResnum())]) + '\n')
        fo.close()    

        for i, ind in enumerate(inds):
            plt.plot(inds, rmsf[inds], 'r*')
            plt.text(ind - np.mod(i,2)*2, 
                     rmsf[ind], resnums[ind])

        return [plotter]

    def _viewAllModes(self, paramName):
        """ visualisation for 2D covariance and cross-correlation matrices"""
        if paramName == 'displayVmd':
            return [createVmdNmwizView(self.protocol)]

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

def createVmdNmwizView(protocol):
    if isinstance(protocol, SetOfNormalModes):
        nmdFile = protocol._nmdFileName
    else:
        nmdFile = protocol.outputModes._nmdFileName
    return VmdView('-e %s' % nmdFile)
