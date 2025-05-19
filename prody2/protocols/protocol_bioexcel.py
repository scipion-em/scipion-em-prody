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
This module will provide the ProDy interface for parsing files from 
the BioExcel CV19 database
"""
from pwem.protocols import EMProtocol
from pyworkflow.protocol import params
from pyworkflow.object import Float, Integer

from prody2 import Plugin
from prody2.constants import (ENSEMBLE_WEIGHTS, N_FRAMES,
                              N_ATOMS, N_RESIDUES, N_CHAINS)
from prody2.objects import SetOfTrajFrames, TrajFrame

from prody2.objects import HAVE_CHEM
if HAVE_CHEM:
    from prody2.objects import DcdMDSystem as ProDyMDSystem
else:
    ProDyMDSystem = SetOfTrajFrames


class ProDyBioExcelCV19(EMProtocol):
    """
    This module will provide the ProDy interface for parsing files from 
    the BioExcel CV19 database
    """
    _label = 'BioExcelCV19'
    _possibleOutputs = {'outputTrajectory': ProDyMDSystem}

    selections = [None, '_C', 'backbone', 'backbone and _C']
    NONE = 0
    _C = 1

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        form.addSection(label='ProDy BioExcelCV19 parsing')

        form.addParam('accession', params.StringParam, label="Simulation accession",
                      important=True, help='Simulation accession or ID number.')

        form.addParam('selection', params.EnumParam, choices=self.selections,
                      label="Atom selection",
                      default=self._C,
                      help='Select carbon atoms or backbone or only backbone carbon atoms')
        
        form.addParam('frames', params.StringParam, label='Frame selection',
                      default='1-5,11-15', help='Please specify frames as ranges like '
                                                '1-5 or 10:20:2 separated by commas')
        
        form.addParam('useMDSystem', params.BooleanParam, label='Output MDSystem',
                      condition=HAVE_CHEM, default=HAVE_CHEM, 
                      help='Select whether to output an MDSystem object. If not, a SetOfTrajFrames will be used.')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('computeStep')
        self._insertFunctionStep('createOutputStep')

    def computeStep(self):
        accession = self.accession.get()
        args = '--accession {0} --folder {1}'.format(accession, 
                                                     self._getExtraPath())
        
        if self.selection.get() != self.NONE:
            selection = self.selections[self.selection.get()]
            args += ' --selection {0}'.format(selection)

        if self.frames.get() != self.NONE:
            args += ' --frames {0}'.format(self.frames.get())

        self.runJob(Plugin.getProgram('bioexcel.py', script=True), args)

        self.PDB_FILENAME = self._getExtraPath('{0}.pdb'.format(accession))
        self.PSF_FILENAME = self._getExtraPath('{0}.psf'.format(accession))
        self.DCD_FILENAME = self._getExtraPath('{0}.dcd'.format(accession))

        args = '{0} -n --pdb {1} > {2}'.format(self.DCD_FILENAME,
                                               self.PDB_FILENAME,
                                               self._getExtraPath('nframes.txt'))
        self.runJob(Plugin.getProgram('catdcd'), args)

    def createOutputStep(self):
        with open(self._getExtraPath('nframes.txt'), 'r') as fi:
            numFrames = int(fi.readlines()[0])

        if self.useMDSystem:
            outputTrajectory = ProDyMDSystem(filename=self.PDB_FILENAME)
        else:
            outputTrajectory = SetOfTrajFrames().create(self._getExtraPath())
            for j in range(numFrames):
                frame = TrajFrame((j+1, self.DCD_FILENAME), objLabel=str(j+1))
                setattr(frame, ENSEMBLE_WEIGHTS, Float(1/numFrames))
                outputTrajectory.append(frame)

        outputTrajectory.setOriStructFile(self.PDB_FILENAME)
        outputTrajectory.setTopologyFile(self.PSF_FILENAME)
        outputTrajectory.setTrajectoryFile(self.DCD_FILENAME)

        with open(self._getExtraPath('pdb_data.txt'), 'r') as fi:
            line = fi.readlines()[0]

        self.pdbFileName, numAtoms, numResidues, numChains = line.split('\t')
        setattr(outputTrajectory, N_FRAMES, Integer(numFrames))
        setattr(outputTrajectory, N_ATOMS, Integer(numAtoms))
        setattr(outputTrajectory, N_RESIDUES, Integer(numResidues))
        setattr(outputTrajectory, N_CHAINS, Integer(numChains))

        self._defineOutputs(outputTrajectory=outputTrajectory)

    def _summary(self):
        if not hasattr(self, 'outputTrajectory'):
            summ = ['Output trajectory not ready yet']
        else:
            summ = ['Trajectory has *{0}* frames with *{1}* atoms including *{2}* residues'.format(
                self.outputTrajectory.getAttributeValue(N_FRAMES),
                self.outputTrajectory.getAttributeValue(N_ATOMS),
                self.outputTrajectory.getAttributeValue(N_RESIDUES))]
        return summ
