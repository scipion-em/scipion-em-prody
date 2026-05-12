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
from pyworkflow.object import Float

import prody

from prody2 import Plugin
from prody2.constants import ENSEMBLE_WEIGHTS
from prody2.objects import SetOfTrajFrames, TrajFrame

from prody2.objects import HAVE_CHEM
if HAVE_CHEM:
    from prody2.objects import DcdMDSystem as ProDyMDSystem
else:
    ProDyMDSystem = SetOfTrajFrames


class ProDyBioExcelCV19(EMProtocol):
    """
    Downloads and parses molecular dynamics trajectories from the
    BioExcel COVID-19 database.

    AI Generated:

    BioExcel CV19 Parser (ProDyBioExcelCV19) — User Manual
        Overview

        The ProDyBioExcelCV19 protocol provides an interface to access
        molecular dynamics simulations deposited in the BioExcel
        COVID-19 database.

        Its main purpose is to automatically retrieve a selected
        simulation, extract topology and trajectory files, and convert
        them into Scipion-compatible trajectory objects for downstream
        structural analysis.

        For structural biology users, this protocol is especially
        useful when working with publicly available molecular dynamics
        simulations of biologically relevant systems such as viral
        proteins, protein-ligand complexes, or conformational
        ensembles.

        Inputs and General Workflow

        The protocol requires a simulation accession identifier,
        corresponding to a BioExcel CV19 simulation entry.

        During execution, the protocol:

            1. Downloads the selected simulation data
            2. Extracts topology and trajectory files
            3. Applies optional atom selection
            4. Applies optional frame selection
            5. Creates an output trajectory object

        The downloaded files include:

            - PDB coordinate file
            - PSF topology file
            - DCD trajectory file

        Atom Selection

        Users may optionally reduce the trajectory to specific subsets
        of atoms.

        Available selections include:

            - Full structure
            - Carbon atoms only
            - Backbone atoms
            - Backbone carbon atoms only

        Biological Interpretation

        This option is useful when users want to reduce trajectory
        size, focus on protein backbone motions, or simplify
        conformational analysis.

        Frame Selection

        The protocol allows selecting only specific trajectory frames.

        Frame ranges can be defined using expressions such as:

            - 1-5
            - 10:20:2
            - 1-5,11-15

        Practical Use

        This is especially useful for:

            - Sampling representative conformations
            - Reducing computational cost
            - Focusing analysis on specific simulation intervals

        Output Options

        The protocol supports two output formats.

        MDSystem Output

        When enabled, the protocol generates a ProDyMDSystem object
        containing:

            - Coordinate file
            - Topology file
            - Trajectory file

        This output is most appropriate for workflows requiring
        direct molecular dynamics handling.

        Frame-Based Output

        Alternatively, the protocol can create a SetOfTrajFrames.

        In this mode:

            - Each frame is stored as an individual trajectory entry
            - Equal statistical weight is assigned to every frame

        This representation is useful for workflows where individual
        conformations are processed independently.

        Summary Information

        After execution, the protocol reports:

            - Total number of atoms in the trajectory
            - Number of protein residues

        Biological Perspective

        The BioExcel CV19 database contains valuable conformational
        ensembles that often capture biologically meaningful motions
        inaccessible to single static structures.

        Using atom and frame selection carefully allows users to focus
        on the most relevant structural regions while keeping
        downstream analyses computationally efficient.

        Final Perspective

        For most molecular dynamics workflows, this protocol serves as
        a convenient bridge between public BioExcel simulation
        repositories and Scipion structural analysis pipelines.

        It allows users to rapidly transform external simulation data
        into analysis-ready trajectory objects suitable for
        conformational exploration, ensemble comparison, and
        biologically meaningful dynamic interpretation.
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

    def createOutputStep(self):
        if self.useMDSystem:
            outputTrajectory = ProDyMDSystem(filename=self.PDB_FILENAME)
            outputTrajectory.setTopologyFile(self.PSF_FILENAME)
            outputTrajectory.setTrajectoryFile(self.DCD_FILENAME)
        else:
            numFrames = prody.DCDFile(self.DCD_FILENAME).numFrames()

            outputTrajectory = SetOfTrajFrames().create(self._getExtraPath())
            outputTrajectory.setTopologyFile(self.PSF_FILENAME)
            outputTrajectory.setTrajectoryFile(self.DCD_FILENAME)
            for j in range(numFrames):
                frame = TrajFrame((j+1, self.DCD_FILENAME), objLabel=str(j+1))
                setattr(frame, ENSEMBLE_WEIGHTS, Float(1/numFrames))
                outputTrajectory.append(frame)

        self._defineOutputs(outputTrajectory=outputTrajectory)

    def _summary(self):
        if not hasattr(self, 'outputTrajectory'):
            summ = ['Output trajectory not ready yet']
        else:
            outputAg = prody.parsePDB(self._getExtraPath(
                '{0}.pdb'.format(self.accession.get())))
            summ = ['Trajectory has *{0}* atoms including *{1}* protein residues'.format(
                outputAg.numAtoms(), outputAg.ca.numAtoms())]
        return summ
