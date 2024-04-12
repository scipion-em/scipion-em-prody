# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     James Krieger (jmkrieger@cnb.csic.es)
# *              Ricardo Serrano Gutiérrez (rserranogut@hotmail.com)                 
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
This module will provide ProDy Dynamical Domain Decomposition using the Gaussian Network Modeling (GNM).
"""

import os

from pwem.objects import SetOfNormalModes, AtomStruct, EMFile, String
from pwem.protocols import EMProtocol

from pyworkflow.utils import glob
from pyworkflow.protocol.params import PointerParam, IntParam

import prody
from prody2 import fixVerbositySecondary, restoreVerbositySecondary


class  ProDyDomainDecomp(EMProtocol):
    """
    This protocol will perform dynamical domain decomposition
    """
    _label = 'Domain Decomposition'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need params to belong to a section:
        form.addSection(label='ProDy DomainDecomp')

        form.addParam('modesGNM', PointerParam, label="Input SetOfNormalModes GNM",
                      important=True,
                      pointerClass='SetOfNormalModes',
                      help='The input SetOfNormalModes can be from an atomic model '
                           '(true PDB) or a pseudoatomic model '
                           '(an EM volume compared into pseudoatoms).\n'
                           'The set must be GNM modes')
        form.addParam('modeNumber', IntParam, default=2,
                label='Number of modes')
    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('computeDecompStep')
        self._insertFunctionStep('createOutputStep')

    def computeDecompStep(self):
        fixVerbositySecondary(self)

        modesPath = os.path.dirname(os.path.dirname(self.modesGNM.get()[1].getModeFile()))

        modes = prody.parseScipionModes(self.modesGNM.get().getFileName(),
                                            pdb=glob(modesPath+"/*atoms.pdb"))

        numModes = self.modeNumber.get()
        
        try:
            mode = modes[:numModes] 
        except IndexError:
            return [self.errorMessage("Invalid number of modes *%d*\n"
                                      "Display the output Normal Modes to see "
                                      "the availables ones." % numModes,
                                      title="Invalid input")] 
        atoms = prody.parsePDB(glob(modesPath+"/*atoms.pdb"))

        domains = prody.calcGNMDomains(mode)

        self.pdbFilename = self._getPath("atoms.pdb")
        prody.writePDB(self.pdbFilename, atoms, beta=domains)
    
        restoreVerbositySecondary(self)

    def createOutputStep(self):
        fhCmd=open(self._getPath("domains.vmd"),'w')
        fhCmd.write("mol new %s\n" % self.pdbFilename)
        fhCmd.write("mol modcolor 0 0 Beta\n")
        fhCmd.write("mol modstyle 0 0 Beads\n")
        fhCmd.close()

        outputPdb = AtomStruct()
        outputPdb.setFileName(self.pdbFilename)

        outputvmd = EMFile()
        outputvmd.setFileName(self._getPath("domains.vmd"))
        
        self._defineOutputs(outputStructure=outputPdb, outputvmd=outputvmd)
