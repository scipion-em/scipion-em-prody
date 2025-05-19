# **************************************************************************
# *
# * Authors:     James Krieger (jmkrieger@cnb.csic.es)
# *
# * Centro Nacional de Biotecnologia, CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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

import os
import prody2

def getProDyEnvName(version):
    return "prody-%s" % version

DEVEL = 'github'
LATEST = 'master'
RELEASE = '2.4.1'
VERSIONS = [DEVEL]
PRODY_DEFAULT_VER_NUM = DEVEL

PRODY_ENV_ACT = "PRODY_ENV_ACT"

PROJ_COEFFS = "_prodyProjCoefficients"
ENSEMBLE_WEIGHTS = "_prodyWeights"
MEASURES = "_prodyMeasures"

N_ATOMS = "_numAtoms"
N_RESIDUES = "_numResidues"
N_CHAINS = "_numChains"

FIRST_RESNUM = "_firstResnum"
LAST_RESNUM = "_lastResnum"
MAX_RESNUM = "_maxResnum"
MIN_RESNUM = "_minResnum"

PRODY_FRACT_VARS = "_prodyFractVars"

PRODY_SCRIPTS = os.path.join(os.path.dirname(prody2.__file__),
                             "protocols", "scripts")
PRODY_TESTFILE = os.path.join(os.path.dirname(prody2.__file__),
                             "protocols", "tests", "pdb4ake_fixed")


# chain matching methods
BEST_MATCH = 0
SAME_CHID = 1
SAME_POS = 2
CUSTOM = 3

# residue mapping methods
NOTHING = 0 # stop trivial mapping if trivial mapping fails
PWALIGN = 1 # biopython pwalign local pairwise sequence alignment after trivial mapping
CEALIGN = 2 # combinatorial extension (CE) as in PyMOL
DEFAULT = 3 # try pwalign then CE
