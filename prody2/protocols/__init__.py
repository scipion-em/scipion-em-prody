# -*- coding: utf-8 -*-
# **************************************************************************
# Module to declare protocols
# Find documentation here: https://scipion-em.github.io/docs/docs/developer/creating-a-protocol
# **************************************************************************
from .protocol_anm import ProDyANM
from .protocol_rtb import ProDyRTB
from .protocol_gnm import ProDyGNM
from .protocol_domdec import ProDyDomainDecomp 

from .protocol_atoms import ProDySelect, ProDyAlign, ProDyBiomol, ProDyAddPDBs
from .protocol_compare import ProDyCompare
from .protocol_edit import ProDyEdit
from .protocol_defvec import ProDyDefvec
from .protocol_import import ProDyImportModes

from .protocol_clustenm import ProDyClustENM

from .protocol_ensemble import ProDyBuildPDBEnsemble
from .protocol_import import ProDyImportEnsemble
from .protocol_pca import ProDyPCA
from .protocol_project import ProDyProject
from .protocol_rmsd import ProDyRmsd
from .protocol_distance import ProDyMeasure

try:
    from prody import addMissingAtoms
    from pdbfixer import PDBFixer
except ImportError:
    pass
else:
    from .protocol_pdbfixer import ProDyPDBFixer
