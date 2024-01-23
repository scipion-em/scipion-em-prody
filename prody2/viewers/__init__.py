# -*- coding: utf-8 -*-
# **************************************************************************
# Module to declare viewers
# Find documentation here: https://scipion-em.github.io/docs/docs/developer/creating-a-viewer
# **************************************************************************
from .viewer_modes import ProDyModeViewer
from .viewer_compare import ProDyComparisonsViewer
from .viewer_gnm import ProDyGNMViewer
from .viewer_domdec import ProDyDomainViewer
from .viewer_project import ProDyProjectionsViewer
from .viewer_ensemble import ProDyEnsembleViewer
from .viewer_lda import ProDyLDAViewer

from pwem.viewers import BasicMDViewer, showj
from pwem.viewers.viewers_data import RegistryViewerConfig
from prody2.objects import SetOfTrajFrames, SetOfAtoms
BasicMDViewer._targets.extend([SetOfTrajFrames, SetOfAtoms])

labels = labels = 'id enabled label _filename '
RegistryViewerConfig.registerConfig(SetOfTrajFrames,
                                    {showj.ORDER: labels,
                                     showj.VISIBLE: labels,
                                     showj.MODE: showj.MODE_MD,
                                     showj.RENDER: 'no'})
