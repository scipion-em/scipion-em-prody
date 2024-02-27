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

from pwem.viewers import DataViewer, showj
from pwem.viewers.viewers_data import RegistryViewerConfig
from prody2.objects import (SetOfTrajFrames, SetOfAtoms,
                            SetOfClassesTraj)
DataViewer._targets.extend([SetOfTrajFrames, SetOfAtoms,
                            SetOfClassesTraj])

from pwem.viewers.views import Classes3DView

RegistryViewerConfig.registerConfig(SetOfTrajFrames,
                                    {showj.ORDER: 'id enabled label _filename ',
                                     showj.VISIBLE: 'id enabled label _filename ',
                                     showj.MODE: showj.MODE_MD,
                                     showj.RENDER: 'no'})
RegistryViewerConfig.registerConfig(SetOfClassesTraj,
                                    {showj.ORDER: 'id enabled label _size ',
                                     showj.VISIBLE: 'id enabled label _size ',
                                     showj.MODE: showj.MODE_MD,
                                     showj.RENDER: 'no'})
