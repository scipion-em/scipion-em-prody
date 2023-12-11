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

from pwem.viewers import DataViewer
from prody2.objects import SetOfTrajFrames
DataViewer._targets.append(SetOfTrajFrames)
