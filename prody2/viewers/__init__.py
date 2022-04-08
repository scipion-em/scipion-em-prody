# -*- coding: utf-8 -*-
# **************************************************************************
# Module to declare viewers
# Find documentation here: https://scipion-em.github.io/docs/docs/developer/creating-a-viewer
# **************************************************************************
import prody
prody.confProDy(auto_secondary=True)

from .viewer_modes import ProDyModeViewer
from .viewer_compare import ProDyComparisonsViewer