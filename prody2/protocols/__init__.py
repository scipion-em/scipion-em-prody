# -*- coding: utf-8 -*-
# **************************************************************************
# Module to declare protocols
# Find documentation here: https://scipion-em.github.io/docs/docs/developer/creating-a-protocol
# **************************************************************************
import prody
prody.confProDy(auto_secondary=True)

from .protocol_anm import ProDyANM
from .protocol_rtb import ProDyRTB
from .protocol_gnm import ProDyGNM

from .protocol_atoms import ProDySelect, ProDyAlign
from .protocol_compare import ProDyCompare
from .protocol_edit import ProDyEdit
from .protocol_defvec import ProDyDefvec
from .protocol_import import ProDyImportModes