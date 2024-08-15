from ndsl.dsl.typing import Float
from pyMoist.aer_activation_constants import MAPL_TICE
import numpy as np

# ToDo:  Exist in RadiationCoupling, merge
MAPL_TICE = 273.16  # K

TMINTBL = Float(150.0)
TMAXTBL = Float(333.0)
TMINLQU = MAPL_TICE - Float(40.0)
DEGSUBS = 100
DELTA_T = Float(1.0) / DEGSUBS
TABLESIZE = np.int32(TMAXTBL - TMINTBL) * DEGSUBS + 1
TMIX = Float(-20.0)

# ToDo:  Exist in RadiationCoupling, merge
MAPL_AIRMW = Float(28.965)
MAPL_H2OMW = Float(18.015)  # kg/Kmole

ESFAC = MAPL_H2OMW / MAPL_AIRMW
ERFAC = DEGSUBS / ESFAC

MAX_MIXING_RATIO = Float(1.0)
