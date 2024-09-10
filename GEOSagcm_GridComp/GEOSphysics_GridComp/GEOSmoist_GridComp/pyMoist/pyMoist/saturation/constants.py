from ndsl.dsl.typing import Float
import numpy as np

# Generic constants to be moved to pyMoist global constant file
MAPL_TICE = Float(273.16)  # K
MAPL_AIRMW = Float(28.965)
MAPL_H2OMW = Float(18.015)  # kg/Kmole

# Saturation specific constants
TMINTBL = Float(150.0)
TMAXTBL = Float(333.0)
TMINLQU = MAPL_TICE - Float(40.0)
DEGSUBS = np.int32(100)
DELTA_T = Float(1.0) / DEGSUBS
TABLESIZE = np.int32(TMAXTBL - TMINTBL) * DEGSUBS + 1
TMIX = Float(-20.0)
ESFAC = MAPL_H2OMW / MAPL_AIRMW
ERFAC = DEGSUBS / ESFAC

MAX_MIXING_RATIO = Float(1.0)
