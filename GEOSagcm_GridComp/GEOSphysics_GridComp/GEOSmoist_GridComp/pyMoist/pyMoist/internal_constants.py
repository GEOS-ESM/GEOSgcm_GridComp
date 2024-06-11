'''
Original Fortran file from 'src/Shared/@MAPL/shared/Constants/Constants.F90'
This is a port so NDSL has constants for math
'''

import numpy as np

MAPL_R8 = np.float64
MAPL_R4 = np.float32
MAPL_RN = float
MAPL_I8 = np.int64
MAPL_I4 = np.int32
MAPL_IN = int

MAPL_UNDEFINED_INTEGER = 1 - np.iinfo(np.int32).max
MAPL_UNDEFINED_REAL = np.finfo(np.float32).max