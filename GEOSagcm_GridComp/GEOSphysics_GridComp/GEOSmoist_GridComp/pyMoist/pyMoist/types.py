import gt4py.cartesian.gtscript as gtscript

import pyMoist.aer_activation_constants as constants
import pyMoist.pyMoist_constants as moist_constants
from ndsl.dsl.typing import Float


# Global space
FloatField_NModes = gtscript.Field[gtscript.IJK, (Float, (constants.n_modes))]
FloatField_NTracers = gtscript.Field[gtscript.IJK, (Float, (moist_constants.ncnst))]
