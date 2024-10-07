from gt4py.cartesian import gtscript
from ndsl.dsl.typing import Float

import pyMoist.aer_activation_constants as constants

# Global space
FloatField_NModes = gtscript.Field[gtscript.IJK, (Float, (constants.n_modes))]
