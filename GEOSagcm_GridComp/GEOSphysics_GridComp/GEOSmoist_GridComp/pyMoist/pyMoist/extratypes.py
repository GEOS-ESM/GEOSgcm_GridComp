import gt4py.cartesian.gtscript as gtscript

from pyMoist.aer_activation_constants import n_modes
from pyMoist.saturation.constants import TABLESIZE
from ndsl.dsl.typing import Float


# Global space
FloatField_NModes = gtscript.Field[gtscript.IJK, (Float, (n_modes))]
FloatField_Extra_Dim = gtscript.Field[gtscript.K, (Float, (int(TABLESIZE)))]