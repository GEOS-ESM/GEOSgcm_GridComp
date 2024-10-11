import gt4py.cartesian.gtscript as gtscript

from ndsl.dsl.typing import Float
from pyMoist.aer_activation_constants import n_modes
from pyMoist.saturation.constants import TABLESIZE


FloatField_NModes = gtscript.Field[gtscript.IJK, (Float, (n_modes))]
FloatField_VaporSaturationTable = gtscript.Field[gtscript.K, (Float, (int(TABLESIZE)))]
