import gt4py.cartesian.gtscript as gtscript

from ndsl.dsl.typing import Float
from pyMoist.constants import N_MODES, NCNST
from pyMoist.saturation.constants import TABLESIZE


FloatField_NModes = gtscript.Field[gtscript.IJK, (Float, (N_MODES))]
FloatField_NTracers = gtscript.Field[gtscript.IJK, (Float, (int(NCNST)))]
FloatField_VaporSaturationTable = gtscript.Field[gtscript.K, (Float, (int(TABLESIZE)))]
