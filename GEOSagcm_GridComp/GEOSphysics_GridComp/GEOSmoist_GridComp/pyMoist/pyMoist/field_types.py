import gt4py.cartesian.gtscript as gtscript
from ndsl.dsl.gt4py import IJK, Field
from ndsl.dsl.typing import Float

from pyMoist.constants import N_MODES, NCNST


FloatField_NModes = Field[IJK, (Float, (N_MODES))]
FloatField_NTracers = gtscript.Field[gtscript.IJK, (Float, (int(NCNST)))]
FloatFieldIJ_NTracers = gtscript.Field[gtscript.IJ, (Float, (int(NCNST)))]
