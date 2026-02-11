from ndsl.dsl.gt4py import IJ, IJK, Field
from ndsl.dsl.typing import Float
from pyMoist.constants import N_MODES, NCNST


FloatField_NModes = Field[IJK, (Float, (N_MODES))]
FloatField_NTracers = Field[IJK, (Float, (int(NCNST)))]
FloatFieldIJ_NTracers = Field[IJ, (Float, (int(NCNST)))]
