from ndsl.dsl.gt4py import IJ, IJK, Field, GlobalTable
from ndsl.dsl.typing import Bool, Float

from pyMoist.constants import NUMBER_OF_TRACERS

FloatField_NTracers = Field[IJK, (Float, (int(NUMBER_OF_TRACERS)))]
FloatFieldIJ_NTracers = Field[IJ, (Float, (int(NUMBER_OF_TRACERS)))]
