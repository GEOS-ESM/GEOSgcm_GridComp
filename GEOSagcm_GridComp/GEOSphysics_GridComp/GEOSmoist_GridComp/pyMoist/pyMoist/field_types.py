from ndsl.dsl.gt4py import IJ, IJK, Field, GlobalTable
from ndsl.dsl.typing import Bool, Float

from pyMoist.constants import NUMBER_OF_TRACERS, N_MODES

FloatField_NModes = Field[IJK, (Float, (N_MODES))]
FloatField_NTracers = Field[IJK, (Float, (int(NUMBER_OF_TRACERS)))]
FloatFieldIJ_NTracers = Field[IJ, (Float, (int(NUMBER_OF_TRACERS)))]
ConvectionTracerMetaDataTable_Float = GlobalTable[(Float, int(NUMBER_OF_TRACERS))]
ConvectionTracerMetaDataTable_Bool = GlobalTable[(Bool, int(NUMBER_OF_TRACERS))]
ConvectionTracerMetaDataTable_x3 = GlobalTable[(Float, (int(NUMBER_OF_TRACERS), 3))]
ConvectionTracerMetaDataTable_x4 = GlobalTable[(Float, (int(NUMBER_OF_TRACERS), 4))]
