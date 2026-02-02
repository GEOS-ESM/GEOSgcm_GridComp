from ndsl.dsl.gt4py import IJK, IJ, Field, GlobalTable
from ndsl.dsl.typing import Float, Int, Bool
from pyMoist.constants import N_MODES, NCNST


FloatField_NModes = Field[IJK, (Float, (N_MODES))]
FloatField_NTracers = Field[IJK, (Float, (int(NCNST)))]
FloatFieldIJ_NTracers = Field[IJ, (Float, (int(NCNST)))]
ConvectionTracerMetaDataTable_Float = GlobalTable[(Float, int(NCNST))]
ConvectionTracerMetaDataTable_Bool = GlobalTable[(Bool, int(NCNST))]
ConvectionTracerMetaDataTable_x3 = GlobalTable[(Float, (int(NCNST), 3))]
ConvectionTracerMetaDataTable_x4 = GlobalTable[(Float, (int(NCNST), 4))]
