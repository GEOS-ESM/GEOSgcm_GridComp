from ndsl.dsl.gt4py import IJK, IJ, Field, GlobalTable
from ndsl.dsl.typing import Float
from pyMoist.constants import N_MODES, NCNST


FloatField_NModes = Field[IJK, (Float, (N_MODES))]
FloatField_NTracers = Field[IJK, (Float, (int(NCNST)))]
FloatFieldIJ_NTracers = Field[IJ, (Float, (int(NCNST)))]
TracerMetaDataTable = GlobalTable[(Float, int(NCNST))]
TracerMetaDataTable_extra_dim_3 = GlobalTable[(Float, (int(NCNST), 3))]
TracerMetaDataTable_extra_dim_4 = GlobalTable[(Float, (int(NCNST), 4))]
