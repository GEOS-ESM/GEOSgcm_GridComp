from ndsl.dsl.gt4py import IJK, Field, GlobalTable, K
from ndsl.dsl.typing import Float
from pyMoist.constants import N_MODES
from pyMoist.saturation_tables.constants import TABLESIZE


FloatField_NModes = Field[IJK, (Float, (N_MODES))]
FloatField_VaporSaturationTable = Field[K, (Float, (int(TABLESIZE)))]
GlobalTable_saturaion_tables = GlobalTable[(Float, (int(TABLESIZE)))]
