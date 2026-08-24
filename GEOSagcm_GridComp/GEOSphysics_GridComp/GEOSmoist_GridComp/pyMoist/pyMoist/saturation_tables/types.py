from ndsl.dsl.gt4py import GlobalTable
from ndsl.dsl.typing import Float

from pyMoist.saturation_tables.constants import TABLESIZE


GlobalTable_saturation_tables = GlobalTable[(Float, (int(TABLESIZE)))]
