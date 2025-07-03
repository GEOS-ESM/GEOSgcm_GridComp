from dataclasses import dataclass

from ndsl.dsl.typing import Float, Int


@dataclass
class PhaseChangeConfiguration:
    DT_MOIST: Float
    PDF_SHAPE: Int
    CCW_EVAP_EFF: Float
    CCI_EVAP_EFF: Float
    TURNRHCRIT_PARAM: Float
    DW_LAND: Float
    DW_OCEAN: Float
    DO_QA: bool
    DO_MELT_FREEZE: bool
    USE_BERGERON: bool
