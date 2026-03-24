from dataclasses import dataclass

from ndsl.dsl.typing import Float, Int


@dataclass
class GF2020Config:
    DT_MOIST: Float
    LHYDROSTATIC: bool
    STOCHASTIC_CNV: bool
    STOCH_TOP: Float
    STOCH_BOT: Float
    GF_MIN_AREA: Float
    GF_ENV_SETTING: Int
    ENTRVERSION: Int
    CONVECTION_TRACER: Int
    C1: Float
    ADV_TRIGGER: Int
    AUTOCONV: Int
    USE_TRACER_TRANSPORT: Int
    SCLM_DEEP: Float
    FIX_CONVECTIVE_CLOUD: bool
    APPLY_SUBSIDENCE_MICROPHYSICS: Int
    NUMBER_OF_TRACERS: Int
    USE_MOMENTUM_TRANSPORT: Int
