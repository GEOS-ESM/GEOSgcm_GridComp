import dataclasses

from ndsl.dsl.typing import Float, Int


@dataclasses.dataclass
class GF2020PlumeDependentConstants:
    DOWNDRAFT_MAX_HEIGHT_LAND: Float = 0.0
    DOWNDRAFT_MAX_HEIGHT_OCEAN: Float = 0.0
    UPDRAFT_MAX_HEIGHT_LAND: Float = 0.0
    UPDRAFT_MAX_HEIGHT_OCEAN: Float = 0.0
    MINIMUM_EVAP_FRACTION_LAND: Float = 0.0
    MINIMUM_EVAP_FRACTION_OCEAN: Float = 0.0
    MAXIMUM_EVAP_FRACTION_LAND: Float = 0.0
    MAXIMUM_EVAP_FRACTION_OCEAN: Float = 0.0
    CLOUD_BASE_MASS_FLUX_FACTOR: Float = 0.0
    USE_EXCESS: Int = 0
    ENABLE_PLUME: Float = 0.0
    CAP_MAX_INC: Float = 0.0
    LAMBDA_DEEP: Float = 0.0
    LAMBDA_DOWN: Float = 0.0
