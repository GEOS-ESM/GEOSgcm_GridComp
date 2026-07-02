import dataclasses

from ndsl.dsl.typing import Float, Int


@dataclasses.dataclass
class GF2020PlumeDependentConstants:
    PLUME_INDEX: Int = Int(0)
    DOWNDRAFT_MAX_HEIGHT_LAND: Float = Float(0.0)
    DOWNDRAFT_MAX_HEIGHT_OCEAN: Float = Float(0.0)
    UPDRAFT_MAX_HEIGHT_LAND: Float = Float(0.0)
    UPDRAFT_MAX_HEIGHT_OCEAN: Float = Float(0.0)
    MINIMUM_EVAP_FRACTION_LAND: Float = Float(0.0)
    MINIMUM_EVAP_FRACTION_OCEAN: Float = Float(0.0)
    MAXIMUM_EVAP_FRACTION_LAND: Float = Float(0.0)
    MAXIMUM_EVAP_FRACTION_OCEAN: Float = Float(0.0)
    CLOUD_BASE_MASS_FLUX_FACTOR: Float = Float(0.0)
    USE_EXCESS: Int = Int(0)
    ENTRAINMENT_RATE: Float = Float(0.0)
    ENABLE_PLUME: Int = Int(0)
    AVERAGE_LAYER_DEPTH = Float(0.0)
    CAP_MAX_INC: Float = Float(0.0)
    LAMBDA_DEEP: Float = Float(0.0)
    LAMBDA_DOWN: Float = Float(0.0)
    MINIMUM_DEPTH: Float = Float(0.0)
    MAX_UPDRAFT_ORIGIN_HEIGHT: Float = Float(0.0)
    MAX_DOWNDRAFT_ORIGIN_HEIGHt: Float = Float(0.0)
    DETRAINMENT_CRITICAL_DEPTH: Float = Float(0.0)
    C0: Float = Float(0.0)
    T_STAR: Float = Float(0.0)
    TAU_CAPE_REMOVAL: Float = Float(0.0)
    CLOSURE_CHOICE: Int = Int(0.0)
