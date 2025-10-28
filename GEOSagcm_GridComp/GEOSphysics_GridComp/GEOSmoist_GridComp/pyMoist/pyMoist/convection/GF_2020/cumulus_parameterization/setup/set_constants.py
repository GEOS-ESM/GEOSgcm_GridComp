from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)
from ndsl.dsl.typing import Float, Int


def set_constants(
    cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
    plume_dependent_constants: GF2020PlumeDependentConstants,
    plume: str,
):
    if plume == "shallow":
        # set a number of plume dependent constants
        plume_dependent_constants.PLUME_INDEX = Int(0)
        plume_dependent_constants.DOWNDRAFT_MAX_HEIGHT_LAND = (
            cumulus_parameterization_config.DOWNDRAFT_MAX_HEIGHT_LAND_SHALLOW
        )
        plume_dependent_constants.DOWNDRAFT_MAX_HEIGHT_OCEAN = (
            cumulus_parameterization_config.DOWNDRAFT_MAX_HEIGHT_OCEAN_SHALLOW
        )
        plume_dependent_constants.UPDRAFT_MAX_HEIGHT_LAND = (
            cumulus_parameterization_config.UPDRAFT_MAX_HEIGHT_LAND_SHALLOW
        )
        plume_dependent_constants.UPDRAFT_MAX_HEIGHT_OCEAN = (
            cumulus_parameterization_config.UPDRAFT_MAX_HEIGHT_OCEAN_SHALLOW
        )
        plume_dependent_constants.MINIMUM_EVAP_FRACTION_LAND = (
            cumulus_parameterization_config.MINIMUM_EVAP_FRACTION_LAND_SHALLOW
        )
        plume_dependent_constants.MINIMUM_EVAP_FRACTION_OCEAN = (
            cumulus_parameterization_config.MINIMUM_EVAP_FRACTION_OCEAN_SHALLOW
        )
        plume_dependent_constants.MAXIMUM_EVAP_FRACTION_LAND = (
            cumulus_parameterization_config.MAXIMUM_EVAP_FRACTION_LAND_SHALLOW
        )
        plume_dependent_constants.MAXIMUM_EVAP_FRACTION_OCEAN = (
            cumulus_parameterization_config.MAXIMUM_EVAP_FRACTION_OCEAN_SHALLOW
        )
        plume_dependent_constants.CLOUD_BASE_MASS_FLUX_FACTOR = (
            cumulus_parameterization_config.CLOUD_BASE_MASS_FLUX_FACTOR_SHALLOW
        )
        plume_dependent_constants.USE_EXCESS = cumulus_parameterization_config.USE_EXCESS_SHALLOW
        plume_dependent_constants.ENTRAINMENT_RATE = cumulus_parameterization_config.ENTRAINMENT_RATE_SHALLOW
        plume_dependent_constants.ENABLE_PLUME = cumulus_parameterization_config.ENABLE_SHALLOW

        # maximum depth (mb) of capping inversion (larger cap = no convection)
        if (
            cumulus_parameterization_config.ZERO_DIFF == 1
            or cumulus_parameterization_config.MOIST_TRIGGER == 0
        ):
            plume_dependent_constants.CAP_MAX_INC = Float(25.0)
        else:
            plume_dependent_constants.CAP_MAX_INC = Float(10.0)

        # lambda_U parameter for momentum transport
        if cumulus_parameterization_config.PRESSURE_GRADIENT_CONSTANT != 0.0:
            plume_dependent_constants.LAMBDA_DEEP = Float(0.0)
            plume_dependent_constants.LAMBDA_DOWN = Float(0.0)
        else:
            plume_dependent_constants.LAMBDA_DEEP = cumulus_parameterization_config.LAMBDA_DEEP
            plume_dependent_constants.LAMBDA_DOWN = cumulus_parameterization_config.LAMBDA_SHALLOW_DOWN

        # minimum depth (m) clouds must have
        plume_dependent_constants.DEPTH_MIN = Float(500.0)

        # max height(m) above ground where updraft air can originate
        plume_dependent_constants.MAX_UPDRAFT_ORIGIN_HEIGHT = Float(2000.0)

        # height(m) above which no downdrafts are allowed to originate
        plume_dependent_constants.MAX_DOWNDRAFT_ORIGIN_HEIGHt = Float(3000.0)

    elif plume == "mid":
        # set a number of plume dependent constants
        plume_dependent_constants.PLUME_INDEX = Int(1)
        plume_dependent_constants.DOWNDRAFT_MAX_HEIGHT_LAND = (
            cumulus_parameterization_config.DOWNDRAFT_MAX_HEIGHT_LAND_MID
        )
        plume_dependent_constants.DOWNDRAFT_MAX_HEIGHT_OCEAN = (
            cumulus_parameterization_config.DOWNDRAFT_MAX_HEIGHT_OCEAN_MID
        )
        plume_dependent_constants.UPDRAFT_MAX_HEIGHT_LAND = (
            cumulus_parameterization_config.UPDRAFT_MAX_HEIGHT_LAND_MID
        )
        plume_dependent_constants.UPDRAFT_MAX_HEIGHT_OCEAN = (
            cumulus_parameterization_config.UPDRAFT_MAX_HEIGHT_OCEAN_MID
        )
        plume_dependent_constants.MINIMUM_EVAP_FRACTION_LAND = (
            cumulus_parameterization_config.MINIMUM_EVAP_FRACTION_LAND_MID
        )
        plume_dependent_constants.MINIMUM_EVAP_FRACTION_OCEAN = (
            cumulus_parameterization_config.MINIMUM_EVAP_FRACTION_OCEAN_MID
        )
        plume_dependent_constants.MAXIMUM_EVAP_FRACTION_LAND = (
            cumulus_parameterization_config.MAXIMUM_EVAP_FRACTION_LAND_MID
        )
        plume_dependent_constants.MAXIMUM_EVAP_FRACTION_OCEAN = (
            cumulus_parameterization_config.MAXIMUM_EVAP_FRACTION_OCEAN_MID
        )
        plume_dependent_constants.CLOUD_BASE_MASS_FLUX_FACTOR = (
            cumulus_parameterization_config.CLOUD_BASE_MASS_FLUX_FACTOR_MID
        )
        plume_dependent_constants.USE_EXCESS = cumulus_parameterization_config.USE_EXCESS_MID
        plume_dependent_constants.ENTRAINMENT_RATE = cumulus_parameterization_config.ENTRAINMENT_RATE_MID
        plume_dependent_constants.ENABLE_PLUME = cumulus_parameterization_config.ENABLE_MID

        # maximum depth (mb) of capping inversion (larger cap = no convection)
        if (
            cumulus_parameterization_config.ZERO_DIFF == 1
            or cumulus_parameterization_config.MOIST_TRIGGER == 0
        ):
            plume_dependent_constants.CAP_MAX_INC = Float(10.0)
        else:
            plume_dependent_constants.CAP_MAX_INC = Float(90.0)

        # lambda_U parameter for momentum transport
        if cumulus_parameterization_config.PRESSURE_GRADIENT_CONSTANT != 0.0:
            plume_dependent_constants.LAMBDA_DEEP = Float(0.0)
            plume_dependent_constants.LAMBDA_DOWN = Float(0.0)
        else:
            plume_dependent_constants.LAMBDA_DEEP = cumulus_parameterization_config.LAMBDA_SHALLOW_DOWN
            plume_dependent_constants.LAMBDA_DOWN = cumulus_parameterization_config.LAMBDA_SHALLOW_DOWN

        # minimum depth (m) clouds must have
        plume_dependent_constants.DEPTH_MIN = Float(1000.0)

        # max height(m) above ground where updraft air can originate
        plume_dependent_constants.MAX_UPDRAFT_ORIGIN_HEIGHT = Float(3000.0)

        # height(m) above which no downdrafts are allowed to originate
        plume_dependent_constants.MAX_DOWNDRAFT_ORIGIN_HEIGHt = Float(3000.0)

    elif plume == "deep":
        # set a number of plume dependent constants
        plume_dependent_constants.PLUME_INDEX = Int(2)
        plume_dependent_constants.DOWNDRAFT_MAX_HEIGHT_LAND = (
            cumulus_parameterization_config.DOWNDRAFT_MAX_HEIGHT_LAND_DEEP
        )
        plume_dependent_constants.DOWNDRAFT_MAX_HEIGHT_OCEAN = (
            cumulus_parameterization_config.DOWNDRAFT_MAX_HEIGHT_OCEAN_DEEP
        )
        plume_dependent_constants.UPDRAFT_MAX_HEIGHT_LAND = (
            cumulus_parameterization_config.UPDRAFT_MAX_HEIGHT_LAND_DEEP
        )
        plume_dependent_constants.UPDRAFT_MAX_HEIGHT_OCEAN = (
            cumulus_parameterization_config.UPDRAFT_MAX_HEIGHT_OCEAN_DEEP
        )
        plume_dependent_constants.MINIMUM_EVAP_FRACTION_LAND = (
            cumulus_parameterization_config.MINIMUM_EVAP_FRACTION_LAND_DEEP
        )
        plume_dependent_constants.MINIMUM_EVAP_FRACTION_OCEAN = (
            cumulus_parameterization_config.MINIMUM_EVAP_FRACTION_OCEAN_DEEP
        )
        plume_dependent_constants.MAXIMUM_EVAP_FRACTION_LAND = (
            cumulus_parameterization_config.MAXIMUM_EVAP_FRACTION_LAND_DEEP
        )
        plume_dependent_constants.MAXIMUM_EVAP_FRACTION_OCEAN = (
            cumulus_parameterization_config.MAXIMUM_EVAP_FRACTION_OCEAN_DEEP
        )
        plume_dependent_constants.CLOUD_BASE_MASS_FLUX_FACTOR = (
            cumulus_parameterization_config.CLOUD_BASE_MASS_FLUX_FACTOR_DEEP
        )
        plume_dependent_constants.USE_EXCESS = cumulus_parameterization_config.USE_EXCESS_DEEP
        plume_dependent_constants.ENTRAINMENT_RATE = cumulus_parameterization_config.ENTRAINMENT_RATE_DEEP
        plume_dependent_constants.ENABLE_PLUME = cumulus_parameterization_config.ENABLE_DEEP

        # maximum depth (mb) of capping inversion (larger cap = no convection)
        if (
            cumulus_parameterization_config.ZERO_DIFF == 1
            or cumulus_parameterization_config.MOIST_TRIGGER == 0
        ):
            plume_dependent_constants.CAP_MAX_INC = Float(20.0)
        else:
            plume_dependent_constants.CAP_MAX_INC = Float(90.0)

        # lambda_U parameter for momentum transport
        if cumulus_parameterization_config.PRESSURE_GRADIENT_CONSTANT != 0.0:
            plume_dependent_constants.LAMBDA_DEEP = Float(0.0)
            plume_dependent_constants.LAMBDA_DOWN = Float(0.0)
        else:
            plume_dependent_constants.LAMBDA_DEEP = cumulus_parameterization_config.LAMBDA_SHALLOW_DOWN
            plume_dependent_constants.LAMBDA_DOWN = cumulus_parameterization_config.LAMBDA_SHALLOW_DOWN

        # minimum depth (m) clouds must have
        plume_dependent_constants.DEPTH_MIN = Float(1000.0)

        # max height(m) above ground where updraft air can originate
        plume_dependent_constants.MAX_UPDRAFT_ORIGIN_HEIGHT = Float(4000.0)

        # height(m) above which no downdrafts are allowed to originate
        plume_dependent_constants.MAX_DOWNDRAFT_ORIGIN_HEIGHt = Float(3000.0)

    else:
        raise NotImplementedError("Unknown plume specified, corresponding constants are unavailable.")

    return plume_dependent_constants
