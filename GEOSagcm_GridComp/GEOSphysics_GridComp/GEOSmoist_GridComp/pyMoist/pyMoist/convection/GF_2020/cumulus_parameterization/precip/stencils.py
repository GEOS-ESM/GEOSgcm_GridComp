from ndsl.dsl.gt4py import PARALLEL, computation, interval, FORWARD, function
from ndsl.dsl.typing import FloatField, FloatFieldIJ, Float, IntField, Int
import pyMoist.constants as constants
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import FloatField_Plume, IntFieldIJ_Plume
from pyMoist.shared_incloud_processes import ice_fraction


@function
def liquid_fraction(
    t,
    convection_fraction,
    surface_type,
    MODIS_FRACTION,
):
    """
    Get the fraction of liquid condensates

    Args:
        t (in): temperature
        convection_fraction (in)
        surface_type (in)
        MODIS_FRACTION (in): use fraction liq/ice content derived from MODIS/CALIPO sensors
    """
    if MODIS_FRACTION == 1:
        liquid_fraction = 1.0 - ice_fraction(t, convection_fraction, surface_type)
    else:
        liquid_fraction = min(
            1.0,
            (
                max(0.0, (t - cumulus_parameterization_constants.T_ICE))
                / (cumulus_parameterization_constants.T_0 - cumulus_parameterization_constants.T_ICE)
            )
            ** 2,
        )

    return liquid_fraction


def partition_liquid_ice(
    t: FloatField,
    p: FloatField_Plume,
    geopotential_height: FloatField,
    topography_height_no_negative: FloatFieldIJ,
    surface_type: FloatFieldIJ,
    convection_fraction: FloatFieldIJ,
    error_code: IntFieldIJ_Plume,
    melting_layer: FloatField,
    part_liquid_ice: FloatField,
    plume: Int,
):
    """
    Partition total condensate into liquid and ice phases

    Args:

    """
    from __externals__ import MELT_ICE, MODIS_FRACTION, k_end

    with computation(PARALLEL), interval(...):
        # constants, set internally because they may differ from global constants
        # and need to only exist inside this stencil
        T1 = 276.16
        Z_meltlayer1 = 4000.0
        Z_meltlayer2 = 6000.0
        delT = 3.0

        # prefill some fields
        part_liquid_ice = 1.0
        melting_layer = 0.0

    with computation(PARALLEL), interval(0, -1):
        if MELT_ICE == 1 and plume == 2:
            if error_code[0, 0][plume] == 0:
                # get function of T for partition of total condensate into liq and ice phases
                part_liquid_ice = liquid_fraction(t, convection_fraction, surface_type, MODIS_FRACTION)

    with computation(PARALLEL), interval(0, -1):
        if MELT_ICE == 1 and plume == 2:
            if error_code[0, 0][plume] == 0:
                # define the melting layer (the layer will be between T_0+1 < TEMP < T_1
                if t <= (cumulus_parameterization_constants.T_0 - delT):
                    melting_layer = 0.0

                elif t < (cumulus_parameterization_constants.T_0 + delT) and t > (
                    cumulus_parameterization_constants.T_0 - delT
                ):
                    melting_layer = (
                        (t - (cumulus_parameterization_constants.T_0 - delT)) / (2.0 * delT)
                    ) ** 2

                else:
                    melting_layer = 1.0

                melting_layer = melting_layer * (1.0 - melting_layer)

    with computation(FORWARD), interval(0, 1):
        if MELT_ICE == 1 and plume == 2:
            # normalize vertical integral of melting_layer to 1
            norm: FloatFieldIJ = 0.0

    with computation(FORWARD), interval(...):
        if MELT_ICE == 1 and plume == 2:
            if error_code[0, 0][plume] == 0:
                # normalize vertical integral of melting_layer to 1
                dp = 100.0 * (p[0, 0, 0][plume] - p[0, 0, 1][plume])
                norm = norm + melting_layer * dp / constants.MAPL_GRAV

    with computation(PARALLEL), interval(...):
        if MELT_ICE == 1 and plume == 2:
            if error_code[0, 0][plume] == 0:
                # normalize vertical integral of melting_layer to 1
                melting_layer = (
                    melting_layer
                    / (norm + 1.0e-6)
                    * (100 * (p[0, 0, 0][plume] - p.at(K=k_end - 1, ddim=[plume])) / constants.MAPL_GRAV)
                )
