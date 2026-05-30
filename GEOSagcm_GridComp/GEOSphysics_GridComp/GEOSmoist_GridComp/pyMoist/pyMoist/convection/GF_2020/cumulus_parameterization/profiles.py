from ndsl.dsl.gt4py import FORWARD, PARALLEL, computation, interval
from ndsl.dsl.typing import FloatField, FloatFieldIJ, Int

import pyMoist.constants as constants
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import FloatField_Plume, IntFieldIJ_Plume


def melting_profile(
    error_code: IntFieldIJ_Plume,
    plume: Int,
    melting_layer: FloatField,
    partition_liquid_ice: FloatField,
    p_cloud_levels_forced: FloatField_Plume,
    condensate_to_fall_forced: FloatField_Plume,
    melting: FloatField,
):
    """Generate the melting profile, given the excess water in the updraft.

    Args:
        error_code (IntFieldIJ_Plume)
        plume (Int)
        melting_layer (FloatField)
        partition_liquid_ice (FloatField)
        p_cloud_levels_forced (FloatField_Plume)
        condensate_to_fall_forced (FloatField_Plume)
        melting (FloatField)
    """
    from __externals__ import k_end

    with computation(FORWARD), interval(...):
        if cumulus_parameterization_constants.MELT_GLAC and plume == cumulus_parameterization_constants.DEEP:
            solid_phase_precipitable_water = 0.0
            effective_precipitable_water = 0.0
            melting = 0.0

            total_solid_phase_precipitable_water: FloatFieldIJ = 0.0

    with computation(FORWARD), interval(0, -2):
        if cumulus_parameterization_constants.MELT_GLAC and plume == cumulus_parameterization_constants.DEEP and error_code[0, 0][plume] == 0:
            dp = 100.0 * (p_cloud_levels_forced[0, 0, 0][plume] - p_cloud_levels_forced[0, 0, 1][plume])

            effective_precipitable_water = 0.5 * (condensate_to_fall_forced[0, 0, 0][plume] + condensate_to_fall_forced[0, 0, 1][plume])

            solid_phase_precipitable_water = (1.0 - partition_liquid_ice) * effective_precipitable_water

            total_solid_phase_precipitable_water = total_solid_phase_precipitable_water + solid_phase_precipitable_water * dp / constants.MAPL_GRAV

    with computation(PARALLEL), interval(0, -1):
        if cumulus_parameterization_constants.MELT_GLAC and plume == cumulus_parameterization_constants.DEEP and error_code[0, 0][plume] == 0:
            melting = melting_layer * (
                total_solid_phase_precipitable_water
                / (100.0 * (p_cloud_levels_forced.at(K=0, ddim=[plume]) - p_cloud_levels_forced.at(K=k_end - 1, ddim=[plume])) / constants.MAPL_GRAV)
            )

    with computation(PARALLEL), interval(...):
        if not (cumulus_parameterization_constants.MELT_GLAC and plume == cumulus_parameterization_constants.DEEP):
            melting = 0.0


class C1DProfile:
    """
    C1D Profile generator. This code is manually disables with a
    "do not enable" note in fortran, so it is not going to be completed until specifically requested.
    """

    def __init__(
        self,
        config: GF2020Config,
        cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
    ):
        if cumulus_parameterization_constants.FIRST_GUESS_W or config.AUTOCONV == 4:
            raise NotImplementedError(
                "[NDSL] GF2020-->CumulusParameterization-->C1DProfile: C1DProfile option has not been"
                "implemented. You should have been caught before getting here by the config checker."
                "Beware, something likely failing in the config checker as well - you may be unknowingly"
                "calling other untested/unimplemented sections."
            )

    def __call__(self, *args, **kwds):
        pass
