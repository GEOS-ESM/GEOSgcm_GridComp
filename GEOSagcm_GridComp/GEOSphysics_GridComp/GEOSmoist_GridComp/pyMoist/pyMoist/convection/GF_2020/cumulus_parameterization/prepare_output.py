from ndsl.dsl.gt4py import (
    PARALLEL,
    computation,
    interval,
    FORWARD,
    function,
    BACKWARD,
    K,
)
from ndsl.dsl.typing import FloatField, FloatFieldIJ, Float, IntField, Int, IntFieldIJ
import pyMoist.constants as constants
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import (
    FloatField_Plume,
    IntFieldIJ_Plume,
    FloatFieldIJ_Plume,
)
from pyMoist.shared_incloud_processes import ice_fraction
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig


def total_evaporation_flux(
    error_code: IntFieldIJ_Plume,
    cloud_top_level: IntFieldIJ_Plume,
    p_cloud_levels_forced: FloatField_Plume,
    evaporation_flux: FloatField,
    evaporation_sublimation_tendency: FloatField,
    plume: Int,
):
    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0:
            if K <= cloud_top_level[0, 0][plume]:
                dp = 100.0 * (p_cloud_levels_forced[0, 0, 0][plume] - p_cloud_levels_forced[0, 0, 1][plume])
                evaporation_sublimation_tendency = (
                    evaporation_sublimation_tendency + evaporation_flux * constants.MAPL_GRAV / dp
                )


def deep_precipitation_output(
    error_code: IntFieldIJ_Plume,
    cloud_top_level: IntFieldIJ_Plume,
    precipitation_flux: FloatField,
    convective_precip_flux: FloatField,
    plume: Int,
):
    with computation(PARALLEL), interval(...):
        if (
            error_code[0, 0][plume] == 0
            and plume == cumulus_parameterization_constants.DEEP
            and K <= cloud_top_level[0, 0][plume] + 1
        ):
            convective_precip_flux = precipitation_flux


def tracer_output(
    error_code: IntFieldIJ_Plume,
    updraft_column_temperature_forced: FloatField,
    t_cloud_levels: FloatField,
    t_updraft: FloatField_Plume,
    plume: Int,
):
    with computation(PARALLEL), interval(0, -1):
        if error_code[0, 0][plume] == 0:
            t_updraft[0, 0, 0][plume] = updraft_column_temperature_forced
    with computation(PARALLEL), interval(-1, None):
        if error_code[0, 0][plume] == 0:
            t_updraft[0, 0, 0][plume] = t_cloud_levels


def prepare_output(
    error_code: IntFieldIJ_Plume,
    cloud_base_mass_flux_modified: FloatFieldIJ_Plume,
    total_normalized_integrated_condensate_forced: FloatFieldIJ_Plume,
    total_normalized_integrated_evaporate_forced: FloatFieldIJ,
    normalized_massflux_updraft_forced: FloatField_Plume,
    normalized_massflux_downdraft_forced: FloatField_Plume,
    condensate_to_fall_forced: FloatField_Plume,
    evaporate_in_downdraft_forced: FloatField_Plume,
    mass_entrainment_updraft_forced: FloatField_Plume,
    mass_detrainment_updraft_forced: FloatField_Plume,
    mass_entrainment_downdraft_forced: FloatField_Plume,
    mass_detrainment_downdraft_forced: FloatField_Plume,
    environment_massflux: FloatField,
    vapor_tendency_from_environmental_subsidence: FloatField,
    moist_static_energy_tendency_from_environmental_subsidence: FloatField,
    t_tendency_from_environmental_subsidence: FloatField,
    plume: Int,
):
    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            total_normalized_integrated_condensate_forced[0, 0][plume] = (
                cloud_base_mass_flux_modified[0, 0][plume]
                * total_normalized_integrated_condensate_forced[0, 0][plume]
            )
            total_normalized_integrated_evaporate_forced = (
                cloud_base_mass_flux_modified[0, 0][plume]
                * total_normalized_integrated_evaporate_forced
            )
    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0:
            normalized_massflux_updraft_forced[0, 0, 0][plume] = (
                cloud_base_mass_flux_modified[0, 0][plume]
                * normalized_massflux_updraft_forced[0, 0, 0][plume]
            )
            normalized_massflux_downdraft_forced[0, 0, 0][plume] = (
                cloud_base_mass_flux_modified[0, 0][plume]
                * normalized_massflux_downdraft_forced[0, 0, 0][plume]
            )
            condensate_to_fall_forced[0, 0, 0][plume] = (
                cloud_base_mass_flux_modified[0, 0][plume] * condensate_to_fall_forced[0, 0, 0][plume]
            )
            evaporate_in_downdraft_forced[0, 0, 0][plume] = (
                cloud_base_mass_flux_modified[0, 0][plume] * evaporate_in_downdraft_forced[0, 0, 0][plume]
            )
            mass_entrainment_updraft_forced[0, 0, 0][plume] = (
                cloud_base_mass_flux_modified[0, 0][plume] * mass_entrainment_updraft_forced[0, 0, 0][plume]
            )
            mass_detrainment_updraft_forced[0, 0, 0][plume] = (
                cloud_base_mass_flux_modified[0, 0][plume] * mass_detrainment_updraft_forced[0, 0, 0][plume]
            )
            mass_entrainment_downdraft_forced[0, 0, 0][plume] = (
                cloud_base_mass_flux_modified[0, 0][plume] * mass_entrainment_downdraft_forced[0, 0, 0][plume]
            )
            mass_detrainment_downdraft_forced[0, 0, 0][plume] = (
                cloud_base_mass_flux_modified[0, 0][plume] * mass_detrainment_downdraft_forced[0, 0, 0][plume]
            )
            environment_massflux = cloud_base_mass_flux_modified[0, 0][plume] * environment_massflux

    with computation(PARALLEL), interval(...):
        vapor_tendency_from_environmental_subsidence = (
            cloud_base_mass_flux_modified[0, 0][plume] * vapor_tendency_from_environmental_subsidence
        )
        moist_static_energy_tendency_from_environmental_subsidence = (
            cloud_base_mass_flux_modified[0, 0][plume]
            * moist_static_energy_tendency_from_environmental_subsidence
        )
        t_tendency_from_environmental_subsidence = (
            cloud_base_mass_flux_modified[0, 0][plume] * t_tendency_from_environmental_subsidence
        )


class LightningFlashDensity:
    def __init__(self, cumulus_parameterization_config: GF2020CumulusParameterizationConfig):
        if cumulus_parameterization_config.LIGHTNING_DIAGNOSTICS:
            raise NotImplementedError(
                "GF2020 lightning output has not been implemented. You should have"
                "been caught before this error - something is wrong with the config checker!"
            )

    def __call__(self, *args, **kwds):
        pass
