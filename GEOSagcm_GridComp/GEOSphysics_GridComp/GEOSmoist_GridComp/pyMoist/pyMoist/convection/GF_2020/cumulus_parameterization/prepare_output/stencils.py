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


def prepare_output(
    error_code: IntFieldIJ_Plume,
    plume: Int,
    cloud_base_mass_flux_modified: FloatFieldIJ_Plume,
    total_normalized_integrated_condensate_forced: FloatFieldIJ_Plume,
    total_normalized_integrated_evaporate_forced: FloatFieldIJ_Plume,
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
):
    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            total_normalized_integrated_condensate_forced[0, 0][plume] = (
                cloud_base_mass_flux_modified[0, 0][plume]
                * total_normalized_integrated_condensate_forced[0, 0][plume]
            )
            total_normalized_integrated_evaporate_forced[0, 0][plume] = (
                cloud_base_mass_flux_modified[0, 0][plume]
                * total_normalized_integrated_evaporate_forced[0, 0][plume]
            )
    with computation(FORWARD), interval(...):
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

    with computation(FORWARD), interval(...):
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
