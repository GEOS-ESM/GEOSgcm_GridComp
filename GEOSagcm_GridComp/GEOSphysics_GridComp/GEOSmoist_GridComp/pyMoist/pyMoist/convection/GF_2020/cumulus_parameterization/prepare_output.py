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
    FloatFieldIJ_Ensemble,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from ndsl.stencils.column_operations import column_max


def ensemble_output_and_feedback(
    error_code: IntFieldIJ_Plume,
    error_code_2: IntFieldIJ,
    error_code_3: IntFieldIJ,
    cloud_top_level: IntFieldIJ_Plume,
    updraft_lfc_level: IntFieldIJ_Plume,
    p_cloud_levels_forced: FloatField_Plume,
    normalized_massflux_updraft_forced: FloatField_Plume,
    precip: FloatFieldIJ_Plume,
    effective_condensate_to_fall_forced: FloatField,
    cloud_base_mass_flux_modified: FloatFieldIJ_Plume,
    scale_dependence_factor: FloatFieldIJ_Plume,
    ocean_fraction: FloatFieldIJ,
    f_dicycle_modified: FloatFieldIJ,
    del_u_cloud_ensemble: FloatField,
    del_v_cloud_ensemble: FloatField,
    del_t_cloud_ensemble: FloatField,
    del_vapor_cloud_ensemble: FloatField,
    del_cloud_liquid_cloud_ensemble: FloatField,
    del_buoyancy_cloud_ensemble: FloatField,
    del_convective_ice_cloud_ensemble: FloatField,
    del_large_scale_ice_cloud_ensemble: FloatField,
    del_convective_liquid_cloud_ensemble: FloatField,
    del_large_scale_liquid_cloud_ensemble: FloatField,
    del_convective_cloud_fraction_cloud_ensemble: FloatField,
    del_large_scale_cloud_fraction_cloud_ensemble: FloatField,
    dtdt: FloatField_Plume,
    dvapordt: FloatField_Plume,
    dcloudicedt: FloatField_Plume,
    dudt: FloatField_Plume,
    dvdt: FloatField_Plume,
    dbuoyancydt: FloatField_Plume,
    dconvectiveicedt: FloatField_Plume,
    dlargescaleicedt: FloatField_Plume,
    dconvectiveliquiddt: FloatField_Plume,
    dlargescaleliquiddt: FloatField_Plume,
    dconvectivecloudfractiondt: FloatField_Plume,
    dlargescalecloudfractiondt: FloatField_Plume,
    mass_flux_ensemble: FloatFieldIJ_Ensemble,
    precipitation_ensemble: FloatFieldIJ_Ensemble,
    xff_mid: FloatFieldIJ_Ensemble,
    CLOSURE_CHOICE: Int,
    CLOUD_BASE_MASS_FLUX_FACTOR: Float,
    plume: Int,
):
    from __externals__ import (
        DTIME,
        MAX_TEMP_VAPOR_TENDENCY,
        APPLY_SUBSIDENCE_MICROPHYSICS,
        USE_SMOOTH_TENDENCIES,
        k_end,
    )

    # ensure outputs are all zero
    with computation(PARALLEL), interval(...):
        dtdt[0, 0, 0][plume] = 0.0
        dvapordt[0, 0, 0][plume] = 0.0
        dcloudicedt[0, 0, 0][plume] = 0.0
        dudt[0, 0, 0][plume] = 0.0
        dvdt[0, 0, 0][plume] = 0.0
        dbuoyancydt[0, 0, 0][plume] = 0.0

    with computation(FORWARD), interval(0, 1):
        precip[0, 0][plume] = 0.0
        cloud_base_mass_flux_modified[0, 0][plume] = 0.0

        if error_code[0, 0][plume] == 0:
            member = 0
            while member < cumulus_parameterization_constants.MAXENS3:
                if precipitation_ensemble[0, 0][member] <= 0.0:
                    mass_flux_ensemble[0, 0][member] = 0.0
                member += 1

    # calculate ensemble average mass fluxes
    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0.0:
            average_mass_flux = 0.0
            if plume == cumulus_parameterization_constants.DEEP:
                count = 0
                member = 0
                while member < cumulus_parameterization_constants.MAXENS3:
                    count = count + 1
                    average_mass_flux = average_mass_flux + mass_flux_ensemble[0, 0][member]
                    member += 1
                # ensemble average mass flux
                average_mass_flux = average_mass_flux / count

            elif plume == cumulus_parameterization_constants.MID:
                if CLOSURE_CHOICE <= 5:
                    if CLOSURE_CHOICE == 0:
                        average_mass_flux = 0.3333 * (xff_mid[0, 0][0] + xff_mid[0, 0][1] + xff_mid[0, 0][2])
                    else:
                        average_mass_flux = xff_mid[0, 0][CLOSURE_CHOICE - 1]

            elif plume == cumulus_parameterization_constants.SHALLOW:
                option_not_implemented = True

    # set the updradt mass flux and do not allow negative values and apply the diurnal cycle closure
    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            # mass flux of updradt at cloud base
            cloud_base_mass_flux_modified[0, 0][plume] = average_mass_flux

            # apply the adjust factor for tunning
            cloud_base_mass_flux_modified[0, 0][plume] = (
                CLOUD_BASE_MASS_FLUX_FACTOR * cloud_base_mass_flux_modified[0, 0][plume]
            )

            # diurnal cycle closure
            cloud_base_mass_flux_modified[0, 0][plume] = (
                cloud_base_mass_flux_modified[0, 0][plume] - f_dicycle_modified
            )

            if cloud_base_mass_flux_modified[0, 0][plume] <= 0.0:
                error_code[0, 0][plume] = 13
                cloud_base_mass_flux_modified[0, 0][plume] = 0.0

    # apply the scale-dependence Arakawa's approach
    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            # scale dependence
            cloud_base_mass_flux_modified[0, 0][plume] = (
                scale_dependence_factor[0, 0][plume] * cloud_base_mass_flux_modified[0, 0][plume]
            )

            if cloud_base_mass_flux_modified[0, 0][plume] == 0.0:
                error_code[0, 0][plume] = 14
            if cloud_base_mass_flux_modified[0, 0][plume] > 100.0:
                error_code[0, 0][plume] = 15

    # sanity check for mass flux
    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            max_mass_flux = (
                100.0
                * (
                    p_cloud_levels_forced.at(K=updraft_lfc_level[0, 0][plume], ddim=[plume])
                    - p_cloud_levels_forced.at(K=updraft_lfc_level[0, 0][plume] + 1, ddim=[plume])
                )
                / (constants.MAPL_GRAV * DTIME)
            )
            cloud_base_mass_flux_modified[0, 0][plume] = min(
                cloud_base_mass_flux_modified[0, 0][plume], max_mass_flux
            )

    with computation(PARALLEL), interval(...):
        # prepare field from which max value is pulled
        abs_del_t_cloud_ensemble = abs(del_t_cloud_ensemble)
        abs_del_vapor_cloud_ensemble = abs(del_vapor_cloud_ensemble)

    # check dtdt and and dvapordt for high values
    # criteria: if abs (dtdt or dvapordt) > MAX_TEMP_VAPOR_TENDENCY (default 100 K/day) => fix mass flux
    with computation(FORWARD), interval(0, 1):
        if MAX_TEMP_VAPOR_TENDENCY > 0.0 and error_code[0, 0][plume] == 0:
            max_val_t, _ = column_max(abs_del_t_cloud_ensemble, 0, cloud_top_level[0, 0][plume])
            max_val_vapor, _ = column_max(abs_del_vapor_cloud_ensemble, 0, cloud_top_level[0, 0][plume])
            fixouts = cloud_base_mass_flux_modified[0, 0][plume] * 86400.0 * max(
                max_val_t,
                (cumulus_parameterization_constants.XLV / cumulus_parameterization_constants.CP)
                * max_val_vapor,
            )

            if fixouts > MAX_TEMP_VAPOR_TENDENCY:  # K/day
                fixouts = MAX_TEMP_VAPOR_TENDENCY / fixouts
                cloud_base_mass_flux_modified[0, 0][plume] = (
                    cloud_base_mass_flux_modified[0, 0][plume] * fixouts
                )
                member = 0
                while member < cumulus_parameterization_constants.MAXENS3:
                    mass_flux_ensemble[0, 0][member] = mass_flux_ensemble[0, 0][member] * fixouts
                    member += 1

    # now do feedback
    with computation(FORWARD), interval(...):
        if error_code[0, 0][plume] == 0:
            if K <= cloud_top_level[0, 0][plume]:
                precip[0, 0][plume] = (
                    precip[0, 0][plume]
                    + effective_condensate_to_fall_forced * cloud_base_mass_flux_modified[0, 0][plume]
                )

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0:
            if K <= cloud_top_level[0, 0][plume]:
                dtdt[0, 0, 0][plume] = del_t_cloud_ensemble * cloud_base_mass_flux_modified[0, 0][plume]
                dvapordt[0, 0, 0][plume] = (
                    del_vapor_cloud_ensemble * cloud_base_mass_flux_modified[0, 0][plume]
                )
                dcloudicedt[0, 0, 0][plume] = (
                    del_cloud_liquid_cloud_ensemble * cloud_base_mass_flux_modified[0, 0][plume]
                )
                dudt[0, 0, 0][plume] = del_u_cloud_ensemble * cloud_base_mass_flux_modified[0, 0][plume]
                dvdt[0, 0, 0][plume] = del_v_cloud_ensemble * cloud_base_mass_flux_modified[0, 0][plume]
                dbuoyancydt[0, 0, 0][plume] = (
                    del_buoyancy_cloud_ensemble * cloud_base_mass_flux_modified[0, 0][plume]
                )

            if APPLY_SUBSIDENCE_MICROPHYSICS == 1:
                if K <= cloud_top_level[0, 0][plume]:
                    dconvectiveicedt[0, 0, 0][plume] = (
                        del_convective_ice_cloud_ensemble * cloud_base_mass_flux_modified[0, 0][plume]
                    )
                    dlargescaleicedt[0, 0, 0][plume] = (
                        del_large_scale_ice_cloud_ensemble * cloud_base_mass_flux_modified[0, 0][plume]
                    )
                    dconvectiveliquiddt[0, 0, 0][plume] = (
                        del_convective_liquid_cloud_ensemble * cloud_base_mass_flux_modified[0, 0][plume]
                    )
                    dlargescaleliquiddt[0, 0, 0][plume] = (
                        del_large_scale_liquid_cloud_ensemble * cloud_base_mass_flux_modified[0, 0][plume]
                    )
                    dconvectivecloudfractiondt[0, 0, 0][plume] = (
                        del_convective_cloud_fraction_cloud_ensemble
                        * cloud_base_mass_flux_modified[0, 0][plume]
                    )
                    dlargescalecloudfractiondt[0, 0, 0][plume] = (
                        del_large_scale_cloud_fraction_cloud_ensemble
                        * cloud_base_mass_flux_modified[0, 0][plume]
                    )

                if K >= cloud_top_level[0, 0][plume] and K < k_end:
                    dconvectiveicedt[0, 0, 0][plume] = 0.0
                    dlargescaleicedt[0, 0, 0][plume] = 0.0
                    dconvectiveliquiddt[0, 0, 0][plume] = 0.0
                    dlargescaleliquiddt[0, 0, 0][plume] = 0.0
                    dconvectivecloudfractiondt[0, 0, 0][plume] = 0.0
                    dlargescalecloudfractiondt[0, 0, 0][plume] = 0.0

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            member = 0
            while (
                member
                < cumulus_parameterization_constants.MAXENS1
                * cumulus_parameterization_constants.MAXENS2
                * cumulus_parameterization_constants.MAXENS3
            ):
                mass_flux_ensemble[0, 0][member] = (
                    scale_dependence_factor[0, 0][plume] * mass_flux_ensemble[0, 0][member]
                )
                member += 1

    # smooth the tendencies (future work: include outbuoy, outmpc* and tracers)
    with computation(PARALLEL), interval(...):
        if USE_SMOOTH_TENDENCIES < 0:
            option_not_implemented = 0.0


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
                cloud_base_mass_flux_modified[0, 0][plume] * total_normalized_integrated_evaporate_forced
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
