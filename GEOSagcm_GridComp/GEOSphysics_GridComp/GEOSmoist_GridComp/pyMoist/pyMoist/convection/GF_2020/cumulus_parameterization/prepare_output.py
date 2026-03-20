import pyMoist.constants as constants
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
from ndsl import Local, NDSLRuntime, Quantity, QuantityFactory, StencilFactory
from ndsl.constants import I_DIM, J_DIM, K_DIM
from ndsl.dsl.gt4py import FORWARD, PARALLEL, K, computation, interval
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ, Int, IntFieldIJ
from ndsl.stencils.column_operations import column_max
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import (
    FloatField_Plume,
    FloatFieldIJ_Ensemble,
    FloatFieldIJ_Plume,
    IntFieldIJ_Plume,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.shared_functions import liquid_fraction
from pyMoist.shared_incloud_processes import (
    G_RATIO,
    RADIATIVE_EFFECTIVE_RADIUS,
    G_RATIO_Table_Type,
    RADIATIVE_EFFECTIVE_RADIUS_Table_Type,
    make_droplet_number,
    make_ice_number,
)


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
    """Compute output tendencies for microphysical properties (ice/liquid/vapor concentrations), wind (u/v),
    buoyancy, and temperature.

    Behavior of this stencil is heavily dependent on the value of the CLOSURE_CHOICE external.

    Args:
        error_code (IntFieldIJ_Plume)
        error_code_2 (IntFieldIJ)
        error_code_3 (IntFieldIJ)
        cloud_top_level (IntFieldIJ_Plume)
        updraft_lfc_level (IntFieldIJ_Plume)
        p_cloud_levels_forced (FloatField_Plume)
        normalized_massflux_updraft_forced (FloatField_Plume)
        precip (FloatFieldIJ_Plume)
        effective_condensate_to_fall_forced (FloatField)
        cloud_base_mass_flux_modified (FloatFieldIJ_Plume)
        scale_dependence_factor (FloatFieldIJ_Plume)
        ocean_fraction (FloatFieldIJ)
        f_dicycle_modified (FloatFieldIJ)
        del_u_cloud_ensemble (FloatField)
        del_v_cloud_ensemble (FloatField)
        del_t_cloud_ensemble (FloatField)
        del_vapor_cloud_ensemble (FloatField)
        del_cloud_liquid_cloud_ensemble (FloatField)
        del_buoyancy_cloud_ensemble (FloatField)
        del_convective_ice_cloud_ensemble (FloatField)
        del_large_scale_ice_cloud_ensemble (FloatField)
        del_convective_liquid_cloud_ensemble (FloatField)
        del_large_scale_liquid_cloud_ensemble (FloatField)
        del_convective_cloud_fraction_cloud_ensemble (FloatField)
        del_large_scale_cloud_fraction_cloud_ensemble (FloatField)
        dtdt (FloatField_Plume)
        dvapordt (FloatField_Plume)
        dcloudicedt (FloatField_Plume)
        dudt (FloatField_Plume)
        dvdt (FloatField_Plume)
        dbuoyancydt (FloatField_Plume)
        dconvectiveicedt (FloatField_Plume)
        dlargescaleicedt (FloatField_Plume)
        dconvectiveliquiddt (FloatField_Plume)
        dlargescaleliquiddt (FloatField_Plume)
        dconvectivecloudfractiondt (FloatField_Plume)
        dlargescalecloudfractiondt (FloatField_Plume)
        mass_flux_ensemble (FloatFieldIJ_Ensemble)
        precipitation_ensemble (FloatFieldIJ_Ensemble)
        xff_mid (FloatFieldIJ_Ensemble)
        CLOSURE_CHOICE (Int)
        CLOUD_BASE_MASS_FLUX_FACTOR (Float)
        plume (Int)
    """
    from __externals__ import (
        APPLY_SUBSIDENCE_MICROPHYSICS,
        DTIME,
        MAX_TEMP_VAPOR_TENDENCY,
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
            fixouts = (
                cloud_base_mass_flux_modified[0, 0][plume]
                * 86400.0
                * max(
                    max_val_t,
                    (cumulus_parameterization_constants.XLV / cumulus_parameterization_constants.CP)
                    * max_val_vapor,
                )
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
    """Compute evaporation of rain below cloud_top_level.

    Args:
        error_code (IntFieldIJ_Plume)
        cloud_top_level (IntFieldIJ_Plume)
        p_cloud_levels_forced (FloatField_Plume)
        evaporation_flux (FloatField)
        evaporation_sublimation_tendency (FloatField)
        plume (Int)
    """
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
    """Compute precipitation specifically from deep convection.

    Args:
        error_code (IntFieldIJ_Plume)
        cloud_top_level (IntFieldIJ_Plume)
        precipitation_flux (FloatField)
        convective_precip_flux (FloatField)
        plume (Int)
    """
    with computation(PARALLEL), interval(...):
        if (
            error_code[0, 0][plume] == 0
            and plume == cumulus_parameterization_constants.DEEP
            and K <= cloud_top_level[0, 0][plume] + 1
        ):
            convective_precip_flux = precipitation_flux


def output_updraft_temperature(
    error_code: IntFieldIJ_Plume,
    updraft_column_temperature_forced: FloatField,
    t_cloud_levels: FloatField,
    t_updraft: FloatField_Plume,
    plume: Int,
):
    """Copy the internal updraft temperature to an output array so that it can be
    fed back to the rest of the model.

    Args:
        error_code (IntFieldIJ_Plume)
        updraft_column_temperature_forced (FloatField)
        t_cloud_levels (FloatField)
        t_updraft (FloatField_Plume)
        plume (Int)
    """
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
    """Scale output mass fluxes based on the cloud base mass flux before output.

    Args:
        error_code (IntFieldIJ_Plume)
        cloud_base_mass_flux_modified (FloatFieldIJ_Plume)
        total_normalized_integrated_condensate_forced (FloatFieldIJ_Plume)
        total_normalized_integrated_evaporate_forced (FloatFieldIJ)
        normalized_massflux_updraft_forced (FloatField_Plume)
        normalized_massflux_downdraft_forced (FloatField_Plume)
        condensate_to_fall_forced (FloatField_Plume)
        evaporate_in_downdraft_forced (FloatField_Plume)
        mass_entrainment_updraft_forced (FloatField_Plume)
        mass_detrainment_updraft_forced (FloatField_Plume)
        mass_entrainment_downdraft_forced (FloatField_Plume)
        mass_detrainment_downdraft_forced (FloatField_Plume)
        environment_massflux (FloatField)
        vapor_tendency_from_environmental_subsidence (FloatField)
        moist_static_energy_tendency_from_environmental_subsidence (FloatField)
        t_tendency_from_environmental_subsidence (FloatField)
        plume (Int)
    """
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


def output_workfunctions_and_precip_concentrations(
    error_code: IntFieldIJ_Plume,
    cloud_top_level: IntFieldIJ_Plume,
    convection_fraction: FloatFieldIJ,
    surface_type: FloatFieldIJ,
    cloud_workfunction_0_output: FloatFieldIJ,
    cloud_workfunction_1_output: FloatFieldIJ,
    cloud_workfunction_0: FloatFieldIJ,
    cloud_workfunction_1: FloatFieldIJ,
    air_density: FloatField,
    updraft_column_temperature_forced: FloatField,
    dcloudicedt: FloatField_Plume,
    dnliquiddt: FloatField_Plume,
    dnicedt: FloatField_Plume,
    G_RATIO: G_RATIO_Table_Type,
    RADIATIVE_EFFECTIVE_RADIUS: RADIATIVE_EFFECTIVE_RADIUS_Table_Type,
    plume: Int,
):
    """Push the internal copy of workfunctions 0 and 1 to the state fields for output, and get precipitate
    concentrations using the functions make_droplet_number and make_ice_number

    Args:
        error_code (IntFieldIJ_Plume)
        cloud_top_level (IntFieldIJ_Plume)
        convection_fraction (FloatFieldIJ)
        surface_type (FloatFieldIJ)
        cloud_workfunction_0_output (FloatFieldIJ)
        cloud_workfunction_1_output (FloatFieldIJ)
        cloud_workfunction_0 (FloatFieldIJ)
        cloud_workfunction_1 (FloatFieldIJ)
        air_density (FloatField)
        updraft_column_temperature_forced (FloatField)
        dcloudicedt (FloatField_Plume)
        dnliquiddt (FloatField_Plume)
        dnicedt (FloatField_Plume)
        G_RATIO (G_RATIO_Table_Type)
        RADIATIVE_EFFECTIVE_RADIUS (RADIATIVE_EFFECTIVE_RADIUS_Table_Type)
        plume (Int)
    """
    from __externals__ import DTIME, FRAC_MODIS

    with computation(PARALLEL), interval(...):
        n_water_friendly_aerosols = 99.0e7  # in the future set this as NCPL

    with computation(FORWARD), interval(0, 1):
        if plume == cumulus_parameterization_constants.DEEP:
            cloud_workfunction_0_output = 0.0
            cloud_workfunction_1_output = 0.0

            if error_code[0, 0][plume] == 0:
                cloud_workfunction_0_output = cloud_workfunction_0
                cloud_workfunction_1_output = cloud_workfunction_1

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0 and K <= cloud_top_level[0, 0][plume] + 1:
            fraction = liquid_fraction(
                updraft_column_temperature_forced, convection_fraction, surface_type, FRAC_MODIS
            )
            cloud_liquid = DTIME * dcloudicedt[0, 0, 0][plume] * air_density * fraction
            cloud_ice = DTIME * dcloudicedt[0, 0, 0][plume] * air_density * (1.0 - fraction)

            dnicedt[0, 0, 0][plume] = max(
                0.0,
                make_ice_number(cloud_ice, updraft_column_temperature_forced, RADIATIVE_EFFECTIVE_RADIUS)
                / air_density,
            )
            dnliquiddt[0, 0, 0][plume] = max(
                0.0, make_droplet_number(cloud_liquid, n_water_friendly_aerosols, G_RATIO) / air_density
            )

            # convert to tendencies
            dnicedt[0, 0, 0][plume] = dnicedt[0, 0, 0][plume] * (1 / DTIME)  # unit [1/s]
            dnliquiddt[0, 0, 0][plume] = dnliquiddt[0, 0, 0][plume] * (1 / DTIME)  # unit [1/s]


class OutputWorkfunctionsAndPrecipConcentrations(NDSLRuntime):
    """Push the internal copy of workfunctions 0 and 1 to the state fields for output, and get precipitate
    concentrations using the functions make_droplet_number and make_ice_number
    """

    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
    ):
        # init NDSLRuntime
        super().__init__(stencil_factory)

        # add dimension to quantityfactory and create classes for constants
        quantity_factory.add_data_dimensions({"G_RATIO_Table": len(G_RATIO)})
        quantity_factory.add_data_dimensions(
            {"RADIATIVE_EFFECTIVE_RADIUS_Table": len(RADIATIVE_EFFECTIVE_RADIUS)}
        )

        self._G_RATIO: Local = quantity_factory.zeros(["G_RATIO_Table"], "n/a")
        self._RADIATIVE_EFFECTIVE_RADIUS: Local = quantity_factory.zeros(
            ["RADIATIVE_EFFECTIVE_RADIUS_Table"], "n/a"
        )

        self._G_RATIO.field[:] = G_RATIO
        self._RADIATIVE_EFFECTIVE_RADIUS.field[:] = RADIATIVE_EFFECTIVE_RADIUS

        # construct stencils
        self._output_workfunctions_and_precip_concentrations = stencil_factory.from_dims_halo(
            func=output_workfunctions_and_precip_concentrations,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "FRAC_MODIS": cumulus_parameterization_config.FRAC_MODIS,
                "DTIME": cumulus_parameterization_config.DTIME,
            },
        )

    def __call__(
        self,
        error_code: Quantity,
        cloud_top_level: Quantity,
        convection_fraction: Quantity,
        surface_type: Quantity,
        cloud_workfunction_0_output: Quantity,
        cloud_workfunction_1_output: Quantity,
        cloud_workfunction_0: Quantity,
        cloud_workfunction_1: Quantity,
        air_density: Quantity,
        updraft_column_temperature_forced: Quantity,
        dcloudicedt: Quantity,
        dnliquiddt: Quantity,
        dnicedt: Quantity,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._output_workfunctions_and_precip_concentrations(
            error_code=error_code,
            cloud_top_level=cloud_top_level,
            convection_fraction=convection_fraction,
            surface_type=surface_type,
            cloud_workfunction_0_output=cloud_workfunction_0_output,
            cloud_workfunction_1_output=cloud_workfunction_1_output,
            cloud_workfunction_0=cloud_workfunction_0,
            cloud_workfunction_1=cloud_workfunction_1,
            air_density=air_density,
            updraft_column_temperature_forced=updraft_column_temperature_forced,
            dcloudicedt=dcloudicedt,
            dnliquiddt=dnliquiddt,
            dnicedt=dnicedt,
            G_RATIO=self._G_RATIO,
            RADIATIVE_EFFECTIVE_RADIUS=self._RADIATIVE_EFFECTIVE_RADIUS,
            plume=plume_dependent_constants.PLUME_INDEX,
        )


class LightningFlashDensity:
    def __init__(self, cumulus_parameterization_config: GF2020CumulusParameterizationConfig):
        if cumulus_parameterization_config.LIGHTNING_DIAGNOSTICS:
            raise NotImplementedError(
                "[NDSL] GF2020-->CumulusParameterization-->LightningFlashDensity: lightning output has not"
                "been implemented. You should have been caught before getting here by the config checker."
                "Beware, something likely failing in the config checker as well - you may be unknowingly"
                "calling other untested/unimplemented sections."
            )

    def __call__(self, *args, **kwds):
        pass
