from ndsl.dsl.typing import FloatField, FloatFieldIJ, Float, IntFieldIJ, Int
from ndsl.dsl.gt4py import computation, PARALLEL, interval, FORWARD
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
import pyMoist.constants as constants
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import (
    IntFieldIJ_Plume,
    FloatFieldIJ_Plume,
    FloatField_Plume,
    FloatFieldIJ_Ensemble,
)
from pyMoist.shared_generic_math import sigma


def set_plume_dependent_fields(
    t_excess: FloatFieldIJ,
    t_excess_local: FloatFieldIJ,
    vapor_excess: FloatFieldIJ,
    vapor_excess_local: FloatFieldIJ,
    ocean_fraction: FloatFieldIJ,
    use_excess: Int,
    t_old: FloatField,
    vapor_old: FloatField,
    grid_scale_forcing_t: FloatField,
    grid_scale_forcing_vapor: FloatField,
    subgrid_scale_forcing_t: FloatField,
    subgrid_scale_forcing_vapor: FloatField,
    t_new: FloatField,
    vapor_forced: FloatField,
    t_new_pbl: FloatField,
    vapor_forced_pbl: FloatField,
    moist_static_energy: FloatField,
):
    """
    fill or modify existing values of temperature and vapor excess before cumulus parameterization

    Args:
        t_excess_cu_param_input: temperature excess input to cumulus parameterization
        t_excess_cu_param_internal: temperature excess internal for cumulus parameterization
        vapor_excess_cu_param_input: water vapor mixing ratio excess input to cumulus parameterization
        vapor_excess_cu_param_internal: excess internal for cumulus parameterization

        use_excess: trigger to control behavior

    Possible behaviors (use_excess options):
        0: fill with zero
        1: no change (retain existing values)
        2: enforce min/max everywhere
        other: enforce min/max only over ocean
    """
    from __externals__ import DT_MOIST

    with computation(FORWARD), interval(0, 1):
        if use_excess == 0:
            t_excess_local = 0.0
            vapor_excess_local = 0.0
        elif use_excess == 2:
            t_excess_local = min(0.5, max(0.2, t_excess))  # Kelvin
            vapor_excess_local = min(5.0e-4, max(1.0e-4, vapor_excess))  # kg kg^-1
        else:
            if ocean_fraction > 0.98:  # ocean
                t_excess_local = min(0.5, max(0.2, t_excess))  # Kelvin
                vapor_excess_local = min(5.0e-4, max(1.0e-4, vapor_excess))  # kg kg^-1

    with computation(PARALLEL), interval(0, -1):
        t_new = t_old + (subgrid_scale_forcing_t + grid_scale_forcing_t) * DT_MOIST
        vapor_forced = vapor_old + (subgrid_scale_forcing_vapor + grid_scale_forcing_vapor) * DT_MOIST
        vapor_forced = max(cumulus_parameterization_constants.smaller_qv, vapor_forced)

        # temp/water vapor modified only by bl processes
        t_new_pbl = t_old + (subgrid_scale_forcing_t) * DT_MOIST
        vapor_forced_pbl = vapor_old + (subgrid_scale_forcing_vapor) * DT_MOIST

        # moist static energy
        moist_static_energy = cumulus_parameterization_constants.CP * (
            subgrid_scale_forcing_t + grid_scale_forcing_t
        ) + cumulus_parameterization_constants.XLV * (subgrid_scale_forcing_vapor + grid_scale_forcing_vapor)


def prefil_internal_fields(
    plume: Int,
    maximum_updraft_origin_level: IntFieldIJ,
    kstabm: IntFieldIJ,
    ocean_fraction: FloatFieldIJ,
    ocean_fraction_local: FloatFieldIJ,
    cap_max: FloatFieldIJ,
    error_code_2: IntFieldIJ,
    error_code_3: IntFieldIJ,
    CAP_MAX_INC: Float,
    cap_max_increment: FloatFieldIJ,
    geopotential_height: FloatField,
    geopotential_height_local: FloatField,
    geopotential_height_modified_local: FloatField,
    cloud_work_function_0: FloatFieldIJ,
    cloud_work_function_1: FloatFieldIJ,
    cloud_work_function_2: FloatFieldIJ,
    cloud_work_function_3: FloatFieldIJ,
    cloud_work_function_0_pbl: FloatFieldIJ,
    cloud_work_function_1_pbl: FloatFieldIJ,
    cloud_work_function_1_fa: FloatFieldIJ,
    cin1: FloatFieldIJ,
    k_x_modified: FloatFieldIJ,
    epsilon_local: FloatFieldIJ,
    pbl_time_scale: FloatFieldIJ,
    t_wetbulb: FloatFieldIJ,
    vapor_wetbulb: FloatFieldIJ,
    tau_ecmwf: FloatFieldIJ,
    f_dicycle_modified: FloatFieldIJ,
    add_buoyancy: FloatFieldIJ,
    scale_dependence_factor_downdraft: FloatFieldIJ,
    hcdo: FloatField,
    cupclw: FloatField,
    qrcdo: FloatField,
    cloud_moist_static_energy_forced_transported: FloatField,
    c1d: FloatField,
    evap_bcb: FloatField,
    mass_flux_ensemble: FloatFieldIJ_Ensemble,
    precipitation_ensemble: FloatFieldIJ_Ensemble,
    epsilon: FloatFieldIJ_Plume,
    precip: FloatFieldIJ_Plume,
    scale_dependence_factor: FloatFieldIJ_Plume,
    lightning_density: FloatFieldIJ,
):
    from __externals__ import k_end, CAP_MAXS, ENSEMBLE_MEMBERS

    # reset to zero manually. cannot use locals.fill(0) because not all fields are reset to zero b/t plumes
    # internal fields
    with computation(FORWARD), interval(0, 1):
        maximum_updraft_origin_level = 1
        kstabm = k_end - 2
        ocean_fraction_local = ocean_fraction
        cap_max = CAP_MAXS
        error_code_2 = 0
        error_code_3 = 0
        cap_max_increment = CAP_MAX_INC
        cloud_work_function_0 = 0.0
        cloud_work_function_1 = 0.0
        cloud_work_function_2 = 0.0
        cloud_work_function_3 = 0.0
        cloud_work_function_0_pbl = 0.0
        cloud_work_function_1_pbl = 0.0
        cloud_work_function_1_fa = 0.0
        cin1 = 0.0
        k_x_modified = 0.0
        epsilon_local = 0.0
        pbl_time_scale = 0.0
        t_wetbulb = 0.0
        vapor_wetbulb = 0.0
        tau_ecmwf = 0.0
        f_dicycle_modified = 0.0
        add_buoyancy = 0.0
        scale_dependence_factor_downdraft = 0.0

    # internal fields
    with computation(PARALLEL), interval(...):
        geopotential_height_local = geopotential_height
        geopotential_height_modified_local = geopotential_height
        hcdo = 0.0
        cupclw = 0.0
        qrcdo = 0.0
        cloud_moist_static_energy_forced_transported = 0.0
        c1d = 0.0
        evap_bcb = 0.0

    # internal fields
    # reset all ensemble members
    with computation(FORWARD), interval(0, 1):
        member = 0
        while member < ENSEMBLE_MEMBERS:
            mass_flux_ensemble[0, 0][member] = 0.0
            precipitation_ensemble[0, 0][member] = 0.0
            member = member + 1

    # external fields (from outside cumulus parameterizaton)
    with computation(FORWARD), interval(0, 1):
        epsilon[0, 0][plume] = 0.0
        precip[0, 0][plume] = 0.0
        scale_dependence_factor[0, 0][plume] = 0.0
        lightning_density = 0.0


def compute_scale_dependence_factor(
    plume: Int,
    scale_dependence_factor: FloatFieldIJ_Plume,
    seed_convection: FloatFieldIJ,
    error_code: IntFieldIJ_Plume,
    grid_length: FloatFieldIJ,
):
    from __externals__ import USE_SCALE_DEP

    # prepare mask to stop loop in next computation
    with computation(FORWARD), interval(0, 1):
        error_at_point = False

    with computation(FORWARD), interval(0, 1):
        if USE_SCALE_DEP == 0:
            scale_dependence_factor[0, 0][plume] = 1.0
        elif USE_SCALE_DEP == 1:
            if plume == 0:
                scale_dependence_factor[0, 0][plume] = 1.0
            else:
                if seed_convection < 0.0:
                    error_code[0, 0][plume] = 1
                    error_at_point = True
                if error_at_point == False:
                    scale_dependence_factor[0, 0][plume] = sigma(grid_length)
                if seed_convection != 1.0:
                    scale_dependence_factor[0, 0][plume] = scale_dependence_factor[0, 0][plume] ** (
                        seed_convection * max(1.0, scale_dependence_factor[0, 0][plume])
                    )
                scale_dependence_factor[0, 0][plume] = max(
                    0.1, min(scale_dependence_factor[0, 0][plume], 1.0)
                )
                if scale_dependence_factor[0, 0][plume] <= 0.1:
                    error_code[0, 0][plume] = 1
                    error_at_point = True


def get_random_number(
    plume: Int,
    random_number: FloatFieldIJ,
):
    from __externals__ import USE_RANDOM_NUMBER

    with computation(FORWARD), interval(0, 1):
        if plume == 2 and USE_RANDOM_NUMBER > 1.0e-6:
            # need to figure out how to get system clock data
            random_number = random_number  # keep input data from fortran for now
        else:
            random_number = 0.0


def initial_entrainment_detrainment(
    plume: Int,
    lateral_entrainment_rate: FloatField,
    current_plume_rate: Float,
    entrainment_rate: FloatField_Plume,
    detrainment_function_updraft: FloatField,
):
    with computation(PARALLEL), interval(0, -1):
        entrainment_rate[0, 0, 0][plume] = lateral_entrainment_rate * current_plume_rate
        detrainment_function_updraft = lateral_entrainment_rate * current_plume_rate


def epsilon_min_max(
    ocean_fraction: FloatFieldIJ,
    epsilon_min: FloatFieldIJ,
    epsilon_max: FloatFieldIJ,
    MINIMUM_EVAP_FRACTION_OCEAN: Float,
    MAXIMUM_EVAP_FRACTION_OCEAN: Float,
    MINIMUM_EVAP_FRACTION_LAND: Float,
    MAXIMUM_EVAP_FRACTION_LAND: Float,
):
    with computation(FORWARD), interval(0, 1):
        if ocean_fraction > 0.99:  # water
            epsilon_min = MINIMUM_EVAP_FRACTION_OCEAN
            epsilon_max = MAXIMUM_EVAP_FRACTION_OCEAN
        else:  # land
            epsilon_min = MINIMUM_EVAP_FRACTION_LAND
            epsilon_max = MAXIMUM_EVAP_FRACTION_LAND


def calculate_arbitrary_numerical_parameter(
    arbitrary_numerical_parameter: FloatFieldIJ,
):
    with computation(FORWARD), interval(0, 1):
        arbitrary_numerical_parameter = 0.1
        # approx xmb * timescale
        # other options (from fortran):
        # 0.1 * dtime*xmb_nm1(i)
        # 100.*(p_cup(i,kbcon(i))-p_cup(i,kbcon(i)+1))/(g*dtime)
        # 0.1*mbdt(i)
