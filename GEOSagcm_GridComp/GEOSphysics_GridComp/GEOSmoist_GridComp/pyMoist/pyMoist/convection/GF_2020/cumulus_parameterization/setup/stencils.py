from ndsl.dsl.typing import FloatField, FloatFieldIJ, Float, IntFieldIJ, Int
from ndsl.dsl.gt4py import computation, PARALLEL, interval, FORWARD
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as GF_2020_constants
import pyMoist.constants as constants
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import (
    FloatFieldIJ_Plumes,
    FloatField_Plumes,
    FloatField_Ensemble,
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
    vapor_new: FloatField,
    t_new_pbl: FloatField,
    vapor_new_pbl: FloatField,
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
        vapor_new = vapor_old + (subgrid_scale_forcing_vapor + grid_scale_forcing_vapor) * DT_MOIST
        vapor_new = max(GF_2020_constants.smaller_qv, vapor_new)

        # temp/water vapor modified only by bl processes
        t_new_pbl = t_old + (subgrid_scale_forcing_t) * DT_MOIST
        vapor_new_pbl = vapor_old + (subgrid_scale_forcing_vapor) * DT_MOIST

        # moist static energy
        moist_static_energy = GF_2020_constants.CP * (
            subgrid_scale_forcing_t + grid_scale_forcing_t
        ) + GF_2020_constants.XLV * (subgrid_scale_forcing_vapor + grid_scale_forcing_vapor)


def prefil_internal_fields(
    kbmax: IntFieldIJ,
    kstamb: IntFieldIJ,
    ocean_fraction: FloatFieldIJ,
    ocean_fraction_local: FloatFieldIJ,
    cap_max: FloatFieldIJ,
    ierr2: IntFieldIJ,
    ierr3: IntFieldIJ,
    ierrc: IntFieldIJ,
    CAP_MAX_INC: Float,
    max_increment: FloatFieldIJ,
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
    plume: Int,
    t_wetbulb: FloatFieldIJ,
    vapor_wetbulb: FloatFieldIJ,
    tau_ecmwf: FloatFieldIJ,
    f_dicycle_modified: FloatFieldIJ,
    add_buoy_modified: FloatFieldIJ,
    scale_dependence_factor_downdraft: FloatFieldIJ,
    hcdo: FloatField,
    cupclw: FloatField,
    qrcdo: FloatField,
    hcot: FloatField,
    c1d: FloatField,
    evap_bcb: FloatField,
    mass_flux_ensemble: FloatField_Ensemble,
    precipitation_ensemble: FloatField_Ensemble,
    epsilon: FloatField_Plumes,
    precip: FloatField_Plumes,
    scale_dependence_factor: FloatField_Plumes,
    lightning_density: FloatFieldIJ,
):
    from __external__ import k_end, CAP_MAXS, ENSEMBLE_MEMEBRS

    # reset to zero manually. cannot use locals.fill(0) because not all fields are reset to zero b/t plumes
    # internal fields
    with computation(FORWARD), interval(0, 1):
        kbmax = 1
        kstabm = k_end - 2
        ocean_fraction_local = ocean_fraction
        cap_max = CAP_MAXS
        ierr2 = 0
        ierr3 = 0
        ierrc = -999
        max_increment = CAP_MAX_INC
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
        add_buoy_modified = 0.0
        scale_dependence_factor_downdraft = 0.0

    # internal fields
    with computation(PARALLEL), interval(...):
        geopotential_height_local = geopotential_height
        geopotential_height_modified_local = geopotential_height
        hcdo = 0.0
        cupclw = 0.0
        qrcdo = 0.0
        hcot = 0.0
        c1d = 0.0
        evap_bcb = 0.0

        member = 0
        while member < ENSEMBLE_MEMEBRS:
            mass_flux_ensemble[0, 0, 0][member] = 0.0
            precipitation_ensemble[0, 0, 0][member] = 0.0
            member = member + 1

    # external fields (from outside cumulus parameterizaton)
    with computation(FORWARD), interval(0, 1):
        epsilon[0, 0, 0][plume] = 0.0
        precip[0, 0, 0][plume] = 0.0
        scale_dependence_factor[0, 0, 0][plume] = 0.0
        lightning_density = 0.0


def compute_scale_dependence_factor(
    plume: Int,
    scale_dependence_factor: FloatFieldIJ,
    seed_convection: FloatFieldIJ,
    ierr: FloatFieldIJ_Plumes,
    ierrc: FloatFieldIJ,
    grid_length: FloatFieldIJ,
):
    from __externals__ import USE_SCALE_DEP

    # prepare mask to stop loop in next computation
    with computation(FORWARD), interval(0, 1):
        error_at_point = False

    with computation(FORWARD), interval(0, 1):
        if USE_SCALE_DEP == 0:
            scale_dependence_factor = 1.0
        elif USE_SCALE_DEP == 1:
            if plume == 1:
                scale_dependence_factor = 0.0
            else:
                if seed_convection < 0.0:
                    ierr[0, 0, 0][plume] = 1
                    ierrc = 1
                    error_at_point = True
                if error_at_point == False:
                    scale_dependence_factor = sigma(grid_length)
                if seed_convection != 1.0:
                    scale_dependence_factor = scale_dependence_factor ** (
                        seed_convection * max(1.0, scale_dependence_factor)
                    )
                scale_dependence_factor = max(0.1, min(scale_dependence_factor, 1.0))
                if scale_dependence_factor <= 0.1:
                    ierr[0, 0, 0][plume] = 1
                    ierrc = 1
                    error_at_point = True


def get_random_number(
    random_number: FloatFieldIJ,
):
    # need to figure out how to get system clock data
    pass
