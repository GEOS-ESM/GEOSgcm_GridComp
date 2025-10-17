from ndsl.dsl.typing import FloatField, FloatFieldIJ, Float, IntFieldIJ, Int
from ndsl.dsl.gt4py import computation, PARALLEL, interval, FORWARD
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as GF_2020_constants
import pyMoist.constants as constants


def set_plume_dependent_fields(
    t_excess_input: FloatFieldIJ,
    t_excess_internal: FloatFieldIJ,
    vapor_excess_input: FloatFieldIJ,
    vapor_excess_internal: FloatFieldIJ,
    ocean_fraction: FloatFieldIJ,
    use_excess: Int,
    t_input: FloatField,
    vapor_input: FloatField,
    grid_scale_forcing_t: FloatField,
    grid_scale_forcing_vapor: FloatField,
    subgrid_scale_forcing_t: FloatField,
    subgrid_scale_forcing_vapor: FloatField,
    t_internal: FloatField,
    vapor_internal: FloatField,
    t_pbl_internal: FloatField,
    vapor_pbl_internal: FloatField,
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
            t_excess_internal = 0.0
            vapor_excess_internal = 0.0
        elif use_excess == 2:
            t_excess_internal = min(0.5, max(0.2, t_excess_input))  # Kelvin
            vapor_excess_internal = min(5.0e-4, max(1.0e-4, vapor_excess_input))  # kg kg^-1
        else:
            if ocean_fraction > 0.98:  # ocean
                t_excess_internal = min(0.5, max(0.2, t_excess_input))  # Kelvin
                vapor_excess_internal = min(5.0e-4, max(1.0e-4, vapor_excess_input))  # kg kg^-1

    with computation(PARALLEL), interval(0, -1):
        t_internal = t_input + (subgrid_scale_forcing_t + grid_scale_forcing_t) * DT_MOIST
        vapor_internal = vapor_input + (subgrid_scale_forcing_vapor + grid_scale_forcing_vapor) * DT_MOIST
        vapor_internal = max(GF_2020_constants.smaller_qv, vapor_internal)

        # temp/water vapor modified only by bl processes
        t_pbl_internal = t_input + (subgrid_scale_forcing_t) * DT_MOIST
        vapor_pbl_internal = vapor_input + (subgrid_scale_forcing_vapor) * DT_MOIST

        # moist static energy
        moist_static_energy = GF_2020_constants.CP * (
            subgrid_scale_forcing_t + grid_scale_forcing_t
        ) + GF_2020_constants.XLV * (subgrid_scale_forcing_vapor + grid_scale_forcing_vapor)


def prefil_internal_fields(
    kbmax: IntFieldIJ,
    kstamb: IntFieldIJ,
    ocean_fraction_input: FloatFieldIJ,
    ocean_fraction_internal: FloatFieldIJ,
    cap_max: FloatFieldIJ,
    ierrc: FloatFieldIJ,
    CAP_MAX_INC: Float,
    cup_max_increment: FloatFieldIJ,
    geopotential_height_cu_param_input: FloatField,
    z: FloatField,
    xz: FloatField,
):
    from __external__ import k_end, CAP_MAXS

    with computation(FORWARD), interval(0, 1):
        kbmax = 1
        kstabm = k_end - 1
        ocean_fraction_internal = ocean_fraction_input
        cap_max = CAP_MAXS
        ierrc = -999
        cup_max_increment = CAP_MAX_INC

    with computation(PARALLEL), interval(...):
        z = geopotential_height_cu_param_input
        xz = geopotential_height_cu_param_input
