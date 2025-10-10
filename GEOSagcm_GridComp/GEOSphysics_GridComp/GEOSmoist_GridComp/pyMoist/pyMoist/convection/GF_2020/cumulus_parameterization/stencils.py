from ndsl.dsl.typing import FloatField, FloatFieldIJ, Float, Int
from ndsl.dsl.gt4py import computation, PARALLEL, interval, FORWARD
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as GF_2020_constants
import pyMoist.constants as constants


def prefil_excess(
    t_excess_cu_param_input: FloatFieldIJ,
    t_excess_cu_param_internal: FloatFieldIJ,
    vapor_excess_cu_param_input: FloatFieldIJ,
    vapor_excess_cu_param_internal: FloatFieldIJ,
    ocean_fraction: FloatFieldIJ,
    use_excess: Int,
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
    with computation(FORWARD), interval(0, 1):
        if use_excess == 0:
            t_excess_cu_param_internal = 0.0
            vapor_excess_cu_param_internal = 0.0
        elif use_excess == 2:
            t_excess_cu_param_internal = min(0.5, max(0.2, t_excess_cu_param_input))  # Kelvin
            vapor_excess_cu_param_internal = min(5.0e-4, max(1.0e-4, vapor_excess_cu_param_input))  # kg kg^-1
        else:
            if ocean_fraction > 0.98:  # ocean
                t_excess_cu_param_internal = min(0.5, max(0.2, t_excess_cu_param_input))  # Kelvin
                vapor_excess_cu_param_internal = min(
                    5.0e-4, max(1.0e-4, vapor_excess_cu_param_input)
                )  # kg kg^-1


def set_new_t_vapor(
    t: FloatField,
    vapor: FloatField,
    grid_scale_forcing_t: FloatField,
    grid_scale_forcing_vapor: FloatField,
    subgrid_scale_forcing_t: FloatField,
    subgrid_scale_forcing_vapor: FloatField,
    t_cu_param_internal: FloatField,
    vapor_cu_param_internal: FloatField,
    t_pbl_cu_param_internal: FloatField,
    vapor_pbl_cu_param_internal: FloatField,
    moist_static_energy: FloatField,
):
    from __externals__ import DT_MOIST

    with computation(PARALLEL), interval(0, -1):
        t_cu_param_internal = t + (subgrid_scale_forcing_t + grid_scale_forcing_t) * DT_MOIST
        vapor_cu_param_internal = vapor + (subgrid_scale_forcing_vapor + grid_scale_forcing_vapor) * DT_MOIST
        vapor_cu_param_internal = max(GF_2020_constants.smaller_qv, vapor_cu_param_internal)

        # temp/water vapor modified only by bl processes
        t_pbl_cu_param_internal = t + (subgrid_scale_forcing_t) * DT_MOIST
        vapor_pbl_cu_param_internal = vapor + (subgrid_scale_forcing_vapor) * DT_MOIST

        # moist static energy
        moist_static_energy = GF_2020_constants.CP * (
            subgrid_scale_forcing_t + grid_scale_forcing_t
        ) + GF_2020_constants.XLV * (subgrid_scale_forcing_vapor + grid_scale_forcing_vapor)
