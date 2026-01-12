from ndsl import StencilFactory, QuantityFactory, ndsl_log
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import (
    GF2020CumulusParameterizationConfig,
)
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from pyMoist.convection.GF_2020.cumulus_parameterization.state import (
    GF2020CumulusParameterizationState,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.locals import (
    GF2020CumulusParameterizationLocals,
)
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)
from ndsl.dsl.typing import FloatField, FloatFieldIJ, Float, Int
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
import pyMoist.constants as constants
from gt4py.cartesian.gtscript import (
    PARALLEL,
    FORWARD,
    BACKWARD,
    computation,
    interval,
    exp,
)
from ndsl.dsl.gt4py import function
from pyMoist.convection.GF_2020.cumulus_parameterization.shared_functions import (
    saturation_vapor_pressure,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import (
    FloatField_Plume,
    FloatFieldIJ_Plume,
    IntFieldIJ_Plume,
)


@function
def saturation_specific_humidity(
    potential_t: Float,
    p: Float,
):
    """
    Compute saturation specific humidity in kg/kg

    Args:
        pt (in): potential temperature in Kelvin
        p (in): pressure in mb
    """
    # local constants, should not be used/visible outside of this function
    # not really a suggestion. don't use them anywhere else, use the global constants
    rd = 287.06
    rv = 461.52
    rtt = 273.16
    retv = rv / rd - 1.0
    r2es = 611.21 * rd / rv
    r3les = 17.502
    r3ies = 22.587
    r4les = 32.19
    r4ies = -0.7
    rtwat = rtt
    rtice = rtt - 23.0
    rticecu = rtt - 23.0
    rtwat_rtice_r = 1.0 / (rtwat - rtice)
    rtwat_rticecu_r = 1.0 / (rtwat - rticecu)

    foealfcu = min(1.0, ((max(rticecu, min(rtwat, potential_t)) - rticecu) * rtwat_rticecu_r) ** 2)
    foeewmcu = r2es * (
        foealfcu * exp(r3les * (potential_t - rtt) / (potential_t - r4les))
        + (1.0 - foealfcu) * exp(r3ies * (potential_t - rtt) / (potential_t - r4ies))
    )

    zew = foeewmcu
    zqs = zew / (100.0 * p)
    if 1.0 - retv * zqs > 0.0:
        zcor = 1.0 / (1.0 - retv * zqs)
        saturation_specific_humidity = zqs * zcor
    else:
        saturation_specific_humidity = cumulus_parameterization_constants.MAX_QSAT

    return saturation_specific_humidity


def environment_conditions(
    p: FloatField,
    p_surface: FloatFieldIJ,
    t: FloatField,
    vapor: FloatField,
    topography_height_no_negative: FloatFieldIJ,
    moist_static_energy: FloatField,
    saturation_moist_static_energy: FloatField,
    saturation_mixing_ratio: FloatField,
    geopotential_height: FloatField,
    error_code: IntFieldIJ_Plume,
    plume: Int,
):
    """
    Compute moist static energy, geopotential height and saturation
    mixing ratio for a set of environmental fields.

    Args:
        p (in): pressure
        p_surface (in): surface pressure
        t (in): temperature
        vapor (in): water vapor mixing ratio
        topography_height_no_negative (in): topography height, without features extending below sea level
        moist_static_energy (out)
        saturation_moist_static_energy (out)
        saturation_mixing_ratio (out)
        geopotential_height (in)
        error_code (out)
        plume (in): specifies the current plume

    """
    from __externals__ import SATURATION_CALCULATION_CHOICE

    with computation(PARALLEL), interval(...):
        itest_internal = -1
        moist_static_energy = 0.0
        saturation_moist_static_energy = 0.0
        saturation_mixing_ratio = 0.0
    with computation(PARALLEL), interval(0, -1):
        if SATURATION_CALCULATION_CHOICE == 0:
            if error_code[0, 0][plume] == 0:
                e = saturation_vapor_pressure(t)
                saturation_mixing_ratio = 0.622 * e / max(1.0e-8, (p - e))

                if saturation_mixing_ratio <= 1.0e-08:
                    saturation_mixing_ratio = 1.0e-08
                if saturation_mixing_ratio > cumulus_parameterization_constants.MAX_QSAT:
                    saturation_mixing_ratio = cumulus_parameterization_constants.MAX_QSAT
                if saturation_mixing_ratio < vapor:
                    saturation_mixing_ratio = vapor

                TV = t + 0.608 * vapor * t
        else:
            if error_code[0, 0][plume] == 0:
                saturation_mixing_ratio = saturation_specific_humidity(t, p)
                saturation_mixing_ratio = min(
                    cumulus_parameterization_constants.MAX_QSAT,
                    max(1.0e-08, saturation_mixing_ratio),
                )
                saturation_mixing_ratio = max(saturation_mixing_ratio, vapor)
                tv = t + 0.608 * vapor * t

    with computation(PARALLEL), interval(0, -1):
        if error_code[0, 0][plume] == 0:
            moist_static_energy = (
                constants.MAPL_GRAV * geopotential_height
                + cumulus_parameterization_constants.CP * t
                + cumulus_parameterization_constants.XLV * vapor
            )

            saturation_moist_static_energy = (
                constants.MAPL_GRAV * geopotential_height
                + cumulus_parameterization_constants.CP * t
                + cumulus_parameterization_constants.XLV * saturation_mixing_ratio
            )

            if moist_static_energy >= saturation_moist_static_energy:
                moist_static_energy = saturation_moist_static_energy


@function
def get_interp(
    p,
    t,
    vapor,
):
    # bunch of interal constants, not to be used outside of this function
    rd = 287.06
    rv = 461.52
    rcpd = 1004.71
    rtt = 273.16
    rhoh2o = 1000.0
    rlvtt = 2.5008e6
    rlstt = 2.8345e6
    retv = rv / rd - 1.0
    rlmlt = rlstt - rlvtt
    rcpv = 4.0 * rv
    r2es = 611.21 * rd / rv
    r3les = 17.502
    r3ies = 22.587
    r4les = 32.19
    r4ies = -0.7
    r5les = r3les * (rtt - r4les)
    r5ies = r3ies * (rtt - r4ies)
    r5alvcp = r5les * rlvtt / rcpd
    r5alscp = r5ies * rlstt / rcpd
    ralvdcp = rlvtt / rcpd
    ralsdcp = rlstt / rcpd
    ralfdcp = rlmlt / rcpd
    rtwat = rtt
    rtber = rtt - 5.0
    rtbercu = rtt - 5.0
    rtice = rtt - 23.0
    rticecu = rtt - 23.0
    rtwat_rtice_r = 1.0 / (rtwat - rtice)
    rtwat_rticecu_r = 1.0 / (rtwat - rticecu)
    rvtmp2 = rcpv / rcpd - 1.0
    zqmax = 0.5

    pt = t  # k
    pq = vapor  # kg/kg
    psp = p * 100.0  # hpa

    zqp = 1.0 / psp

    iteration = 0
    while iteration <= 1:
        ptare = pt

        foealfcu = min(1.0, ((max(rticecu, min(rtwat, ptare)) - rticecu) * rtwat_rticecu_r) ** 2)
        foeewmcu = r2es * (
            foealfcu * exp(r3les * (ptare - rtt) / (ptare - r4les))
            + (1.0 - foealfcu) * exp(r3ies * (ptare - rtt) / (ptare - r4ies))
        )
        zqsat = foeewmcu * zqp

        zcor = 1.0 / (1.0 - retv * zqsat)
        zqsat = zqsat * zcor

        foedemcu = foealfcu * r5alvcp * (1.0 / (ptare - r4les) ** 2) + (1.0 - foealfcu) * r5alscp * (
            1.0 / (ptare - r4ies) ** 2
        )

        zcond1 = (pq - zqsat) / (1.0 + zqsat * zcor * foedemcu)

        foeldcpmcu = foealfcu * ralvdcp + (1.0 - foealfcu) * ralsdcp
        pt = pt + foeldcpmcu * zcond1
        pq = pq - zcond1

        iteration += 1

    t_new = pt
    vapor_forced = pq

    return t_new, vapor_forced


def environment_cloud_levels(
    p: FloatField,
    p_surface: FloatFieldIJ,
    p_cloud_levels: FloatField,
    topography_height_no_negative: FloatFieldIJ,
    geopotential_height: FloatField,
    geopotential_height_cloud_levels: FloatField,
    t: FloatField,
    t_surface: FloatFieldIJ,
    t_cloud_levels: FloatField,
    vapor: FloatField,
    vapor_cloud_levels: FloatField,
    u: FloatField,
    v: FloatField,
    u_cloud_levels: FloatField,
    v_cloud_levels: FloatField,
    environment_saturation_mixing_ratio: FloatField,
    environment_saturation_mixing_ratio_cloud_levels: FloatField,
    environment_moist_static_energy: FloatField,
    environment_moist_static_energy_cloud_levels: FloatField,
    environment_saturation_moist_static_energy: FloatField,
    environment_saturation_moist_static_energy_cloud_levels: FloatField,
    gamma_cloud_levels: FloatField,
    error_code: IntFieldIJ_Plume,
    plume: Int,
):
    """
    Compute environmental fields at cloud levels.

    Args:
        p (in): pressure
        p_surface (in): surface pressure
        p_cloud_levels (out): pressure at cloud levels
        topography_height_no_negative (in): topography height, without features extending below sea level
        geopotential_height (in)
        geopotential_height_cloud_levels (out): geopotential height at cloud levels
        t (in): temperature
        t_surface (in): surface temperature
        t_cloud_levels: temperature at cloud levels
        vapor (in): water vapor mixing ratio
        vapor_cloud_levels (out): water vapor mixing ratio at cloud levels
        u (in): zonal wind
        v (in): meridional wind
        u_cloud_levels (out): zonal wind at cloud levels
        v_cloud_levels (out): meridional wind at cloud levels
        environment_saturation_mixing_ratio (in)
        environment_saturation_mixing_ratio_cloud_levels (out)
        environment_moist_static_energy (in)
        environment_moist_static_energy_cloud_levels (out)
        environment_saturation_moist_static_energy (in)
        environment_saturation_moist_static_energy_cloud_levels (out)
        gamma_cloud_levels (out)
        error_code (in)
        plume (in): specifies the current plume
    """
    from __externals__ import CLOUD_LEVEL_GRID

    with computation(PARALLEL), interval(...):
        # ensure no data from previous call slips through untouched
        environment_saturation_mixing_ratio_cloud_levels = 0.0
        vapor_cloud_levels = 0.0
        environment_saturation_moist_static_energy_cloud_levels = 0.0
        environment_moist_static_energy_cloud_levels = 0.0
        geopotential_height_cloud_levels = 0.0
        p_cloud_levels = 0.0
        t_cloud_levels = 0.0
        gamma_cloud_levels = 0.0
        u_cloud_levels = 0.0
        v_cloud_levels = 0.0

    with computation(PARALLEL), interval(1, -1):
        if CLOUD_LEVEL_GRID == 2:
            # original formulation
            if error_code[0, 0][plume] == 0:
                environment_saturation_mixing_ratio_cloud_levels = 0.5 * (
                    environment_saturation_mixing_ratio[0, 0, -1] + environment_saturation_mixing_ratio
                )
                vapor_cloud_levels = 0.5 * (vapor[0, 0, -1] + vapor)
                environment_saturation_moist_static_energy_cloud_levels = 0.5 * (
                    environment_saturation_moist_static_energy[0, 0, -1]
                    + environment_saturation_moist_static_energy
                )
                environment_moist_static_energy_cloud_levels = 0.5 * (
                    environment_moist_static_energy[0, 0, -1] + environment_moist_static_energy
                )
                if (
                    environment_moist_static_energy_cloud_levels
                    > environment_saturation_moist_static_energy_cloud_levels
                ):
                    environment_moist_static_energy_cloud_levels = (
                        environment_saturation_moist_static_energy_cloud_levels
                    )
                geopotential_height_cloud_levels = 0.5 * (geopotential_height[0, 0, -1] + geopotential_height)
                p_cloud_levels = 0.5 * (p[0, 0, -1] + p)
                t_cloud_levels = 0.5 * (t[0, 0, -1] + t)
                gamma_cloud_levels = (
                    (cumulus_parameterization_constants.XLV / cumulus_parameterization_constants.CP)
                    * (
                        cumulus_parameterization_constants.XLV
                        / (cumulus_parameterization_constants.RV * t_cloud_levels * t_cloud_levels)
                    )
                    * environment_saturation_mixing_ratio_cloud_levels
                )
                u_cloud_levels = 0.5 * (u[0, 0, -1] + u)
                v_cloud_levels = 0.5 * (v[0, 0, -1] + v)

    with computation(FORWARD), interval(0, 1):
        if CLOUD_LEVEL_GRID == 2:
            # original formulation
            if error_code[0, 0][plume] == 0:
                environment_saturation_mixing_ratio_cloud_levels = environment_saturation_mixing_ratio
                vapor_cloud_levels = vapor
                environment_saturation_moist_static_energy_cloud_levels = (
                    constants.MAPL_GRAV * topography_height_no_negative
                    + cumulus_parameterization_constants.CP * t
                    + cumulus_parameterization_constants.XLV * environment_saturation_mixing_ratio
                )
                environment_moist_static_energy_cloud_levels = (
                    constants.MAPL_GRAV * topography_height_no_negative
                    + cumulus_parameterization_constants.CP * t
                    + cumulus_parameterization_constants.XLV * vapor
                )
                geopotential_height_cloud_levels = topography_height_no_negative
                p_cloud_levels = p_surface
                t_cloud_levels = t
                gamma_cloud_levels = (
                    cumulus_parameterization_constants.XLV
                    / cumulus_parameterization_constants.CP
                    * (
                        cumulus_parameterization_constants.XLV
                        / (cumulus_parameterization_constants.RV * t_cloud_levels * t_cloud_levels)
                    )
                    * environment_saturation_mixing_ratio_cloud_levels
                )
                u_cloud_levels = u
                v_cloud_levels = v

    with computation(FORWARD), interval(0, 1):
        if CLOUD_LEVEL_GRID == 0:
            # weighted mean
            if error_code[0, 0][plume] == 0:
                p_cloud_levels = p_surface
                geopotential_height_cloud_levels = topography_height_no_negative

    with computation(FORWARD), interval(0, -2):
        if CLOUD_LEVEL_GRID == 0:
            # weighted mean
            if error_code[0, 0][plume] == 0:
                p_cloud_levels[0, 0, 1] = 2.0 * p - p_cloud_levels
                geopotential_height_cloud_levels[0, 0, 1] = (
                    2.0 * geopotential_height - geopotential_height_cloud_levels
                )

    with computation(FORWARD), interval(0, -2):
        if CLOUD_LEVEL_GRID == 0:
            # weighted mean
            if error_code[0, 0][plume] == 0:
                p1 = abs((p[0, 0, 1] - p_cloud_levels[0, 0, 1]) / (p[0, 0, 1] - p))
                p2 = abs((p_cloud_levels[0, 0, 1] - p) / (p[0, 0, 1] - p))
                t_cloud_levels[0, 0, 1] = p1 * t + p2 * t[0, 0, 1]
                u_cloud_levels[0, 0, 1] = p1 * u + p2 * u[0, 0, 1]
                v_cloud_levels[0, 0, 1] = p1 * v + p2 * v[0, 0, 1]
                vapor_cloud_levels[0, 0, 1] = p1 * vapor + p2 * vapor[0, 0, 1]
                environment_moist_static_energy_cloud_levels[0, 0, 1] = (
                    p1 * environment_moist_static_energy + p2 * environment_moist_static_energy[0, 0, 1]
                )
                environment_saturation_mixing_ratio_cloud_levels[0, 0, 1] = (
                    p1 * environment_saturation_mixing_ratio
                    + p2 * environment_saturation_mixing_ratio[0, 0, 1]
                )
                environment_saturation_moist_static_energy_cloud_levels[0, 0, 1] = (
                    p1 * environment_saturation_moist_static_energy
                    + p2 * environment_saturation_moist_static_energy[0, 0, 1]
                )
                if (
                    environment_moist_static_energy_cloud_levels[0, 0, 1]
                    > environment_saturation_moist_static_energy_cloud_levels[0, 0, 1]
                ):
                    environment_moist_static_energy_cloud_levels[0, 0, 1] = (
                        environment_saturation_moist_static_energy_cloud_levels[0, 0, 1]
                    )

                gamma_cloud_levels[0, 0, 1] = (
                    (cumulus_parameterization_constants.XLV / cumulus_parameterization_constants.CP)
                    * (
                        cumulus_parameterization_constants.XLV
                        / (
                            cumulus_parameterization_constants.RV
                            * t_cloud_levels[0, 0, 1]
                            * t_cloud_levels[0, 0, 1]
                        )
                    )
                    * environment_saturation_mixing_ratio_cloud_levels[0, 0, 1]
                )

    with computation(FORWARD), interval(0, 1):
        if CLOUD_LEVEL_GRID == 0:
            # weighted mean
            # surface level: using level 0 and 1 to determine level 0
            if error_code[0, 0][plume] == 0:
                p1 = abs(p - p_cloud_levels)
                p2 = abs(p_cloud_levels[0, 0, 1] - p_cloud_levels)

                ct1 = (p1 + p2) / p2
                ct2 = p1 / p2

                t_cloud_levels = ct1 * t - ct2 * t_cloud_levels[0, 0, 1]
                vapor_cloud_levels = ct1 * vapor - ct2 * vapor_cloud_levels[0, 0, 1]

                u_cloud_levels = ct1 * u - ct2 * u_cloud_levels[0, 0, 1]
                v_cloud_levels = ct1 * v - ct2 * v_cloud_levels[0, 0, 1]
                environment_saturation_mixing_ratio_cloud_levels = (
                    ct1 * environment_saturation_mixing_ratio
                    - ct2 * environment_saturation_mixing_ratio_cloud_levels[0, 0, 1]
                )

                environment_saturation_moist_static_energy_cloud_levels = (
                    constants.MAPL_GRAV * geopotential_height_cloud_levels
                    + cumulus_parameterization_constants.CP * t_cloud_levels
                    + cumulus_parameterization_constants.XLV
                    * environment_saturation_mixing_ratio_cloud_levels
                )
                environment_moist_static_energy_cloud_levels = (
                    constants.MAPL_GRAV * geopotential_height_cloud_levels
                    + cumulus_parameterization_constants.CP * t_cloud_levels
                    + cumulus_parameterization_constants.XLV * vapor_cloud_levels
                )

                if (
                    environment_moist_static_energy_cloud_levels
                    > environment_saturation_moist_static_energy_cloud_levels
                ):
                    environment_moist_static_energy_cloud_levels = (
                        environment_saturation_moist_static_energy_cloud_levels
                    )

                gamma_cloud_levels = (
                    cumulus_parameterization_constants.XLV
                    / cumulus_parameterization_constants.CP
                    * (
                        cumulus_parameterization_constants.XLV
                        / (cumulus_parameterization_constants.RV * t_cloud_levels * t_cloud_levels)
                    )
                    * environment_saturation_mixing_ratio_cloud_levels
                )

    with computation(BACKWARD), interval(1, -1):
        if CLOUD_LEVEL_GRID == 1:
            # based on Tiedke (1989)
            if error_code[0, 0][plume] == 0:
                environment_saturation_mixing_ratio_cloud_levels = environment_saturation_mixing_ratio
                vapor_cloud_levels = vapor
                p_cloud_levels = 0.5 * (p[0, 0, -1] + p)
                geopotential_height_cloud_levels = 0.5 * (geopotential_height[0, 0, -1] + geopotential_height)
                t_cloud_levels = (
                    max(
                        cumulus_parameterization_constants.CP * t[0, 0, -1]
                        + constants.MAPL_GRAV * geopotential_height[0, 0, -1],
                        cumulus_parameterization_constants.CP * t + constants.MAPL_GRAV * geopotential_height,
                    )
                    - constants.MAPL_GRAV * geopotential_height_cloud_levels
                ) / cumulus_parameterization_constants.CP

                if environment_saturation_mixing_ratio < cumulus_parameterization_constants.MAX_QSAT:
                    t_cloud_levels, environment_saturation_mixing_ratio_cloud_levels = get_interp(
                        p=p_cloud_levels,
                        t=t_cloud_levels,
                        vapor=environment_saturation_mixing_ratio_cloud_levels,
                    )

                vapor_cloud_levels = (
                    min(vapor, environment_saturation_mixing_ratio)
                    + environment_saturation_mixing_ratio_cloud_levels
                    - environment_saturation_mixing_ratio
                )
                vapor_cloud_levels = max(vapor_cloud_levels, 0.0)

    with computation(FORWARD), interval(0, 1):
        if CLOUD_LEVEL_GRID == 1:
            # based on Tiedke (1989)
            if error_code[0, 0][plume] == 0:
                # surface
                environment_saturation_mixing_ratio_cloud_levels = environment_saturation_mixing_ratio
                vapor_cloud_levels = vapor
                geopotential_height_cloud_levels = topography_height_no_negative
                p_cloud_levels = p_surface

                t_cloud_levels = (
                    cumulus_parameterization_constants.CP * t
                    + constants.MAPL_GRAV * geopotential_height
                    - constants.MAPL_GRAV * geopotential_height_cloud_levels
                ) / cumulus_parameterization_constants.CP

                environment_saturation_moist_static_energy_cloud_levels = (
                    constants.MAPL_GRAV * geopotential_height_cloud_levels
                    + cumulus_parameterization_constants.CP * t_cloud_levels
                    + cumulus_parameterization_constants.XLV
                    * environment_saturation_mixing_ratio_cloud_levels
                )
                environment_moist_static_energy_cloud_levels = (
                    constants.MAPL_GRAV * geopotential_height_cloud_levels
                    + cumulus_parameterization_constants.CP * t_cloud_levels
                    + cumulus_parameterization_constants.XLV * vapor_cloud_levels
                )

                gamma_cloud_levels = (
                    cumulus_parameterization_constants.XLV
                    / cumulus_parameterization_constants.CP
                    * (
                        cumulus_parameterization_constants.XLV
                        / (cumulus_parameterization_constants.RV * t_cloud_levels * t_cloud_levels)
                    )
                    * environment_saturation_mixing_ratio_cloud_levels
                )
                u_cloud_levels = u
                v_cloud_levels = v

    with computation(BACKWARD), interval(1, -1):
        if CLOUD_LEVEL_GRID == 1:
            # based on Tiedke (1989)
            if error_code[0, 0][plume] == 0:
                p1 = max(
                    cumulus_parameterization_constants.CP * t_cloud_levels
                    + constants.MAPL_GRAV * geopotential_height_cloud_levels,
                    cumulus_parameterization_constants.CP * t_cloud_levels[0, 0, -1]
                    + constants.MAPL_GRAV * geopotential_height_cloud_levels[0, 0, -1],
                )
                t_cloud_levels = (
                    p1 - constants.MAPL_GRAV * geopotential_height_cloud_levels
                ) / cumulus_parameterization_constants.CP

                environment_saturation_moist_static_energy_cloud_levels = (
                    cumulus_parameterization_constants.CP * t_cloud_levels
                    + cumulus_parameterization_constants.XLV
                    * environment_saturation_mixing_ratio_cloud_levels
                    + constants.MAPL_GRAV * geopotential_height_cloud_levels
                )
                environment_moist_static_energy_cloud_levels = (
                    cumulus_parameterization_constants.CP * t_cloud_levels
                    + cumulus_parameterization_constants.XLV * vapor_cloud_levels
                    + constants.MAPL_GRAV * geopotential_height_cloud_levels
                )
                if (
                    environment_moist_static_energy_cloud_levels
                    > environment_saturation_moist_static_energy_cloud_levels
                ):
                    environment_moist_static_energy_cloud_levels = (
                        environment_saturation_moist_static_energy_cloud_levels
                    )

                gamma_cloud_levels = (
                    (cumulus_parameterization_constants.XLV / cumulus_parameterization_constants.CP)
                    * (
                        cumulus_parameterization_constants.XLV
                        / (cumulus_parameterization_constants.RV * t_cloud_levels * t_cloud_levels)
                    )
                    * environment_saturation_mixing_ratio_cloud_levels
                )
                u_cloud_levels = u
                v_cloud_levels = v


def environment_mass_flux(
    error_code: IntFieldIJ_Plume,
    epsilon_forced: FloatFieldIJ_Plume,
    normalized_massflux_updraft_forced: FloatField_Plume,
    normalized_massflux_downdraft_forced: FloatField_Plume,
    environment_massflux: FloatField,
    plume: Int,
):
    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0:
            environment_massflux = (
                normalized_massflux_updraft_forced[0, 0, 0][plume]
                - epsilon_forced[0, 0][plume] * normalized_massflux_downdraft_forced[0, 0, 0][plume]
            )
        else:
            environment_massflux = 0.0


def modify_environment_profiles(
    
):
class EnvironmentalSubsidence:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GF2020Config,
        cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
    ):
        self.config = config
        self.cumulus_parameterization_config = cumulus_parameterization_config

    def __call__(self, *args, **kwds):
        if self.cumulus_parameterization_config.APPLY_SUB_MP != 0:
            raise NotImplementedError(
                "[NDSL] GF2020-->CumulusParameterization-->EnvironmentalSubsidence this code"
                "has not been impemented. You should have been caught before getting here, but here we are."
                "Please choose another option for APPLY_SUB_MP or implement to continue."
            )
