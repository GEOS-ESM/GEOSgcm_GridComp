import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import (
    FORWARD,
    PARALLEL,
    computation,
    i32,
    interval,
    trunc,
)

from pyMoist.GFDL_1M.driver.constants import constants
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ
from pyMoist.GFDL_1M.driver.sat_tables import GlobalTable_driver_qsat


@gtscript.function
def wqs2(
    ta: Float,
    den: Float,
    table2: GlobalTable_driver_qsat,
    des2: GlobalTable_driver_qsat,
):
    """
    compute the saturated specific humidity for table2
    with additional calculation of gradient (dq/dt)

    pure water phase; universal dry / moist formular using air density
    input "den" can be either dry or moist air density

    reference Fortran: gfdl_cloud_microphys.F90: function wqs2
    """
    tmin = constants.TABLE_ICE - 160.0

    if ta - tmin > 0:
        ans = ta - tmin
    else:
        ans = 0
    ap1 = 10.0 * ans + 1.0
    ap1 = min(2621.0, ap1)
    it = i32(trunc(ap1))
    es = table2.A[it - 1] + (ap1 - it) * des2.A[it - 1]
    qsat = es / (constants.RVGAS * ta * den)
    it = i32(
        trunc(ap1 - 0.5)
    )  # check if this rounds or truncates. need truncation here
    # finite diff, del_t = 0.1:
    dqdt = (
        10.0
        * (des2.A[it - 1] + (ap1 - it) * (des2.A[it] - des2.A[it - 1]))
        / (constants.RVGAS * ta * den)
    )

    return qsat, dqdt


def implicit_fall(
    q: FloatField,
    vt: FloatField,
    ze: FloatField,
    dp1: FloatField,
    m1: FloatField,
    m1_sol: FloatField,
    precip: FloatFieldIJ,
    precip_fall: FloatFieldIJ,
):
    """
    reference Fortran: gfdl_cloud_microphys.F90: subroutine implicit_fall
    this code computes the time-implicit monotonic scheme
    Fortran author: Shian-Jiann Lin, 2016
    """
    from __externals__ import dts, use_ppm

    with computation(PARALLEL), interval(...):
        if precip_fall == 1:
            if use_ppm == False:  # noqa
                height_diff = ze - ze[0, 0, 1]
                dd = dts * vt
                q = q * dp1

    # -----------------------------------------------------------------------
    # sedimentation: non - vectorizable loop
    # -----------------------------------------------------------------------
    with computation(FORWARD), interval(0, 1):
        if precip_fall == 1:
            if use_ppm == False:  # noqa
                qm = q / (height_diff + dd)

    with computation(FORWARD), interval(1, None):
        if precip_fall == 1:
            if use_ppm == False:  # noqa
                qm = (q + dd[0, 0, -1] * qm[0, 0, -1]) / (height_diff + dd)

    # -----------------------------------------------------------------------
    # qm is density at this stage
    # -----------------------------------------------------------------------
    with computation(PARALLEL), interval(...):
        if precip_fall == 1:
            if use_ppm == False:  # noqa
                qm = qm * height_diff

    # -----------------------------------------------------------------------
    # output mass fluxes: non - vectorizable loop
    # -----------------------------------------------------------------------
    with computation(FORWARD), interval(0, 1):
        if precip_fall == 1:
            if use_ppm == False:  # noqa
                m1 = q - qm

    with computation(FORWARD), interval(1, None):
        if precip_fall == 1:
            if use_ppm == False:  # noqa
                m1 = m1[0, 0, -1] + q - qm

    with computation(FORWARD), interval(-1, None):
        if precip_fall == 1:
            if use_ppm == False:  # noqa
                precip = m1

    # -----------------------------------------------------------------------
    # update:
    # -----------------------------------------------------------------------
    with computation(PARALLEL), interval(...):
        if precip_fall == 1:
            if use_ppm == False:  # noqa
                q = qm / dp1

    with computation(PARALLEL), interval(...):
        if precip_fall == 1:
            m1_sol = m1_sol + m1
