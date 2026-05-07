from ndsl.dsl.gt4py import FORWARD, PARALLEL, computation, function, int32, interval, trunc
from ndsl.dsl.typing import BoolFieldIJ, Float, FloatField, FloatFieldIJ

from pyMoist.microphysics.GFDL_1M.driver.constants import constants
from pyMoist.microphysics.GFDL_1M.driver.sat_tables import GlobalTable_driver_qsat


@function
def wqs2(
    ta: Float,
    den: Float,
    table2: GlobalTable_driver_qsat,
    des2: GlobalTable_driver_qsat,
):
    """
    Compute the saturated specific humidity for table2
    with additional calculation of gradient (dq/dt)

    pure water phase; universal dry / moist formula using air density
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
    it = int32(trunc(ap1))
    es = table2.A[it - 1] + (ap1 - it) * des2.A[it - 1]
    qsat = es / (constants.RVGAS * ta * den)
    it = int32(trunc(ap1 - 0.5))  # check if this rounds or truncates. need truncation here
    # finite diff, del_t = 0.1:
    dqdt = 10.0 * (des2.A[it - 1] + (ap1 - it) * (des2.A[it] - des2.A[it - 1])) / (constants.RVGAS * ta * den)

    return qsat, dqdt


def implicit_fall(
    mixing_ratio: FloatField,
    terminal_speed: FloatField,
    z_interface: FloatField,
    dp: FloatField,
    mass: FloatField,
    precip_flux: FloatField,
    precip: FloatFieldIJ,
    precip_fall: BoolFieldIJ,
):
    """
    Compute the time-implicit monotonic scheme
    """
    from __externals__ import dts

    with computation(PARALLEL), interval(...):
        if precip_fall == True:  # noqa
            height_diff = z_interface - z_interface[0, 0, 1]
            dd = dts * terminal_speed
            mixing_ratio = mixing_ratio * dp

    # sedimentation: non - vectorizable loop
    with computation(FORWARD), interval(0, 1):
        if precip_fall == True:  # noqa
            qm = mixing_ratio / (height_diff + dd)

    with computation(FORWARD), interval(1, None):
        if precip_fall == True:  # noqa
            qm = (mixing_ratio + dd[0, 0, -1] * qm[0, 0, -1]) / (height_diff + dd)

    # qm is density at this stage
    with computation(PARALLEL), interval(...):
        if precip_fall == True:  # noqa
            qm = qm * height_diff

    # output mass fluxes: non - vectorizable loop
    with computation(FORWARD), interval(0, 1):
        if precip_fall == True:  # noqa
            mass = mixing_ratio - qm

    with computation(FORWARD), interval(1, None):
        if precip_fall == True:  # noqa
            mass = mass[0, 0, -1] + mixing_ratio - qm

    with computation(FORWARD), interval(-1, None):
        if precip_fall == True:  # noqa
            precip = mass
        else:
            precip = 0

    # update:
    with computation(PARALLEL), interval(...):
        if precip_fall == True:  # noqa
            mixing_ratio = qm / dp

    with computation(PARALLEL), interval(...):
        if precip_fall == True:  # noqa
            precip_flux = precip_flux + mass
