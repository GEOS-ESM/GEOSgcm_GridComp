from gt4py.cartesian.gtscript import (
    FORWARD,
    PARALLEL,
    computation,
    interval,
    sqrt,
)

import pyMoist.GFDL_1M.driver.constants as constants
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ


def update_tendencies(
    qv0: FloatField,
    ql0: FloatField,
    qr0: FloatField,
    qi0: FloatField,
    qs0: FloatField,
    qg0: FloatField,
    qa0: FloatField,
    qv1: FloatField,
    ql1: FloatField,
    qr1: FloatField,
    qi1: FloatField,
    qs1: FloatField,
    qg1: FloatField,
    qv_dt: FloatField,
    ql_dt: FloatField,
    qr_dt: FloatField,
    qi_dt: FloatField,
    qs_dt: FloatField,
    qg_dt: FloatField,
    qa_dt: FloatField,
    t: FloatField,
    t1: FloatField,
    t_dt: FloatField,
    w: FloatField,
    w1: FloatField,
    uin: FloatField,
    u1: FloatField,
    udt: FloatField,
    vin: FloatField,
    v1: FloatField,
    vdt: FloatField,
    dp: FloatField,
    dp1: FloatField,
    m1: FloatField,
    rain: FloatFieldIJ,  # strict output
    snow: FloatFieldIJ,  # strict output
    ice: FloatFieldIJ,  # strict output
    graupel: FloatFieldIJ,  # strict output
):
    """
    compute output tendencies of the microphysics driver

    reference Fortran: gfdl_cloud_microphys.F90:
    subroutine mpdrv, subroutine gfdl_cloud_microphys_driver
    """
    from __externals__ import c_air, c_vap, do_qa, do_sedi_w, rdt, sedi_transport

    # -----------------------------------------------------------------------
    # momentum transportation during sedimentation
    # note: dp1 is dry mass; dp0 is the old moist (total) mass
    # -----------------------------------------------------------------------

    with computation(PARALLEL), interval(1, None):
        if sedi_transport == True:  # noqa
            u1 = (dp * u1 + m1[0, 0, -1] * u1[0, 0, -1]) / (dp + m1[0, 0, -1])
            v1 = (dp * v1 + m1[0, 0, -1] * v1[0, 0, -1]) / (dp + m1[0, 0, -1])
            udt = udt + (u1 - uin) * rdt
            vdt = vdt + (v1 - vin) * rdt

    with computation(PARALLEL), interval(...):
        if do_sedi_w:
            w = w1

    # -----------------------------------------------------------------------
    # update moist air mass (actually hydrostatic pressure)
    # convert to dry mixing ratios
    # -----------------------------------------------------------------------

    with computation(PARALLEL), interval(...):
        omq = dp1 / dp
        qv_dt = qv_dt + rdt * (qv1 - qv0) * omq
        ql_dt = ql_dt + rdt * (ql1 - ql0) * omq
        qr_dt = qr_dt + rdt * (qr1 - qr0) * omq
        qi_dt = qi_dt + rdt * (qi1 - qi0) * omq
        qs_dt = qs_dt + rdt * (qs1 - qs0) * omq
        qg_dt = qg_dt + rdt * (qg1 - qg0) * omq
        cvm = (
            c_air
            + qv1 * c_vap
            + (qr1 + ql1) * constants.C_LIQ
            + (qi1 + qs1 + qg1) * constants.C_ICE
        )
        t_dt = t_dt + rdt * (t1 - t) * cvm / constants.CP_AIR

        # -----------------------------------------------------------------------
        # update cloud fraction tendency
        # -----------------------------------------------------------------------
        if do_qa == False:  # noqa
            qa_dt = qa_dt + rdt * (
                qa0 * sqrt((qi1 + ql1) / max(qi0 + ql0, constants.QCMIN)) - qa0
            )  # New Cloud - Old CloudCloud

    with computation(FORWARD), interval(0, 1):
        # convert to mm / day
        conversion_factor = 86400.0 * rdt * constants.RGRAV
        rain = rain * conversion_factor
        snow = snow * conversion_factor
        ice = ice * conversion_factor
        graupel = graupel * conversion_factor
