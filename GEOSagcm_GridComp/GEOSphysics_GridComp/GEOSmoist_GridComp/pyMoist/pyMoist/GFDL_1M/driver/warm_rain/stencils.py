import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import (
    BACKWARD,
    FORWARD,
    PARALLEL,
    computation,
    exp,
    interval,
    log,
    max,
    sqrt,
)

from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ
from pyMoist.GFDL_1M.driver.constants import constants
from pyMoist.GFDL_1M.driver.sat_tables import GlobalTable_driver_qsat
from pyMoist.GFDL_1M.driver.stencils import wqs2


@gtscript.function
def revap_racc(
    t1: Float,
    qv1: Float,
    ql1: Float,
    qr1: Float,
    qi1: Float,
    qs1: Float,
    qg1: Float,
    qa1: Float,
    den1: Float,
    denfac: Float,
    rh_limited: Float,
    table1: GlobalTable_driver_qsat,
    table2: GlobalTable_driver_qsat,
    table3: GlobalTable_driver_qsat,
    table4: GlobalTable_driver_qsat,
    des1: GlobalTable_driver_qsat,
    des2: GlobalTable_driver_qsat,
    des3: GlobalTable_driver_qsat,
    des4: GlobalTable_driver_qsat,
    c_air: Float,
    c_vap: Float,
    cracw: Float,
    crevp_0: Float,
    crevp_1: Float,
    crevp_2: Float,
    crevp_3: Float,
    crevp_4: Float,
    d0_vap: Float,
    lv00: Float,
    tau_revp: Float,
    dts: Float,
):
    """
    Evaporate rain

    reference Fortran: gfdl_cloud_microphys.F90: subroutine revap_racc
    """

    revap = 0.0

    if t1 > constants.T_WFR and qr1 > constants.QPMIN:
        # area and timescale efficiency on revap
        fac_revp = 1.0 - exp(-(0.5 * dts) / tau_revp)

        # -----------------------------------------------------------------------
        # define heat capacity and latent heat coefficient
        # -----------------------------------------------------------------------

        lhl = lv00 + d0_vap * t1
        q_liq = ql1 + qr1
        q_sol = qi1 + qs1 + qg1
        cvm = c_air + qv1 * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
        lcpk = lhl / cvm

        tin = t1 - lcpk * ql1  # presence of clouds suppresses the rain evap
        qpz = qv1 + ql1
        qsat, dqsdt = wqs2(tin, den1, table2, des2)
        dqh = max(ql1, rh_limited * max(qpz, constants.QCMIN))
        dqh = min(dqh, 0.2 * qpz)  # new limiter
        dqv = qsat - qv1  # use this to prevent super - sat the gird box
        q_minus = qpz - dqh
        q_plus = qpz + dqh

        # -----------------------------------------------------------------------
        # qsat must be > q_minus to activate evaporation
        # qsat must be < q_plus to activate accretion
        # -----------------------------------------------------------------------

        # -----------------------------------------------------------------------
        # rain evaporation
        # -----------------------------------------------------------------------

        if dqv > constants.QVMIN and qsat > q_minus:
            if qsat > q_plus:
                dq = qsat - qpz
            else:
                # -----------------------------------------------------------------------
                # q_minus < qsat < q_plus
                # dq == dqh if qsat == q_minus
                # -----------------------------------------------------------------------
                dq = 0.25 * (q_minus - qsat) ** 2 / dqh
            qden = qr1 * den1
            t2 = tin * tin
            evap = (
                crevp_0
                * t2
                * dq
                * (crevp_1 * sqrt(qden) + crevp_2 * exp(0.725 * log(qden)))
                / (crevp_3 * t2 + crevp_4 * qsat * den1)
            )
            evap = min(
                qr1, min(0.5 * dts * fac_revp * evap, dqv / (1.0 + lcpk * dqsdt))
            )
            qr1 = qr1 - evap
            qv1 = qv1 + evap
            q_liq = q_liq - evap
            cvm = (
                c_air + qv1 * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
            )
            t1 = t1 - evap * lhl / cvm
            revap = evap / (0.5 * dts)

        # -----------------------------------------------------------------------
        # accretion: pracc
        # -----------------------------------------------------------------------

        if qr1 > constants.QPMIN and ql1 > constants.QCMIN and qsat < q_minus:
            sink = 0.5 * dts * denfac * cracw * exp(0.95 * log(qr1 * den1))
            sink = sink / (1.0 + sink) * ql1
            ql1 = ql1 - sink
            qr1 = qr1 + sink

    return (
        t1,
        qv1,
        qr1,
        ql1,
        qi1,
        qs1,
        qg1,
        qa1,
        revap,
    )


def warm_rain_step_1(
    dp1: FloatField,
    dz1: FloatField,
    t1: FloatField,
    qv1: FloatField,
    ql1: FloatField,
    qr1: FloatField,
    qi1: FloatField,
    qs1: FloatField,
    qg1: FloatField,
    qa1: FloatField,
    ccn: FloatField,
    den1: FloatField,
    denfac: FloatField,
    c_praut: FloatField,
    vtr: FloatField,
    evap1: FloatField,
    rh_limited: FloatField,
    eis: FloatFieldIJ,
    onemsig: FloatFieldIJ,
    ze: FloatField,
    dm: FloatField,
    precip_fall: FloatFieldIJ,
    table1: GlobalTable_driver_qsat,
    table2: GlobalTable_driver_qsat,
    table3: GlobalTable_driver_qsat,
    table4: GlobalTable_driver_qsat,
    des1: GlobalTable_driver_qsat,
    des2: GlobalTable_driver_qsat,
    des3: GlobalTable_driver_qsat,
    des4: GlobalTable_driver_qsat,
):
    """
    Warm rain cloud microphysics: evaporation, accretion

    first half of the timestep
    """
    from __externals__ import (
        c_air,
        c_vap,
        const_vr,
        cracw,
        crevp_0,
        crevp_1,
        crevp_2,
        crevp_3,
        crevp_4,
        d0_vap,
        do_qa,
        do_sedi_w,
        dts,
        irain_f,
        lv00,
        ql0_max,
        rthreshs,
        rthreshu,
        tau_revp,
        vr_fac,
        vr_max,
        z_slope_liq,
    )

    # reference Fortran: gfdl_cloud_microphys.F90: subroutine check_column
    # determine if any precip falls in the column
    # if it falls anywhere in the column, the entire column becomes true
    # initalized to 0 (false), potentially changed to 1 (true)
    with computation(FORWARD), interval(...):
        if qr1 > constants.QPMIN:
            precip_fall = 1
    # end reference Fortran: gfdl_cloud_microphys.F90: subroutine check_column

    # -----------------------------------------------------------------------
    # auto - conversion
    # assuming linear subgrid vertical distribution of cloud water
    # following lin et al. 1994, mwr
    # -----------------------------------------------------------------------

    with computation(PARALLEL), interval(...):
        # Use In-Cloud condensates
        if do_qa == False:  # noqa
            qadum = max(qa1, constants.QCMIN)
        else:
            qadum = 1.0
        ql1 = ql1 / qadum
        qi1 = qi1 / qadum

        fac_rc = (
            min(1.0, eis / 15.0) ** 2
        )  # Estimated inversion strength determine stable regime
        fac_rc = constants.RC * (rthreshs * fac_rc + rthreshu * (1.0 - fac_rc)) ** 3

    with computation(PARALLEL), interval(...):
        if irain_f != 0:
            # -----------------------------------------------------------------------
            # no subgrid varaibility
            # -----------------------------------------------------------------------
            if qadum > onemsig:
                if t1 > constants.T_WFR:
                    qc = fac_rc * ccn / den1
                    dq = ql1 - qc
                    if dq > 0.0:
                        sink = min(
                            dq,
                            dts * c_praut * den1 * exp(constants.SO3 * log(ql1)),
                        )
                        sink = min(ql0_max, min(ql1, max(0.0, sink)))
                        ql1 = ql1 - sink
                        qr1 = qr1 + sink * qadum

    with computation(FORWARD), interval(1, None):
        if irain_f == 0:
            # -----------------------------------------------------------------------
            # with subgrid variability
            # -----------------------------------------------------------------------

            # begin reference Fortran: gfdl_cloud_microphys.F90: subroutine linear_prof
            # definition of vertical subgrid variability
            # used for cloud ice and cloud water autoconversion
            # qi -- > ql & ql -- > qr
            # edges: qe == qbar + / - dm
            if z_slope_liq == True:  # noqa
                dql = 0.5 * (ql1 - ql1[0, 0, -1])

    # -----------------------------------------------------------------------
    # use twice the strength of the positive definiteness limiter (lin et al 1994)
    # -----------------------------------------------------------------------

    with computation(FORWARD), interval(1, -1):
        if irain_f == 0:
            if z_slope_liq == True:  # noqa
                dl = 0.5 * min(abs(dql + dql[0, 0, 1]), 0.5 * ql1)
                if dql * dql[0, 0, 1] <= 0.0:
                    if dql > 0.0:  # local max
                        dl = min(dl, min(dql, -dql[0, 0, 1]))
                    else:
                        dl = 0.0

    with computation(FORWARD), interval(0, 1):
        if irain_f == 0:
            if z_slope_liq == True:  # noqa
                dl = 0

    with computation(FORWARD), interval(-1, None):
        if irain_f == 0:
            if z_slope_liq == True:  # noqa
                dl = 0

    # -----------------------------------------------------------------------
    # impose a presumed background horizontal variability
    # that is proportional to the value itself
    # -----------------------------------------------------------------------
    with computation(PARALLEL), interval(...):
        if irain_f == 0:
            if z_slope_liq == True:  # noqa
                dl = max(dl, max(constants.QVMIN, rh_limited * ql1))
            if z_slope_liq == False:  # noqa
                dl = max(constants.QVMIN, rh_limited * ql1)

            # end reference Fortran: gfdl_cloud_microphys.F90: subroutine linear_prof

            if qadum > onemsig:
                if t1 > constants.T_WFR + constants.DT_FR:
                    dl = min(max(constants.QCMIN, dl), 0.5 * ql1)
                    # --------------------------------------------------------------------
                    # as in klein's gfdl am2 stratiform scheme (with subgrid variations)
                    # --------------------------------------------------------------------
                    qc = fac_rc * ccn / den1
                    dq = 0.5 * (ql1 + dl - qc)
                    # --------------------------------------------------------------------
                    # dq = dl if qc == q_minus = ql - dl
                    # dq = 0 if qc == q_plus = ql + dl
                    # --------------------------------------------------------------------
                    if dq > 0.0:  # q_plus > qc
                        # --------------------------------------------------------------------
                        # revised continuous form: linearly decays
                        # (with subgrid dl) to zero at qc == ql + dl
                        # --------------------------------------------------------------------
                        sink = (
                            min(1.0, dq / dl)
                            * dts
                            * c_praut
                            * den1
                            * exp(constants.SO3 * log(ql1))
                        )
                        sink = min(ql0_max, min(ql1, max(0.0, sink)))
                        ql1 = ql1 - sink
                        qr1 = qr1 + sink * qadum

        # Revert In-Cloud condensate
        ql1 = ql1 * qadum
        qi1 = qi1 * qadum

        # -----------------------------------------------------------------------
        # fall speed of rain
        # -----------------------------------------------------------------------

        if precip_fall == 0:
            vtr = constants.VF_MIN
        elif const_vr == True:  # noqa
            vtr = vr_fac  # ifs_2016: 4.0
        else:
            qden = qr1 * den1
            if qr1 < constants.THR:
                vtr = constants.VR_MIN
            else:
                vtr = (
                    vr_fac
                    * constants.VCONR
                    * sqrt(min(10.0, constants.SFCRHO / den1))
                    * exp(0.2 * log(qden / constants.NORMR))
                )
                vtr = min(vr_max, max(constants.VR_MIN, vtr))

    with computation(FORWARD), interval(-1, None):
        ze[0, 0, 1] = constants.ZS

    with computation(BACKWARD), interval(...):
        ze = ze[0, 0, 1] - dz1  # dz < 0

    # -----------------------------------------------------------------------
    # evaporation and accretion of rain for the first 1 / 2 time step
    # -----------------------------------------------------------------------
    with computation(PARALLEL), interval(...):
        (t1, qv1, qr1, ql1, qi1, qs1, qg1, qa1, revap,) = revap_racc(
            t1,
            qv1,
            ql1,
            qr1,
            qi1,
            qs1,
            qg1,
            qa1,
            den1,
            denfac,
            rh_limited,
            table1,
            table2,
            table3,
            table4,
            des1,
            des2,
            des3,
            des4,
            c_air,
            c_vap,
            cracw,
            crevp_0,
            crevp_1,
            crevp_2,
            crevp_3,
            crevp_4,
            d0_vap,
            lv00,
            tau_revp,
            dts,
        )

        evap1 = revap

    with computation(PARALLEL), interval(...):
        if do_sedi_w == True:  # noqa
            dm = dp1 * (1.0 + qv1 + ql1 + qr1 + qi1 + qs1 + qg1)


def warm_rain_step_2(
    t1: FloatField,
    qv1: FloatField,
    ql1: FloatField,
    qr1: FloatField,
    qi1: FloatField,
    qs1: FloatField,
    qg1: FloatField,
    qa1: FloatField,
    den1: FloatField,
    denfac: FloatField,
    vtr: FloatField,
    evap1: FloatField,
    m1_rain: FloatField,
    w1: FloatField,
    rh_limited: FloatField,
    dm: FloatField,
    table1: GlobalTable_driver_qsat,
    table2: GlobalTable_driver_qsat,
    table3: GlobalTable_driver_qsat,
    table4: GlobalTable_driver_qsat,
    des1: GlobalTable_driver_qsat,
    des2: GlobalTable_driver_qsat,
    des3: GlobalTable_driver_qsat,
    des4: GlobalTable_driver_qsat,
):
    """
    Warm rain cloud microphysics: evaporation, accretion

    first half of the timestep
    """
    from __externals__ import (
        c_air,
        c_vap,
        cracw,
        crevp_0,
        crevp_1,
        crevp_2,
        crevp_3,
        crevp_4,
        d0_vap,
        do_sedi_w,
        dts,
        lv00,
        tau_revp,
    )

    with computation(PARALLEL), interval(...):
        # -----------------------------------------------------------------------
        # terminal speed of rain
        # -----------------------------------------------------------------------
        m1_rain = 0.0

    # -----------------------------------------------------------------------
    # vertical velocity transportation during sedimentation
    # -----------------------------------------------------------------------
    with computation(FORWARD), interval(0, 1):
        if do_sedi_w == True:  # noqa
            w1 = (dm * w1 + m1_rain * vtr) / (dm - m1_rain)

    with computation(FORWARD), interval(1, None):
        if do_sedi_w == True:  # noqa
            w1 = (dm * w1 - m1_rain[0, 0, -1] * vtr[0, 0, -1] + m1_rain * vtr) / (
                dm + m1_rain[0, 0, -1] - m1_rain
            )

    # -----------------------------------------------------------------------
    # evaporation and accretion of rain for the remaing 1 / 2 time step
    # -----------------------------------------------------------------------
    with computation(PARALLEL), interval(...):
        (t1, qv1, qr1, ql1, qi1, qs1, qg1, qa1, revap,) = revap_racc(
            t1,
            qv1,
            ql1,
            qr1,
            qi1,
            qs1,
            qg1,
            qa1,
            den1,
            denfac,
            rh_limited,
            table1,
            table2,
            table3,
            table4,
            des1,
            des2,
            des3,
            des4,
            c_air,
            c_vap,
            cracw,
            crevp_0,
            crevp_1,
            crevp_2,
            crevp_3,
            crevp_4,
            d0_vap,
            lv00,
            tau_revp,
            dts,
        )

        evap1 = evap1 + revap


def update_outputs(
    m1_rain: FloatField,
    m1_sol: FloatField,
    rain1: FloatFieldIJ,
    evap1: FloatField,
    revap: FloatField,
    m2_rain: FloatField,
    m2_sol: FloatField,
    m1: FloatField,
    rain: FloatFieldIJ,
):
    """
    Update precipitation totals with results of warm_rain stencil

    reference Fortran: gfdl_cloud_microphys.F90: subroutine mpdrv
    """
    with computation(PARALLEL), interval(...):
        revap = revap + evap1
        m2_rain = m2_rain + m1_rain
        m2_sol = m2_sol + m1_sol
        m1 = m1 + m1_rain + m1_sol

        evap1 = 0
        m1_rain = 0
        m1_sol = 0

    with computation(FORWARD), interval(0, 1):
        rain = rain + rain1

        rain1 = 0
