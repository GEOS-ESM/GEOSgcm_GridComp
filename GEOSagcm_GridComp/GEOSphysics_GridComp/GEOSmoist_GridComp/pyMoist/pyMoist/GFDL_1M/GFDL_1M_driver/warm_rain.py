import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import (
    BACKWARD,
    FORWARD,
    PARALLEL,
    computation,
    exp,
    i32,
    interval,
    log,
    max,
    sqrt,
    trunc,
)

import pyMoist.GFDL_1M.GFDL_1M_driver.GFDL_1M_driver_constants as driver_constants
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ


GlobalTable_driver_qsat = gtscript.GlobalTable[(Float, (int(driver_constants.length)))]


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
    tmin = driver_constants.table_ice - 160.0

    if ta - tmin > 0:
        ans = ta - tmin
    else:
        ans = 0
    ap1 = 10.0 * ans + 1.0
    ap1 = min(2621.0, ap1)
    it = i32(trunc(ap1))
    es = table2.A[it - 1] + (ap1 - it) * des2.A[it - 1]
    qsat = es / (driver_constants.rvgas * ta * den)
    it = i32(
        trunc(ap1 - 0.5)
    )  # check if this rounds or truncates. need truncation here
    # finite diff, del_t = 0.1:
    dqdt = (
        10.0
        * (des2.A[it - 1] + (ap1 - it) * (des2.A[it] - des2.A[it - 1]))
        / (driver_constants.rvgas * ta * den)
    )

    return qsat, dqdt


@gtscript.function
def revap_racc(
    half_dt: Float,
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
):
    """
    evaporate rain

    reference Fortran: gfdl_cloud_microphys.F90: subroutine revap_racc
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
        lv00,
        tau_revp,
    )

    prec_ls = (qr1 + qs1 + qg1) * den1
    area_ls_prc = qa1 * (qr1 + qs1 + qg1) * den1
    revap = 0.0

    if t1 > driver_constants.t_wfr and qr1 > driver_constants.qpmin:
        # area and timescale efficiency on revap
        fac_revp = 1.0 - exp(-half_dt / tau_revp)

        # -----------------------------------------------------------------------
        # define heat capacity and latent heat coefficient
        # -----------------------------------------------------------------------

        lhl = lv00 + d0_vap * t1
        q_liq = ql1 + qr1
        q_sol = qi1 + qs1 + qg1
        cvm = (
            c_air
            + qv1 * c_vap
            + q_liq * driver_constants.c_liq
            + q_sol * driver_constants.c_ice
        )
        lcpk = lhl / cvm

        tin = t1 - lcpk * ql1  # presence of clouds suppresses the rain evap
        qpz = qv1 + ql1
        qsat, dqsdt = wqs2(tin, den1, table2, des2)
        dqh = max(ql1, rh_limited * max(qpz, driver_constants.qcmin))
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

        if dqv > driver_constants.qvmin and qsat > q_minus:
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
            evap = min(qr1, min(half_dt * fac_revp * evap, dqv / (1.0 + lcpk * dqsdt)))
            qr1 = qr1 - evap
            qv1 = qv1 + evap
            q_liq = q_liq - evap
            cvm = (
                c_air
                + qv1 * c_vap
                + q_liq * driver_constants.c_liq
                + q_sol * driver_constants.c_ice
            )
            t1 = t1 - evap * lhl / cvm
            revap = evap / half_dt

        # -----------------------------------------------------------------------
        # accretion: pracc
        # -----------------------------------------------------------------------

        if (
            qr1 > driver_constants.qpmin
            and ql1 > driver_constants.qcmin
            and qsat < q_minus
        ):
            sink = half_dt * denfac * cracw * exp(0.95 * log(qr1 * den1))
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
        prec_ls,
        area_ls_prc,
    )


def warm_rain(
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
    m1_rain: FloatField,
    w1: FloatField,
    rh_limited: FloatField,
    eis: FloatFieldIJ,
    onemsig: FloatFieldIJ,
    precip_rain: FloatFieldIJ,
    ze: FloatField,
    zt: FloatField,
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
    warm rain cloud microphysics: evaporation, accretion

    reference Fortran: gfdl_cloud_microphys.F90: subroutine warm_rain
    """
    from __externals__ import (
        const_vr,
        do_qa,
        do_sedi_w,
        dts,
        irain_f,
        k_end,
        ql0_max,
        rthreshs,
        rthreshu,
        use_ppm,
        vr_fac,
        vr_max,
        z_slope_liq,
    )

    # warm rain cloud microphysics
    with computation(PARALLEL), interval(...):
        half_dt = 0.5 * dts

        # -----------------------------------------------------------------------
        # terminal speed of rain
        # -----------------------------------------------------------------------
        m1_rain = 0.0

    # reference Fortran: gfdl_cloud_microphys.F90: subroutine check_column
    # determine if any precip falls in the column
    # if it falls anywhere in the column, the entire column becomes true
    # initalized to 0 (false), potentially changed to 1 (true)
    with computation(FORWARD), interval(...):
        if qr1 > driver_constants.qpmin:
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
            qadum = max(qa1, driver_constants.qcmin)
        else:
            qadum = 1.0
        ql1 = ql1 / qadum
        qi1 = qi1 / qadum

        fac_rc = (
            min(1.0, eis / 15.0) ** 2
        )  # Estimated inversion strength determine stable regime
        fac_rc = (
            driver_constants.rc * (rthreshs * fac_rc + rthreshu * (1.0 - fac_rc)) ** 3
        )

    with computation(PARALLEL), interval(...):
        if irain_f != 0:
            # -----------------------------------------------------------------------
            # no subgrid varaibility
            # -----------------------------------------------------------------------
            if qadum > onemsig:
                if t1 > driver_constants.t_wfr:
                    qc = fac_rc * ccn / den1
                    dq = ql1 - qc
                    if dq > 0.0:
                        sink = min(
                            dq,
                            dts * c_praut * den1 * exp(driver_constants.so3 * log(ql1)),
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
                dl = max(dl, max(driver_constants.qvmin, rh_limited * ql1))
            if z_slope_liq == False:  # noqa
                dl = max(driver_constants.qvmin, rh_limited * ql1)

            # end reference Fortran: gfdl_cloud_microphys.F90: subroutine linear_prof

            if qadum > onemsig:
                if t1 > driver_constants.t_wfr + driver_constants.dt_fr:
                    dl = min(max(driver_constants.qcmin, dl), 0.5 * ql1)
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
                            * exp(driver_constants.so3 * log(ql1))
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
            vtr = driver_constants.vf_min
        elif const_vr == True:  # noqa
            vtr = vr_fac  # ifs_2016: 4.0
        else:
            qden = qr1 * den1
            if qr1 < driver_constants.thr:
                vtr = driver_constants.vr_min
            else:
                vtr = (
                    vr_fac
                    * driver_constants.vconr
                    * sqrt(min(10.0, driver_constants.sfcrho / den1))
                    * exp(0.2 * log(qden / driver_constants.normr))
                )
                vtr = min(vr_max, max(driver_constants.vr_min, vtr))

    with computation(FORWARD), interval(-1, None):
        ze[0, 0, 1] = driver_constants.zs

    with computation(BACKWARD), interval(...):
        ze = ze[0, 0, 1] - dz1  # dz < 0

    # -----------------------------------------------------------------------
    # evaporation and accretion of rain for the first 1 / 2 time step
    # -----------------------------------------------------------------------
    with computation(PARALLEL), interval(...):
        (
            t1,
            qv1,
            qr1,
            ql1,
            qi1,
            qs1,
            qg1,
            qa1,
            revap,
            prec_ls,
            area_ls_prc,
        ) = revap_racc(
            half_dt,
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
        )

        evap1 = revap

    with computation(FORWARD), interval(1, None):
        prec_ls = prec_ls + prec_ls[0, 0, -1]
        area_ls_prc = area_ls_prc + area_ls_prc[0, 0, -1]

    with computation(FORWARD), interval(0, 1):
        total_prec_ls = prec_ls.at(K=k_end)
        total_area_ls_prc = area_ls_prc.at(K=k_end)

    with computation(PARALLEL), interval(...):
        if do_sedi_w == True:  # noqa
            dm = dp1 * (1.0 + qv1 + ql1 + qr1 + qi1 + qs1 + qg1)

    # -----------------------------------------------------------------------
    # mass flux induced by falling rain
    # -----------------------------------------------------------------------

    with computation(FORWARD), interval(0, 1):
        if precip_fall == 0:
            precip_rain = 0.0

    # begin reference Fortran: gfdl_cloud_microphys.F90: subroutine implicit_fall
    # this code computes the time-implicit monotonic scheme
    # Fortran author: Shian-Jiann Lin, 2016
    # set up inputs to the "function"
    with computation(PARALLEL), interval(...):
        if precip_fall == 1:
            if use_ppm == False:  # noqa
                q_implicit_fall = qr1
                vt_implicit_fall = vtr

    with computation(PARALLEL), interval(...):
        if precip_fall == 1:
            if use_ppm == False:  # noqa
                hold_data = ze - ze[0, 0, 1]
                dd = dts * vt_implicit_fall
                q_implicit_fall = q_implicit_fall * dp1

    # -----------------------------------------------------------------------
    # sedimentation: non - vectorizable loop
    # -----------------------------------------------------------------------
    with computation(FORWARD), interval(0, 1):
        if precip_fall == 1:
            if use_ppm == False:  # noqa
                qm = q_implicit_fall / (hold_data + dd)

    with computation(FORWARD), interval(1, None):
        if precip_fall == 1:
            if use_ppm == False:  # noqa
                qm = (q_implicit_fall + dd[0, 0, -1] * qm[0, 0, -1]) / (hold_data + dd)

    # -----------------------------------------------------------------------
    # qm is density at this stage
    # -----------------------------------------------------------------------
    with computation(PARALLEL), interval(...):
        if precip_fall == 1:
            if use_ppm == False:  # noqa
                qm = qm * hold_data

    # -----------------------------------------------------------------------
    # output mass fluxes: non - vectorizable loop
    # -----------------------------------------------------------------------
    with computation(FORWARD), interval(0, 1):
        if precip_fall == 1:
            if use_ppm == False:  # noqa
                m1 = q_implicit_fall - qm

    with computation(FORWARD), interval(1, None):
        if precip_fall == 1:
            if use_ppm == False:  # noqa
                m1 = m1[0, 0, -1] + q_implicit_fall - qm

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
                q_implicit_fall = qm / dp1

    # update "outputs" after "function"
    with computation(PARALLEL), interval(...):
        if precip_fall == 1:
            if use_ppm == False:  # noqa
                qr1 = q_implicit_fall
                m1_rain = (
                    m1_rain + m1
                )  # NOTE: setting this to just m1_rain = m1 gives WILD values (1e31)

    with computation(FORWARD), interval(-1, None):
        if precip_fall == 1:
            if use_ppm == False:  # noqa
                precip_rain = precip
    # end reference Fortran: gfdl_cloud_microphys.F90: subroutine implicit_fall

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
        (
            t1,
            qv1,
            qr1,
            ql1,
            qi1,
            qs1,
            qg1,
            qa1,
            revap,
            prec_ls,
            area_ls_prc,
        ) = revap_racc(
            half_dt,
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
        )

        evap1 = evap1 + revap
