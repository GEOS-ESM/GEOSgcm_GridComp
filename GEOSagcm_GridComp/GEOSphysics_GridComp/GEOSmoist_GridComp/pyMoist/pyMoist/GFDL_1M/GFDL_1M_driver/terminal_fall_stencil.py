import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import (
    PARALLEL,
    FORWARD,
    BACKWARD,
    computation,
    interval,
    exp,
    max,
    i32,
    f32,
)
from ndsl.dsl.typing import Float, FloatFieldIJ, FloatField, IntField, BoolField
import pyMoist.GFDL_1M.GFDL_1M_driver.GFDL_1M_driver_constants as driver_constants


@gtscript.function
def terminal_fall_speed_component(
    t: Float,
    qv: Float,
    ql: Float,
    qr: Float,
    qg: Float,
    qs: Float,
    qi: Float,
    m1_sol: Float,
    is_frozen: bool,
):
    from __externals__ import tau_imlt, ql_mlt

    """
    Component of the reference Fortran: gfdl_cloud_microphys.F90: subroutine terminal_fall
    """
    from __externals__ import c_air, c_vap, d0_vap, lv00, dts

    fac_imlt = 1.0 - exp(-dts / tau_imlt)

    # -----------------------------------------------------------------------
    # define heat capacity and latent heat coefficient
    # -----------------------------------------------------------------------

    m1_sol = 0.0
    lhl = lv00 + d0_vap * t
    lhi = driver_constants.li00 + driver_constants.dc_ice * t
    q_liq = ql + qr
    q_sol = qi + qs + qg
    cvm = (
        c_air
        + qv * c_vap
        + q_liq * driver_constants.c_liq
        + q_sol * driver_constants.c_ice
    )
    lcpk = lhl / cvm
    icpk = lhi / cvm

    # -----------------------------------------------------------------------
    # melting of cloud_ice (before fall) :
    # -----------------------------------------------------------------------

    tc = t - driver_constants.tice
    if is_frozen == False and (qi > driver_constants.qcmin and tc > 0.0):
        sink = min(qi, fac_imlt * tc / icpk)
        if ql_mlt - ql > 0:
            ans = ql_mlt - ql
        else:
            ans = 0
        tmp = min(sink, ans)
        ql = ql + tmp
        qr = qr + sink - tmp
        qi = qi - sink
        q_liq = q_liq + sink
        q_sol = q_sol - sink
        cvm = (
            c_air
            + qv * c_vap
            + q_liq * driver_constants.c_liq
            + q_sol * driver_constants.c_ice
        )
        t = t - sink * lhi / cvm

    return t, ql, qr, qi, cvm, lhi, icpk, m1_sol


def terminal_fall_stencil(
    t1: FloatField,
    qv1: FloatField,
    ql1: FloatField,
    qr1: FloatField,
    qg1: FloatField,
    qs1: FloatField,
    qi1: FloatField,
    dz1: FloatField,
    dp1: FloatField,
    den1: FloatField,
    vtg: FloatField,
    vts: FloatField,
    vti: FloatField,
    precip_rain: FloatFieldIJ,
    precip_graupel: FloatFieldIJ,
    precip_snow: FloatFieldIJ,
    precip_ice: FloatFieldIJ,
    m1_sol: FloatField,
    w1: FloatField,
    ze: FloatField,
    zt: FloatField,
    is_frozen: BoolField,
    precip_fall: FloatFieldIJ,
    melting_mask_1: BoolField,
    melting_mask_2: BoolField,
    current_k_level: IntField,
    TESTVAR: FloatField,
):
    from __externals__ import (
        dts,
        tau_imlt,
        ql_mlt,
        vi_fac,
        do_sedi_w,
        use_ppm,
        k_end,
        tau_smlt,
        tau_g2r,
    )

    # begin reference Fortran: gfdl_cloud_microphys.F90: subroutine terminal_fall

    # determine frozen levels
    # later operations will only be executed if frozen/melted
    # initalized to is_frozen = False, True = frozen, False = melted
    with computation(PARALLEL), interval(...):
        if t1 <= driver_constants.tice:
            is_frozen = True

    # we only want the melting layer closest to the surface
    with computation(BACKWARD), interval(0, -1):
        if is_frozen[0, 0, 1] == True and is_frozen[0, 0, 0] == False:
            is_frozen = True

    # force surface to "melt" for later calculations
    with computation(PARALLEL), interval(-1, None):
        is_frozen = False

    with computation(PARALLEL), interval(...):
        t1, ql1, qr1, qi1, cvm, lhi, icpk, m1_sol = terminal_fall_speed_component(
            t1,
            qv1,
            ql1,
            qr1,
            qg1,
            qs1,
            qi1,
            m1_sol,
            is_frozen,
        )

    # if timestep is too small turn off  melting
    with computation(PARALLEL), interval(0, -1):
        if dts < 300.0:
            is_frozen = True
            disable_melt = True
        else:
            disable_melt = False

    with computation(PARALLEL), interval(-1, None):
        if dts < 300.0:
            is_frozen = False
            disable_melt = True
        else:
            disable_melt = False

    with computation(BACKWARD), interval(...):
        ze = ze[0, 0, 1] - dz1  # dz < 0

    with computation(FORWARD), interval(0, 1):
        zt = ze

    # -----------------------------------------------------------------------
    # update capacity heat and latend heat coefficient
    # -----------------------------------------------------------------------

    with computation(PARALLEL), interval(...):
        if is_frozen == False:
            lhi = driver_constants.li00 + driver_constants.dc_ice * t1
            icpk = lhi / cvm

    with computation(FORWARD), interval(-1, None):
        lhi = driver_constants.li00 + driver_constants.dc_ice * t1
        icpk = lhi / cvm

    # -----------------------------------------------------------------------
    # -----------------------------------------------------------------------
    # melting of falling cloud ice into rain
    # -----------------------------------------------------------------------
    # -----------------------------------------------------------------------

    # reference Fortran: gfdl_cloud_microphys.F90: subroutine check_column
    # determine if any precip falls in the column
    # if it falls anywhere in the column, the entire column becomes true
    # initalized to 0 (false), potentially changed to 1 (true)
    with computation(FORWARD), interval(...):
        if qi1 > driver_constants.qpmin:
            precip_fall = 1
    # end reference Fortran: gfdl_cloud_microphys.F90: subroutine check_column

    with computation(FORWARD), interval(0, 1):
        if vi_fac < 1.0e-5 or precip_fall == 0:
            precip_ice = 0

    with computation(FORWARD), interval(1, None):
        if vi_fac >= 1.0e-5 and precip_fall == 1:
            zt = ze - dts * (vti[0, 0, -1] + vti) / 2.0

    with computation(FORWARD), interval(-1, None):
        if vi_fac >= 1.0e-5 and precip_fall == 1:
            zt[0, 0, 1] = driver_constants.zs - dts * vti

    with computation(FORWARD), interval(...):
        if vi_fac >= 1.0e-5 and precip_fall == 1:
            if zt[0, 0, 1] >= zt:
                zt[0, 0, 1] = zt - driver_constants.dz_min

    with computation(PARALLEL), interval(...):
        if vi_fac >= 1.0e-5 and precip_fall == 1 and disable_melt == False:
            # copy frozen mask to make modifyable version
            melting_mask_1 = is_frozen
            # flip logic for clarity
            if melting_mask_1 == True:
                melting_mask_1 = False
            else:
                melting_mask_1 = True

    with computation(BACKWARD), interval(-1, None):
        if vi_fac >= 1.0e-5 and precip_fall == 1 and disable_melt == False:
            # ensure no operations are performed on the surface
            melting_mask_1 = False

    with computation(BACKWARD), interval(...):
        # THIS HAS NEVER BEEN TESTED B/C DTS WAS LESS THAN 300 IN THE TEST CASE
        if vi_fac >= 1.0e-5 and precip_fall == 1 and disable_melt == False:
            # only operate on melted layers
            if melting_mask_1 == True and qi1 > driver_constants.qcmin:
                # initalize exit trigger
                stop_melting = False
                m = 0
                while m < k_end and stop_melting == False:
                    mplus1 = m + 1  # TODO remove this
                    # only opterate on previously iterated k-levels
                    if melting_mask_2 == True:
                        if zt[0, 0, 1] >= ze.at(K=m):
                            stop_melting = (
                                True  # if true exit early for ONLY this k level
                            )
                        if stop_melting == False:
                            dtime = min(
                                1.0,
                                (ze.at(K=m) - ze.at(K=mplus1))
                                / (max(driver_constants.vr_min, vti) * tau_imlt),
                            )
                            sink = min(
                                qi1 * dp1 / dp1.at(K=m),
                                dtime
                                * (t1.at(K=m) - driver_constants.tice)
                                / icpk.at(K=m),
                            )
                        if ql_mlt - ql1.at(K=m) > 0:
                            ans = ql_mlt - ql1.at(K=m)
                        else:
                            ans = 0
                        tmp = min(sink, ans)
                        offset = i32(m - current_k_level)
                        hold_ql1 = ql1.at(K=m)
                        ql1[0, 0, offset] = hold_ql1 + tmp
                        hold_qr1 = qr1.at(K=m)
                        qr1[0, 0, offset] = hold_qr1 - tmp + sink
                        hold_t1 = t1.at(K=m)
                        t1[0, 0, offset] = hold_t1 - sink * icpk.at(K=m)
                        qi1 = qi1 - sink * dp1.at(K=m) / dp1
                    m = m + 1
            # set current layer as iterated layer
            melting_mask_2 = True

    with computation(PARALLEL), interval(...):
        if vi_fac >= 1.0e-5 and precip_fall == 1:
            if do_sedi_w == True:
                dm = dp1 * (1.0 + qv1 + ql1 + qr1 + qi1 + qs1 + qg1)

    # begin reference Fortran: gfdl_cloud_microphys.F90: subroutine implicit_fall
    # this code computes the time-implicit monotonic scheme
    # Fortran author: Shian-Jiann Lin, 2016
    with computation(PARALLEL), interval(...):
        if vi_fac >= 1.0e-5 and precip_fall == 1:
            if use_ppm == False:
                hold_data = ze - ze[0, 0, 1]
                dd = dts * vti
                qi1 = qi1 * dp1

    # -----------------------------------------------------------------------
    # sedimentation: non - vectorizable loop
    # -----------------------------------------------------------------------
    with computation(FORWARD), interval(0, 1):
        if vi_fac >= 1.0e-5 and precip_fall == 1:
            if use_ppm == False:
                qm = qi1 / (hold_data + dd)

    with computation(FORWARD), interval(1, None):
        if vi_fac >= 1.0e-5 and precip_fall == 1:
            if use_ppm == False:
                qm = (qi1 + dd[0, 0, -1] * qm[0, 0, -1]) / (hold_data + dd)

    # -----------------------------------------------------------------------
    # qm is density at this stage
    # -----------------------------------------------------------------------
    with computation(PARALLEL), interval(...):
        if vi_fac >= 1.0e-5 and precip_fall == 1:
            if use_ppm == False:
                qm = qm * hold_data

    # -----------------------------------------------------------------------
    # output mass fluxes: non - vectorizable loop
    # -----------------------------------------------------------------------
    with computation(FORWARD), interval(0, 1):
        if vi_fac >= 1.0e-5 and precip_fall == 1:
            if use_ppm == False:
                m1_sol = qi1 - qm

    with computation(FORWARD), interval(1, None):
        if vi_fac >= 1.0e-5 and precip_fall == 1:
            if use_ppm == False:
                m1_sol = m1_sol[0, 0, -1] + qi1 - qm

    with computation(FORWARD), interval(-1, None):
        if vi_fac >= 1.0e-5 and precip_fall == 1:
            if use_ppm == False:
                precip_ice = m1_sol

    # -----------------------------------------------------------------------
    # update:
    # -----------------------------------------------------------------------
    with computation(PARALLEL), interval(...):
        if vi_fac >= 1.0e-5 and precip_fall == 1:
            if use_ppm == False:
                qi1 = qm / dp1
    # end reference Fortran: gfdl_cloud_microphys.F90: subroutine implicit_fall

    with computation(FORWARD), interval(0, 1):
        if vi_fac >= 1.0e-5 and precip_fall == 1:
            if do_sedi_w:
                w1 = (dm * w1 + m1_sol * vti) / (dm - m1_sol)

    with computation(FORWARD), interval(1, None):
        if vi_fac >= 1.0e-5 and precip_fall == 1:
            if do_sedi_w:
                w1 = (dm * w1 - m1_sol[0, 0, -1] * vti[0, 0, -1] + m1_sol * vti) / (
                    dm + m1_sol[0, 0, -1] - m1_sol
                )

    # reset masks
    with computation(FORWARD), interval(0, 1):
        precip_fall = 0

    with computation(PARALLEL), interval(...):
        melting_mask_1 = False
        melting_mask_2 = False

    # -----------------------------------------------------------------------
    # -----------------------------------------------------------------------
    # melting of falling snow into rain
    # -----------------------------------------------------------------------
    # -----------------------------------------------------------------------

    # reference Fortran: gfdl_cloud_microphys.F90: subroutine check_column
    # determine if any precip falls in the column
    # if it falls anywhere in the column, the entire column becomes true
    # initalized to 0 (false), potentially changed to 1 (true)
    with computation(FORWARD), interval(...):
        if qs1 > driver_constants.qpmin:
            precip_fall = 1
    # end reference Fortran: gfdl_cloud_microphys.F90: subroutine check_column

    with computation(FORWARD), interval(0, 1):
        if precip_fall == 0:
            precip_snow = 0

    with computation(FORWARD), interval(1, None):
        if precip_fall == 1:
            zt = ze - dts * (vts[0, 0, -1] + vts) / 2.0

    with computation(FORWARD), interval(-1, None):
        if precip_fall == 1:
            zt[0, 0, 1] = driver_constants.zs - dts * vts

    with computation(FORWARD), interval(...):
        if precip_fall == 1:
            if zt[0, 0, 1] >= zt:
                zt[0, 0, 1] = zt - driver_constants.dz_min

    with computation(PARALLEL), interval(...):
        if precip_fall == 1 and disable_melt == False:
            # copy frozen mask to make modifyable version
            melting_mask_1 = is_frozen
            # flip logic for clarity
            if melting_mask_1 == True:
                melting_mask_1 = False
            else:
                melting_mask_1 = True

    with computation(BACKWARD), interval(-1, None):
        if precip_fall == 1 and disable_melt == False:
            # ensure no operations are performed on the surface
            melting_mask_1 = False

    with computation(BACKWARD), interval(...):
        # THIS HAS NEVER BEEN TESTED B/C DTS WAS LESS THAN 300 IN THE TEST CASE
        if precip_fall == 1 and disable_melt == False:
            # only operate on melted layers
            if melting_mask_1 == True and qs1 > driver_constants.qpmin:
                # initalize exit trigger
                stop_melting = False
                m = 0
                while m < k_end and stop_melting == False:
                    mplus1 = m + 1  # TODO remove this
                    # only opterate on previously iterated k-levels
                    if melting_mask_2 == True:
                        if zt[0, 0, 1] >= ze.at(K=m):
                            stop_melting = (
                                True  # if true exit early for ONLY this k level
                            )
                        if stop_melting == False:
                            dtime = min(
                                dts,
                                (ze.at(K=m) - ze.at(K=mplus1))
                                / (driver_constants.vr_min + vts),
                            )
                            if (
                                zt < ze.at(K=mplus1)
                                and t1.at(K=m) > driver_constants.tice
                            ):
                                dtime = min(1.0, dtime / tau_smlt)
                                sink = min(
                                    qs1 * dp1 / dp1.at(K=m),
                                    dtime
                                    * (t1.at(K=m) - driver_constants.tice)
                                    / icpk.at(K=m),
                                )
                                offset = i32(m - current_k_level)
                                hold_t1 = t1.at(K=m)
                                t1[0, 0, offset] = hold_t1 - sink * icpk.at(K=m)
                                qs1 = qs1 - sink * dp1.at(K=m) / dp1
                                if zt < 0:
                                    precip_rain = precip_rain + sink * dp1.at(
                                        K=m
                                    )  # precip as rain
                                else:
                                    # qr source here will fall next time step (therefore, can evap)
                                    hold_qr1 = qr1.at(K=m)
                                    qr1[0, 0, offset] = hold_qr1 + sink
                        if stop_melting == False:
                            if qs1 < driver_constants.qpmin:
                                stop_melting = True
                    m = m + 1
            # set current layer as iterated layer
            melting_mask_2 = True

    with computation(PARALLEL), interval(...):
        if precip_fall == 1:
            if do_sedi_w == True:
                dm = dp1 * (1.0 + qv1 + ql1 + qr1 + qi1 + qs1 + qg1)

    # begin reference Fortran: gfdl_cloud_microphys.F90: subroutine implicit_fall
    # this code computes the time-implicit monotonic scheme
    # Fortran author: Shian-Jiann Lin, 2016
    with computation(PARALLEL), interval(...):
        if precip_fall == 1:
            if use_ppm == False:
                hold_data = ze - ze[0, 0, 1]
                dd = dts * vts
                qs1 = qs1 * dp1

    # -----------------------------------------------------------------------
    # sedimentation: non - vectorizable loop
    # -----------------------------------------------------------------------
    with computation(FORWARD), interval(0, 1):
        if precip_fall == 1:
            if use_ppm == False:
                qm = qs1 / (hold_data + dd)

    with computation(FORWARD), interval(1, None):
        if precip_fall == 1:
            if use_ppm == False:
                qm = (qs1 + dd[0, 0, -1] * qm[0, 0, -1]) / (hold_data + dd)

    # -----------------------------------------------------------------------
    # qm is density at this stage
    # -----------------------------------------------------------------------
    with computation(PARALLEL), interval(...):
        if precip_fall == 1:
            if use_ppm == False:
                qm = qm * hold_data

    # -----------------------------------------------------------------------
    # output mass fluxes: non - vectorizable loop
    # -----------------------------------------------------------------------
    with computation(FORWARD), interval(0, 1):
        if precip_fall == 1:
            if use_ppm == False:
                m1 = qs1 - qm

    with computation(FORWARD), interval(1, None):
        if precip_fall == 1:
            if use_ppm == False:
                m1 = m1[0, 0, -1] + qs1 - qm

    with computation(FORWARD), interval(-1, None):
        if precip_fall == 1:
            if use_ppm == False:
                precip_snow = m1

    # -----------------------------------------------------------------------
    # update:
    # -----------------------------------------------------------------------
    with computation(PARALLEL), interval(...):
        if precip_fall == 1:
            if use_ppm == False:
                qs1 = qm / dp1
    # end reference Fortran: gfdl_cloud_microphys.F90: subroutine implicit_fall

    with computation(FORWARD), interval(0, 1):
        if precip_fall == 1:
            if do_sedi_w:
                w1 = (dm * w1 + m1 * vts) / (dm - m1)

    with computation(FORWARD), interval(1, None):
        if precip_fall == 1:
            if do_sedi_w:
                w1 = (dm * w1 - m1[0, 0, -1] * vts[0, 0, -1] + m1 * vts) / (
                    dm + m1[0, 0, -1] - m1
                )

    with computation(PARALLEL), interval(...):
        m1_sol = m1_sol + m1

    # reset masks and temporaries
    with computation(FORWARD), interval(0, 1):
        precip_fall = 0

    with computation(PARALLEL), interval(...):
        melting_mask_1 = False
        melting_mask_2 = False
        # must be reset to zero everywhere in case there is no graupel
        # in a column and no calculations are performed
        m1 = 0

    # -----------------------------------------------------------------------
    # -----------------------------------------------------------------------
    # melting of falling graupel into rain
    # -----------------------------------------------------------------------
    # -----------------------------------------------------------------------

    # reference Fortran: gfdl_cloud_microphys.F90: subroutine check_column
    # determine if any precip falls in the column
    # if it falls anywhere in the column, the entire column becomes true
    # initalized to 0 (false), potentially changed to 1 (true)
    with computation(FORWARD), interval(...):
        if qg1 > driver_constants.qpmin:
            precip_fall = 1
    # end reference Fortran: gfdl_cloud_microphys.F90: subroutine check_column

    with computation(FORWARD), interval(0, 1):
        if precip_fall == 0:
            precip_graupel = 0

    with computation(FORWARD), interval(1, None):
        if precip_fall == 1:
            zt = ze - dts * (vtg[0, 0, -1] + vtg) / 2.0

    with computation(FORWARD), interval(-1, None):
        if precip_fall == 1:
            zt[0, 0, 1] = driver_constants.zs - dts * vtg

    with computation(FORWARD), interval(...):
        if precip_fall == 1:
            if zt[0, 0, 1] >= zt:
                zt[0, 0, 1] = zt - driver_constants.dz_min

    with computation(PARALLEL), interval(...):
        if precip_fall == 1 and disable_melt == False:
            # copy frozen mask to make modifyable version
            melting_mask_1 = is_frozen
            # flip logic for clarity
            if melting_mask_1 == True:
                melting_mask_1 = False
            else:
                melting_mask_1 = True

    with computation(BACKWARD), interval(-1, None):
        if precip_fall == 1 and disable_melt == False:
            # ensure no operations are performed on the surface
            melting_mask_1 = False

    with computation(BACKWARD), interval(...):
        # THIS HAS NEVER BEEN TESTED B/C DTS WAS LESS THAN 300 IN THE TEST CASE
        if precip_fall == 1 and disable_melt == False:
            # only operate on melted layers
            if melting_mask_1 == True and qg1 > driver_constants.qpmin:
                # initalize exit trigger
                stop_melting = False
                m = 0
                while m < k_end and stop_melting == False:
                    mplus1 = m + 1  # TODO remove this
                    # only opterate on previously iterated k-levels
                    if melting_mask_2 == True:
                        if zt[0, 0, 1] >= ze.at(K=m):
                            stop_melting = (
                                True  # if true exit early for ONLY this k level
                            )
                        if stop_melting == False:
                            dtime = min(dts, (ze.at(K=m) - ze.at(K=mplus1)) / vtg)
                            if (
                                zt < ze.at(K=mplus1)
                                and t1.at(K=m) > driver_constants.tice
                            ):
                                dtime = min(1.0, dtime / tau_g2r)
                                sink = min(
                                    qg1 * dp1 / dp1.at(K=m),
                                    dtime
                                    * (t1.at(K=m) - driver_constants.tice)
                                    / icpk.at(K=m),
                                )
                                offset = i32(m - current_k_level)
                                hold_t1 = t1.at(K=m)
                                t1[0, 0, offset] = hold_t1 - sink * icpk.at(K=m)
                                qg1 = qg1 - sink * dp1.at(K=m) / dp1
                                if zt < 0:
                                    precip_rain = precip_rain + sink * dp1.at(
                                        K=m
                                    )  # precip as rain
                                else:
                                    # qr source here will fall next time step (therefore, can evap)
                                    hold_qr1 = qr1.at(K=m)
                                    qr1[0, 0, offset] = hold_qr1 + sink
                        if stop_melting == False:
                            if qg1 < driver_constants.qpmin:
                                stop_melting = True
                    m = m + 1
            # set current layer as iterated layer
            melting_mask_2 = True

    with computation(PARALLEL), interval(...):
        if precip_fall == 1:
            if do_sedi_w == True:
                dm = dp1 * (1.0 + qv1 + ql1 + qr1 + qi1 + qs1 + qg1)

    # begin reference Fortran: gfdl_cloud_microphys.F90: subroutine implicit_fall
    # this code computes the time-implicit monotonic scheme
    # Fortran author: Shian-Jiann Lin, 2016
    with computation(PARALLEL), interval(...):
        if precip_fall == 1:
            if use_ppm == False:
                hold_data = ze - ze[0, 0, 1]
                dd = dts * vtg
                qg1 = qg1 * dp1

    # -----------------------------------------------------------------------
    # sedimentation: non - vectorizable loop
    # -----------------------------------------------------------------------
    with computation(FORWARD), interval(0, 1):
        if precip_fall == 1:
            if use_ppm == False:
                qm = qg1 / (hold_data + dd)

    with computation(FORWARD), interval(1, None):
        if precip_fall == 1:
            if use_ppm == False:
                qm = (qg1 + dd[0, 0, -1] * qm[0, 0, -1]) / (hold_data + dd)

    # -----------------------------------------------------------------------
    # qm is density at this stage
    # -----------------------------------------------------------------------
    with computation(PARALLEL), interval(...):
        if precip_fall == 1:
            if use_ppm == False:
                qm = qm * hold_data

    # -----------------------------------------------------------------------
    # output mass fluxes: non - vectorizable loop
    # -----------------------------------------------------------------------
    with computation(FORWARD), interval(0, 1):
        if precip_fall == 1:
            if use_ppm == False:
                m1 = qg1 - qm

    with computation(FORWARD), interval(1, None):
        if precip_fall == 1:
            if use_ppm == False:
                m1 = m1[0, 0, -1] + qg1 - qm

    with computation(FORWARD), interval(-1, None):
        if precip_fall == 1:
            if use_ppm == False:
                precip_graupel = m1

    # -----------------------------------------------------------------------
    # update:
    # -----------------------------------------------------------------------
    with computation(PARALLEL), interval(...):
        if precip_fall == 1:
            if use_ppm == False:
                qg1 = qm / dp1
    # end reference Fortran: gfdl_cloud_microphys.F90: subroutine implicit_fall

    with computation(FORWARD), interval(0, 1):
        if precip_fall == 1:
            if do_sedi_w:
                w1 = (dm * w1 + m1 * vtg) / (dm - m1)

    with computation(FORWARD), interval(1, None):
        if precip_fall == 1:
            if do_sedi_w:
                w1 = (dm * w1 - m1[0, 0, -1] * vtg[0, 0, -1] + m1 * vtg) / (
                    dm + m1[0, 0, -1] - m1
                )

    with computation(PARALLEL), interval(...):
        m1_sol = m1_sol + m1

    # end reference Fortran: gfdl_cloud_microphys.F90: subroutine terminal_fall
