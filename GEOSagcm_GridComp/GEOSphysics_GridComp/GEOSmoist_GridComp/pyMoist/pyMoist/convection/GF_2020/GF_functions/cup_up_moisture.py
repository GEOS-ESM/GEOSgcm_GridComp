import copy
from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField, IntField, Float, Int, FloatFieldIJ

import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
import pyMoist.constants as constants
from gt4py.cartesian.gtscript import (
    BACKWARD,
    PARALLEL,
    FORWARD,
    K,
    computation,
    interval,
    function,
    int32,
    stencil,
    float32,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.shared_functions import (
    get_cloud_boundary_conditions,
)

# Parameters needed for cup_up_moisture
bdispm = Float(0.366)  # berry--size dispersion (maritime)
bdispc = Float(0.146)  # berry--size dispersion (continental)
T_BF = Float(268.16)
T_ice_BF = Float(235.16)
rk = Int(3)
xexp = Float(2.0)
frac_ave_layer_ocean = Float(0.3)
cum_ave_layer_deep = Float(50.0)
cum_ave_layer_shal = Float(30.0)
cum_ave_layer_mid = Float(50.0)
AUTOCONV = Int(1)
ZERO_DIFF = 0


# @function
# def get_cloud_bc(
#     cumulus: Int,
#     xland: Float,
#     po: FloatField,
#     array: FloatField,
#     k22: Int,
#     # add: Float, Optional input
#     # Tpert: Float, Optional input
# ):

#     if bc_meth == 0:

#         order_aver = 3
#         local_order_aver = min(k22, order_aver)

#         x_aver = 0.0
#         idx = 0
#         while idx <= local_order_aver - 1:
#             x_aver = x_aver + array.at(K=k22 - idx)
#             idx += 1

#         x_aver = x_aver / float(local_order_aver - 1)

#     # elif bc_meth == 1:
#     #     effec_frac = (1.0 - xland) + xland * frac_ave_layer_ocean
#     #     if cumulus == cumulus_parameterization_constants.deep:
#     #         x_ave_layer = cum_ave_layer_deep * effec_frac
#     #     if cumulus == cumulus_parameterization_constants.shallow:
#     #         x_ave_layer = cum_ave_layer_shal * effec_frac
#     #     if cumulus == cumulus_parameterization_constants.mid:
#     #         x_ave_layer = cum_ave_layer_mid * effec_frac

#     #     # i_beg = minloc(abs(po(kts:ktf)-(po(k22)+0.5*x_ave_layer)),1)
#     #     # i_end = minloc(abs(po(kts:ktf)-(po(k22)-0.5*x_ave_layer)),1)
#     #     i_beg = min(70, max(i_beg, 0))
#     #     i_end = min(70, max(i_end, 0))

#     #     if i_beg >= i_end:
#     #         x_aver = array.at(K=k22 - 1)
#     #         dp_layer = 0.0
#     #         ic = i_beg

#     # else:
#     #     dp_layer = 1.0e-06
#     #     x_aver = 0.0
#     #     ic = 0
#     #     if K >= i_beg and K <= 70:
#     #         dp = -(po[0, 0, 1] - po)
#     #         if dp_layer + dp <= x_ave_layer:
#     #             dp_layer = dp_layer + dp
#     #             x_aver = x_aver + array * dp

#     #         else:
#     #             dp = x_ave_layer - dp_layer
#     #             dp_layer = dp_layer + dp
#     #             x_aver = x_aver + array * dp

#     #             # exit ?

#     #     x_aver = x_aver / dp_layer
#     # ic  = max(i_beg,i)

#     # if(present(Tpert)) x_aver = x_aver + cp*maxval(Tpert(i_beg:ic))

#     # IF(present(add)) x_aver = x_aver + add

#     return x_aver


def cup_up_moisture(
    # In
    c1d: FloatField,
    ccn: FloatField,
    cd: FloatField,
    cnvfrc: FloatField,
    name: IntField,
    dby: FloatField,
    entr_rate: FloatField,
    gamma_cup: FloatField,
    hc: FloatField,
    k22: IntField,
    kbcon: IntField,
    klcl: IntField,
    ktop: IntField,
    p_cup: IntField,
    po: FloatField,
    q: FloatField,
    qe_cup: FloatField,
    qes_cup: FloatField,
    rho: FloatField,
    srftype: FloatField,
    start_level: IntField,
    t_cup: FloatField,
    up_massdetr: FloatField,
    up_massentr: FloatField,
    use_linear_subcl_mf: IntField,
    vvel1d: FloatField,
    vvel2d: FloatField,
    x_add_buoy: FloatField,
    xland: FloatField,
    z_cup: FloatField,
    zqexec: FloatField,
    zu: FloatField,
    zws: FloatField,
    c0_shal: FloatField,
    c0_mid: FloatField,
    c0_deep: FloatField,
    qrc_crit_lnd: FloatField,
    qrc_crit_ocn: FloatField,
    bc_meth: IntField,
    perturbation_field: FloatField,
    # Out
    clw_all: FloatField,
    ierr: IntField,
    psum: FloatField,
    psumh: FloatField,
    pw: FloatField,
    pwav: FloatField,
    qc: FloatField,
    qrc: FloatField,
    tempc: FloatField,
):
    from __externals__ import k_end

    with computation(PARALLEL), interval(...):
        # No precip for small clouds
        if name == cumulus_parameterization_constants.SHALLOW:
            c0 = c0_shal
            AVERAGE_LAYER_DEPTH = cum_ave_layer_shal
        if name == cumulus_parameterization_constants.MID:
            c0 = c0_mid
            AVERAGE_LAYER_DEPTH = cum_ave_layer_mid
        if name == cumulus_parameterization_constants.DEEP:
            c0 = c0_deep
            AVERAGE_LAYER_DEPTH = cum_ave_layer_deep

        pwav = 0.0
        psum = 0.0
        psumh = 0.0

    with computation(PARALLEL), interval(0, -1):
        pw = 0.0
        clw_all = 0.0
        tempc = t_cup
        qrc = 0.0
        qc = qe_cup
    with computation(FORWARD), interval(0, 1):
        if ierr == 0:
            qaver: FloatFieldIJ = get_cloud_boundary_conditions(
                qe_cup,
                0,
                po,
                k22 - 1,
                xland,
                bc_meth,
                AVERAGE_LAYER_DEPTH,
                k_end,
                False,
                perturbation_field,
            )

    with computation(PARALLEL), interval(0, -1):
        if ierr == 0:
            if K <= start_level - 1:
                qc = qaver + zqexec + 0.5 * x_add_buoy / cumulus_parameterization_constants.XLV
                qrc = 0.0

    # with computation(PARALLEL), interval(...):
    #    if name == cumulus_parameterization_constants.shallow and use_linear_subcl_mf == 1:
    #         if ierr == 0:
    #             call get_delmix(name,kts,kte,ktf,xland(i),start_level(i),po(i,kts:kte) &
    #                            ,qe_cup(i,kts:kte), qc(i,kts:kte))

    # with computation(PARALLEL), interval(...):
    #     if ierr == 0:
    #         if K >= start_level and K <= ktop:
    #             dz = z_cup - z_cup[0, 0, -1]

    #             qrch = (
    #                 qes_cup
    #                 + (1.0 / cumulus_parameterization_constants.XLV)
    #                 * (gamma_cup / (1.0 + gamma_cup))
    #                 * dby
    #             )

    #             denom = (
    #                 zu[0, 0, -1] - 0.5 * up_massdetr[0, 0, -1] + up_massentr[0, 0, -1]
    #             )
    #             if denom > 0.0:

    #                 qc = (
    #                     qc[0, 0, -1] * zu[0, 0, -1]
    #                     - 0.5 * up_massdetr[0, 0, -1] * qc[0, 0, -1]
    #                     + up_massentr[0, 0, -1] * q[0, 0, -1]
    #                 ) / denom

    #                 if K == start_level:
    #                     qc = qc + zqexec * up_massentr[0, 0, -1] / denom

    #                 qrc = (
    #                     qrc[0, 0, -1] * zu[0, 0, -1]
    #                     - 0.5 * up_massdetr[0, 0, -1] * qrc[0, 0, -1]
    #                 ) / denom

    #             else:
    #                 qc = qc[0, 0, -1]
    #                 qrc = qrc[0, 0, -1]

    #             tempc = (1.0 / cumulus_parameterization_constants.CP) * (
    #                 hc
    #                 - constants.MAPL_GRAV * z_cup
    #                 - cumulus_parameterization_constants.XLV * qrch
    #             )

    #             clw_all = max(0.0, qc - qrch)

    #             qrc = min(clw_all, qrc)

    #             cup = max(0.0, qc - qrch - qrc) / dz

    #             if c0 < 1.0e-6:
    #                 qrc = clw_all
    #                 qc = qrc + min(qc, qrch)
    #                 pwav = 0.0
    #                 psum = psum + clw_all * zu * dz

    # if AUTOCONV == 1:
    #     min_liq = xland * qrc_crit_ocn + (1.0 - xland) * qrc_crit_lnd
    #     cx0 = (c1d + c0) * dz
    #     qrc = clw_all / (1.0 + cx0)
    #     pw = cx0 * max(0.0, qrc - min_liq)

    #     pw = pw * zu

    # qc = qrc + min(qc, qrch)

    # pwav = pwav + pw
    # psum = psum + clw_all * zu * dz

    # with computation(PARALLEL), interval(...):
    #     if ZERO_DIFF == 0:
    #         if pwav < 0.0:
    #             ierr = 66
    #             # ierrc="pwav negative"

    #     if ierr == 0:
    #         if K <= ktop:
    #             qc = qc - qrc


class CupUpMoisture:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
    ) -> None:

        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory

        self._cup_up_moisture = self.stencil_factory.from_dims_halo(
            func=cup_up_moisture,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self.perturbation_field = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a"
        )

    def __call__(
        self,
        # In
        c1d: FloatField,
        ccn: FloatField,
        cd: FloatField,
        cnvfrc: FloatField,
        name: IntField,
        dby: FloatField,
        entr_rate: FloatField,
        gamma_cup: FloatField,
        hc: FloatField,
        k22: IntField,
        kbcon: IntField,
        klcl: IntField,
        ktop: IntField,
        p_cup: IntField,
        po: FloatField,
        q: FloatField,
        qe_cup: FloatField,
        qes_cup: FloatField,
        rho: FloatField,
        srftype: FloatField,
        start_level: IntField,
        t_cup: FloatField,
        up_massdetr: FloatField,
        up_massentr: FloatField,
        use_linear_subcl_mf: IntField,
        vvel1d: FloatField,
        vvel2d: FloatField,
        x_add_buoy: FloatField,
        xland: FloatField,
        z_cup: FloatField,
        zqexec: FloatField,
        zu: FloatField,
        zws: FloatField,
        c0_shal: FloatField,
        c0_mid: FloatField,
        c0_deep: FloatField,
        qrc_crit_lnd: FloatField,
        qrc_crit_ocn: FloatField,
        bc_meth: IntField,
        # Out
        clw_all: FloatField,
        ierr: IntField,
        psum: FloatField,
        psumh: FloatField,
        pw: FloatField,
        pwav: FloatField,
        qc: FloatField,
        qrc: FloatField,
        tempc: FloatField,
    ):
        # Raise an error if bc_meth is not zero
        # if bc_meth.view[:].all() == Int(0):
        #     raise NotImplementedError(
        #         f"Warning: This code has not been ported!! Expected bc_meth == 0"
        #     )

        # Raise error if cumulus is not 'deep'
        # if name.view[:].all() == Int(1):
        #     raise NotImplementedError(
        #         f"Warning: This code has not been ported!! Expected cumulus == 1"
        #     )

        self._cup_up_moisture(
            # In
            c1d=c1d,
            ccn=ccn,
            cd=cd,
            cnvfrc=cnvfrc,
            name=name,
            dby=dby,
            entr_rate=entr_rate,
            gamma_cup=gamma_cup,
            hc=hc,
            k22=k22,
            kbcon=kbcon,
            klcl=klcl,
            ktop=ktop,
            p_cup=p_cup,
            po=po,
            q=q,
            qe_cup=qe_cup,
            qes_cup=qes_cup,
            rho=rho,
            srftype=srftype,
            start_level=start_level,
            t_cup=t_cup,
            up_massdetr=up_massdetr,
            up_massentr=up_massentr,
            use_linear_subcl_mf=use_linear_subcl_mf,
            vvel1d=vvel1d,
            vvel2d=vvel2d,
            x_add_buoy=x_add_buoy,
            xland=xland,
            z_cup=z_cup,
            zqexec=zqexec,
            zu=zu,
            zws=zws,
            c0_shal=c0_shal,
            c0_mid=c0_mid,
            c0_deep=c0_deep,
            qrc_crit_lnd=qrc_crit_lnd,
            qrc_crit_ocn=qrc_crit_ocn,
            bc_meth=bc_meth,
            perturbation_field=self.perturbation_field,
            # Out
            clw_all=clw_all,
            ierr=ierr,
            psum=psum,
            psumh=psumh,
            pw=pw,
            pwav=pwav,
            qc=qc,
            qrc=qrc,
            tempc=tempc,
        )
