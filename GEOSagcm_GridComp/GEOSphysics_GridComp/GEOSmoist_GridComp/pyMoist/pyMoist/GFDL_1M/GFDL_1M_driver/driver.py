"""GFDL_1M driver"""

from ndsl import QuantityFactory, StencilFactory, orchestrate
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ
from pyMoist.GFDL_1M.GFDL_1M_driver.config import config
from pyMoist.GFDL_1M.GFDL_1M_driver.connections import (
    fall_speed,
    fix_negative_values,
    icloud_update,
    init_temporaries,
    update_tendencies,
    warm_rain_update,
)
from pyMoist.GFDL_1M.GFDL_1M_driver.icloud import icloud
from pyMoist.GFDL_1M.GFDL_1M_driver.sat_tables import get_tables
from pyMoist.GFDL_1M.GFDL_1M_driver.support import (
    check_flags,
    ConfigConstants,
    Masks,
    Outputs,
    Temporaries,
)
from pyMoist.GFDL_1M.GFDL_1M_driver.terminal_fall.terminal_fall import TerminalFall
from pyMoist.GFDL_1M.GFDL_1M_driver.warm_rain import warm_rain


class driver:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        GFDL_1M_config: config,
    ):

        self.config_dependent_constants = ConfigConstants(GFDL_1M_config)

        # Check values for untested code paths
        check_flags(
            GFDL_1M_config,
            self.config_dependent_constants.DTS,
        )

        # -----------------------------------------------------------------------
        # initialize precipitation outputs
        # -----------------------------------------------------------------------

        self.outputs = Outputs(quantity_factory)

        # -----------------------------------------------------------------------
        # initialize temporaries
        # -----------------------------------------------------------------------

        self.temporaries = Temporaries(quantity_factory)

        # -----------------------------------------------------------------------
        # initialize masks
        # -----------------------------------------------------------------------

        self.masks = Masks(quantity_factory)

        # -----------------------------------------------------------------------
        # generate saturation specific humidity tables
        # -----------------------------------------------------------------------

        self.sat_tables = get_tables(stencil_factory.backend)

        # -----------------------------------------------------------------------
        # initialize stencils
        # -----------------------------------------------------------------------

        orchestrate(obj=self, config=stencil_factory.config.dace_config)
        self._init_temporaries = stencil_factory.from_dims_halo(
            func=init_temporaries,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "cpaut": self.config_dependent_constants.CPAUT,
            },
        )

        self._gfdl_1m_driver_preloop = stencil_factory.from_dims_halo(
            func=fix_negative_values,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "c_air": self.config_dependent_constants.C_AIR,
                "c_vap": self.config_dependent_constants.C_VAP,
                "p_nonhydro": self.config_dependent_constants.P_NONHYDRO,
                "d0_vap": self.config_dependent_constants.D0_VAP,
                "lv00": self.config_dependent_constants.LV00,
                "latv": self.config_dependent_constants.LATV,
                "lati": self.config_dependent_constants.LATI,
                "lats": self.config_dependent_constants.LATS,
                "lat2": self.config_dependent_constants.LAT2,
                "lcp": self.config_dependent_constants.LCP,
                "icp": self.config_dependent_constants.ICP,
                "tcp": self.config_dependent_constants.TCP,
                "mpdt": self.config_dependent_constants.MPDT,
                "rdt": self.config_dependent_constants.RDT,
                "ntimes": self.config_dependent_constants.NTIMES,
                "dts": self.config_dependent_constants.DTS,
                "do_sedi_w": GFDL_1M_config.do_sedi_w,
                "cpaut": self.config_dependent_constants.CPAUT,
                "hydrostatic": GFDL_1M_config.hydrostatic,
                "phys_hydrostatic": GFDL_1M_config.phys_hydrostatic,
                "fix_negative": GFDL_1M_config.fix_negative,
                "sedi_transport": GFDL_1M_config.sedi_transport,
                "const_vi": GFDL_1M_config.const_vi,
                "const_vs": GFDL_1M_config.const_vs,
                "const_vg": GFDL_1M_config.const_vg,
            },
        )

        self._fall_speed = stencil_factory.from_dims_halo(
            func=fall_speed,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "p_nonhydro": self.config_dependent_constants.P_NONHYDRO,
                "const_vi": GFDL_1M_config.const_vi,
                "const_vs": GFDL_1M_config.const_vs,
                "const_vg": GFDL_1M_config.const_vg,
                "vi_fac": GFDL_1M_config.vi_fac,
                "vi_max": GFDL_1M_config.vi_max,
                "vs_fac": GFDL_1M_config.vs_fac,
                "vs_max": GFDL_1M_config.vs_max,
                "vg_fac": GFDL_1M_config.vg_fac,
                "vg_max": GFDL_1M_config.vg_max,
            },
        )

        self._terminal_fall = TerminalFall(
            stencil_factory,
            quantity_factory,
            GFDL_1M_config,
            self.config_dependent_constants,
        )

        self._warm_rain = stencil_factory.from_dims_halo(
            func=warm_rain,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "dts": self.config_dependent_constants.DTS,
                "do_qa": GFDL_1M_config.do_qa,
                "rthreshs": GFDL_1M_config.rthreshs,
                "rthreshu": GFDL_1M_config.rthreshu,
                "irain_f": GFDL_1M_config.irain_f,
                "ql0_max": GFDL_1M_config.ql0_max,
                "z_slope_liq": GFDL_1M_config.z_slope_liq,
                "vr_fac": GFDL_1M_config.vr_fac,
                "const_vr": GFDL_1M_config.const_vr,
                "vr_max": GFDL_1M_config.vr_max,
                "tau_revp": GFDL_1M_config.tau_revp,
                "lv00": self.config_dependent_constants.LV00,
                "d0_vap": self.config_dependent_constants.D0_VAP,
                "c_air": self.config_dependent_constants.C_AIR,
                "c_vap": self.config_dependent_constants.C_VAP,
                "crevp_0": self.config_dependent_constants.CREVP_0,
                "crevp_1": self.config_dependent_constants.CREVP_1,
                "crevp_2": self.config_dependent_constants.CREVP_2,
                "crevp_3": self.config_dependent_constants.CREVP_3,
                "crevp_4": self.config_dependent_constants.CREVP_4,
                "cracw": self.config_dependent_constants.CRACW,
                "do_sedi_w": GFDL_1M_config.do_sedi_w,
                "use_ppm": GFDL_1M_config.use_ppm,
            },
        )

        self._warm_rain_update = stencil_factory.from_dims_halo(
            func=warm_rain_update,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._icloud = stencil_factory.from_dims_halo(
            func=icloud,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "c_air": self.config_dependent_constants.C_AIR,
                "c_vap": self.config_dependent_constants.C_VAP,
                "dts": self.config_dependent_constants.DTS,
                "rdts": self.config_dependent_constants.RDTS,
                "const_vi": GFDL_1M_config.const_vi,
                "fac_g2v": self.config_dependent_constants.FAC_G2V,
                "fac_i2s": self.config_dependent_constants.FAC_I2S,
                "fac_imlt": self.config_dependent_constants.FAC_IMLT,
                "fac_frz": self.config_dependent_constants.FAC_FRZ,
                "fac_l2v": self.config_dependent_constants.FAC_L2V,
                "fac_s2v": self.config_dependent_constants.FAC_S2V,
                "fac_v2s": self.config_dependent_constants.FAC_V2S,
                "fac_v2g": self.config_dependent_constants.FAC_V2G,
                "cgacs": self.config_dependent_constants.CGACS,
                "csacw": self.config_dependent_constants.CSACW,
                "csaci": self.config_dependent_constants.CSACI,
                "cgacw": self.config_dependent_constants.CGACW,
                "cgaci": self.config_dependent_constants.CGACI,
                "cgfr_0": self.config_dependent_constants.CGFR_0,
                "cgfr_1": self.config_dependent_constants.CGFR_1,
                "csmlt_0": self.config_dependent_constants.CSMLT_0,
                "csmlt_1": self.config_dependent_constants.CSMLT_1,
                "csmlt_2": self.config_dependent_constants.CSMLT_2,
                "csmlt_3": self.config_dependent_constants.CSMLT_3,
                "csmlt_4": self.config_dependent_constants.CSMLT_4,
                "cgmlt_0": self.config_dependent_constants.CGMLT_0,
                "cgmlt_1": self.config_dependent_constants.CGMLT_1,
                "cgmlt_2": self.config_dependent_constants.CGMLT_2,
                "cgmlt_3": self.config_dependent_constants.CGMLT_3,
                "cgmlt_4": self.config_dependent_constants.CGMLT_4,
                "cssub_0": self.config_dependent_constants.CSSUB_0,
                "cssub_1": self.config_dependent_constants.CSSUB_1,
                "cssub_2": self.config_dependent_constants.CSSUB_2,
                "cssub_3": self.config_dependent_constants.CSSUB_3,
                "cssub_4": self.config_dependent_constants.CSSUB_4,
                "qi0_crt": GFDL_1M_config.qi0_crt,
                "qs0_crt": GFDL_1M_config.qs0_crt,
                "qs_mlt": GFDL_1M_config.qs_mlt,
                "ql_mlt": GFDL_1M_config.ql_mlt,
                "z_slope_ice": GFDL_1M_config.z_slope_ice,
                "lv00": self.config_dependent_constants.LV00,
                "d0_vap": self.config_dependent_constants.D0_VAP,
                "lat2": self.config_dependent_constants.LAT2,
                "do_qa": GFDL_1M_config.do_qa,
                "do_evap": GFDL_1M_config.do_evap,
                "do_bigg": GFDL_1M_config.do_bigg,
                "qc_crt": GFDL_1M_config.qc_crt,
                "qi_lim": GFDL_1M_config.qi_lim,
                "rh_inc": GFDL_1M_config.rh_inc,
                "rh_inr": GFDL_1M_config.rh_inr,
                "t_min": GFDL_1M_config.t_min,
                "t_sub": GFDL_1M_config.t_sub,
                "preciprad": GFDL_1M_config.preciprad,
                "icloud_f": GFDL_1M_config.icloud_f,
            },
        )

        self._icloud_update = stencil_factory.from_dims_halo(
            func=icloud_update,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._update_tendencies = stencil_factory.from_dims_halo(
            func=update_tendencies,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "c_air": self.config_dependent_constants.C_AIR,
                "c_vap": self.config_dependent_constants.C_VAP,
                "rdt": self.config_dependent_constants.RDT,
                "do_sedi_w": GFDL_1M_config.do_sedi_w,
                "sedi_transport": GFDL_1M_config.sedi_transport,
                "do_qa": GFDL_1M_config.do_qa,
            },
        )

    def __call__(
        self,
        GFDL_1M_config: config,
        qv: FloatField,
        ql: FloatField,
        qr: FloatField,
        qi: FloatField,
        qs: FloatField,
        qg: FloatField,
        qa: FloatField,
        qn: FloatField,  # NACTL + NACTI
        qv_dt: FloatField,
        ql_dt: FloatField,
        qr_dt: FloatField,
        qi_dt: FloatField,
        qs_dt: FloatField,
        qg_dt: FloatField,
        qa_dt: FloatField,
        t_dt: FloatField,
        t: FloatField,
        w: FloatField,
        uin: FloatField,
        vin: FloatField,
        udt: FloatField,
        vdt: FloatField,
        dz: FloatField,
        dp: FloatField,
        area: FloatFieldIJ,
        fr_land: FloatFieldIJ,
        cnv_frc: FloatFieldIJ,
        srf_type: FloatFieldIJ,
        eis: FloatFieldIJ,
        rhcrit3d: FloatField,
        anv_icefall: Float,
        ls_icefall: Float,
    ):
        # The driver modifies a number of variables (t, p, qX) but does not pass
        # the changes back to the rest of the model. To replicate this behavior,
        # temporary copies of these variables are used throughout the driver.
        self._init_temporaries(
            t,
            dp,
            rhcrit3d,
            qv,
            ql,
            qi,
            qr,
            qs,
            qg,
            qa,
            qn,
            self.temporaries.qv0,
            self.temporaries.ql0,
            self.temporaries.qr0,
            self.temporaries.qi0,
            self.temporaries.qs0,
            self.temporaries.qg0,
            self.temporaries.qa0,
            self.temporaries.qv1,
            self.temporaries.ql1,
            self.temporaries.qr1,
            self.temporaries.qi1,
            self.temporaries.qs1,
            self.temporaries.qg1,
            self.temporaries.qa1,
            dz,
            uin,
            vin,
            w,
            area,
            self.temporaries.t1,
            self.temporaries.dp1,
            self.temporaries.omq,
            self.temporaries.den,
            self.temporaries.p_dry,
            self.temporaries.m1,
            self.temporaries.u1,
            self.temporaries.v1,
            self.temporaries.w1,
            self.temporaries.onemsig,
            self.temporaries.ccn,
            self.temporaries.c_praut,
            self.temporaries.rh_limited,
            self.outputs.rain,
            self.outputs.snow,
            self.outputs.graupel,
            self.outputs.ice,
            self.outputs.m2_rain,
            self.outputs.m2_sol,
            self.outputs.revap,
            self.outputs.isubl,
        )

        self._gfdl_1m_driver_preloop(
            self.temporaries.t1,
            self.temporaries.qv1,
            self.temporaries.ql1,
            self.temporaries.qr1,
            self.temporaries.qi1,
            self.temporaries.qs1,
            self.temporaries.qg1,
            self.temporaries.dp1,
        )

        for n in range(self.config_dependent_constants.NTIMES):
            self._fall_speed(
                self.temporaries.ql1,
                self.temporaries.qi1,
                self.temporaries.qs1,
                self.temporaries.qg1,
                t,
                self.temporaries.t1,
                dz,
                self.temporaries.dz1,
                self.temporaries.den,
                self.temporaries.den1,
                self.temporaries.denfac,
                self.temporaries.p_dry,
                self.temporaries.vti,
                self.temporaries.vts,
                self.temporaries.vtg,
                cnv_frc,
                anv_icefall,
                ls_icefall,
            )

            self._terminal_fall(
                self.temporaries.t1,
                self.temporaries.qv1,
                self.temporaries.ql1,
                self.temporaries.qr1,
                self.temporaries.qg1,
                self.temporaries.qs1,
                self.temporaries.qi1,
                self.temporaries.dz1,
                self.temporaries.dp1,
                self.temporaries.den1,
                self.temporaries.vtg,
                self.temporaries.vts,
                self.temporaries.vti,
                self.outputs.rain,
                self.outputs.graupel,
                self.outputs.snow,
                self.outputs.ice,
                self.temporaries.rain1,
                self.temporaries.graupel1,
                self.temporaries.snow1,
                self.temporaries.ice1,
                self.temporaries.m1_sol,
                self.temporaries.w1,
                self.temporaries.ze,
                self.temporaries.zt,
                self.masks.is_frozen,
                self.masks.precip_fall,
            )

            self._warm_rain(
                self.temporaries.dp1,
                self.temporaries.dz1,
                self.temporaries.t1,
                self.temporaries.qv1,
                self.temporaries.ql1,
                self.temporaries.qr1,
                self.temporaries.qi1,
                self.temporaries.qs1,
                self.temporaries.qg1,
                self.temporaries.qa1,
                self.temporaries.ccn,
                self.temporaries.den,
                self.temporaries.denfac,
                self.temporaries.c_praut,
                self.temporaries.vtr,
                self.temporaries.evap1,
                self.temporaries.m1_rain,
                self.temporaries.w1,
                self.temporaries.rh_limited,
                eis,
                self.temporaries.onemsig,
                self.temporaries.rain1,
                self.temporaries.ze,
                self.temporaries.zt,
                self.masks.precip_fall,
                self.sat_tables.table1,
                self.sat_tables.table2,
                self.sat_tables.table3,
                self.sat_tables.table4,
                self.sat_tables.des1,
                self.sat_tables.des2,
                self.sat_tables.des3,
                self.sat_tables.des4,
            )

            self._warm_rain_update(
                self.outputs.rain,
                self.temporaries.rain1,
                self.temporaries.evap1,
                self.outputs.revap,
                self.temporaries.m1_rain,
                self.outputs.m2_rain,
                self.temporaries.m1_sol,
                self.outputs.m2_sol,
                self.temporaries.m1,
            )

            self._icloud(
                self.temporaries.t1,
                self.temporaries.p_dry,
                self.temporaries.dp1,
                self.temporaries.qv1,
                self.temporaries.ql1,
                self.temporaries.qr1,
                self.temporaries.qi1,
                self.temporaries.qs1,
                self.temporaries.qg1,
                self.temporaries.qa1,
                self.temporaries.den1,
                self.temporaries.denfac,
                self.temporaries.vts,
                self.temporaries.vtg,
                self.temporaries.vtr,
                self.temporaries.subl1,
                self.temporaries.rh_limited,
                self.temporaries.ccn,
                cnv_frc,
                srf_type,
                self.sat_tables.table1,
                self.sat_tables.table2,
                self.sat_tables.table3,
                self.sat_tables.table4,
                self.sat_tables.des1,
                self.sat_tables.des2,
                self.sat_tables.des3,
                self.sat_tables.des4,
            )

            self._icloud_update(
                self.outputs.isubl,
                self.temporaries.subl1,
            )

        self._update_tendencies(
            self.temporaries.qv0,
            self.temporaries.ql0,
            self.temporaries.qr0,
            self.temporaries.qi0,
            self.temporaries.qs0,
            self.temporaries.qg0,
            self.temporaries.qa0,
            self.temporaries.qv1,
            self.temporaries.ql1,
            self.temporaries.qr1,
            self.temporaries.qi1,
            self.temporaries.qs1,
            self.temporaries.qg1,
            self.temporaries.qa1,
            self.temporaries.ccn,
            qv_dt,
            ql_dt,
            qr_dt,
            qi_dt,
            qs_dt,
            qg_dt,
            qa_dt,
            t,
            self.temporaries.t1,
            t_dt,
            w,
            self.temporaries.w1,
            uin,
            self.temporaries.u1,
            udt,
            vin,
            self.temporaries.v1,
            vdt,
            dz,
            dp,
            self.temporaries.dp1,
            self.temporaries.den,
            self.temporaries.p_dry,
            area,
            GFDL_1M_config.dt_moist,
            fr_land,
            cnv_frc,
            srf_type,
            eis,
            self.temporaries.rh_limited,
            self.temporaries.m1,
            anv_icefall,
            ls_icefall,
            self.outputs.revap,
            self.outputs.isubl,
            self.outputs.rain,
            self.outputs.snow,
            self.outputs.ice,
            self.outputs.graupel,
            self.outputs.m2_rain,
            self.outputs.m2_sol,
        )
