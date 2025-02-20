from ndsl import QuantityFactory, StencilFactory, orchestrate
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from pyMoist.GFDL_1M.driver.config import config
from pyMoist.GFDL_1M.driver.icecloud.stencils import (
    icloud_core,
    update_output,
)
from pyMoist.GFDL_1M.driver.support import ConfigConstants


class IceCloud:
    """
    ice cloud microphysics processes
    bulk cloud micro - physics; processes splitting
    with some un - split sub - grouping
    time implicit (when possible) accretion and autoconversion
    """

    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        GFDL_1M_config: config,
        config_dependent_constants: ConfigConstants,
    ):

        # Initalize stencils
        orchestrate(obj=self, config=stencil_factory.config.dace_config)
        self._icloud_core = stencil_factory.from_dims_halo(
            func=icloud_core,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "c_air": config_dependent_constants.C_AIR,
                "c_vap": config_dependent_constants.C_VAP,
                "dts": config_dependent_constants.DTS,
                "rdts": config_dependent_constants.RDTS,
                "const_vi": GFDL_1M_config.const_vi,
                "fac_g2v": config_dependent_constants.FAC_G2V,
                "fac_i2s": config_dependent_constants.FAC_I2S,
                "fac_imlt": config_dependent_constants.FAC_IMLT,
                "fac_frz": config_dependent_constants.FAC_FRZ,
                "fac_l2v": config_dependent_constants.FAC_L2V,
                "fac_s2v": config_dependent_constants.FAC_S2V,
                "fac_v2s": config_dependent_constants.FAC_V2S,
                "fac_v2g": config_dependent_constants.FAC_V2G,
                "cgacs": config_dependent_constants.CGACS,
                "csacw": config_dependent_constants.CSACW,
                "csaci": config_dependent_constants.CSACI,
                "cgacw": config_dependent_constants.CGACW,
                "cgaci": config_dependent_constants.CGACI,
                "cgfr_0": config_dependent_constants.CGFR_0,
                "cgfr_1": config_dependent_constants.CGFR_1,
                "csmlt_0": config_dependent_constants.CSMLT_0,
                "csmlt_1": config_dependent_constants.CSMLT_1,
                "csmlt_2": config_dependent_constants.CSMLT_2,
                "csmlt_3": config_dependent_constants.CSMLT_3,
                "csmlt_4": config_dependent_constants.CSMLT_4,
                "cgmlt_0": config_dependent_constants.CGMLT_0,
                "cgmlt_1": config_dependent_constants.CGMLT_1,
                "cgmlt_2": config_dependent_constants.CGMLT_2,
                "cgmlt_3": config_dependent_constants.CGMLT_3,
                "cgmlt_4": config_dependent_constants.CGMLT_4,
                "cssub_0": config_dependent_constants.CSSUB_0,
                "cssub_1": config_dependent_constants.CSSUB_1,
                "cssub_2": config_dependent_constants.CSSUB_2,
                "cssub_3": config_dependent_constants.CSSUB_3,
                "cssub_4": config_dependent_constants.CSSUB_4,
                "qi0_crt": GFDL_1M_config.qi0_crt,
                "qs0_crt": GFDL_1M_config.qs0_crt,
                "qs_mlt": GFDL_1M_config.qs_mlt,
                "ql_mlt": GFDL_1M_config.ql_mlt,
                "z_slope_ice": GFDL_1M_config.z_slope_ice,
                "lv00": config_dependent_constants.LV00,
                "d0_vap": config_dependent_constants.D0_VAP,
                "lat2": config_dependent_constants.LAT2,
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

        self._update_output = stencil_factory.from_dims_halo(
            func=update_output,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        t1,
        p_dry,
        dp1,
        qv1,
        ql1,
        qr1,
        qi1,
        qs1,
        qg1,
        qa1,
        den1,
        denfac,
        vts,
        vtg,
        vtr,
        subl1,
        isubl,
        rh_limited,
        ccn,
        cnv_frc,
        srf_type,
        table1,
        table2,
        table3,
        table4,
        des1,
        des2,
        des3,
        des4,
    ):
        self._icloud_core(
            t1,
            p_dry,
            dp1,
            qv1,
            ql1,
            qr1,
            qi1,
            qs1,
            qg1,
            qa1,
            den1,
            denfac,
            vts,
            vtg,
            vtr,
            subl1,
            rh_limited,
            ccn,
            cnv_frc,
            srf_type,
            table1,
            table2,
            table3,
            table4,
            des1,
            des2,
            des3,
            des4,
        )

        self._update_output(
            isubl,
            subl1,
        )
