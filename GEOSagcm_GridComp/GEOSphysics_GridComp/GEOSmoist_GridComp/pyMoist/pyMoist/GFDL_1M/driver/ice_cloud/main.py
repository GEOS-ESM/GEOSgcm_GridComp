from ndsl import StencilFactory, orchestrate
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField, FloatFieldIJ
from pyMoist.GFDL_1M.driver.config import MicrophysicsConfiguration
from pyMoist.GFDL_1M.driver.config_constants import ConfigConstants
from pyMoist.GFDL_1M.driver.ice_cloud.stencils import icloud_core, update_precip_total
from pyMoist.GFDL_1M.driver.sat_tables import GlobalTable_driver_qsat


class IceCloud:
    """
    Ice cloud microphysics processes
    bulk cloud micro - physics; processes splitting
    with some un - split sub - grouping
    time implicit (when possible) accretion and autoconversion
    """

    def __init__(
        self,
        stencil_factory: StencilFactory,
        GFDL_1M_config: MicrophysicsConfiguration,
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
                "const_vi": GFDL_1M_config.CONST_VI,
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
                "qi0_crt": GFDL_1M_config.QI0_CRT,
                "qs0_crt": GFDL_1M_config.QS0_CRT,
                "qs_mlt": GFDL_1M_config.QS_MLT,
                "ql_mlt": GFDL_1M_config.QL_MLT,
                "z_slope_ice": GFDL_1M_config.Z_SLOPE_ICE,
                "lv00": config_dependent_constants.LV00,
                "d0_vap": config_dependent_constants.D0_VAP,
                "lat2": config_dependent_constants.LAT2,
                "do_qa": GFDL_1M_config.DO_QA,
                "do_evap": GFDL_1M_config.DO_EVAP,
                "do_bigg": GFDL_1M_config.DO_BIGG,
                "qc_crt": GFDL_1M_config.QC_CRT,
                "qi_lim": GFDL_1M_config.QI_LIM,
                "rh_inc": GFDL_1M_config.RH_INC,
                "rh_inr": GFDL_1M_config.RH_INR,
                "t_min": GFDL_1M_config.T_MIN,
                "t_sub": GFDL_1M_config.T_SUB,
                "preciprad": GFDL_1M_config.PRECIPRAD,
                "icloud_f": GFDL_1M_config.ICLOUD_F,
            },
        )

        self._update_output = stencil_factory.from_dims_halo(
            func=update_precip_total,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        t1: FloatField,
        p_dry: FloatField,
        dp1: FloatField,
        qv1: FloatField,
        ql1: FloatField,
        qr1: FloatField,
        qi1: FloatField,
        qs1: FloatField,
        qg1: FloatField,
        qa1: FloatField,
        den1: FloatField,
        denfac: FloatField,
        vts: FloatField,
        vtg: FloatField,
        vtr: FloatField,
        subl1: FloatField,
        isubl: FloatField,
        rh_limited: FloatField,
        ccn: FloatField,
        cnv_frc: FloatFieldIJ,
        srf_type: FloatFieldIJ,
        table2: GlobalTable_driver_qsat,
        table3: GlobalTable_driver_qsat,
        des2: GlobalTable_driver_qsat,
        des3: GlobalTable_driver_qsat,
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
            table2,
            table3,
            des2,
            des3,
        )

        self._update_output(
            isubl,
            subl1,
        )
