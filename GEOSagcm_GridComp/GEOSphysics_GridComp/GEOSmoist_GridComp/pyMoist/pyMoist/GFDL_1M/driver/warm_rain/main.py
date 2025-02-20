from ndsl import QuantityFactory, StencilFactory, orchestrate
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from pyMoist.GFDL_1M.driver.config import config
from pyMoist.GFDL_1M.driver.warm_rain.stencils import (
    warm_rain_core,
    update_outputs,
)
from pyMoist.GFDL_1M.driver.terminal_fall.temporaries import Temporaries
from pyMoist.GFDL_1M.driver.support import ConfigConstants


class WarmRain:
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

        self._warm_rain_core = stencil_factory.from_dims_halo(
            func=warm_rain_core,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "dts": config_dependent_constants.DTS,
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
                "lv00": config_dependent_constants.LV00,
                "d0_vap": config_dependent_constants.D0_VAP,
                "c_air": config_dependent_constants.C_AIR,
                "c_vap": config_dependent_constants.C_VAP,
                "crevp_0": config_dependent_constants.CREVP_0,
                "crevp_1": config_dependent_constants.CREVP_1,
                "crevp_2": config_dependent_constants.CREVP_2,
                "crevp_3": config_dependent_constants.CREVP_3,
                "crevp_4": config_dependent_constants.CREVP_4,
                "cracw": config_dependent_constants.CRACW,
                "do_sedi_w": GFDL_1M_config.do_sedi_w,
                "use_ppm": GFDL_1M_config.use_ppm,
            },
        )

        self.update_outputs = stencil_factory.from_dims_halo(
            func=update_outputs,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        dz1,
        t1,
        qv1,
        dp1,
        ql1,
        qr1,
        qi1,
        qs1,
        qg1,
        qa1,
        ccn,
        den,
        denfac,
        c_praut,
        vtr,
        evap1,
        m1_rain,
        w1,
        rh_limited,
        eis,
        onemsig,
        rain1,
        ze,
        zt,
        m1,
        m1_sol,
        rain,
        revap,
        m2_rain,
        m2_sol,
        precip_fall,
        table1,
        table2,
        table3,
        table4,
        des1,
        des2,
        des3,
        des4,
    ):
        self._warm_rain_core(
            dp1,
            dz1,
            t1,
            qv1,
            ql1,
            qr1,
            qi1,
            qs1,
            qg1,
            qa1,
            ccn,
            den,
            denfac,
            c_praut,
            vtr,
            evap1,
            m1_rain,
            w1,
            rh_limited,
            eis,
            onemsig,
            rain1,
            ze,
            zt,
            precip_fall,
            table1,
            table2,
            table3,
            table4,
            des1,
            des2,
            des3,
            des4,
        )

        self.update_outputs(
            m1_rain,
            m1_sol,
            rain1,
            evap1,
            revap,
            m2_rain,
            m2_sol,
            m1,
            rain,
        )
