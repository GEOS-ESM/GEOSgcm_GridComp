from ndsl import QuantityFactory, StencilFactory, orchestrate
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import BoolFieldIJ, FloatField, FloatFieldIJ
from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.GFDL_1M.driver.config_constants import ConfigConstants
from pyMoist.GFDL_1M.driver.sat_tables import GlobalTable_driver_qsat
from pyMoist.GFDL_1M.driver.stencils import implicit_fall
from pyMoist.GFDL_1M.driver.warm_rain.stencils import (
    update_outputs,
    warm_rain_step_1,
    warm_rain_step_2,
)
from pyMoist.GFDL_1M.driver.warm_rain.temporaries import Temporaries


class WarmRain:
    """
    Warm rain cloud microphysics: evaporation, accretion
    """

    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        GFDL_1M_config: GFDL1MConfig,
        config_dependent_constants: ConfigConstants,
    ):
        self.GFDL_1M_config = GFDL_1M_config

        # Initalize temporaries
        self.temporaries = Temporaries.make(quantity_factory)

        # Initalize stencils
        orchestrate(obj=self, config=stencil_factory.config.dace_config)

        self._warm_rain_step_1 = stencil_factory.from_dims_halo(
            func=warm_rain_step_1,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "dts": config_dependent_constants.DTS,
                "do_qa": GFDL_1M_config.DO_QA,
                "rthreshs": GFDL_1M_config.RTHRESHS,
                "rthreshu": GFDL_1M_config.RTHRESHU,
                "irain_f": GFDL_1M_config.IRAIN_F,
                "ql0_max": GFDL_1M_config.QL0_MAX,
                "z_slope_liq": GFDL_1M_config.Z_SLOPE_LIQ,
                "vr_fac": GFDL_1M_config.VR_FAC,
                "const_vr": GFDL_1M_config.CONST_VR,
                "vr_max": GFDL_1M_config.VR_MAX,
                "tau_revp": GFDL_1M_config.TAU_REVP,
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
                "do_sedi_w": GFDL_1M_config.DO_SEDI_W,
                "use_ppm": GFDL_1M_config.USE_PPM,
            },
        )

        self._implicit_fall = stencil_factory.from_dims_halo(
            func=implicit_fall,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "dts": config_dependent_constants.DTS,
                "use_ppm": GFDL_1M_config.USE_PPM,
            },
        )

        self._warm_rain_step_2 = stencil_factory.from_dims_halo(
            func=warm_rain_step_2,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "dts": config_dependent_constants.DTS,
                "do_qa": GFDL_1M_config.DO_QA,
                "rthreshs": GFDL_1M_config.RTHRESHS,
                "rthreshu": GFDL_1M_config.RTHRESHU,
                "irain_f": GFDL_1M_config.IRAIN_F,
                "ql0_max": GFDL_1M_config.QL0_MAX,
                "z_slope_liq": GFDL_1M_config.Z_SLOPE_LIQ,
                "vr_fac": GFDL_1M_config.VR_FAC,
                "const_vr": GFDL_1M_config.CONST_VR,
                "vr_max": GFDL_1M_config.VR_MAX,
                "tau_revp": GFDL_1M_config.TAU_REVP,
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
                "do_sedi_w": GFDL_1M_config.DO_SEDI_W,
                "use_ppm": GFDL_1M_config.USE_PPM,
            },
        )

        self.update_outputs = stencil_factory.from_dims_halo(
            func=update_outputs,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        dz1: FloatField,
        t1: FloatField,
        qv1: FloatField,
        dp1: FloatField,
        ql1: FloatField,
        qr1: FloatField,
        qi1: FloatField,
        qs1: FloatField,
        qg1: FloatField,
        qa1: FloatField,
        ccn: FloatField,
        den: FloatField,
        denfac: FloatField,
        c_praut: FloatField,
        vtr: FloatField,
        evap1: FloatField,
        m1_rain: FloatField,
        w1: FloatField,
        rh_limited: FloatField,
        eis: FloatFieldIJ,
        onemsig: FloatFieldIJ,
        rain1: FloatFieldIJ,
        ze: FloatField,
        zt: FloatField,
        m1: FloatField,
        m1_sol: FloatField,
        rain: FloatFieldIJ,
        revap: FloatField,
        m2_rain: FloatField,
        m2_sol: FloatField,
        precip_fall: BoolFieldIJ,
        table1: GlobalTable_driver_qsat,
        table2: GlobalTable_driver_qsat,
        table3: GlobalTable_driver_qsat,
        table4: GlobalTable_driver_qsat,
        des1: GlobalTable_driver_qsat,
        des2: GlobalTable_driver_qsat,
        des3: GlobalTable_driver_qsat,
        des4: GlobalTable_driver_qsat,
    ):
        self._warm_rain_step_1(
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
            rh_limited,
            eis,
            onemsig,
            ze,
            self.temporaries.dm,
            precip_fall,
            table1,
            table2,
            table3,
            table4,
            des1,
            des2,
            des3,
            des4,
            self.temporaries.test_var_1,
        )
        if self.GFDL_1M_config.USE_PPM == False:
            # NOTE: somehow errors pop up in rain1 and m1_sol within implicit fall, despite all of the
            # imputs being correct and implicit_fall verifying at three separate calls
            # within the terminal_fall module. May be a similar issue to the warm_rain_part_1 error
            # (different result despite inputs being identical, possible registry issue??).
            self._implicit_fall(
                qr1,
                vtr,
                ze,
                dp1,
                m1,
                m1_sol,
                rain1,
                precip_fall,
            )

        self._warm_rain_step_2(
            t1,
            qv1,
            ql1,
            qr1,
            qi1,
            qs1,
            qg1,
            qa1,
            den,
            denfac,
            vtr,
            evap1,
            m1_rain,
            w1,
            rh_limited,
            self.temporaries.dm,
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
