from ndsl import QuantityFactory, StencilFactory, orchestrate
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from pyMoist.GFDL_1M.driver.config import config
from pyMoist.GFDL_1M.driver.terminal_fall.stencils import (
    check_precip_get_zt,
    melting_loop,
    reset,
    setup,
    update_dm,
    update_w1,
    update_outputs,
)
from pyMoist.GFDL_1M.driver.stencils import implicit_fall
from pyMoist.GFDL_1M.driver.terminal_fall.temporaries import Temporaries
from pyMoist.GFDL_1M.driver.support import ConfigConstants


class TerminalFall:
    """
    calculate terminal fall speed, accounting for
    melting of ice, snow, and graupel during fall

    reference Fortran: gfdl_cloud_microphys.F90: subroutine terminal_fall
    """

    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        GFDL_1M_config: config,
        config_dependent_constants: ConfigConstants,
    ):

        self.GFDL_1M_config = GFDL_1M_config

        # Initalize temporaries
        self.temporaries = Temporaries(quantity_factory)

        # Initalize stencils
        orchestrate(obj=self, config=stencil_factory.config.dace_config)
        self._setup = stencil_factory.from_dims_halo(
            func=setup,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "dts": config_dependent_constants.DTS,
                "tau_imlt": GFDL_1M_config.TAU_IMLT,
                "ql_mlt": GFDL_1M_config.QL_MLT,
                "c_air": config_dependent_constants.C_AIR,
                "c_vap": config_dependent_constants.C_VAP,
                "d0_vap": config_dependent_constants.D0_VAP,
                "lv00": config_dependent_constants.LV00,
            },
        )

        self._check_precip_get_zt = stencil_factory.from_dims_halo(
            func=check_precip_get_zt,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "dts": config_dependent_constants.DTS,
            },
        )

        self._melting_loop = stencil_factory.from_dims_halo(
            func=melting_loop,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._update_dm = stencil_factory.from_dims_halo(
            func=update_dm,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "do_sedi_w": GFDL_1M_config.DO_SEDI_W,
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

        self._update_w1 = stencil_factory.from_dims_halo(
            func=update_w1,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "do_sedi_w": GFDL_1M_config.DO_SEDI_W,
            },
        )

        self._reset = stencil_factory.from_dims_halo(
            func=reset,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._update_outputs = stencil_factory.from_dims_halo(
            func=update_outputs,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        t1,
        qv1,
        ql1,
        qr1,
        qg1,
        qs1,
        qi1,
        dz1,
        dp1,
        den1,
        vtg,
        vts,
        vti,
        rain,
        graupel,
        snow,
        ice,
        precip_rain,
        precip_graupel,
        precip_snow,
        precip_ice,
        m1_sol,
        w1,
        ze,
        zt,
        is_frozen,
        precip_fall,
    ):
        self._setup(
            t1,
            qv1,
            ql1,
            qr1,
            qg1,
            qs1,
            qi1,
            dz1,
            m1_sol,
            ze,
            zt,
            self.temporaries.lhi,
            self.temporaries.icpk,
            self.temporaries.cvm,
            is_frozen,
        )

        # -----------------------------------------------------------------------
        # -----------------------------------------------------------------------
        # melting of falling cloud ice into rain
        # -----------------------------------------------------------------------
        # -----------------------------------------------------------------------
        if self.GFDL_1M_config.VI_FAC < 1.0e-5:
            precip_ice.view[:] = 0
        else:
            self._check_precip_get_zt(
                qi1,
                vti,
                ze,
                zt,
                precip_ice,
                precip_fall,
            )

            # placeholder function. need the listed features to implement
            self._melting_loop(
                self.temporaries.need_2d_temporaries_feature,
                self.temporaries.need_double_k_loop_feature,
            )

            self._update_dm(
                self.temporaries.dm,
                dp1,
                qv1,
                ql1,
                qr1,
                qi1,
                qs1,
                qg1,
                precip_fall,
            )

            self._implicit_fall(
                qi1,
                vti,
                ze,
                dp1,
                self.temporaries.m1,
                m1_sol,
                precip_ice,
                precip_fall,
            )

            self._update_w1(
                w1,
                self.temporaries.dm,
                self.temporaries.m1,
                vti,
                precip_fall,
            )

            self._reset(
                self.temporaries.m1,
                precip_fall,
            )

        # -----------------------------------------------------------------------
        # -----------------------------------------------------------------------
        # melting of falling snow into rain
        # -----------------------------------------------------------------------
        # -----------------------------------------------------------------------

        self._check_precip_get_zt(
            qs1,
            vts,
            ze,
            zt,
            precip_snow,
            precip_fall,
        )

        # placeholder function. need the listed features to implement
        self._melting_loop(
            self.temporaries.need_2d_temporaries_feature,
            self.temporaries.need_double_k_loop_feature,
        )

        self._update_dm(
            self.temporaries.dm,
            dp1,
            qv1,
            ql1,
            qr1,
            qi1,
            qs1,
            qg1,
            precip_fall,
        )

        self._implicit_fall(
            qs1,
            vts,
            ze,
            dp1,
            self.temporaries.m1,
            m1_sol,
            precip_snow,
            precip_fall,
        )

        self._update_w1(
            w1,
            self.temporaries.dm,
            self.temporaries.m1,
            vts,
            precip_fall,
        )

        self._reset(
            self.temporaries.m1,
            precip_fall,
        )

        # -----------------------------------------------------------------------
        # -----------------------------------------------------------------------
        # melting of falling graupel into rain
        # -----------------------------------------------------------------------
        # -----------------------------------------------------------------------

        self._check_precip_get_zt(
            qg1,
            vtg,
            ze,
            zt,
            precip_graupel,
            precip_fall,
        )

        # placeholder function. need the listed features to implement
        self._melting_loop(
            self.temporaries.need_2d_temporaries_feature,
            self.temporaries.need_double_k_loop_feature,
        )

        self._update_dm(
            self.temporaries.dm,
            dp1,
            qv1,
            ql1,
            qr1,
            qi1,
            qs1,
            qg1,
            precip_fall,
        )

        self._implicit_fall(
            qg1,
            vtg,
            ze,
            dp1,
            self.temporaries.m1,
            m1_sol,
            precip_graupel,
            precip_fall,
        )

        self._update_w1(
            w1,
            self.temporaries.dm,
            self.temporaries.m1,
            vtg,
            precip_fall,
        )

        self._reset(
            self.temporaries.m1,
            precip_fall,
        )

        self._update_outputs(
            rain,
            graupel,
            snow,
            ice,
            precip_rain,
            precip_graupel,
            precip_snow,
            precip_ice,
        )
