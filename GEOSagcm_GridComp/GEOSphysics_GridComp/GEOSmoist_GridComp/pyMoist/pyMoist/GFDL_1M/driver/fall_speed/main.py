from ndsl import QuantityFactory, StencilFactory, orchestrate
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from pyMoist.GFDL_1M.driver.config import config
from pyMoist.GFDL_1M.driver.fall_speed.stencils import (
    fall_speed_core,
)
from pyMoist.GFDL_1M.driver.support import ConfigConstants


class FallSpeed:
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

        self._fall_speed_core = stencil_factory.from_dims_halo(
            func=fall_speed_core,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "p_nonhydro": config_dependent_constants.P_NONHYDRO,
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

    def __call__(
        self,
        ql1,
        qi1,
        qs1,
        qg1,
        t,
        t1,
        dz,
        dz1,
        den,
        den1,
        denfac,
        p_dry,
        vti,
        vts,
        vtg,
        cnv_frc,
        anv_icefall,
        ls_icefall,
    ):
        self._fall_speed_core(
            ql1,
            qi1,
            qs1,
            qg1,
            t,
            t1,
            dz,
            dz1,
            den,
            den1,
            denfac,
            p_dry,
            vti,
            vts,
            vtg,
            cnv_frc,
            anv_icefall,
            ls_icefall,
        )
