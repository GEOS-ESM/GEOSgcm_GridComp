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
                "const_vi": GFDL_1M_config.CONST_VI,
                "const_vs": GFDL_1M_config.CONST_VS,
                "const_vg": GFDL_1M_config.CONST_VG,
                "vi_fac": GFDL_1M_config.VI_FAC,
                "vi_max": GFDL_1M_config.VI_MAX,
                "vs_fac": GFDL_1M_config.VS_FAC,
                "vs_max": GFDL_1M_config.VS_MAX,
                "vg_fac": GFDL_1M_config.VG_FAC,
                "vg_max": GFDL_1M_config.VG_MAX,
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
