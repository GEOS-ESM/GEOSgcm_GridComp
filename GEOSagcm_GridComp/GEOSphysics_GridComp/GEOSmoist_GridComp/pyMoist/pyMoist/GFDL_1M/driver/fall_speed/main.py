from ndsl import QuantityFactory, StencilFactory, orchestrate
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ
from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.GFDL_1M.driver.config_constants import ConfigConstants
from pyMoist.GFDL_1M.driver.fall_speed.stencils import fall_speed


class FallSpeed:
    """
    Ice cloud microphysics processes
    bulk cloud micro - physics; processes splitting
    with some un - split sub - grouping
    time implicit (when possible) accretion and autoconversion
    """

    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        GFDL_1M_config: GFDL1MConfig,
        config_dependent_constants: ConfigConstants,
    ):

        orchestrate(obj=self, config=stencil_factory.config.dace_config)
        # Initalize stencils
        self._fall_speed_core = stencil_factory.from_dims_halo(
            func=fall_speed,
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
        ql1: FloatField,
        qi1: FloatField,
        qs1: FloatField,
        qg1: FloatField,
        t: FloatField,
        t1: FloatField,
        dz: FloatField,
        dz1: FloatField,
        den: FloatField,
        den1: FloatField,
        denfac: FloatField,
        p_dry: FloatField,
        vti: FloatField,
        vts: FloatField,
        vtg: FloatField,
        cnv_frc: FloatFieldIJ,
        anv_icefall: Float,
        ls_icefall: Float,
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
            vti,
            vts,
            vtg,
            cnv_frc,
            anv_icefall,
            ls_icefall,
        )
