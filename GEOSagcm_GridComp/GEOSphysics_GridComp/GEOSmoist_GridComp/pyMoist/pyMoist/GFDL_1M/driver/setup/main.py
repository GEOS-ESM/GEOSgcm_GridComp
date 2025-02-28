from ndsl import QuantityFactory, StencilFactory, orchestrate
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatFieldIJ, FloatField
from pyMoist.GFDL_1M.driver.config import config
from pyMoist.GFDL_1M.driver.setup.stencils import (
    init_temporaries,
    fix_negative_values,
)
from pyMoist.GFDL_1M.driver.config_constants import ConfigConstants


class Setup:

    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        GFDL_1M_config: config,
        config_dependent_constants: ConfigConstants,
    ):
        """
        The driver modifies a number of variables (t, p, qX) but does not pass
        the changes back to the rest of the model. To replicate this behavior,
        temporary copies of these variables are used throughout the driver.
        """
        # Initalize stencils
        orchestrate(obj=self, config=stencil_factory.config.dace_config)

        self._init_temporaries = stencil_factory.from_dims_halo(
            func=init_temporaries,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "cpaut": config_dependent_constants.CPAUT,
            },
        )

        self._fix_negative_values = stencil_factory.from_dims_halo(
            func=fix_negative_values,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "c_air": config_dependent_constants.C_AIR,
                "c_vap": config_dependent_constants.C_VAP,
                "p_nonhydro": config_dependent_constants.P_NONHYDRO,
                "d0_vap": config_dependent_constants.D0_VAP,
                "lv00": config_dependent_constants.LV00,
                "latv": config_dependent_constants.LATV,
                "lati": config_dependent_constants.LATI,
                "lats": config_dependent_constants.LATS,
                "lat2": config_dependent_constants.LAT2,
                "lcp": config_dependent_constants.LCP,
                "icp": config_dependent_constants.ICP,
                "tcp": config_dependent_constants.TCP,
                "mpdt": config_dependent_constants.MPDT,
                "rdt": config_dependent_constants.RDT,
                "ntimes": config_dependent_constants.NTIMES,
                "dts": config_dependent_constants.DTS,
                "do_sedi_w": GFDL_1M_config.DO_SEDI_W,
                "cpaut": config_dependent_constants.CPAUT,
                "hydrostatic": GFDL_1M_config.HYDROSTATIC,
                "phys_hydrostatic": GFDL_1M_config.PHYS_HYDROSTATIC,
                "fix_negative": GFDL_1M_config.FIX_NEGATIVE,
                "sedi_transport": GFDL_1M_config.SEDI_TRANSPORT,
                "const_vi": GFDL_1M_config.CONST_VI,
                "const_vs": GFDL_1M_config.CONST_VS,
                "const_vg": GFDL_1M_config.CONST_VG,
            },
        )

    def __call__(
        self,
        t: FloatField,
        dp: FloatField,
        rhcrit3d: FloatField,
        qv: FloatField,
        ql: FloatField,
        qi: FloatField,
        qr: FloatField,
        qs: FloatField,
        qg: FloatField,
        qa: FloatField,
        qn: FloatField,
        qv0: FloatField,
        ql0: FloatField,
        qr0: FloatField,
        qi0: FloatField,
        qs0: FloatField,
        qg0: FloatField,
        qa0: FloatField,
        qv1: FloatField,
        ql1: FloatField,
        qr1: FloatField,
        qi1: FloatField,
        qs1: FloatField,
        qg1: FloatField,
        qa1: FloatField,
        dz: FloatField,
        u: FloatField,
        v: FloatField,
        w: FloatField,
        area: FloatFieldIJ,
        t1: FloatField,
        dp1: FloatField,
        omq: FloatField,
        den: FloatField,
        p_dry: FloatField,
        m1: FloatField,
        u1: FloatField,
        v1: FloatField,
        w1: FloatField,
        onemsig: FloatFieldIJ,
        ccn: FloatField,
        c_praut: FloatField,
        rh_limited: FloatField,
        rain: FloatFieldIJ,
        snow: FloatFieldIJ,
        graupel: FloatFieldIJ,
        ice: FloatFieldIJ,
        m2_rain: FloatField,
        m2_sol: FloatField,
        revap: FloatField,
        isubl: FloatField,
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
            qv0,
            ql0,
            qr0,
            qi0,
            qs0,
            qg0,
            qa0,
            qv1,
            ql1,
            qr1,
            qi1,
            qs1,
            qg1,
            qa1,
            dz,
            u,
            v,
            w,
            area,
            t1,
            dp1,
            omq,
            den,
            p_dry,
            m1,
            u1,
            v1,
            w1,
            onemsig,
            ccn,
            c_praut,
            rh_limited,
            rain,
            snow,
            graupel,
            ice,
            m2_rain,
            m2_sol,
            revap,
            isubl,
        )

        self._fix_negative_values(
            t1,
            qv1,
            ql1,
            qr1,
            qi1,
            qs1,
            qg1,
            dp1,
        )
