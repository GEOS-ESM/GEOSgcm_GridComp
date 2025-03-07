from ndsl import StencilFactory, orchestrate
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField, FloatFieldIJ
from pyMoist.GFDL_1M.driver.config import MicrophysicsConfiguration
from pyMoist.GFDL_1M.driver.config_constants import ConfigConstants
from pyMoist.GFDL_1M.driver.finish.stencils import update_tendencies


class Finish:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        GFDL_1M_config: MicrophysicsConfiguration,
        config_dependent_constants: ConfigConstants,
    ):
        """
        Save driver outputs: tentencies & precipitation
        """
        # Initalize stencils
        orchestrate(obj=self, config=stencil_factory.config.dace_config)

        self._update_tendencies = stencil_factory.from_dims_halo(
            func=update_tendencies,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "c_air": config_dependent_constants.C_AIR,
                "c_vap": config_dependent_constants.C_VAP,
                "rdt": config_dependent_constants.RDT,
                "do_sedi_w": GFDL_1M_config.DO_SEDI_W,
                "sedi_transport": GFDL_1M_config.SEDI_TRANSPORT,
                "do_qa": GFDL_1M_config.DO_QA,
            },
        )

    def __call__(
        self,
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
        qv_dt: FloatField,
        ql_dt: FloatField,
        qr_dt: FloatField,
        qi_dt: FloatField,
        qs_dt: FloatField,
        qg_dt: FloatField,
        qa_dt: FloatField,
        t: FloatField,
        t1: FloatField,
        t_dt: FloatField,
        w: FloatField,
        w1: FloatField,
        u: FloatField,
        u1: FloatField,
        udt: FloatField,
        v: FloatField,
        v1: FloatField,
        vdt: FloatField,
        dp: FloatField,
        dp1: FloatField,
        m1: FloatField,
        rain: FloatFieldIJ,
        snow: FloatFieldIJ,
        ice: FloatFieldIJ,
        graupel: FloatFieldIJ,
    ):
        self._update_tendencies(
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
            qv_dt,
            ql_dt,
            qr_dt,
            qi_dt,
            qs_dt,
            qg_dt,
            qa_dt,
            t,
            t1,
            t_dt,
            w,
            w1,
            u,
            u1,
            udt,
            v,
            v1,
            vdt,
            dp,
            dp1,
            m1,
            rain,
            snow,
            ice,
            graupel,
        )
