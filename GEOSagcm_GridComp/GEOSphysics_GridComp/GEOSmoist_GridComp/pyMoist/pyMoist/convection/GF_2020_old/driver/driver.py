from ndsl import StencilFactory, QuantityFactory
from pyMoist.convective_parameterization.GF_2020.driver.setup import setup
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from pyMoist.convective_parameterization.GF_2020.config import GF2020Config
from pyMoist.convective_parameterization.GF_2020.driver.temporaries import GF2020DriverTemporaries


class GF2020Driver:
    def __init__(
        self, stencil_factory: StencilFactory, quantity_factory: QuantityFactory, GF_2020_config: GF2020Config
    ):
        # Construct stencils
        self.setup = stencil_factory.from_dims_halo(
            func=setup,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "C1": GF_2020_config.C1,
                "ADV_TRIGGER": GF_2020_config.ADV_TRIGGER,
                "AUTOCONV": GF_2020_config.AUTOCONV,
                "USE_TRACER_TRANSP": GF_2020_config.USE_TRACER_TRANSP,
            },
        )

        self.driver_temporaries = GF2020DriverTemporaries.make(quantity_factory)

    def setup(self, *args, **kwds):
        pass

    def finalize(self, *args, **kwds):
        pass
