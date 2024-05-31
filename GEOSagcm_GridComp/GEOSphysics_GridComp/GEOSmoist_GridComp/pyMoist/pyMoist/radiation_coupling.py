from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField
from gt4py.cartesian.gtscript import computation, interval, PARALLEL


def _fix_up_clouds_stencil(Q: FloatField):
    # Stencil code goes here
    with computation(PARALLEL), interval(...):
        Q = 0


def _radcouple_stencil(Q: FloatField):
    # Stencil code goes here
    with computation(PARALLEL), interval(...):
        Q = 0


class RadiationCoupling:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        do_qa: bool,
    ) -> None:
        self._fix_up_clouds = stencil_factory.from_dims_halo(
            func=_fix_up_clouds_stencil, compute_dims=[X_DIM, Y_DIM, Z_DIM]
        )
        self._radcouple = stencil_factory.from_dims_halo(
            func=_radcouple_stencil, compute_dims=[X_DIM, Y_DIM, Z_DIM]
        )
        self.do_qa = do_qa

    def __call__(
        self,
        # QILS,
        # QRAIN,
        # T,
        # MAX_RL,
        # QSNOW,
        # QGRAUPEL,
        # FAC_RL,
        # MIN_RL,
        # CLLS,
        # NACTL,
        # MIN_RI,
        # NACTI,
        # CLCN,
        # QICN,
        # MAX_RI,
        # FAC_RI,
        # PLmb,
        # QLCN,
        Q,
        # QLLS,
    ) -> None:
        self._fix_up_clouds(Q)
        self._radcouple(Q)
        if self.do_qa:
            pass
