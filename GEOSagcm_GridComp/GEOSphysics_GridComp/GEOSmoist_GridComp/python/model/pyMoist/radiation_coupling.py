from ndsl import (
    StencilFactory,
    QuantityFactory
)
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from gt4py.cartesian.gtscript import (
    PARALLEL,
    interval    
)

def _fix_up_clouds_stencil():
    pass

def _radcouple_stencil():
    pass

class RadiationCoupling:
    def __init__(
            self,
            stencil_factory: StencilFactory,
            quantity_factory: QuantityFactory,
            do_qa: bool
    ) -> None:
        self._fix_up_clouds = stencil_factory.from_dims_halo(
            func=_fix_up_clouds_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM]
        )
        self._radcouple = stencil_factory.from_dims_halo(
            func=_radcouple_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM]
        )
        self.do_qa = do_qa

    def __call__(self) -> None:
        self._fix_up_clouds()
        self._radcouple()
        if self.do_qa:
            pass
