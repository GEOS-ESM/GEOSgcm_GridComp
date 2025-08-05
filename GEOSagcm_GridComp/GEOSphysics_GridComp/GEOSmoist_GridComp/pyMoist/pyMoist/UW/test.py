from ndsl.dsl.gt4py import (
    computation,
    interval,
    PARALLEL,
)
from ndsl.boilerplate import get_factories_single_tile
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField, Float, IntField, Int
from ndsl import StencilFactory, orchestrate
import numpy as np
from gt4py.cartesian.gtscript import K


domain = (2, 2, 5)
backend = "dace:cpu"

stcil_fctry, qty_factry = get_factories_single_tile(
    domain[0], domain[1], domain[2], 0, backend
)


def the_stencil(in_field: IntField, out_field: FloatField):
    with computation(PARALLEL), interval(...):
        out_field = in_field.at(K=K + 1)


class Code:
    def __init__(self, stencil_factory: StencilFactory):
        orchestrate(obj=self, config=stencil_factory.config.dace_config)
        self._the_stencil = stcil_fctry.from_dims_halo(
            func=the_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(self, in_field: IntField, out_field: FloatField):
        self._the_stencil(in_field, out_field)


if __name__ == "__main__":
    in_arr = qty_factry.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="inputs", dtype=Int)
    sz = domain[0] * domain[1] * domain[2]
    in_arr.view[:] = np.arange(sz, dtype=Int).reshape(domain)

    out_arr = qty_factry.zeros([X_DIM, Y_DIM, Z_DIM], units="outputs")

    code = Code(stcil_fctry)
    code(in_arr, out_arr)

    print("Done 🚀")
    print(out_arr.field)
