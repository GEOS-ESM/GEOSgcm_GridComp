import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import (
    computation,
    interval,
    PARALLEL
)
from ndsl.boilerplate import get_factories_single_tile_numpy
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField 
from ndsl import Quantity, QuantityFactory, StencilFactory, orchestrate

stcil_fctry, qty_fctry = get_factories_single_tile_numpy(5, 5, 1, 1)

out_q = qty_fctry.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
in_q = qty_fctry.ones([X_DIM, Y_DIM, Z_DIM], "n/a")

@gtscript.function
def add_one(field: FloatField) -> FloatField:
    from __externals__ import UNDEF_CST
    return field + UNDEF_CST

def copy(in_field: FloatField, out_field: FloatField):
    with computation(PARALLEL), interval(...):
        out_field = in_field + add_one(in_field)

copy_stcil = stcil_fctry.from_dims_halo(
    func=copy,
    compute_dims=[X_DIM, Y_DIM, Z_DIM],
    externals={
         "UNDEF_CST": 5
    }
)

copy_stcil(in_q, out_q)

print(in_q)
print(out_q)