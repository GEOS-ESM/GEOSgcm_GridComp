from ndsl import QuantityFactory, StencilFactory  # type: ignore
from ndsl.constants import X_DIM, Y_DIM, Z_DIM  # type: ignore
import ndsl.constants as constants
from ndsl.comm.communicator import Communicator  # type: ignore
from ndsl.dsl.typing import FloatField, Int # type: ignore
import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import computation, interval, PARALLEL, BACKWARD, FORWARD, log10, exp  # type: ignore
import numpy as np 

@gtscript.function
def add_one(field: FloatField) -> FloatField:
    from __externals__ import UNDEF_CST
    return field + UNDEF_CST

def copy(in_field: FloatField, out_field: FloatField):
    out_field = in_field + add_one(in_field)
    print(out_field)

copy_stencil = stencil_factory.from_dims_halo(
    func=copy,
    compute_dims=[X_DIM, Y_DIM, Z_DIM],
    externals = {
        "UNDEF_CST": 5
    }
)