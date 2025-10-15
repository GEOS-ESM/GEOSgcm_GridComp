from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import (
    Float,
    FloatField,
    Int,
    FloatFieldIJ,
    IntFieldIJ,
    IntField,
)
import gt4py.cartesian.gtscript as gtscript

@gtscript.function
def cup_minimi(
    array: Float,
    ks: Int,
    kend: Int,
    ierr: Int,
):

    kt = ks

    if ierr == 0:
        x = array.at(K=ks)
        kstop = max(ks+1,kend)
        
        idx=ks+1
        while idx >= ks+1 and idx <= kstop:
            if array.at(K=idx) < x:
                x = array.at(K=idx)
                kt=idx
            idx+=1

    return kt

def test_cup_minimi(
    # In
    array: FloatField,
    ks: IntField,
    kend: IntField,
    ierr: IntField,
    # Out
    kt: IntField,
):
    with computation(PARALLEL), interval(...):
        kt = cup_minimi(array, ks-1, kend-1, ierr)


class CupMinimi:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
    ) -> None:

     
        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory

        self._test_cup_minimi = self.stencil_factory.from_dims_halo(
            func=test_cup_minimi,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        # In
        array: FloatField,
        ks: IntField,
        kend: IntField,
        ierr: IntField,
        # Out
        kt: IntField,
    ):

        self._test_cup_minimi(
            # In
            array=array,
            ks=ks,
            kend=kend,
            ierr=ierr,
            # Out
            kt=kt,
        )

           
