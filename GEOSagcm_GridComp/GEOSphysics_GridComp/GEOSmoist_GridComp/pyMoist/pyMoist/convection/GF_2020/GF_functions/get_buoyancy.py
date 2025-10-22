from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import (
    Float,
    FloatField,
    Int,
    IntField,
)
import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import (
    PARALLEL,
    K,
    computation,
    interval,
)


@gtscript.function
def get_buoyancy(
    hc: Float,
    he_cup: Float,
    hes_cup: Float,
    ierr: Int,
    kbcon: Int,
    klcl: Int,
    ktop: Int,
):

    dby = 0.0
    if ierr == 0:
        if K <= klcl:
            dby = hc - he_cup

        if K >= klcl + 1 and K <= ktop + 1:
            dby = hc - hes_cup

    return dby


def test_get_buoyancy(
    # In
    hc: FloatField,
    he_cup: FloatField,
    hes_cup: FloatField,
    ierr: IntField,
    kbcon: IntField,
    klcl: IntField,
    ktop: IntField,
    # Out
    dby: FloatField,
):

    with computation(PARALLEL), interval(...):
        # Subtract 1 from k-level arrays to account for Python starting at zero
        dby = get_buoyancy(hc, he_cup, hes_cup, ierr, kbcon - 1, klcl - 1, ktop - 1)


class GetBuoyancy:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
    ) -> None:

        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory

        self._test_get_buoyancy = self.stencil_factory.from_dims_halo(
            func=test_get_buoyancy,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        # In
        hc: FloatField,
        he_cup: FloatField,
        hes_cup: FloatField,
        ierr: IntField,
        kbcon: IntField,
        klcl: IntField,
        ktop: IntField,
        # Out
        dby: FloatField,
    ):

        self._test_get_buoyancy(
            # In
            hc=hc,
            he_cup=he_cup,
            hes_cup=hes_cup,
            ierr=ierr,
            kbcon=kbcon,
            klcl=klcl,
            ktop=ktop,
            # Out
            dby=dby,
        )
