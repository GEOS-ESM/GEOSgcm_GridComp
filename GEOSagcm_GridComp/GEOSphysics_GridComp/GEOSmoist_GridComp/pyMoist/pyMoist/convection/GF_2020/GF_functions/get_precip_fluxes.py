from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import (
    FloatField,
    IntField,
)

from gt4py.cartesian.gtscript import (
    BACKWARD,
    K,
    computation,
    interval,
)


def get_precip_fluxes(
    # In
    edto: FloatField,
    ierr: IntField,
    ktop: IntField,
    pwdo: FloatField,
    pwo: FloatField,
    xmb: FloatField,
    # Out
    prec_flx: FloatField,
    evap_flx: FloatField,
):
    with computation(BACKWARD), interval(...):
        if ierr == 0:
            if K < ktop:
                prec_flx = prec_flx[0, 0, 1] + xmb * (pwo + edto * pwdo)
                prec_flx = max(0.0, prec_flx)

                evap_flx = evap_flx[0, 0, 1] - xmb * edto * pwdo
                evap_flx = max(0.0, evap_flx)


class GetPrecipFluxes:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
    ) -> None:

        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory

        self._get_precip_fluxes = self.stencil_factory.from_dims_halo(
            func=get_precip_fluxes,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        # In
        edto: FloatField,
        ierr: IntField,
        ktop: IntField,
        pwdo: FloatField,
        pwo: FloatField,
        xmb: FloatField,
        # Out
        prec_flx: FloatField,
        evap_flx: FloatField,
    ):

        self._get_precip_fluxes(
            # In
            edto=edto,
            ierr=ierr,
            ktop=ktop,
            pwdo=pwdo,
            pwo=pwo,
            xmb=xmb,
            # Out
            prec_flx=prec_flx,
            evap_flx=evap_flx,
        )
