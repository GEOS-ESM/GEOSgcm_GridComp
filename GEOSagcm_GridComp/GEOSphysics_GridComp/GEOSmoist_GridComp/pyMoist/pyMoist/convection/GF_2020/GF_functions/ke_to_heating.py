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
import pyMoist.convection.shared_constants as constants
from gt4py.cartesian.gtscript import (
    sqrt,
)


def ke_to_heating(
    # In
    dellu: FloatField,
    dellv: FloatField,
    ierr: IntField,
    ktop: IntField,
    po_cup: FloatField,
    us: FloatField,
    vs: FloatField,
    dts: FloatFieldIJ,
    fpi: FloatFieldIJ,
    # Out
    dellat: FloatField,
):
    with computation(FORWARD), interval(...):
        if ierr == 0:
            dts=0.
            fpi=0.
    with computation(FORWARD), interval(...):
        if ierr == 0:
            if K <= ktop-1:
                dp=(po_cup-po_cup[0,0,1])*100.

                dts = dts - (((dellu*us)+(dellv*vs))*dp/constants.MAPL_GRAV)
        
                fpi = fpi + (sqrt((dellu*dellu) + (dellv*dellv))*dp)

    with computation(FORWARD), interval(...):
        if ierr == 0:
            if fpi > 0.:
                if K <= ktop-1:
                    fp= sqrt((dellu*dellu+dellv*dellv))/fpi

                    dellat = dellat + (fp*dts*constants.MAPL_GRAV/constants.CP)

    

class KeToHeating:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
    ) -> None:

     
        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory

        self._ke_to_heating = self.stencil_factory.from_dims_halo(
            func=ke_to_heating,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self.dts = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM], units="n/a"
        )
        self.fpi = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM], units="n/a"
        )
    

    def __call__(
        self,
        # In
        dellu: FloatField,
        dellv: FloatField,
        ierr: IntField,
        ktop: IntField,
        po_cup: FloatField,
        us: FloatField,
        vs: FloatField,
        # Out
        dellat: FloatField,
    ):
     
        self._ke_to_heating(
            # In
            dellu=dellu,
            dellv=dellv,
            ierr=ierr,
            ktop=ktop,
            po_cup=po_cup,
            us=us,
            vs=vs,
            dts=self.dts,
            fpi=self.fpi,
            # Out
            dellat=dellat,
        )

           
