from ndsl import Namelist, Quantity, StencilFactory
from ndsl.dsl.typing import FloatField, Float
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.saturation import QSat
import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import computation, PARALLEL, interval, exp, FORWARD, sqrt
import gt4py.cartesian.gtscript as gtscript
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import (
    FloatFieldIJ,
    Float,
    FloatField,
    IntField,
    Int,
)
from ndsl import StencilFactory, QuantityFactory
from pyMoist.UW.roots import roots,sign

def harness_stencil(
    aquad: FloatField,
    bquad: FloatField,
    cquad: FloatField,
    xc1: FloatField,
    xc2: FloatField,
    status: IntField,
    ):

    with computation(PARALLEL), interval(...):
        xc1, xc2, status = roots(aquad,bquad,cquad)


class TranslateRoots(TranslateFortranData2Py):
    def __init__(self, grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory
        self._grid = grid

        # FloatField Inputs
        self.in_vars["data_vars"] = {"aquad": {}, "bquad": {}, "cquad": {}}

        # Float/Int Inputs
        self.in_vars["parameters"] = []

        # FloatField Outputs
        self.out_vars = {
            "xc1": (),
            "xc2": (),
            #"status": (),
        }

    def compute(self, inputs):
        domain = (1, 1, inputs["aquad"].shape[0])
        stencil = self.stencil_factory.from_origin_domain(
            func=harness_stencil,
            origin=(0, 0, 0),
            domain=domain,
        )

        aquad = self.quantity_factory._numpy.empty(domain, dtype=Float)
        aquad[0, 0, :] = inputs["aquad"]
        bquad = self.quantity_factory._numpy.empty(domain, dtype=Float)
        bquad[0, 0, :] = inputs["bquad"]
        cquad = self.quantity_factory._numpy.empty(domain, dtype=Float)
        cquad[0, 0, :] = inputs["cquad"]

        xc1 = self.quantity_factory._numpy.empty(domain, dtype=Float)
        xc2 = self.quantity_factory._numpy.empty(domain, dtype=Float)
        status = self.quantity_factory._numpy.empty(domain, dtype=Int)

        stencil(aquad, bquad, cquad, xc1 ,xc2, status)

        return {"xc1": xc1[0, 0, :],
                "xc2": xc2[0, 0, :]}
                #"status": status[0, 0, :]}