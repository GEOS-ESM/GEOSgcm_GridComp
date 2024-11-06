from ndsl import Namelist,StencilFactory
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from gt4py.cartesian.gtscript import computation, PARALLEL, interval
import gt4py.cartesian.gtscript as gtscript
from ndsl.dsl.typing import (
    Float,
    FloatField
)
import pyMoist.constants as constants
from pyMoist.UW.compute_mumin2 import compute_mumin2

def harness_stencil(
    mulcl: FloatField,
    rmaxfrax: FloatField,
    mulow: FloatField,
    mumin2: FloatField,
    ):

    with computation(PARALLEL), interval(...):
        mumin2 = compute_mumin2(mulcl, rmaxfrax, mulow)


class TranslateComputeMumin2(TranslateFortranData2Py):
    def __init__(self, grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory
        self._grid = grid

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "mulcl": {}, 
            "rmaxfrax": {},
            "mu": {},
            }

        # Float Inputs
        self.in_vars["parameters"] = []

        # FloatField Outputs
        self.out_vars = {
            "mumin2": (),
        }

        self.stencil = stencil_factory.from_origin_domain(
            func=harness_stencil,
            origin=(0, 0, 0),
            domain=(1, 1, 10),
        )

    def compute(self, inputs):
        mulcl = self.quantity_factory._numpy.empty((1, 1, 10), dtype=Float)
        mulcl[0, 0, :] = inputs["mulcl"]
        rmaxfrax = self.quantity_factory._numpy.empty((1, 1, 10), dtype=Float)
        rmaxfrax[0, 0, :] = inputs["rmaxfrax"]
        mulow = self.quantity_factory._numpy.empty((1, 1, 10), dtype=Float)
        mulow[0, 0, :] = inputs["mu"]

        mumin2 = self.quantity_factory._numpy.empty((1, 1, 10), dtype=Float)

        self.stencil(mulcl, rmaxfrax, mulow, mumin2)

        return {"mumin2": mumin2[0, 0, :]}