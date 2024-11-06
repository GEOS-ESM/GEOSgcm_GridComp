from ndsl import Namelist, Quantity
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from gt4py.cartesian.gtscript import computation, PARALLEL, interval
import gt4py.cartesian.gtscript as gtscript
from ndsl.dsl.typing import (
    FloatFieldIJ,
    Float,
    FloatField
)
from ndsl import StencilFactory, QuantityFactory
from pyMoist.UW.compute_ppen import compute_ppen

def harness_stencil(
    wtwb: FloatField,
    D: FloatField,
    bogbot: FloatField,
    bogtop: FloatField,
    rh0j: FloatField,
    dpen: FloatField,
    ppen: FloatField,
    ):

    with computation(PARALLEL), interval(...):
        ppen = compute_ppen(wtwb, D, bogbot, bogtop, rh0j, dpen)


class TranslateComputePpen(TranslateFortranData2Py):
    def __init__(self, grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory
        self._grid = grid

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "wtwb": {}, 
            "drage": {},
            "bogbot": {},
            "bogtop": {},
            "rhomid0j": {},
            "dp0": {},
            }

        # Float Inputs
        self.in_vars["parameters"] = []

        # FloatField Outputs
        self.out_vars = {
            "ppen": (),
        }

        self.stencil = stencil_factory.from_origin_domain(
            func=harness_stencil,
            origin=(0, 0, 0),
            domain=(1, 1, 120),
        )

    def compute(self, inputs):
        wtwb = self.quantity_factory._numpy.empty((1, 1, 120), dtype=Float)
        wtwb[0, 0, :] = inputs["wtwb"]
        drage = self.quantity_factory._numpy.empty((1, 1, 120), dtype=Float)
        drage[0, 0, :] = inputs["drage"]
        bogbot = self.quantity_factory._numpy.empty((1, 1, 120), dtype=Float)
        bogbot[0, 0, :] = inputs["bogbot"]
        bogtop = self.quantity_factory._numpy.empty((1, 1, 120), dtype=Float)
        bogtop[0, 0, :] = inputs["bogtop"]
        rhomid0j = self.quantity_factory._numpy.empty((1, 1, 120), dtype=Float)
        rhomid0j[0, 0, :] = inputs["rhomid0j"]
        dp0 = self.quantity_factory._numpy.empty((1, 1, 120), dtype=Float)
        dp0[0, 0, :] = inputs["dp0"]


        ppen = self.quantity_factory._numpy.empty((1, 1, 120), dtype=Float)

        self.stencil(wtwb, drage, bogbot, bogtop, rhomid0j, dp0, ppen)

        return {"ppen": ppen[0, 0, :]}