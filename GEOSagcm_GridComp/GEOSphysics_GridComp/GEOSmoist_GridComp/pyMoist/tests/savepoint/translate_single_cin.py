from ndsl import Namelist, StencilFactory
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.saturation import QSat
import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import computation, PARALLEL, interval
import gt4py.cartesian.gtscript as gtscript
import pyMoist.constants as constants
from ndsl.dsl.typing import (
    Float,
    FloatField,
)
from pyMoist.UW.single_cin import single_cin
from pyMoist.UW.compute_uwshcu import exnerfn


def harness_stencil(
    pbot: FloatField,
    thv0bot: FloatField,
    ptop: FloatField,
    thv0top: FloatField,
    thvubot: FloatField,
    thvutop: FloatField,
    cin1: FloatField,
    cin2: FloatField,
    ):

    with computation(PARALLEL), interval(...):
        cin2 = cin1 + single_cin(pbot, thv0bot,ptop,thv0top,thvubot,thvutop)


class TranslateSingleCin(TranslateFortranData2Py):
    def __init__(self, grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory
        self._grid = grid

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "pifc0_bot": {}, 
            "thv0bot": {}, 
            "pifc0_top": {},
            "thv0top": {},
            "thvubot": {},
            "thvutop": {},
            "cin1": {},
            }

        # Float Inputs
        self.in_vars["parameters"] = []

        # FloatField Outputs
        self.out_vars = {
            "cin2": (),
        }

        self.stencil = stencil_factory.from_origin_domain(
            func=harness_stencil,
            origin=(0, 0, 0),
            domain=(1, 1, 54),
        )

    def compute(self, inputs):
        pbot = self.quantity_factory._numpy.empty((1, 1, 54), dtype=Float)
        pbot[0, 0, :] = inputs["pifc0_bot"]
        thv0bot = self.quantity_factory._numpy.empty((1, 1, 54), dtype=Float)
        thv0bot[0, 0, :] = inputs["thv0bot"]
        ptop = self.quantity_factory._numpy.empty((1, 1, 54), dtype=Float)
        ptop[0, 0, :] = inputs["pifc0_top"]
        thv0top = self.quantity_factory._numpy.empty((1, 1, 54), dtype=Float)
        thv0top[0, 0, :] = inputs["thv0top"]
        thvubot = self.quantity_factory._numpy.empty((1, 1, 54), dtype=Float)
        thvubot[0, 0, :] = inputs["thvubot"]
        thvutop = self.quantity_factory._numpy.empty((1, 1, 54), dtype=Float)
        thvutop[0, 0, :] = inputs["thvutop"]
        cin1 = self.quantity_factory._numpy.empty((1, 1, 54), dtype=Float)
        cin1[0, 0, :] = inputs["cin1"]


        cin2 = self.quantity_factory._numpy.empty((1, 1, 54), dtype=Float)

        self.stencil(pbot, thv0bot, ptop, thv0top, thvubot, thvutop, cin1, cin2)


        return {
            "cin2": cin2[0, 0, :],
        }