from ndsl import Namelist,StencilFactory
from ndsl.dsl.typing import FloatField, Float
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from gt4py.cartesian.gtscript import computation, PARALLEL, interval
import gt4py.cartesian.gtscript as gtscript
from ndsl.dsl.typing import (
    Float,
    FloatField
)
from pyMoist.UW.getbuoy import getbuoy


def harness_stencil(
    pbot: FloatField,
    thv0bot: FloatField,
    ptop: FloatField,
    thv0top: FloatField,
    thvubot: FloatField,
    thvutop: FloatField,
    cin_in: FloatField,
    plfc_in: FloatField,
    plfc: FloatField,
    cin: FloatField,
    ):

    with computation(PARALLEL), interval(...):
        if thvubot > thv0bot and thvutop > thv0top:
            plfc, _ = getbuoy(pbot, thv0bot, ptop, thv0top, thvubot,thvutop,cin_in,plfc_in)

        elif thvubot < thv0bot and thvutop < thv0top:
            _, cin = getbuoy(pbot, thv0bot, ptop, thv0top, thvubot,thvutop,cin_in,plfc_in)

        elif thvubot > thv0bot and thvutop < thv0top:
            _, cin = getbuoy(pbot, thv0bot, ptop, thv0top, thvubot,thvutop,cin_in,plfc_in)

        else:
            plfc, cin = getbuoy(pbot, thv0bot, ptop, thv0top, thvubot,thvutop,cin_in,plfc_in)
        
        


class TranslateGetbuoy(TranslateFortranData2Py):
    def __init__(self, grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory
        self._grid = grid

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "pifc0_below": {}, 
            "thv0bot_buoy": {},
            "pifc0_above": {},
            "thv0top_buoy": {},
            "thvubot_buoy": {},
            "thvutop_buoy": {},
            "cin_in": {},
            "plfc_in": {}
            }

        # Float Inputs
        self.in_vars["parameters"] = []

        # FloatField Outputs
        self.out_vars = {
            "plfc_out": (),
            "cin_out": (),
        }

        self.stencil = stencil_factory.from_origin_domain(
            func=harness_stencil,
            origin=(0, 0, 0),
            domain=(1, 1, 94),
        )

    def compute(self, inputs):
        pifc0_below = self.quantity_factory._numpy.empty((1, 1, 94), dtype=Float)
        pifc0_below[0, 0, :] = inputs["pifc0_below"]
        thv0bot = self.quantity_factory._numpy.empty((1, 1, 94), dtype=Float)
        thv0bot[0, 0, :] = inputs["thv0bot_buoy"]
        pifc0_above = self.quantity_factory._numpy.empty((1, 1, 94), dtype=Float)
        pifc0_above[0, 0, :] = inputs["pifc0_above"]
        thv0top = self.quantity_factory._numpy.empty((1, 1, 94), dtype=Float)
        thv0top[0, 0, :] = inputs["thv0top_buoy"]
        thvubot = self.quantity_factory._numpy.empty((1, 1, 94), dtype=Float)
        thvubot[0, 0, :] = inputs["thvubot_buoy"]
        thvutop = self.quantity_factory._numpy.empty((1, 1, 94), dtype=Float)
        thvutop[0, 0, :] = inputs["thvutop_buoy"]
        cin_in = self.quantity_factory._numpy.empty((1, 1, 94), dtype=Float)
        cin_in[0, 0, :] = inputs["cin_in"]
        plfc_in = self.quantity_factory._numpy.empty((1, 1, 94), dtype=Float)
        plfc_in[0, 0, :] = inputs["plfc_in"]
        


        plfc = self.quantity_factory._numpy.empty((1, 1, 94), dtype=Float)
        cin = self.quantity_factory._numpy.empty((1, 1, 94), dtype=Float)

        self.stencil(pifc0_below, thv0bot, pifc0_above, thv0top, thvubot, thvutop, cin_in, plfc_in,plfc, cin)


        return {
            "plfc_out": plfc[0, 0, :],
            "cin_out": cin[0, 0, :],
        }
