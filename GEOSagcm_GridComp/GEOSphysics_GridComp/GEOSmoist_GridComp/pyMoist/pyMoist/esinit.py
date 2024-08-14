import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import computation, interval, PARALLEL, FORWARD, atan, sin, tan, sqrt, tanh, exp, log10
from ndsl.boilerplate import get_factories_single_tile_numpy
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField, FloatFieldIJ, Float, IntField, IntFieldIJ
from ndsl import StencilFactory, QuantityFactory, orchestrate
import numpy as np
import pyMoist.pyMoist_constants as constants
import fileinput
import sys



def replaceAll(file,searchExp,replaceExp):
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp,replaceExp)
        sys.stdout.write(line)

class esinit:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
    ):
        orchestrate(obj=self, config=stencil_factory.config.dace_config)
        self.ESTBLE = np.zeros(constants.TABLESIZE)
        self.ESTBLW = np.zeros(constants.TABLESIZE)
        self.ESTBLX = np.zeros(constants.TABLESIZE)

        # need to turn this into a stencil
        # likely doable, but unsure how to create non-standard 1D stencil dimentions
        for i in range(len(self.ESTBLE)):
            T = i * constants.DELTA_T + constants.TMINTBL
            self.ESTBLW(i) = 
        
        


if __name__ == "__main__":
    domain = (3, 3, 4)

    stcil_fctry, ijk_qty_fctry = get_factories_single_tile_numpy(
        domain[0], domain[1], domain[2], 0
    )

    code = esinit(stcil_fctry, ijk_qty_fctry)
    code()