from ndsl import Namelist,StencilFactory
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from gt4py.cartesian.gtscript import computation, PARALLEL, interval
import gt4py.cartesian.gtscript as gtscript
from ndsl.dsl.typing import (
    Float,
    FloatField
)
from pyMoist.UW.compute_alpha import compute_alpha


def harness_stencil(
    del_CIN: FloatField,
    ke: FloatField,
    alpha: FloatField,
):
    with computation(PARALLEL), interval(...):
        alpha = compute_alpha(del_CIN, ke)


class TranslateComputeAlpha(TranslateFortranData2Py):
    def __init__(self, grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory
        self._grid = grid

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "del_CIN": {},
            "ke": {},
        }

        # Float/Int Inputs
        self.in_vars["parameters"] = []

        # FloatField Outputs
        self.out_vars = {
            "alpha": (),
        }

    def compute(self, inputs):
        domain = (1, 1, inputs["del_CIN"].shape[0])
        stencil = self.stencil_factory.from_origin_domain(
            func=harness_stencil,
            origin=(0, 0, 0),
            domain=domain,
        )

        del_CIN = self.quantity_factory._numpy.empty(domain, dtype=Float)
        del_CIN[0, 0, :] = inputs["del_CIN"]
        ke = self.quantity_factory._numpy.empty(domain, dtype=Float)
        ke[0, 0, :] = inputs["ke"]

        alpha = self.quantity_factory._numpy.empty(domain, dtype=Float)

        stencil(del_CIN, ke, alpha)

        return {"alpha": alpha[0, 0, :]}