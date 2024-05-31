from ndsl import Namelist, StencilFactory
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.radiation_coupling import RadiationCoupling


class TranslateRadCouple(TranslateFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.compute_func = RadiationCoupling(  # type: ignore
            self.stencil_factory,
            self.grid.quantity_factory,
            do_qa=namelist.do_qa,
        )

        # ADAPT BELOW TO INPUTS
        #
        # fillq_info = self.grid.compute_dict()
        # fillq_info["serialname"] = "fq"
        self.in_vars["data_vars"] = {
            "Q": self.grid.compute_dict(),
        }
        self.out_vars = {
            "Q": self.grid.compute_dict(),
        }

    def compute_from_storage(self, inputs):
        self.compute_func(**inputs)
        return inputs
