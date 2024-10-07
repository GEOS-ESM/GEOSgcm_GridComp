from ndsl import Namelist, StencilFactory
from ndsl.stencils.testing.translate import (
    TranslateFortranData2Py,
)

from pyMoist.redistribute_clouds import RedistributeClouds


class TranslateRedistributeClouds(TranslateFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.compute_func = RedistributeClouds(
            self.stencil_factory,
        )

        self.max_error = 1e-9

        print(grid.__dict__)
        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "RAD_CF": {},
            "RAD_QL": {},
            "RAD_QI": {},
            "CLCN": {},
            "CLLS": {},
            "QLCN": {},
            "QLLS": {},
            "QICN": {},
            "QILS": {},
            "RAD_QV": {},
            "T": {},
        }

        # Float Inputs
        # self.in_vars["parameters"] = ["alhlbcp", "alhsbcp"]

        # FloatField Outputs
        self.out_vars = {
            "RAD_CF": self.grid.compute_dict(),
            "RAD_QL": self.grid.compute_dict(),
            "RAD_QI": self.grid.compute_dict(),
            "CLCN": self.grid.compute_dict(),
            "CLLS": self.grid.compute_dict(),
            "QLCN": self.grid.compute_dict(),
            "QLLS": self.grid.compute_dict(),
            "QICN": self.grid.compute_dict(),
            "QILS": self.grid.compute_dict(),
            "RAD_QV": self.grid.compute_dict(),
            "T": self.grid.compute_dict(),
        }

    # Calculated Outputs
    def compute_from_storage(self, inputs):
        self.compute_func(**inputs)
        return {
            "RAD_CF": inputs["RAD_CF"],
            "RAD_QL": inputs["RAD_QL"],
            "RAD_QI": inputs["RAD_QI"],
            "CLCN": inputs["CLCN"],
            "CLLS": inputs["CLLS"],
            "QLCN": inputs["QLCN"],
            "QLLS": inputs["QLLS"],
            "QICN": inputs["QICN"],
            "QILS": inputs["QILS"],
            "RAD_QV": inputs["RAD_QV"],
            "T": inputs["T"],
        }
