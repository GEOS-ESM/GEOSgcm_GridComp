from ndsl import Namelist, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.UW.roots import Get_Roots
import numpy as np


class TranslateRoots(TranslateFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.compute_func = Get_Roots(
            self.stencil_factory,
            self.grid.quantity_factory,
        )

        self.max_error = 1e-9

        print(grid.__dict__)

        # FloatField Inputs
        self.in_vars["data_vars"] = {
        }

      
        # Float/Int Inputs
        self.in_vars["parameters"] = [
            "aquad",
            "bquad",
            "cquad",
        ]

        # FloatField Outputs
        self.out_vars = {
            "xc1": self.grid.compute_dict(),
            "xc2": self.grid.compute_dict(),
            "status": self.grid.compute_dict(),
        }
        

    # Calculated Outputs
    def compute(self, inputs):
        roots = Get_Roots(
            self.stencil_factory,
            self.grid.quantity_factory,
        )

        aquad = inputs["aquad"]
        bquad = inputs["bquad"]
        cquad = inputs["cquad"]

        xc1 = inputs["aquad"]
        xc2 = inputs["aquad"]
        status = 0
     
        xc1, xc2, status = roots(
            aquad=aquad,
            bquad=bquad,
            cquad=cquad,
            xc1=xc1,
            xc2=xc2,
            status=status,
        )

        return {
            "xc1": xc1,
            "xc2": xc2,
            "status": status,
            }