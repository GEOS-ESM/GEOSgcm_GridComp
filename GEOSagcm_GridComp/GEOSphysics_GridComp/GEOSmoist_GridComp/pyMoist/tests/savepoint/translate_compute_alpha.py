from ndsl import Namelist, StencilFactory, Quantity
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.UW.compute_alpha import Compute_Alpha
import numpy as np
from ndsl.dsl.typing import (
    FloatFieldIJ,
    FloatField,
)


class TranslateComputeAlpha(TranslateFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.quantity_factory = grid.quantity_factory
        self.compute_func = Compute_Alpha(
            self.stencil_factory,
            self.grid.quantity_factory,
        )

        self.max_error = 1e-9

        print(grid.__dict__)

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "del_CIN": {},
            "ke": {},
        }
      
        # Float/Int Inputs
        self.in_vars["parameters"] = [
        ]

        # FloatField Outputs
        self.out_vars = {
            "alpha": self.grid.compute_dict(),
        }

    def reshape_before(self, inputs):
        # Reshape input fields to the necessary shape
        i, j, k = self.grid.nic, self.grid.njc, self.grid.npz
        reshaped_inputs = {}
        for key, array in inputs.items():
            reshaped_inputs[key] = np.reshape(array[:,1,0,0], (i, j)).astype(np.float32)
      
        return reshaped_inputs


    def reshape_after(self, outputs):
        # Reshape output fields back to original shape
        i, j, k = self.grid.nic, self.grid.njc, self.grid.npz
        reshaped_outputs = np.reshape(outputs, (i * j)).astype(np.float64)

        return reshaped_outputs

        
    def make_ij_field(self, data, dtype=FloatField) -> Quantity:
        qty = self.quantity_factory.empty([X_DIM, Y_DIM], "n/a", dtype=dtype)
        qty.view[:, :] = qty.np.asarray(data[:, :])
        return qty

    # Calculated Outputs
    def compute(self, inputs):
        compute_alpha = Compute_Alpha(
            self.stencil_factory,
            self.grid.quantity_factory,
        )

        inputs_reshaped = self.reshape_before(inputs)
        print("Reshaped inputs...")

        # Inputs
        del_CIN = self.make_ij_field(inputs_reshaped["del_CIN"], dtype=FloatFieldIJ)
        ke = self.make_ij_field(inputs_reshaped["ke"], dtype=FloatFieldIJ)

        # Outputs
        alpha = self.make_ij_field(inputs_reshaped["ke"], dtype=FloatFieldIJ)

        
        compute_alpha(
            del_CIN=del_CIN,
            ke=ke,
            alpha=alpha
        )
        
        print("Executed stencil computation...")

        # Reshape output variables back to original shape
        alpha_out = self.reshape_after(alpha.view[:,:])
        print("Reshaped outputs back to original shape...")

        alpha_4D = inputs["del_CIN"]
        alpha_4D[:,1,0,0] = alpha_out

        return {"alpha": alpha_4D
            }