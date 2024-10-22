from ndsl import Namelist, StencilFactory, Quantity
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.UW.qsinvert import QsInvert
import numpy as np
from ndsl.dsl.typing import (
    FloatFieldIJ,
    FloatField,
)


class TranslateQsInvert(TranslateFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.quantity_factory = grid.quantity_factory
        self.compute_func = QsInvert(
            self.stencil_factory,
            self.grid.quantity_factory,
        )

        self.max_error = 1e-9

        print(grid.__dict__)

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "qtsrc": {},
            "thlsrc": {},
            "pifc0": {},
        }
      
        # Float/Int Inputs
        self.in_vars["parameters"] = [
        ]

        # FloatField Outputs
        self.out_vars = {
            "plcl": self.grid.compute_dict(),
        }

    def reshape_before(self, inputs):
        # Reshape input fields to the necessary shape
        i, j, k = self.grid.nic, self.grid.njc, self.grid.npz
        reshaped_inputs = {}
        for key, array in inputs.items():
            reshaped_inputs[key] = np.reshape(array[:,0,0,0], (i, j)).astype(np.float32)
      
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
        qsinvert = QsInvert(
            self.stencil_factory,
            self.grid.quantity_factory,
        )

        inputs_reshaped = self.reshape_before(inputs)
        print("Reshaped inputs...")

        # Inputs
        qtsrc = self.make_ij_field(inputs_reshaped["qtsrc"], dtype=FloatFieldIJ)
        thlsrc = self.make_ij_field(inputs_reshaped["thlsrc"], dtype=FloatFieldIJ)
        pifc0 = self.make_ij_field(inputs_reshaped["pifc0"], dtype=FloatFieldIJ)

        # Outputs
        plcl = self.make_ij_field(inputs_reshaped["pifc0"], dtype=FloatFieldIJ)

        
        qsinvert(
            qtsrc=qtsrc,
            thlsrc=thlsrc,
            pifc0=pifc0,
            plcl=plcl
        )
        
        print("Executed stencil computation...")

        # Reshape output variables back to original shape
        plcl_out = self.reshape_after(plcl.view[:,:])
        print("Reshaped outputs back to original shape...")

        plcl_4D = inputs["pifc0"]
        plcl_4D[:,0,0,0] = plcl_out

        return {"plcl": plcl_4D
            }