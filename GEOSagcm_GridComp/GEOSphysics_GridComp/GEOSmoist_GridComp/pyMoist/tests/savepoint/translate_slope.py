from ndsl import Namelist, StencilFactory, QuantityFactory 
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.slope import Slope
import numpy as np


class TranslateSlope(TranslateFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.compute_func = Slope(
            self.stencil_factory,
            self.grid.quantity_factory,
        )

        self.max_error = 1e-9

        print(grid.__dict__)

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "thl0_test": {},
            "pmid0_test": {},
        }

        print(self.in_vars["data_vars"])

        # FloatField Outputs
        self.out_vars = {
            "ssthl0_test": self.grid.compute_dict(),
        }

    def reshape_before(self,inputs):
        # Reshape input fields from (i*j, k) to (i, j, k)
        i, j, k = self.grid.nic, self.grid.njc, self.grid.npz 
        reshaped_inputs = {key: np.reshape(inputs[key], (i, j, k)).astype(np.float32) for key in inputs}
    
        # Add halo of 3 in i and j dimensions
        halo_arrays = {}
        for key, array in reshaped_inputs.items():
            new_array = np.zeros((i+6, j+6, k), dtype=np.float32)
            new_array[3:-3, 3:-3, :] = array

            halo_arrays[key] = new_array

        return halo_arrays

    def reshape_after(self,outputs):
        # Reshape output fields from (i, j, k) back to (i*j, k)
        i, j, k = self.grid.nic, self.grid.njc, self.grid.npz 
        reshaped_outputs = [np.reshape(outputs[key], (i, j, k)).astype(np.float32) for key in outputs]
        return reshaped_outputs
        

    # Calculated Outputs
    def compute(self, inputs):

        inputs_2D = self.reshape_before(inputs)
        print("Reshaped 2D fields to 3D")

        outputs = {
            "ssthl0_test": self.grid.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="unknown"),
        }
     
        self.compute_func(**inputs_2D,**outputs)
        print("Calculated slope on 3D fields")

        outputs_2D = self.reshape_after(outputs["ssthl0_test"].view[:,:,:])
        print("Reshaped back to 2D")

        return {"ssthl0_test": outputs_2D
            }
