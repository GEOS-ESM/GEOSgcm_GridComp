from ndsl import Namelist, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.compute_uwshcu import ComputeUwshcu
import numpy as np
from ndsl.dsl.typing import IntField 
from pyMoist.pyMoist_constants import ncnst


class TranslateComputeUwshcu(TranslateFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.compute_func = ComputeUwshcu(
            self.stencil_factory,
            self.grid.quantity_factory,
        )

        self.max_error = 1e-9

        print(grid.__dict__)

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "dotransport": {},
            "exnifc0_in": {},
            "pmid0_in": {},
            "zmid0_in": {},
            "exnmid0_in": {},
            "u0_in": {},
            "v0_in": {},
            "qv0_in": {},
            "ql0_in": {},
            "qi0_in": {},
            "th0_in": {},
            "tr0_inout": {},
        }

        print(self.in_vars["data_vars"])

        # FloatField Outputs
        self.out_vars = {
            "tr0_test": self.grid.compute_dict(),
            "ssthl0_test": self.grid.compute_dict(),
            "ssqt0_test": self.grid.compute_dict(),
            "ssu0_test": self.grid.compute_dict(),
            "ssv0_test": self.grid.compute_dict(),
            "sstr0_test": self.grid.compute_dict(),
        }

    def reshape_before(self,inputs):
        # Reshape input fields from (i*j, k) to (i, j, k)
        i, j, k = self.grid.nic, self.grid.njc, self.grid.npz
        n = ncnst # Number of tracers

        reshaped_inputs = {}

        for key, array in inputs.items():
            print(key,array)
            # Check if the input is a numpy array and has 2 dimensions
            if isinstance(array, np.ndarray) and array.ndim == 2:
                reshaped_inputs[key] = np.reshape(array, (i, j, k)).astype(np.float32)
            elif isinstance(array, np.ndarray) and array.ndim == 3:
                reshaped_inputs[key] = np.reshape(array, (i, j, k, n)).astype(np.float32)
            else:
                # If not a 2D or 3D array, keep as is
                reshaped_inputs[key] = array

        #reshaped_inputs = {key: np.reshape(inputs[key], (i, j, k)).astype(np.float32) for key in inputs}
    
        # Add halo of 3 in i and j dimensions
        halo_arrays = {}
        for key, array in reshaped_inputs.items():
            # Add conditional statements to deal wit 2D and 3D arrays and Floats
            new_array = np.zeros((i+6, j+6, k), dtype=np.float64)
            new_array[3:-3, 3:-3, :] = array
            halo_arrays[key] = new_array

        return halo_arrays

    def reshape_after(self,outputs):
        # Reshape output fields from (i, j, k) back to (i*j, k)
        i, j, k = self.grid.nic, self.grid.njc, self.grid.npz
        if len(outputs) > 1:
            reshaped_outputs = {}
            for key, array in outputs.items():
                reshaped_array = np.reshape(array.view[:,:,:], (i*j, k))
                reshaped_outputs[key] = reshaped_array
        else:
            reshaped_outputs = np.reshape(outputs, (i*j, k))

        return reshaped_outputs
        

    # Calculated Outputs
    def compute(self, inputs):

        inputs_2D = self.reshape_before(inputs)
        print("Reshaped 2D fields to 3D")

        outputs = {
            "tr0_test": self.grid.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="unknown"),
            "ssthl0_test": self.grid.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="unknown"),
            "ssqt0_test": self.grid.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="unknown"),
            "ssu0_test": self.grid.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="unknown"),
            "ssv0_test": self.grid.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="unknown"),
            "sstr0_test": self.grid.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], dtype=IntField, units="unknown"),
        }

     
        self.compute_func(**inputs_2D,**outputs)
        print("Calculated conden on 3D fields")

        outputs_2D = self.reshape_after(outputs)
        print("Reshaped back to 2D")

        return {"tr0_test": outputs_2D["tr0_test"],
                "ssthl0_test": outputs_2D["ssthl0_test"],
                "ssqt0_test": outputs_2D["ssqt0_test"],
                "ssu0_test": outputs_2D["ssu0_test"],
                "ssv0_test": outputs_2D["ssv0_test"],
                "sstr0_test": outputs_2D["sstr0_test"],
            }