from ndsl import Namelist, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.conden import Conden
import numpy as np
from ndsl.dsl.typing import IntField 


class TranslateConden(TranslateFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.compute_func = Conden(
            self.stencil_factory,
            self.grid.quantity_factory,
        )

        self.max_error = 1e-9

        print(grid.__dict__)

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "pifc0_test": {},
            "thl0bot_test": {},
            "qt0bot_test": {},
        }

        print(self.in_vars["data_vars"])

        # FloatField Outputs
        self.out_vars = {
            "thj_test": self.grid.compute_dict(),
            "qvj_test": self.grid.compute_dict(),
            "qlj_test": self.grid.compute_dict(),
            "qij_test": self.grid.compute_dict(),
            "qse_test": self.grid.compute_dict(),
            "id_check_test": self.grid.compute_dict(),
        }

    def reshape_before(self,inputs):
        # Reshape input fields from (i*j, k) to (i, j, k)
        i, j, k = self.grid.nic, self.grid.njc, self.grid.npz 
        reshaped_inputs = {key: np.reshape(inputs[key], (i, j, k)).astype(np.float32) for key in inputs}
    
        # Add halo of 3 in i and j dimensions
        halo_arrays = {}
        for key, array in reshaped_inputs.items():
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
            "thj_test": self.grid.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="unknown"),
            "qvj_test": self.grid.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="unknown"),
            "qlj_test": self.grid.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="unknown"),
            "qij_test": self.grid.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="unknown"),
            "qse_test": self.grid.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="unknown"),
            "id_check_test": self.grid.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], dtype=IntField, units="unknown"),
        }
     
        self.compute_func(**inputs_2D,**outputs)
        print("Calculated conden on 3D fields")

        outputs_2D = self.reshape_after(outputs)
        print("Reshaped back to 2D")

        return {"thj_test": outputs_2D["thj_test"],
                "qvj_test": outputs_2D["qvj_test"],
                "qlj_test": outputs_2D["qlj_test"],
                "qij_test": outputs_2D["qij_test"],
                "qse_test": outputs_2D["qse_test"],
                "id_check_test": outputs_2D["id_check_test"],
            }
