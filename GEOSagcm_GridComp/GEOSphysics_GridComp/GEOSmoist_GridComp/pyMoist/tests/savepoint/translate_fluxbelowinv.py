from ndsl import Namelist, StencilFactory, Quantity
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.UW.fluxbelowinv import FluxBelowInv
import numpy as np
from ndsl.dsl.typing import (
    FloatFieldIJ,
    FloatField,
    IntFieldIJ,
    IntField
)
import xarray as xr


class TranslateFluxBelowInv(TranslateFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.quantity_factory = grid.quantity_factory
        self.compute_func = FluxBelowInv(
            self.stencil_factory,
            self.grid.quantity_factory,
        )

        self.max_error = 1e-9

        print(grid.__dict__)

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "cbmf": {},
            "ps0": {},
            "kinv_in": {},
            "dt": {},
            "xsrc": {},
            "xmean": {},
            "xtop": {},
            "xbot": {},
            "xflx": {},
        }
      
        # Float/Int Inputs
        self.in_vars["parameters"] = [
        ]

        # FloatField Outputs
        self.out_vars = {
            "xflx_out": self.grid.compute_dict(),
        }

    def reshape_before(self, inputs):
        # Reshape input fields to the necessary shape
        i, j, k = self.grid.nic, self.grid.njc, self.grid.npz
        reshaped_inputs = {}
        for key, array in inputs.items():
            if key == "ps0" or key == "xflx":
                reshaped_inputs[key] = np.reshape(array[:,0,:,0], (i, j,k+1)).astype(np.float32)
            else:
                reshaped_inputs[key] = np.reshape(array[:,0,72,0], (i, j)).astype(np.float32)

        return reshaped_inputs


    def reshape_after(self, outputs):
        # Reshape output fields back to original shape
        i, j, k = self.grid.nic, self.grid.njc, self.grid.npz
        reshaped_outputs = np.reshape(outputs, (i * j, k+1)).astype(np.float32)

        return reshaped_outputs

    def make_ijk_field(self, data, dtype=FloatField) -> Quantity:
        qty = self.quantity_factory.empty([X_DIM, Y_DIM, Z_DIM], "n/a", dtype=dtype)
        qty.view[:, :, :] = qty.np.asarray(data[:, :, :])
        return qty

    def make_ij_field(self, data, dtype=FloatField) -> Quantity:
        qty = self.quantity_factory.empty([X_DIM, Y_DIM], "n/a", dtype=dtype)
        qty.view[:, :] = qty.np.asarray(data[:, :])
        return qty

    def make_zinterface_field(self, data, dtype=FloatField) -> Quantity:
        qty = self.quantity_factory.empty(
            [X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a", dtype=dtype
        )
        qty.view[:, :, :] = qty.np.asarray(data[:, :, :])
        return qty


    # Calculated Outputs
    def compute(self, inputs):
        fluxbelowinv = FluxBelowInv(
            self.stencil_factory,
            self.grid.quantity_factory,
        )
      
        inputs_reshaped = self.reshape_before(inputs)
        print("Reshaped inputs...")
        
        # Inputs
        cbmf = self.make_ij_field(inputs_reshaped["cbmf"], dtype=FloatFieldIJ)
        ps0 = self.make_zinterface_field(inputs_reshaped["ps0"], dtype=FloatField)
        kinv = self.make_ij_field(inputs_reshaped["kinv_in"], dtype=IntFieldIJ)
        dt = self.make_ij_field(inputs_reshaped["dt"], dtype=FloatFieldIJ)
        xsrc = self.make_ij_field(inputs_reshaped["xsrc"], dtype=FloatFieldIJ)
        xmean = self.make_ij_field(inputs_reshaped["xmean"], dtype=FloatFieldIJ)
        xtopin = self.make_ij_field(inputs_reshaped["xtop"], dtype=FloatFieldIJ)
        xbotin = self.make_ij_field(inputs_reshaped["xbot"], dtype=FloatFieldIJ)

        # Outputs
        xflx = self.make_zinterface_field(inputs_reshaped["xflx"], dtype=FloatField)

        fluxbelowinv(
            cbmf=cbmf,
            ps0=ps0,
            kinv=kinv,
            dt=dt,
            xsrc=xsrc,
            xmean=xmean,
            xtopin=xtopin,
            xbotin=xbotin,
            xflx=xflx,
        )

        print("Executed stencil computation...")
        # Reshape output variables back to original shape

        #xflx_out = self.reshape_after(xflx.view[:,:,:])
        #print("Reshaped outputs back to original shape...")
        #print(xflx_out)

        #xflx4D=inputs["xflx"]
        #xflx4D[:,0,:,0] = xflx_out

        # Replace nans with zeros
        #xflx_out_zeros = self.subset_output(self, self.out_vars["xflx_out"])

        with xr.open_dataset("/Users/kfandric/netcdf/FluxBelowInv-Out.nc") as ds:
            # Load in netcdf test var
            xflx_nan = ds.variables["xflx_out"].data[:,0,:,0]
            # Replace nans with zero
            xflx_zeros= np.nan_to_num(xflx_nan,nan=0)
            # Reshape and make testvar quantity
            i, j, k = self.grid.nic, self.grid.njc, self.grid.npz
            xflx_out_reshaped =  np.reshape(xflx_zeros, (i, j,k+1)).astype(np.float32)
            xflx_out = self.make_zinterface_field(
                xflx_out_reshaped
            )

        # Run translate test by hand
        print("TESTING: xflx_out")
        tot_indicies = 0
        failed_indicies = 0
        testing_variable = xflx
        reference_variable = xflx_out
        for i in range(testing_variable.view[:].shape[0]):
            for j in range(testing_variable.view[:].shape[1]):
                for k in range(testing_variable.view[:].shape[2]):
                    tot_indicies = tot_indicies + 1
                    if (
                        testing_variable.view[i, j, k]
                        != reference_variable.view[i, j, k]
                    ):
                        failed_indicies = failed_indicies + 1
                        print(
                            "DIFF: ",
                            i,
                            j,
                            k,
                            "computed: ",
                            testing_variable.view[i, j, k],
                            "reference: ",
                            reference_variable.view[i, j, k],
                            "difference: ",
                            testing_variable.view[i, j, k]
                            - reference_variable.view[i, j, k],
                        )
        print("failures: ", failed_indicies, "/", tot_indicies, "(",(failed_indicies/tot_indicies)*100,"%)")
        
        #return {"xflx_out": xflx4D
        #    }
            