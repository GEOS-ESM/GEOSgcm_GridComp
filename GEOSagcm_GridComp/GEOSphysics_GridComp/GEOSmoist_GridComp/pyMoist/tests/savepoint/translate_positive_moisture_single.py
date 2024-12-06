from ndsl import Namelist, StencilFactory, Quantity
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.UW.positive_moisture_single import PositiveMoistureSingle
import numpy as np
from ndsl.dsl.typing import (
    FloatFieldIJ,
    FloatField,
)
import xarray as xr


class TranslatePositiveMoistureSingle(TranslateFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.quantity_factory = grid.quantity_factory
        self.compute_func = PositiveMoistureSingle(
            self.stencil_factory,
            self.grid.quantity_factory,
        )

        self.max_error = 1e-9

        print(grid.__dict__)

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "qlmin": {},
            "xlv": {},
            "ql0_star": {},
            "dt": {},
            "qvmin": {},
            "qiten": {},
            "qv0_star": {},
            "dp0_in": {},
            "qlten": {},
            "qi0_star": {},
            "qvten": {},
            "qimin": {},
            "xls": {},
            "sten": {},
            "s0_star": {},
        }
      
        # Float/Int Inputs
        self.in_vars["parameters"] = [
            "mkx",
        ]

        # FloatField Outputs
        self.out_vars = {
            "ql0_star": self.grid.compute_dict(),
            "qi0_star": self.grid.compute_dict(),
            "qiten": self.grid.compute_dict(),
            "qv0_star": self.grid.compute_dict(),
            "qlten": self.grid.compute_dict(),
            "qvten": self.grid.compute_dict(),
            "sten": self.grid.compute_dict(),
            "s0_star": self.grid.compute_dict(),
        }

    def reshape_before(self, inputs):
        # Reshape input fields to the necessary shape
        i, j, k = self.grid.nic, self.grid.njc, self.grid.npz
        reshaped_inputs = {}
        for key, array in inputs.items():
            if key == "dp0_in" or key == "qv0_star" or key == "ql0_star" or key == "qi0_star" or key == "s0_star" or key == "qvten" or key == "qlten" or key == "qiten" or key == "sten":
                reshaped_inputs[key] = np.reshape(array[:,0,0:72,0], (i, j,k)).astype(np.float32)
            elif key == "xlv" or key =="xls" or key=="dt" or key=="qvmin" or key=="qlmin" or key=="qimin":
                reshaped_inputs[key] = np.reshape(array[:,0,72,0], (i, j)).astype(np.float32)
            else:
                reshaped_inputs[key] = np.int64(array)

        return reshaped_inputs


    def reshape_after(self, outputs):
        # Reshape output fields back to original shape
        i, j, k = self.grid.nic, self.grid.njc, self.grid.npz
        reshaped_outputs = np.reshape(outputs, (i * j, k)).astype(np.float32)

        return reshaped_outputs

    def make_ijk_field(self, data, dtype=FloatField) -> Quantity:
        qty = self.quantity_factory.empty([X_DIM, Y_DIM, Z_DIM], "n/a", dtype=dtype)
        qty.view[:, :, :] = qty.np.asarray(data[:, :, :])
        return qty

    def make_ij_field(self, data, dtype=FloatField) -> Quantity:
        qty = self.quantity_factory.empty([X_DIM, Y_DIM], "n/a", dtype=dtype)
        qty.view[:, :] = qty.np.asarray(data[:, :])
        return qty


    # Calculated Outputs
    def compute(self, inputs):
        positive_moisture_single = PositiveMoistureSingle(
            self.stencil_factory,
            self.grid.quantity_factory,
        )
      
        inputs_reshaped = self.reshape_before(inputs)
        print("Reshaped inputs...")
        
        # Inputs
        xlv = self.make_ij_field(inputs_reshaped["xlv"], dtype=FloatFieldIJ)
        xls = self.make_ij_field(inputs_reshaped["xls"], dtype=FloatFieldIJ)
        mkx = inputs_reshaped["mkx"][0]
        dt = self.make_ij_field(inputs_reshaped["dt"], dtype=FloatFieldIJ)
        qvmin = self.make_ij_field(inputs_reshaped["qvmin"], dtype=FloatFieldIJ)
        qlmin = self.make_ij_field(inputs_reshaped["qlmin"], dtype=FloatFieldIJ)
        qimin = self.make_ij_field(inputs_reshaped["qimin"], dtype=FloatFieldIJ)
        dp0 = self.make_ijk_field(inputs_reshaped["dp0_in"], dtype=FloatField)

        # In/Outs
        qv0_star = self.make_ijk_field(inputs_reshaped["qv0_star"], dtype=FloatField)
        ql0_star = self.make_ijk_field(inputs_reshaped["ql0_star"], dtype=FloatField)
        qi0_star = self.make_ijk_field(inputs_reshaped["qi0_star"], dtype=FloatField)
        s0_star = self.make_ijk_field(inputs_reshaped["s0_star"], dtype=FloatField)
        qvten = self.make_ijk_field(inputs_reshaped["qvten"], dtype=FloatField)
        qlten = self.make_ijk_field(inputs_reshaped["qlten"], dtype=FloatField)
        qiten = self.make_ijk_field(inputs_reshaped["qiten"], dtype=FloatField)
        sten = self.make_ijk_field(inputs_reshaped["sten"], dtype=FloatField)

     
        positive_moisture_single(
            xlv=xlv,
            xls=xls,
            mkx=mkx,
            dt=dt,
            qvmin=qvmin,
            qlmin=qlmin,
            qimin=qimin,
            dp=dp0,
            qv=qv0_star,
            ql=ql0_star,
            qi=qi0_star,
            s=s0_star,
            qvten=qvten,
            qlten=qlten,
            qiten=qiten,
            sten=sten,
        )

        print("Executed stencil computation...")

        # Reshape output variables back to original shape
        qv_out = self.reshape_after(qv0_star.view[:,:,:])
        ql_out = self.reshape_after(ql0_star.view[:,:,:])
        qi_out = self.reshape_after(qi0_star.view[:,:,:])
        s_out = self.reshape_after(s0_star.view[:,:,:])
        qvten_out = self.reshape_after(qvten.view[:,:,:])
        qlten_out = self.reshape_after(qlten.view[:,:,:])
        qiten_out = self.reshape_after(qiten.view[:,:,:])
        sten_out = self.reshape_after(sten.view[:,:,:])
        print("Reshaped outputs back to original shape...")

        # Adjust the shapes of outputs before testing
        qvout4D=inputs["qv0_star"]
        qvout4D[:,0,0:72,0] = qv_out
        qlout4D=inputs["ql0_star"]
        qlout4D[:,0,0:72,0] = ql_out
        qiout4D=inputs["qi0_star"]
        qiout4D[:,0,0:72,0] = qi_out
        sout4D=inputs["s0_star"]
        sout4D[:,0,0:72,0] = s_out
        qvtenout4D=inputs["qvten"]
        qvtenout4D[:,0,0:72,0] = qvten_out
        qltenout4D=inputs["qlten"]
        qltenout4D[:,0,0:72,0] = qlten_out
        qitenout4D=inputs["qiten"]
        qitenout4D[:,0,0:72,0] = qiten_out
        stenout4D=inputs["sten"]
        stenout4D[:,0,0:72,0] = sten_out

        return {
            "qv0_star": qvout4D,
            "ql0_star": qlout4D,
            "qi0_star": qiout4D,
            "s0_star": sout4D,
            "qvten": qvtenout4D,
            "qiten": qitenout4D,
            "qlten": qltenout4D,
            "sten": stenout4D,
            }