from ndsl import Namelist, StencilFactory, Quantity
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.UW.compute_uwshcu import ComputeUwshcu
import numpy as np
from pyMoist.pyMoist_constants import ncnst
from ndsl.dsl.typing import FloatField, Float, Int, IntField


class TranslateComputeUwshcu(TranslateFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        self.ntracers_quantity_factory = ComputeUwshcu.make_ntracers_quantity_factory(
            self.quantity_factory
        )

        self.max_error = 1e-9

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "pifc0_in": {},
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

        # Float/Int Inputs
        self.in_vars["parameters"] = [
            "dotransport",
            "ncnst",
        ]

        # FloatField Outputs
        self.out_vars = {
            "tr0_test": self.grid.compute_dict(),
            "ssthl0_test": self.grid.compute_dict(),
            "ssqt0_test": self.grid.compute_dict(),
            "ssu0_test": self.grid.compute_dict(),
            "ssv0_test": self.grid.compute_dict(),
            "sstr0_test": self.grid.compute_dict(),
            "thj_test": self.grid.compute_dict(),
            "qvj_test": self.grid.compute_dict(),
            "qlj_test": self.grid.compute_dict(),
            "qij_test": self.grid.compute_dict(),
            "qse_test": self.grid.compute_dict(),
            "id_check_test": self.grid.compute_dict(),
        }

    def reshape_before(self,inputs):
        # Reshape input fields to the necessary shape
        i, j, k = self.grid.nic, self.grid.njc, self.grid.npz
        reshaped_inputs = {}
        
        for key, array in inputs.items():
            if isinstance(array, np.ndarray) and array.shape[1] == k and array.ndim == 2:
                reshaped_inputs[key] = np.reshape(array, (i, j, k)).astype(np.float32)
            elif isinstance(array, np.ndarray) and array.shape[1] == k+1 and array.ndim == 2:
                reshaped_inputs[key] = np.reshape(array, (i, j, k+1)).astype(np.float32)
            elif isinstance(array, np.ndarray) and array.ndim == 3:
                reshaped_inputs[key] = np.reshape(array, (i, j, k, ncnst)).astype(np.float32)
            else:
                # If not a 2D or 3D array, keep as is
                reshaped_inputs[key] = np.float32(array)

        return reshaped_inputs
    

    def reshape_after(self,outputs):
        # Reshape output fields back to original shape
        i, j, k = self.grid.nic, self.grid.njc, self.grid.npz
        
        if isinstance(outputs, dict):
            reshaped_outputs = {}
            for key, array in outputs.items():
                if array.ndim == 4: 
                    reshaped_array = np.reshape(array, (i * j, k, ncnst)).astype(np.float64)
                else:  
                    reshaped_array = np.reshape(array, (i * j, k)).astype(np.float64)
                reshaped_outputs[key] = reshaped_array
        else:
            # If outputs is not a dictionary, handle single array
            if outputs.ndim == 4: 
                reshaped_outputs = np.reshape(outputs, (i * j, k, ncnst)).astype(np.float64)
            else:
                reshaped_outputs = np.reshape(outputs, (i * j, k)).astype(np.float64)

        return reshaped_outputs
      
     
    def make_ijk_field(self, data, dtype=FloatField) -> Quantity:
        qty = self.quantity_factory.empty(
            [X_DIM, Y_DIM, Z_DIM],
            "n/a", dtype=dtype
        )
        qty.view[:, :, :] = qty.np.asarray(data[:, :, :])
        return qty
    
    def make_zinterface_field(self, data, dtype=FloatField) -> Quantity:
        qty = self.quantity_factory.empty(
            [X_DIM, Y_DIM, Z_INTERFACE_DIM],
            "n/a", dtype=dtype
        )
        qty.view[:, :, :] = qty.np.asarray(data[:, :, :])
        return qty

    def make_ntracers_ijk_field(self, data) -> Quantity:
        qty = self.ntracers_quantity_factory.empty(
            [X_DIM, Y_DIM, Z_DIM, "ntracers"],
            "n/a",
        )
        qty.view[:, :, :, :] = qty.np.asarray(data[:, :, :, :])
        return qty
        

    # Perform stencil computation
    def compute(self, inputs):
        compute_uwshcu = ComputeUwshcu(
            self.stencil_factory,
            self.grid.quantity_factory,
            int(inputs["ncnst"]),
        )

        # Reshape input variables and add halo
        inputs_reshaped = self.reshape_before(inputs)
        print("Reshaped inputs")
        
        # Inputs
        dotransport = inputs_reshaped["dotransport"]
        pifc0_in = self.make_zinterface_field(inputs_reshaped["pifc0_in"])
        pmid0_in = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        zmid0_in = self.make_ijk_field(inputs_reshaped["zmid0_in"])
        exnmid0_in = self.make_ijk_field(inputs_reshaped["exnmid0_in"])
        u0_in = self.make_ijk_field(inputs_reshaped["u0_in"])
        v0_in = self.make_ijk_field(inputs_reshaped["v0_in"])
        qv0_in = self.make_ijk_field(inputs_reshaped["qv0_in"])
        ql0_in = self.make_ijk_field(inputs_reshaped["ql0_in"])
        qi0_in = self.make_ijk_field(inputs_reshaped["qi0_in"])
        th0_in = self.make_ijk_field(inputs_reshaped["th0_in"])
        tr0_inout = self.make_ntracers_ijk_field(inputs_reshaped["tr0_inout"])


        # Outputs
        tr0_test = self.make_ntracers_ijk_field(inputs_reshaped["tr0_inout"])
        ssthl0_test = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        ssqt0_test = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        ssu0_test = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        ssv0_test = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        sstr0_test = self.make_ntracers_ijk_field(inputs_reshaped["tr0_inout"])
        thj_test = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        qvj_test = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        qlj_test = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        qij_test = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        qse_test = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        id_check_test = self.make_ijk_field(inputs_reshaped["pmid0_in"], dtype=Int)


        compute_uwshcu(
            dotransport=dotransport,
            pifc0_in=pifc0_in,
            pmid0_in=pmid0_in,
            zmid0_in=zmid0_in,
            exnmid0_in=exnmid0_in,
            u0_in=u0_in,
            v0_in=v0_in,
            qv0_in=qv0_in,
            ql0_in=ql0_in,
            qi0_in=qi0_in,
            th0_in=th0_in,
            tr0_inout=tr0_inout,
            tr0_test=tr0_test,
            ssthl0_test=ssthl0_test,
            ssqt0_test=ssqt0_test,
            ssu0_test=ssu0_test,
            ssv0_test=ssv0_test,
            sstr0_test=sstr0_test,
            thj_test=thj_test,
            qvj_test=qvj_test,
            qlj_test=qlj_test,
            qij_test=qij_test,
            qse_test=qse_test,
            id_check_test=id_check_test,
        )
        print("Performed compute_uwshcu on reshaped inputs")

        # Reshape output variables back to original shape
        tr0_test_3D = self.reshape_after(tr0_test.view[:,:,:,:])
        ssthl0_test_2D = self.reshape_after(ssthl0_test.view[:,:,:])
        ssqt0_test_2D = self.reshape_after(ssqt0_test.view[:,:,:])
        ssu0_test_2D = self.reshape_after(ssu0_test.view[:,:,:])
        ssv0_test_2D = self.reshape_after(ssv0_test.view[:,:,:])
        sstr0_test_3D = self.reshape_after(sstr0_test.view[:,:,:,:])
        thj_test_3D =  self.reshape_after(thj_test.view[:,:,:])
        qvj_test_3D = self.reshape_after(qvj_test.view[:,:,:])
        qlj_test_3D = self.reshape_after(qlj_test.view[:,:,:])
        qij_test_3D = self.reshape_after(qij_test.view[:,:,:])
        qse_test_3D = self.reshape_after(qse_test.view[:,:,:])
        id_check_test_3D = self.reshape_after(id_check_test.view[:,:,:])
        print("Reshaped outputs back to original shape")

        return {"tr0_test": tr0_test_3D,
                "ssthl0_test": ssthl0_test_2D,
                "ssqt0_test": ssqt0_test_2D,
                "ssu0_test": ssu0_test_2D,
                "ssv0_test": ssv0_test_2D,
                "sstr0_test": sstr0_test_3D,
                "thj_test": thj_test_3D,
                "qvj_test": qvj_test_3D,
                "qlj_test": qlj_test_3D,
                "qij_test": qij_test_3D,
                "qse_test": qse_test_3D,
                "id_check_test": id_check_test_3D,
            }