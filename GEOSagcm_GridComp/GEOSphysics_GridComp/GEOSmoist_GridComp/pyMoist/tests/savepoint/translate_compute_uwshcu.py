from ndsl import Namelist, StencilFactory, Quantity
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.UW.compute_uwshcu import ComputeUwshcu
import numpy as np
from pyMoist.pyMoist_constants import ncnst
from ndsl.dsl.typing import (
    FloatField,
    Float,
    Int,
    IntField,
    IntFieldIJ,
    FloatFieldIJ,
    BoolField,
)


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
            "kpbl_in": {},
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
            "k0",
        ]

        # FloatField Outputs
        self.out_vars = {
            "umf_out": self.grid.compute_dict(),
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
            "dcm_out": self.grid.compute_dict(),
            "qvten_out": self.grid.compute_dict(),
            "qlten_out": self.grid.compute_dict(),
            "qiten_out": self.grid.compute_dict(),
            "sten_out": self.grid.compute_dict(),
            "uten_out": self.grid.compute_dict(),
            "vten_out": self.grid.compute_dict(),
            "qrten_out": self.grid.compute_dict(),
            "qsten_out": self.grid.compute_dict(),
            "cufrc_out": self.grid.compute_dict(),
            "fer_out": self.grid.compute_dict(),
            "fdr_out": self.grid.compute_dict(),
            "qldet_out": self.grid.compute_dict(),
            "qidet_out": self.grid.compute_dict(),
            "qlsub_out": self.grid.compute_dict(),
            "qisub_out": self.grid.compute_dict(),
            "ndrop_out": self.grid.compute_dict(),
            "nice_out": self.grid.compute_dict(),
            "tpert_out": self.grid.compute_dict(),
            "qpert_out": self.grid.compute_dict(),
            "qtflx_out": self.grid.compute_dict(),
            "slflx_out": self.grid.compute_dict(),
            "uflx_out": self.grid.compute_dict(),
            "vflx_out": self.grid.compute_dict(),
        }

    def reshape_before(self, inputs):
        # Reshape input fields to the necessary shape
        i, j, k = self.grid.nic, self.grid.njc, self.grid.npz
        reshaped_inputs = {}

        for key, array in inputs.items():
            if (
                isinstance(array, np.ndarray)
                and array.ndim == 2
                and array.shape[1] == k
            ):
                reshaped_inputs[key] = np.reshape(array, (i, j, k)).astype(np.float32)
            elif (
                isinstance(array, np.ndarray)
                and array.ndim == 2
                and array.shape[1] == k + 1
            ):
                reshaped_inputs[key] = np.reshape(array, (i, j, k + 1)).astype(
                    np.float32
                )
            elif (
                isinstance(array, np.ndarray)
                and array.ndim == 1
                and array.shape[0] == 576
            ):
                reshaped_inputs[key] = np.reshape(array, (i, j)).astype(np.int64)
            elif isinstance(array, np.ndarray) and array.ndim == 3:
                reshaped_inputs[key] = np.reshape(array, (i, j, k, ncnst)).astype(
                    np.float32
                )
            else:
                # If not a 2D or 3D array, keep as is
                reshaped_inputs[key] = np.float32(array)

        return reshaped_inputs

    def reshape_after(self, outputs):
        # Reshape output fields back to original shape
        i, j, k = self.grid.nic, self.grid.njc, self.grid.npz

        if isinstance(outputs, dict):
            reshaped_outputs = {}
            for key, array in outputs.items():
                if array.ndim == 4:
                    reshaped_array = np.reshape(array, (i * j, k, ncnst)).astype(
                        np.float64
                    )
                else:
                    reshaped_array = np.reshape(array, (i * j, k)).astype(np.float64)
                reshaped_outputs[key] = reshaped_array
        else:
            # If outputs is not a dictionary, handle single array
            if outputs.ndim == 4:
                reshaped_outputs = np.reshape(outputs, (i * j, k, ncnst)).astype(
                    np.float64
                )
            elif outputs.ndim == 3 and outputs.shape[2] == 73:
                reshaped_outputs = np.reshape(outputs, (i * j, k + 1)).astype(
                    np.float64
                )
            elif outputs.ndim == 2 and outputs.shape[0] == 24:
                reshaped_outputs = np.reshape(outputs, (i * j)).astype(np.float64)
            else:
                reshaped_outputs = np.reshape(outputs, (i * j, k)).astype(np.float64)

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
        )

        # Reshape input variables and add halo
        inputs_reshaped = self.reshape_before(inputs)
        print("Reshaped inputs")

        # Inputs
        dotransport = inputs_reshaped["dotransport"]
        ncnst = np.int64(inputs_reshaped["ncnst"])
        k0 = np.int64(inputs_reshaped["k0"])
        kpbl_in = self.make_ij_field(inputs_reshaped["kpbl_in"], dtype=IntFieldIJ)
        pifc0_in = self.make_zinterface_field(inputs_reshaped["pifc0_in"])
        pmid0_in = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        exnmid0_in = self.make_ijk_field(inputs_reshaped["exnmid0_in"])
        u0_in = self.make_ijk_field(inputs_reshaped["u0_in"])
        v0_in = self.make_ijk_field(inputs_reshaped["v0_in"])
        qv0_in = self.make_ijk_field(inputs_reshaped["qv0_in"])
        ql0_in = self.make_ijk_field(inputs_reshaped["ql0_in"])
        qi0_in = self.make_ijk_field(inputs_reshaped["qi0_in"])
        th0_in = self.make_ijk_field(inputs_reshaped["th0_in"])
        tr0_inout = self.make_ntracers_ijk_field(inputs_reshaped["tr0_inout"])

        # Outputs
        umf_out = self.make_zinterface_field(inputs_reshaped["pifc0_in"])
        dcm_out = self.make_ijk_field(inputs_reshaped["v0_in"])
        qvten_out = self.make_ijk_field(inputs_reshaped["v0_in"])
        qlten_out = self.make_ijk_field(inputs_reshaped["v0_in"])
        qiten_out = self.make_ijk_field(inputs_reshaped["v0_in"])
        sten_out = self.make_ijk_field(inputs_reshaped["v0_in"])
        uten_out = self.make_ijk_field(inputs_reshaped["v0_in"])
        vten_out = self.make_ijk_field(inputs_reshaped["v0_in"])
        qrten_out = self.make_ijk_field(inputs_reshaped["v0_in"])
        qsten_out = self.make_ijk_field(inputs_reshaped["v0_in"])
        cufrc_out = self.make_ijk_field(inputs_reshaped["v0_in"])
        fer_out = self.make_ijk_field(inputs_reshaped["v0_in"])
        fdr_out = self.make_ijk_field(inputs_reshaped["v0_in"])
        qldet_out = self.make_ijk_field(inputs_reshaped["v0_in"])
        qidet_out = self.make_ijk_field(inputs_reshaped["v0_in"])
        qlsub_out = self.make_ijk_field(inputs_reshaped["v0_in"])
        qisub_out = self.make_ijk_field(inputs_reshaped["v0_in"])
        ndrop_out = self.make_ijk_field(inputs_reshaped["v0_in"])
        nice_out = self.make_ijk_field(inputs_reshaped["v0_in"])
        tpert_out = self.make_ij_field(inputs_reshaped["kpbl_in"], dtype=FloatFieldIJ)
        qpert_out = self.make_ij_field(inputs_reshaped["kpbl_in"], dtype=FloatFieldIJ)
        qtflx_out = self.make_zinterface_field(inputs_reshaped["pifc0_in"])
        slflx_out = self.make_zinterface_field(inputs_reshaped["pifc0_in"])
        uflx_out = self.make_zinterface_field(inputs_reshaped["pifc0_in"])
        vflx_out = self.make_zinterface_field(inputs_reshaped["pifc0_in"])
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
            ncnst=ncnst,
            k0=k0,
            kpbl_in=kpbl_in,
            pifc0_in=pifc0_in,
            pmid0_in=pmid0_in,
            exnmid0_in=exnmid0_in,
            u0_in=u0_in,
            v0_in=v0_in,
            qv0_in=qv0_in,
            ql0_in=ql0_in,
            qi0_in=qi0_in,
            th0_in=th0_in,
            tr0_inout=tr0_inout,
            umf_out=umf_out,
            dcm_out=dcm_out,
            qvten_out=qvten_out,
            qlten_out=qlten_out,
            qiten_out=qiten_out,
            sten_out=sten_out,
            uten_out=uten_out,
            vten_out=vten_out,
            qrten_out=qrten_out,
            qsten_out=qsten_out,
            cufrc_out=cufrc_out,
            fer_out=fer_out,
            fdr_out=fdr_out,
            qldet_out=qldet_out,
            qidet_out=qidet_out,
            qlsub_out=qlsub_out,
            qisub_out=qisub_out,
            ndrop_out=ndrop_out,
            nice_out=nice_out,
            tpert_out=tpert_out,
            qpert_out=qpert_out,
            qtflx_out=qtflx_out,
            slflx_out=slflx_out,
            uflx_out=uflx_out,
            vflx_out=vflx_out,
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
        umf_out_2D = self.reshape_after(umf_out.view[:, :, :])
        dcm_out_2D = self.reshape_after(dcm_out.view[:, :, :])
        qvten_out_2D = self.reshape_after(qvten_out.view[:, :, :])
        qlten_out_2D = self.reshape_after(qlten_out.view[:, :, :])
        qiten_out_2D = self.reshape_after(qiten_out.view[:, :, :])
        sten_out_2D = self.reshape_after(sten_out.view[:, :, :])
        uten_out_2D = self.reshape_after(uten_out.view[:, :, :])
        vten_out_2D = self.reshape_after(vten_out.view[:, :, :])
        qrten_out_2D = self.reshape_after(qrten_out.view[:, :, :])
        qsten_out_2D = self.reshape_after(qsten_out.view[:, :, :])
        cufrc_out_2D = self.reshape_after(cufrc_out.view[:, :, :])
        fer_out_2D = self.reshape_after(fer_out.view[:, :, :])
        fdr_out_2D = self.reshape_after(fdr_out.view[:, :, :])
        qldet_out_2D = self.reshape_after(qldet_out.view[:, :, :])
        qidet_out_2D = self.reshape_after(qidet_out.view[:, :, :])
        qlsub_out_2D = self.reshape_after(qlsub_out.view[:, :, :])
        qisub_out_2D = self.reshape_after(qisub_out.view[:, :, :])
        ndrop_out_2D = self.reshape_after(ndrop_out.view[:, :, :])
        nice_out_2D = self.reshape_after(nice_out.view[:, :, :])
        tpert_out_2D = self.reshape_after(tpert_out.view[:, :])
        qpert_out_2D = self.reshape_after(qpert_out.view[:, :])
        qtflx_out_2D = self.reshape_after(qtflx_out.view[:, :, :])
        slflx_out_2D = self.reshape_after(slflx_out.view[:, :, :])
        uflx_out_2D = self.reshape_after(uflx_out.view[:, :, :])
        vflx_out_2D = self.reshape_after(vflx_out.view[:, :, :])
        tr0_test_3D = self.reshape_after(tr0_test.view[:, :, :, :])
        ssthl0_test_2D = self.reshape_after(ssthl0_test.view[:, :, :])
        ssqt0_test_2D = self.reshape_after(ssqt0_test.view[:, :, :])
        ssu0_test_2D = self.reshape_after(ssu0_test.view[:, :, :])
        ssv0_test_2D = self.reshape_after(ssv0_test.view[:, :, :])
        sstr0_test_3D = self.reshape_after(sstr0_test.view[:, :, :, :])
        thj_test_3D = self.reshape_after(thj_test.view[:, :, :])
        qvj_test_3D = self.reshape_after(qvj_test.view[:, :, :])
        qlj_test_3D = self.reshape_after(qlj_test.view[:, :, :])
        qij_test_3D = self.reshape_after(qij_test.view[:, :, :])
        qse_test_3D = self.reshape_after(qse_test.view[:, :, :])
        id_check_test_3D = self.reshape_after(id_check_test.view[:, :, :])
        print("Reshaped outputs back to original shape")

        return {
            "umf_out": umf_out_2D,
            "dcm_out": dcm_out_2D,
            "qvten_out": qvten_out_2D,
            "qlten_out": qlten_out_2D,
            "qiten_out": qiten_out_2D,
            "sten_out": sten_out_2D,
            "uten_out": uten_out_2D,
            "vten_out": vten_out_2D,
            "qrten_out": qrten_out_2D,
            "qsten_out": qsten_out_2D,
            "cufrc_out": cufrc_out_2D,
            "fer_out": fer_out_2D,
            "fdr_out": fdr_out_2D,
            "qldet_out": qldet_out_2D,
            "qidet_out": qidet_out_2D,
            "qlsub_out": qlsub_out_2D,
            "qisub_out": qisub_out_2D,
            "ndrop_out": ndrop_out_2D,
            "nice_out": nice_out_2D,
            "tpert_out": tpert_out_2D,
            "qpert_out": qpert_out_2D,
            "qtflx_out": qtflx_out_2D,
            "slflx_out": slflx_out_2D,
            "uflx_out": uflx_out_2D,
            "vflx_out": vflx_out_2D,
            "tr0_test": tr0_test_3D,
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
