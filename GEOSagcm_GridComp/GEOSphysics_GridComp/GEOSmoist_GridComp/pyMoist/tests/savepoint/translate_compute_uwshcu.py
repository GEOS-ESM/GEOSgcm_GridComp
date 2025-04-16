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
import xarray as xr


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
            "zifc0_in": {},
            "pmid0_in": {},
            "zmid0_in": {},
            "kpbl_in": {},
            "exnmid0_in": {},
            "exnifc0_in": {},
            "dp0_in": {},
            "u0_in": {},
            "v0_in": {},
            "qv0_in": {},
            "ql0_in": {},
            "qi0_in": {},
            "th0_in": {},
            "tr0_inout": {},
            "frland_in": {},
            "tke_in": {},
            "rkfre": {},
            "cush_inout": {},
            "shfx": {},
            "evap": {},
        }

        # Float/Int Inputs
        self.in_vars["parameters"] = [
            "dotransport",
            "ncnst",
            "k0",
            "windsrcavg",
            "qtsrchgt",
            "qtsrc_fac",
            "thlsrc_fac",
            "frc_rasn",
            "rbuoy",
            "epsvarw",
            "use_CINcin",
            "mumin1",
            "rmaxfrac",
            "PGFc",
            "dt",
            "niter_xc",
            "criqc",
            "rle",
            "cridist_opt",
            "mixscale",
            "rkm",
            "detrhgt",
            "rdrag",
            "use_self_detrain",
            "detrhgt",
            "use_cumpenent",
            "rpen",
            "use_momenflx",
            "rdrop",
        ]

        # FloatField Outputs
        self.out_vars = {
            "ssthl0": self.grid.compute_dict(),
            "ssqt0": self.grid.compute_dict(),
            "ssu0": self.grid.compute_dict(),
            "ssv0": self.grid.compute_dict(),
            "sstr0": self.grid.compute_dict(),
            "tr0": self.grid.compute_dict(),
            "thj": self.grid.compute_dict(),
            "qlj": self.grid.compute_dict(),
            "qij": self.grid.compute_dict(),
            "qvj": self.grid.compute_dict(),
            "qse": self.grid.compute_dict(),
            "id_check": self.grid.compute_dict(),
            "thvl0top": self.grid.compute_dict(),
            "tr0_o": self.grid.compute_dict(),
            "sstr0_o": self.grid.compute_dict(),
            "trflx": self.grid.compute_dict(),
            "trten": self.grid.compute_dict(),
            "tru": self.grid.compute_dict(),
            "tru_emf": self.grid.compute_dict(),
            "kinv": self.grid.compute_dict(),
            "umf_out": self.grid.compute_dict(),
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
            "cushinout": self.grid.compute_dict(),
            "qtflx_out": self.grid.compute_dict(),
            "slflx_out": self.grid.compute_dict(),
            "uflx_out": self.grid.compute_dict(),
            "vflx_out": self.grid.compute_dict(),
            "fer_out": self.grid.compute_dict(),
            "fdr_out": self.grid.compute_dict(),
            "thvlavg": self.grid.compute_dict(),
            "tkeavg": self.grid.compute_dict(),
            "uavg": self.grid.compute_dict(),
            "vavg": self.grid.compute_dict(),
            "thvlmin": self.grid.compute_dict(),
            "qtavg": self.grid.compute_dict(),
            "dpi": self.grid.compute_dict(),
            "thlsrc": self.grid.compute_dict(),
            "usrc": self.grid.compute_dict(),
            "vsrc": self.grid.compute_dict(),
            "plcl": self.grid.compute_dict(),
            "klcl_test": self.grid.compute_dict(),
            "thl0lcl": self.grid.compute_dict(),
            "qt0lcl": self.grid.compute_dict(),
            "thv0lcl": self.grid.compute_dict(),
            "plfc": self.grid.compute_dict(),
            "cin": self.grid.compute_dict(),
            "thvubot_test": self.grid.compute_dict(),
            "thvlsrc": self.grid.compute_dict(),
            "thv0bot_test": self.grid.compute_dict(),
            "thl0top": self.grid.compute_dict(),
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
                reshaped_inputs[key] = np.reshape(array, (i, j)).astype(np.float32)
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

    def make_ntracers_ij_field(self, data) -> Quantity:
        qty = self.ntracers_quantity_factory.empty(
            [X_DIM, Y_DIM, "ntracers"],
            "n/a",
        )
        qty.view[:, :, :] = qty.np.asarray(data[:, :, :])
        return qty

    def make_ntracers_zdim_field(self, data) -> Quantity:
        qty = self.ntracers_quantity_factory.empty(
            [X_DIM, Y_DIM, Z_INTERFACE_DIM, "ntracers"],
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

        # Float/Int Inputs
        dotransport = inputs_reshaped["dotransport"]
        ncnst = Int(inputs_reshaped["ncnst"])
        k0 = Int(inputs_reshaped["k0"])
        windsrcavg = Int(inputs_reshaped["windsrcavg"])
        qtsrchgt = inputs_reshaped["qtsrchgt"]
        qtsrc_fac = inputs_reshaped["qtsrc_fac"]
        thlsrc_fac = inputs_reshaped["thlsrc_fac"]
        frc_rasn = inputs_reshaped["frc_rasn"]
        rbuoy = inputs_reshaped["rbuoy"]
        epsvarw = inputs_reshaped["epsvarw"]
        use_CINcin = np.int32(inputs_reshaped["use_CINcin"])
        mumin1 = inputs_reshaped["mumin1"]
        rmaxfrac = inputs_reshaped["rmaxfrac"]
        PGFc = inputs_reshaped["PGFc"]
        dt = inputs_reshaped["dt"]
        niter_xc = Int(inputs_reshaped["niter_xc"])
        criqc = inputs_reshaped["criqc"]
        rle = inputs_reshaped["rle"]
        cridist_opt = np.int32(inputs_reshaped["cridist_opt"])
        mixscale = inputs_reshaped["mixscale"]
        rdrag = inputs_reshaped["rdrag"]
        rkm = inputs_reshaped["rkm"]
        use_self_detrain = np.int32(inputs_reshaped["use_self_detrain"])
        detrhgt = inputs_reshaped["detrhgt"]
        use_cumpenent = np.int32(inputs_reshaped["use_cumpenent"])
        rpen = inputs_reshaped["rpen"]
        use_momenflx = np.int32(inputs_reshaped["use_momenflx"])
        rdrop = inputs_reshaped["rdrop"]

        # Field inputs
        pifc0_in = self.make_zinterface_field(inputs_reshaped["pifc0_in"])
        zifc0_in = self.make_zinterface_field(inputs_reshaped["zifc0_in"])
        pmid0_in = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        zmid0_in = self.make_ijk_field(inputs_reshaped["zmid0_in"])
        kpbl_in = self.make_ij_field(inputs_reshaped["kpbl_in"], dtype=IntFieldIJ)
        exnmid0_in = self.make_ijk_field(inputs_reshaped["exnmid0_in"])
        exnifc0_in = self.make_zinterface_field(inputs_reshaped["exnifc0_in"])
        dp0_in = self.make_ijk_field(inputs_reshaped["dp0_in"])
        u0_in = self.make_ijk_field(inputs_reshaped["u0_in"])
        v0_in = self.make_ijk_field(inputs_reshaped["v0_in"])
        qv0_in = self.make_ijk_field(inputs_reshaped["qv0_in"])
        ql0_in = self.make_ijk_field(inputs_reshaped["ql0_in"])
        qi0_in = self.make_ijk_field(inputs_reshaped["qi0_in"])
        th0_in = self.make_ijk_field(inputs_reshaped["th0_in"])
        tr0_inout = self.make_ntracers_ijk_field(inputs_reshaped["tr0_inout"])
        frland_in = self.make_ij_field(inputs_reshaped["frland_in"])
        tke_in = self.make_zinterface_field(inputs_reshaped["tke_in"])
        rkfre = self.make_ij_field(inputs_reshaped["rkfre"])
        cush_inout = self.make_ij_field(inputs_reshaped["cush_inout"])
        shfx = self.make_ij_field(inputs_reshaped["shfx"])
        evap = self.make_ij_field(inputs_reshaped["evap"])

        # Outputs
        ssthl0 = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        ssqt0 = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        ssu0 = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        ssv0 = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        sstr0 = self.make_ntracers_ijk_field(inputs_reshaped["tr0_inout"])
        tr0 = self.make_ntracers_ijk_field(inputs_reshaped["tr0_inout"])
        thj = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        qlj = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        qij = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        qvj = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        qse = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        id_check = self.make_ijk_field(inputs_reshaped["pmid0_in"], dtype=Int)
        thv0top = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        thv0bot = self.make_ijk_field(np.zeros(shape=(24, 24, 72)))
        thvl0top = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        tr0_o = self.make_ntracers_ijk_field(np.zeros(shape=[24, 24, 72, 23]))
        sstr0_o = self.make_ntracers_ijk_field(np.zeros(shape=[24, 24, 72, 23]))
        trflx = self.make_ntracers_zdim_field(np.zeros(shape=[24, 24, 73, 23]))
        trsrc = self.make_ntracers_ijk_field(np.zeros(shape=[24, 24, 72, 23]))
        trten = self.make_ntracers_ijk_field(np.zeros(shape=[24, 24, 72, 23]))
        tru = self.make_ntracers_zdim_field(np.zeros(shape=[24, 24, 73, 23]))
        tru_emf = self.make_ntracers_zdim_field(np.zeros(shape=[24, 24, 73, 23]))
        qtu_emf = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        kinv = self.make_ijk_field(np.zeros(shape=[24, 24, 72]), dtype=Int)
        umf_out = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        dcm_out = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qvten_out = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qlten_out = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qiten_out = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        sten_out = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        uten_out = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        vten_out = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qrten_out = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qsten_out = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        cufrc_out = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        cush_inout = self.make_ij_field(inputs_reshaped["cush_inout"])
        qtflx_out = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        slflx_out = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        slflx = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        thlu_emf = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        qty_emf = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        uu_emf = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        vu_emf = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        uemf = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        uflx_out = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        vflx_out = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        fer_out = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        fdr_out = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        thvlavg = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        tkeavg = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        uavg = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        vavg = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        thvlmin = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qtavg = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        zmid0 = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qt0 = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        thvl0 = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        thvl0bot = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        dpi = self.make_ij_field(np.zeros(shape=[24, 24]))
        t0 = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qv0 = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        pmid0 = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        thl0 = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        thlsrc = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        usrc = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        vsrc = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        trsrc_o = self.make_ntracers_ijk_field(np.zeros(shape=[24, 24, 72, 23]))
        plcl = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        klcl = self.make_ijk_field(np.zeros(shape=[24, 24, 72]), dtype=Int)
        thl0lcl = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qt0lcl = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        thv0lcl = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        plfc = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        klfc = self.make_ijk_field(np.zeros(shape=[24, 24, 72]), dtype=Int)
        cin = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        thvubot = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        thvutop = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        thvlsrc = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        thl0top = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qt0top = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        thl0bot = self.make_ij_field(np.zeros(shape=[24, 24]))
        cin_IJ = self.make_ij_field(np.zeros(shape=[24, 24]))
        plfc_IJ = self.make_ij_field(np.zeros(shape=[24, 24]))
        klfc_IJ = self.make_ij_field(np.zeros(shape=[24, 24]), dtype=Int)
        cinlcl_IJ = self.make_ij_field(np.zeros(shape=[24, 24]))
        test_var_3D = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        test_var_2D = self.make_ij_field(np.zeros(shape=[24, 24]))
        qtsrc = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        ufrc = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        umf = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        wu = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        emf = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        thlu = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        qtu = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        thvu = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        uplus = self.make_ij_field(np.zeros(shape=[24, 24]))
        vplus = self.make_ij_field(np.zeros(shape=[24, 24]))
        uu = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        vu = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        tre = self.make_ntracers_ijk_field(np.zeros(shape=[24, 24, 72, 23]))
        uplus_3D = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        vplus_3D = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        cin_i = self.make_ij_field(np.zeros(shape=[24, 24]))
        cinlcl_i = self.make_ij_field(np.zeros(shape=[24, 24]))
        ke = self.make_ij_field(np.zeros(shape=[24, 24]))
        krel = self.make_ijk_field(np.zeros(shape=[24, 24, 72]), dtype=Int)
        prel = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        thv0rel = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        winv = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        cbmf = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        rho0inv = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        ufrcinv = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        tscaleh = self.make_ij_field(np.zeros(shape=[24, 24]))
        wlcl = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qsat_pe = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        thlue = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qtue = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        wue = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        wtwb = self.make_ij_field(np.zeros(shape=[24, 24]))
        umf_zint = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        emf = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        thvu = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        cush = self.make_ij_field(np.zeros(shape=[24, 24]))
        rei = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        pe = self.make_ij_field(np.zeros(shape=[24, 24]))
        thle = self.make_ij_field(np.zeros(shape=[24, 24]))
        qte = self.make_ij_field(np.zeros(shape=[24, 24]))
        dpe = self.make_ij_field(np.zeros(shape=[24, 24]))
        exne = self.make_ij_field(np.zeros(shape=[24, 24]))
        thvebot = self.make_ij_field(np.zeros(shape=[24, 24]))
        ue = self.make_ij_field(np.zeros(shape=[24, 24]))
        ve = self.make_ij_field(np.zeros(shape=[24, 24]))
        drage = self.make_ij_field(np.zeros(shape=[24, 24]))
        bogbot = self.make_ij_field(np.zeros(shape=[24, 24]))
        bogtop = self.make_ij_field(np.zeros(shape=[24, 24]))
        kpen_IJ = self.make_ij_field(np.zeros(shape=[24, 24]), dtype=Int)
        rhomid0j = self.make_ij_field(np.zeros(shape=[24, 24]))
        kpen = self.make_ijk_field(np.zeros(shape=[24, 24, 72]), dtype=Int)
        fer = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        ppen = self.make_ij_field(np.zeros(shape=[24, 24]))
        kbup_IJ = self.make_ij_field(np.zeros(shape=[24, 24]), dtype=Int)
        kbup = self.make_ijk_field(np.zeros(shape=[24, 24, 72]), dtype=Int)
        dwten = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        diten = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        thlu_top = self.make_ij_field(np.zeros(shape=[24, 24]))
        qtu_top = self.make_ij_field(np.zeros(shape=[24, 24]))
        cldhgt = self.make_ij_field(np.zeros(shape=[24, 24]))
        xflx = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        qtflx = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        uflx = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        ql0 = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qi0 = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        uten = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        vten = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        uf = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        vf = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        vflx = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        dwten_temp = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        diten_temp = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        umf_temp = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        qlubelow = self.make_ij_field(np.zeros(shape=[24, 24]))
        qiubelow = self.make_ij_field(np.zeros(shape=[24, 24]))
        qlj_2D = self.make_ij_field(np.zeros(shape=[24, 24]))
        qij_2D = self.make_ij_field(np.zeros(shape=[24, 24]))
        fdr = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qlten_sink = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qiten_sink = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qrten = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qsten = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        s0 = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qvten = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qlten = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        sten = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qiten = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qmin = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        trflx_d = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        trflx_u = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        ufrclcl = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qcubelow = self.make_ij_field(np.zeros(shape=[24, 24]))
        rcwp = self.make_ij_field(np.zeros(shape=[24, 24]))
        rlwp = self.make_ij_field(np.zeros(shape=[24, 24]))
        riwp = self.make_ij_field(np.zeros(shape=[24, 24]))
        qcu = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qlu = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qiu = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        cufrc = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        tr0_s = self.make_ntracers_ijk_field(np.zeros(shape=[24, 24, 72, 23]))
        umf_s = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        dcm = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        slflx_s = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        qtflx_s = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        uflx_s = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        vflx_s = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        xco = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qc = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qlten_det = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qiten_det = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        ufrc_s = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        qv0_s = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        ql0_s = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qi0_s = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        s0_s = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        t0_s = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        slten = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qv0_o = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        kinv_o = self.make_ijk_field(np.zeros(shape=[24, 24, 72]), dtype=Int)
        klcl_o = self.make_ijk_field(np.zeros(shape=[24, 24, 72]), dtype=Int)
        klfc_o = self.make_ijk_field(np.zeros(shape=[24, 24, 72]), dtype=Int)
        plcl_o = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        plfc_o = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        tkeavg_o = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        thvlmin_o = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        thvlmin_IJ = self.make_ij_field(np.zeros(shape=[24, 24]))
        qtsrc_o = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        thvlsrc_o = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        thlsrc_o = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        usrc_o = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        vsrc_o = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        thv0lcl_o = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        ql0_o = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qi0_o = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        t0_o = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        s0_o = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        u0_o = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        v0_o = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qt0_o = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        thl0_o = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        thvl0_o = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        ssthl0_o = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        ssqt0_o = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        thv0bot_o = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        thv0top_o = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        thvl0bot_o = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        thvl0top_o = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        ssu0_o = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        ssv0_o = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        dcm_s = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qvten_s = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qlten_s = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qiten_s = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        sten_s = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        uten_s = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        vten_s = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qrten_s = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qsten_s = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qldet_s = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qidet_s = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qlsub_s = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qisub_s = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        cush_s = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        cufrc_s = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        fer_s = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        fdr_s = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        wcrit = self.make_ij_field(np.zeros(shape=[24, 24]))
        alpha = self.make_ij_field(np.zeros(shape=[24, 24]))
        del_CIN = self.make_ij_field(np.zeros(shape=[24, 24]))
        qldet_out = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qidet_out = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qlsub_out = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qisub_out = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        ndrop_out = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        nice_out = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        umf_outvar = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        dcm_outvar = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qvten_outvar = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qlten_outvar = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qiten_outvar = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        sten_outvar = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        uten_outvar = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        vten_outvar = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qrten_outvar = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qsten_outvar = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        cufrc_outvar = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        cush_inoutvar = self.make_ij_field(np.zeros(shape=[24, 24]))
        qldet_outvar = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qidet_outvar = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qlsub_outvar = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qisub_outvar = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qtflx_outvar = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        slflx_outvar = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        uflx_outvar = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        vflx_outvar = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        fer_outvar = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        fdr_outvar = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        test_tracers = self.make_ntracers_ijk_field(np.zeros(shape=[24, 24, 72, 23]))
        # test_tracersIJ = self.make_ntracers_ij_field(np.zeros(shape=[24, 24, 23]))

        compute_uwshcu(
            windsrcavg=windsrcavg,
            qtsrchgt=qtsrchgt,
            qtsrc_fac=qtsrc_fac,
            thlsrc_fac=thlsrc_fac,
            frc_rasn=frc_rasn,
            rbuoy=rbuoy,
            epsvarw=epsvarw,
            use_CINcin=use_CINcin,
            mumin1=mumin1,
            rmaxfrac=rmaxfrac,
            PGFc=PGFc,
            k0=k0,
            dt=dt,
            ncnst=ncnst,
            pifc0_in=pifc0_in,
            zifc0_in=zifc0_in,
            exnifc0_in=exnifc0_in,
            pmid0_in=pmid0_in,
            zmid0_in=zmid0_in,
            exnmid0_in=exnmid0_in,
            dp0_in=dp0_in,
            u0_in=u0_in,
            v0_in=v0_in,
            qv0_in=qv0_in,
            ql0_in=ql0_in,
            qi0_in=qi0_in,
            th0_in=th0_in,
            tr0_inout=tr0_inout,
            kpbl_in=kpbl_in,
            frland_in=frland_in,
            tke_in=tke_in,
            rkfre=rkfre,
            cush_inout=cush_inout,
            shfx=shfx,
            evap=evap,
            dotransport=dotransport,
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
            qtflx_out=qtflx_out,
            slflx_out=slflx_out,
            uflx_out=uflx_out,
            vflx_out=vflx_out,
            # qtflx=qtflx_out,
            slflx=slflx,
            uflx=uflx,
            vflx=vflx,
            thlu_emf=thlu_emf,
            qtu_emf=qtu_emf,
            uu_emf=uu_emf,
            vu_emf=vu_emf,
            uemf=uemf,
            fer_out=fer_out,
            fdr_out=fdr_out,
            zmid0=zmid0,
            qt0=qt0,
            thvl0=thvl0,
            thvl0bot=thvl0bot,
            # Outputs for testing
            ssthl0=ssthl0,
            ssqt0=ssqt0,
            ssu0=ssu0,
            ssv0=ssv0,
            sstr0=sstr0,
            tr0=tr0,
            thj=thj,
            qij=qij,
            qlj=qlj,
            qse=qse,
            qvj=qvj,
            id_check=id_check,
            thv0top=thv0top,
            thvl0top=thvl0top,
            tr0_o=tr0_o,
            sstr0_o=sstr0_o,
            trflx=trflx,
            trten=trten,
            tru=tru,
            tru_emf=tru_emf,
            kinv=kinv,
            thvlavg=thvlavg,
            tkeavg=tkeavg,
            uavg=uavg,
            vavg=vavg,
            thvlmin=thvlmin,
            qtavg=qtavg,
            dpi=dpi,
            t0=t0,
            qv0=qv0,
            pmid0=pmid0,
            thl0=thl0,
            thlsrc=thlsrc,
            usrc=usrc,
            vsrc=vsrc,
            qtsrc=qtsrc,
            trsrc=trsrc,
            trsrc_o=trsrc_o,
            plcl=plcl,
            klcl=klcl,
            thl0lcl=thl0lcl,
            qt0lcl=qt0lcl,
            thv0lcl=thv0lcl,
            thv0bot=thv0bot,
            plfc=plfc,
            klfc=klfc,
            cin=cin,
            thvubot=thvubot,
            thvlsrc=thvlsrc,
            thl0top=thl0top,
            qt0top=qt0top,
            thvutop=thvutop,
            test_var3D=test_var_3D,
            test_var2D=test_var_2D,
            cin_IJ=cin_IJ,
            plfc_IJ=plfc_IJ,
            klfc_IJ=klfc_IJ,
            cinlcl_IJ=cinlcl_IJ,
            ufrc=ufrc,
            umf_zint=umf_zint,
            wu=wu,
            emf=emf,
            thlu=thlu,
            qtu=qtu,
            uplus=uplus,
            vplus=vplus,
            uu=uu,
            vu=vu,
            tre=tre,
            uplus_3D=uplus_3D,
            vplus_3D=vplus_3D,
            cin_i=cin_i,
            cinlcl_i=cinlcl_i,
            ke=ke,
            krel=krel,
            prel=prel,
            thv0rel=thv0rel,
            winv=winv,
            cbmf=cbmf,
            rho0inv=rho0inv,
            ufrcinv=ufrcinv,
            tscaleh=tscaleh,
            wlcl=wlcl,
            qsat_pe=qsat_pe,
            niter_xc=niter_xc,
            criqc=criqc,
            rle=rle,
            cridist_opt=cridist_opt,
            thlue=thlue,
            qtue=qtue,
            wue=wue,
            wtwb=wtwb,
            mixscale=mixscale,
            rkm=rkm,
            rdrag=rdrag,
            use_self_detrain=use_self_detrain,
            detrhgt=detrhgt,
            thvu=thvu,
            cush=cush,
            rei=rei,
            pe=pe,
            thle=thle,
            qte=qte,
            thvebot=thvebot,
            dpe=dpe,
            ue=ue,
            ve=ve,
            exne=exne,
            drage=drage,
            bogbot=bogbot,
            bogtop=bogtop,
            kpen_IJ=kpen_IJ,
            rhomid0j=rhomid0j,
            kpen=kpen,
            fer=fer,
            ppen=ppen,
            kbup_IJ=kbup_IJ,
            kbup=kbup,
            dwten=dwten,
            diten=diten,
            thlu_top=thlu_top,
            qtu_top=qtu_top,
            cldhgt=cldhgt,
            rpen=rpen,
            use_cumpenent=use_cumpenent,
            xflx=xflx,
            qtflx=qtflx,
            use_momenflx=use_momenflx,
            ql0=ql0,
            qi0=qi0,
            uten=uten,
            vten=vten,
            uf=uf,
            vf=vf,
            dwten_temp=dwten_temp,
            diten_temp=diten_temp,
            umf_temp=umf_temp,
            qlubelow=qlubelow,
            qiubelow=qiubelow,
            qlj_2D=qlj_2D,
            qij_2D=qij_2D,
            fdr=fdr,
            qlten_sink=qlten_sink,
            qiten_sink=qiten_sink,
            qrten=qrten,
            qsten=qsten,
            s0=s0,
            qvten=qvten,
            qlten=qlten,
            sten=sten,
            qiten=qiten,
            qmin=qmin,
            trflx_d=trflx_d,
            trflx_u=trflx_u,
            ufrclcl=ufrclcl,
            qcubelow=qcubelow,
            rcwp=rcwp,
            rlwp=rlwp,
            riwp=riwp,
            qcu=qcu,
            qlu=qlu,
            qiu=qiu,
            cufrc=cufrc,
            tr0_s=tr0_s,
            umf_s=umf_s,
            dcm=dcm,
            slflx_s=slflx_s,
            qtflx_s=qtflx_s,
            uflx_s=uflx_s,
            vflx_s=vflx_s,
            xco=xco,
            qc=qc,
            qlten_det=qlten_det,
            qiten_det=qiten_det,
            ufrc_s=ufrc_s,
            qv0_s=qv0_s,
            ql0_s=ql0_s,
            qi0_s=qi0_s,
            s0_s=s0_s,
            t0_s=t0_s,
            qv0_o=qv0_o,
            ql0_o=ql0_o,
            qi0_o=qi0_o,
            t0_o=t0_o,
            s0_o=s0_o,
            u0_o=u0_o,
            v0_o=v0_o,
            qt0_o=qt0_o,
            thl0_o=thl0_o,
            thvl0_o=thvl0_o,
            ssthl0_o=ssthl0_o,
            ssqt0_o=ssqt0_o,
            thv0bot_o=thv0bot_o,
            thv0top_o=thv0top_o,
            thvl0bot_o=thvl0bot_o,
            thvl0top_o=thvl0top_o,
            ssu0_o=ssu0_o,
            ssv0_o=ssv0_o,
            slten=slten,
            kinv_o=kinv_o,
            klcl_o=klcl_o,
            klfc_o=klfc_o,
            plcl_o=plcl_o,
            plfc_o=plfc_o,
            tkeavg_o=tkeavg_o,
            thvlmin_o=thvlmin_o,
            qtsrc_o=qtsrc_o,
            thvlsrc_o=thvlsrc_o,
            thlsrc_o=thlsrc_o,
            usrc_o=usrc_o,
            vsrc_o=vsrc_o,
            thv0lcl_o=thv0lcl_o,
            thvlmin_IJ=thvlmin_IJ,
            dcm_s=dcm_s,
            qvten_s=qvten_s,
            qlten_s=qlten_s,
            qiten_s=qiten_s,
            sten_s=sten_s,
            uten_s=uten_s,
            vten_s=vten_s,
            qrten_s=qrten_s,
            qsten_s=qsten_s,
            qldet_s=qldet_s,
            qidet_s=qidet_s,
            qlsub_s=qlsub_s,
            qisub_s=qisub_s,
            cush_s=cush_s,
            cufrc_s=cufrc_s,
            fer_s=fer_s,
            fdr_s=fdr_s,
            wcrit=wcrit,
            alpha=alpha,
            del_CIN=del_CIN,
            rdrop=rdrop,
            ndrop_out=ndrop_out,
            nice_out=nice_out,
            umf_outvar=umf_outvar,
            dcm_outvar=dcm_outvar,
            qvten_outvar=qvten_outvar,
            qlten_outvar=qlten_outvar,
            qiten_outvar=qiten_outvar,
            sten_outvar=sten_outvar,
            uten_outvar=uten_outvar,
            vten_outvar=vten_outvar,
            qrten_outvar=qrten_outvar,
            qsten_outvar=qsten_outvar,
            cufrc_outvar=cufrc_outvar,
            cush_inoutvar=cush_inoutvar,
            qldet_outvar=qldet_outvar,
            qldet_out=qldet_out,
            qidet_out=qidet_out,
            qidet_outvar=qidet_outvar,
            qlsub_out=qlsub_out,
            qlsub_outvar=qlsub_outvar,
            qisub_out=qisub_out,
            qisub_outvar=qisub_outvar,
            qtflx_outvar=qtflx_outvar,
            slflx_outvar=slflx_outvar,
            uflx_outvar=uflx_outvar,
            vflx_outvar=vflx_outvar,
            fer_outvar=fer_outvar,
            fdr_outvar=fdr_outvar,
            test_tracers=test_tracers,
        )
        print("Performed compute_uwshcu on reshaped inputs")
        print(test_tracers.view[23, 20, :, 22])
        # print("Reshaped outputs back to original shape")

        with xr.open_dataset("/Users/kfandric/netcdf/ComputeUwshcu-Out.nc") as ds:
            # Load in netcdf test var
            testvar = "sstr02"
            var = test_tracers
            testvar_nan = ds.variables[testvar].data[0, 0, :, :-1, :23, 0]
            # Replace nans with zero
            testvar_zeros = np.nan_to_num(testvar_nan, nan=0)

            # Reshape and make testvar quantity
            i, j, k = self.grid.nic, self.grid.njc, self.grid.npz
            testvar_reshaped = np.reshape(testvar_zeros, (i, j, k, 23)).astype(
                np.float32
            )
            # testvar_out = self.make_zinterface_field(testvar_reshaped)
            # testvar_out = self.make_ijk_field(testvar_reshaped)
            testvar_out = self.make_ntracers_ijk_field(testvar_reshaped)
            # testvar_out = self.make_ntracers_zdim_field(testvar_reshaped)

        # Run translate test by hand
        print("TESTING: " + testvar)
        tot_indicies = 0
        failed_indicies = 0
        testing_variable = var
        reference_variable = testvar_out

        for i in range(testing_variable.view[:].shape[0]):
            for j in range(testing_variable.view[:].shape[1]):
                for k in range(testing_variable.view[:].shape[2]):
                    for n in range(testing_variable.view[:].shape[3]):
                        diff = (
                            testing_variable.view[i, j, k, n]
                            - reference_variable.view[i, j, k, n]
                        )

                        rel_error = abs(diff / reference_variable.view[i, j, k, n])

                        tot_indicies = tot_indicies + 1
                        if (
                            testing_variable.view[i, j, k, n]
                            != reference_variable.view[i, j, k, n]
                        ):
                            if abs(diff) != 0.0:
                                failed_indicies = failed_indicies + 1
                                print(
                                    "DIFF: ",
                                    i,
                                    j,
                                    k,
                                    n,
                                    "computed: ",
                                    testing_variable.view[i, j, k, n],
                                    "reference: ",
                                    reference_variable.view[i, j, k, n],
                                    "abs difference: ",
                                    testing_variable.view[i, j, k, n]
                                    - reference_variable.view[i, j, k, n],
                                    "rel difference: ",
                                    diff / reference_variable.view[i, j, k, n],
                                )
        print(
            "failures: ",
            failed_indicies,
            "/",
            tot_indicies,
            "(",
            (failed_indicies / tot_indicies) * 100,
            "%)",
        )
