from ndsl import Namelist, StencilFactory, Quantity
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.UW.compute_uwshcu_inv import ComputeUwshcuInv
import numpy as np
from pyMoist.pyMoist_constants import ncnst
from ndsl.dsl.typing import (
    FloatField,
    Int,
)
import xarray as xr


class TranslateComputeUwshcuInv(TranslateFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        self.ntracers_quantity_factory = (
            ComputeUwshcuInv.make_ntracers_quantity_factory(self.quantity_factory)
        )

        self.max_error = 1e-9

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "pifc0_inv": {},
            "zifc0_inv": {},
            "pmid0_inv": {},
            "zmid0_inv": {},
            "kpbl_inv": {},
            "exnmid0_inv": {},
            "exnifc0_inv": {},
            "dp0_inv": {},
            "u0_inv": {},
            "v0_inv": {},
            "qv0_inv": {},
            "ql0_inv": {},
            "qi0_inv": {},
            "t0_inv": {},
            "frland": {},
            "tke_inv": {},
            "rkfre": {},
            "cush": {},
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
            "vten_inv": self.grid.compute_dict(),
            "cufrc_inv": self.grid.compute_dict(),
            "qlten_inv": self.grid.compute_dict(),
            "fer_inv": self.grid.compute_dict(),
            "qrten_inv": self.grid.compute_dict(),
            "qisub_inv": self.grid.compute_dict(),
            "tten_inv": self.grid.compute_dict(),
            "qiten_inv": self.grid.compute_dict(),
            "uten_inv": self.grid.compute_dict(),
            "dcm_inv": self.grid.compute_dict(),
            "vflx_inv": self.grid.compute_dict(),
            "qldet_inv": self.grid.compute_dict(),
            "umf_inv": self.grid.compute_dict(),
            "nice_inv": self.grid.compute_dict(),
            "fdr_inv": self.grid.compute_dict(),
            "qidet_inv": self.grid.compute_dict(),
            "qlsub_inv": self.grid.compute_dict(),
            "qtflx_inv": self.grid.compute_dict(),
            "uflx_inv": self.grid.compute_dict(),
            "slflx_inv": self.grid.compute_dict(),
            "dotransport": self.grid.compute_dict(),
            "ndrop_inv": self.grid.compute_dict(),
            "qvten_inv": self.grid.compute_dict(),
            "qsten_inv": self.grid.compute_dict(),
            # "tpert_inv": self.grid.compute_dict(),
            # "qpert_inv": self.grid.compute_dict(),
        }

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
        compute_uwshcu = ComputeUwshcuInv(
            self.stencil_factory,
            self.grid.quantity_factory,
        )

        # Float/Int Inputs
        dotransport = np.float32(inputs["dotransport"])
        ncnst = Int(inputs["ncnst"])
        k0 = Int(inputs["k0"])
        windsrcavg = Int(inputs["windsrcavg"])
        qtsrchgt = np.float32(inputs["qtsrchgt"])
        qtsrc_fac = np.float32(inputs["qtsrc_fac"])
        thlsrc_fac = np.float32(inputs["thlsrc_fac"])
        frc_rasn = np.float32(inputs["frc_rasn"])
        rbuoy = np.float32(inputs["rbuoy"])
        epsvarw = np.float32(inputs["epsvarw"])
        use_CINcin = np.int32(inputs["use_CINcin"])
        mumin1 = np.float32(inputs["mumin1"])
        rmaxfrac = np.float32(inputs["rmaxfrac"])
        PGFc = np.float32(inputs["PGFc"])
        dt = np.float32(inputs["dt"])
        niter_xc = Int(inputs["niter_xc"])
        criqc = np.float32(inputs["criqc"])
        rle = np.float32(inputs["rle"])
        cridist_opt = np.int32(inputs["cridist_opt"])
        mixscale = np.float32(inputs["mixscale"])
        rdrag = np.float32(inputs["rdrag"])
        rkm = np.float32(inputs["rkm"])
        use_self_detrain = np.int32(inputs["use_self_detrain"])
        detrhgt = np.float32(inputs["detrhgt"])
        use_cumpenent = np.int32(inputs["use_cumpenent"])
        rpen = np.float32(inputs["rpen"])
        use_momenflx = np.int32(inputs["use_momenflx"])
        rdrop = np.float32(inputs["rdrop"])

        # Field inputs
        pifc0_inv = self.make_zinterface_field(inputs["pifc0_inv"])
        zifc0_inv = self.make_zinterface_field(inputs["zifc0_inv"])
        pmid0_inv = self.make_ijk_field(inputs["pmid0_inv"])
        zmid0_inv = self.make_ijk_field(inputs["zmid0_inv"])
        kpbl_inv = self.make_ij_field(inputs["kpbl_inv"], dtype=Int)
        exnmid0_inv = self.make_ijk_field(inputs["exnmid0_inv"])
        exnifc0_inv = self.make_zinterface_field(inputs["exnifc0_inv"])
        dp0_inv = self.make_ijk_field(inputs["dp0_inv"])
        u0_inv = self.make_ijk_field(inputs["u0_inv"])
        v0_inv = self.make_ijk_field(inputs["v0_inv"])
        qv0_inv = self.make_ijk_field(inputs["qv0_inv"])
        ql0_inv = self.make_ijk_field(inputs["ql0_inv"])
        qi0_inv = self.make_ijk_field(inputs["qi0_inv"])
        t0_inv = self.make_ijk_field(inputs["t0_inv"])
        frland = self.make_ij_field(inputs["frland"])
        tke_inv = self.make_zinterface_field(inputs["tke_inv"])
        rkfre = self.make_ij_field(inputs["rkfre"])
        cush = self.make_ij_field(inputs["cush"])
        shfx = self.make_ij_field(inputs["shfx"])
        evap = self.make_ij_field(inputs["evap"])
        dotransport = np.float32(inputs["dotransport"])

        # Outputs
        # Z_interface fields
        umf_inv = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        qtflx_inv = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        slflx_inv = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        uflx_inv = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))
        vflx_inv = self.make_zinterface_field(np.zeros(shape=[24, 24, 73]))

        # FloatFields
        dcm_inv = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qvten_inv = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qlten_inv = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qiten_inv = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        tten_inv = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        uten_inv = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        vten_inv = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qrten_inv = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qsten_inv = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        cufrc_inv = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        fer_inv = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        fdr_inv = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        ndrop_inv = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        nice_inv = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qldet_inv = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qlsub_inv = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qidet_inv = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qisub_inv = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))

        # Testing variables
        test_var_3D = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        test_var_2D = self.make_ij_field(np.zeros(shape=[24, 24]))
        test_tracers = self.make_ntracers_ijk_field(np.zeros(shape=[24, 24, 72, 23]))
        test_tracersIJ = self.make_ntracers_ij_field(np.zeros(shape=[24, 24, 23]))

        compute_uwshcu(
            # Field inputs
            pifc0_inv=pifc0_inv,
            zifc0_inv=zifc0_inv,
            pmid0_inv=pmid0_inv,
            zmid0_inv=zmid0_inv,
            kpbl_inv=kpbl_inv,
            exnmid0_inv=exnmid0_inv,
            exnifc0_inv=exnifc0_inv,
            dp0_inv=dp0_inv,
            u0_inv=u0_inv,
            v0_inv=v0_inv,
            qv0_inv=qv0_inv,
            ql0_inv=ql0_inv,
            qi0_inv=qi0_inv,
            t0_inv=t0_inv,
            frland=frland,
            tke_inv=tke_inv,
            rkfre=rkfre,
            cush=cush,
            shfx=shfx,
            evap=evap,
            # Float/Int inputs
            dotransport=dotransport,
            ncnst=ncnst,
            k0=k0,
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
            dt=dt,
            niter_xc=niter_xc,
            criqc=criqc,
            rle=rle,
            cridist_opt=cridist_opt,
            mixscale=mixscale,
            rdrag=rdrag,
            rkm=rkm,
            use_self_detrain=use_self_detrain,
            detrhgt=detrhgt,
            use_cumpenent=use_cumpenent,
            rpen=rpen,
            use_momenflx=use_momenflx,
            rdrop=rdrop,
            # Outputs
            umf_inv=umf_inv,
            dcm_inv=dcm_inv,
            qtflx_inv=qtflx_inv,
            slflx_inv=slflx_inv,
            uflx_inv=uflx_inv,
            vflx_inv=vflx_inv,
            qvten_inv=qvten_inv,
            qlten_inv=qlten_inv,
            qiten_inv=qiten_inv,
            tten_inv=tten_inv,
            uten_inv=uten_inv,
            vten_inv=vten_inv,
            qrten_inv=qrten_inv,
            qsten_inv=qsten_inv,
            cufrc_inv=cufrc_inv,
            fer_inv=fer_inv,
            fdr_inv=fdr_inv,
            ndrop_inv=ndrop_inv,
            nice_inv=nice_inv,
            qldet_inv=qldet_inv,
            qlsub_inv=qlsub_inv,
            qidet_inv=qidet_inv,
            qisub_inv=qisub_inv,
            # Test vars
            test_var3D=test_var_3D,
            test_var2D=test_var_2D,
            test_tracers=test_tracers,
            test_tracersIJ=test_tracersIJ,
        )

        # with xr.open_dataset("/Users/kfandric/netcdf/ComputeUwshcu-Out2.nc") as ds:
        #     # Load in netcdf test var
        #     testvar = "dcm_out1"
        #     var = test_var_3D
        #     testvar_nan = ds.variables[testvar].data[0, 0, :, :]
        #     # Replace nans with zero
        #     testvar_zeros = np.nan_to_num(testvar_nan, nan=0)

        #     # Reshape and make testvar quantity
        #     i, j, k = self.grid.nic, self.grid.njc, self.grid.npz
        #     testvar_reshaped = np.reshape(testvar_zeros, (i, j, k)).astype(np.float32)
        #     testvar_trans = testvar_reshaped.transpose(1, 0, 2)
        #     # testvar_out = self.make_zinterface_field(testvar_reshaped)
        #     testvar_out = self.make_ijk_field(testvar_trans)
        #     # testvar_out = self.make_ntracers_ijk_field(testvar_reshaped)
        #     # testvar_out = self.make_ntracers_zdim_field(testvar_reshaped)

        # # Run translate test by hand
        # print("TESTING: " + testvar)
        # tot_indicies = 0
        # failed_indicies = 0
        # testing_variable = var
        # reference_variable = testvar_out

        # for i in range(testing_variable.view[:].shape[0]):
        #     for j in range(testing_variable.view[:].shape[1]):
        #         for k in range(testing_variable.view[:].shape[2]):
        #             # for n in range(testing_variable.view[:].shape[3]):
        #             diff = (
        #                 testing_variable.view[i, j, k]
        #                 - reference_variable.view[i, j, k]
        #             )

        #             rel_error = abs(diff / reference_variable.view[i, j, k])

        #             tot_indicies = tot_indicies + 1
        #             if (
        #                 testing_variable.view[i, j, k]
        #                 != reference_variable.view[i, j, k]
        #             ):
        #                 if abs(diff) != 0.0:
        #                     failed_indicies = failed_indicies + 1
        #                     print(
        #                         "DIFF: ",
        #                         i,
        #                         j,
        #                         k,
        #                         "computed: ",
        #                         testing_variable.view[i, j, k],
        #                         "reference: ",
        #                         reference_variable.view[i, j, k],
        #                         "abs difference: ",
        #                         testing_variable.view[i, j, k]
        #                         - reference_variable.view[i, j, k],
        #                         "rel difference: ",
        #                         diff / reference_variable.view[i, j, k],
        #                     )
        # print(
        #     "failures: ",
        #     failed_indicies,
        #     "/",
        #     tot_indicies,
        #     "(",
        #     (failed_indicies / tot_indicies) * 100,
        #     "%)",
        # )

        return {
            "umf_inv": umf_inv.view[:],
            "dcm_inv": dcm_inv.view[:],
            "qtflx_inv": qtflx_inv.view[:],
            "slflx_inv": slflx_inv.view[:],
            "uflx_inv": uflx_inv.view[:],
            "vflx_inv": vflx_inv.view[:],
            "qvten_inv": qvten_inv.view[:],
            "qlten_inv": qlten_inv.view[:],
            "qiten_inv": qiten_inv.view[:],
            "tten_inv": tten_inv.view[:],
            "uten_inv": uten_inv.view[:],
            "vten_inv": vten_inv.view[:],
            "qrten_inv": qrten_inv.view[:],
            "qsten_inv": qsten_inv.view[:],
            "cufrc_inv": cufrc_inv.view[:],
            "fer_inv": fer_inv.view[:],
            "fdr_inv": fdr_inv.view[:],
            "ndrop_inv": ndrop_inv.view[:],
            "nice_inv": nice_inv.view[:],
            "qldet_inv": qldet_inv.view[:],
            "qlsub_inv": qlsub_inv.view[:],
            "qidet_inv": qidet_inv.view[:],
            "qisub_inv": qisub_inv.view[:],
            "dotransport": dotransport,
        }
