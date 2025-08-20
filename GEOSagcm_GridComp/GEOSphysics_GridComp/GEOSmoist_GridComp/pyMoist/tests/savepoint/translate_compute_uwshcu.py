from ndsl import Namelist, Quantity, QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Float, Int, FloatField
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array
from pyMoist.UW.compute_uwshcu import ComputeUwshcuInv
from pyMoist.UW.config import UWConfiguration


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

        self.ntracers_quantity_factory = ComputeUwshcuInv.make_ntracers_quantity_factory(
            self.quantity_factory,
        )

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
            "cnvtr": {},
            "CNV_Tracers": {},
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
            "iter_cin",
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
            "tpert_out": self.grid.compute_dict(),
            "qpert_out": self.grid.compute_dict(),
        }

    def make_ntracers_ijk_field(self, data) -> Quantity:
        qty = self.ntracers_quantity_factory.empty(
            [X_DIM, Y_DIM, Z_DIM, "ntracers"],
            "n/a",
        )
        qty.view[:, :, :, :] = qty.np.asarray(data[:, :, :, :])
        return qty

    def make_ijk_field(self, data, dtype=FloatField) -> Quantity:
        qty = self.quantity_factory.empty([X_DIM, Y_DIM, Z_DIM], "n/a", dtype=dtype)
        qty.view[:, :, :] = qty.np.asarray(data[:, :, :])
        return qty

    def make_ij_field(self, data, dtype=FloatField) -> Quantity:
        qty = self.quantity_factory.empty([X_DIM, Y_DIM], "n/a", dtype=dtype)
        qty.view[:, :] = qty.np.asarray(data[:, :])
        return qty

    def compute(self, inputs):
        self.UW_config = UWConfiguration(Int(inputs["ncnst"]), Int(inputs["k0"]), Int(inputs["windsrcavg"]))

        compute_uwshcu = ComputeUwshcuInv(
            self.stencil_factory,
            self.grid.quantity_factory,
            self.UW_config,
        )

        # Float/Int Inputs
        dotransport = Int(inputs["dotransport"])
        k0 = Int(inputs["k0"])
        windsrcavg = Int(inputs["windsrcavg"])
        qtsrchgt = Float(inputs["qtsrchgt"])
        qtsrc_fac = Float(inputs["qtsrc_fac"])
        thlsrc_fac = Float(inputs["thlsrc_fac"])
        frc_rasn = Float(inputs["frc_rasn"])
        rbuoy = Float(inputs["rbuoy"])
        epsvarw = Float(inputs["epsvarw"])
        use_CINcin = Int(inputs["use_CINcin"])
        mumin1 = Float(inputs["mumin1"])
        rmaxfrac = Float(inputs["rmaxfrac"])
        PGFc = Float(inputs["PGFc"])
        dt = Float(inputs["dt"])
        niter_xc = Int(inputs["niter_xc"])
        criqc = Float(inputs["criqc"])
        rle = Float(inputs["rle"])
        cridist_opt = Int(inputs["cridist_opt"])
        mixscale = Float(inputs["mixscale"])
        rdrag = Float(inputs["rdrag"])
        rkm = Float(inputs["rkm"])
        use_self_detrain = Int(inputs["use_self_detrain"])
        detrhgt = Float(inputs["detrhgt"])
        use_cumpenent = Int(inputs["use_cumpenent"])
        rpen = Float(inputs["rpen"])
        use_momenflx = Int(inputs["use_momenflx"])
        rdrop = Float(inputs["rdrop"])
        iter_cin = Int(inputs["iter_cin"])

        # Field inputs
        pifc0_inv = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a"
        )
        safe_assign_array(pifc0_inv.view[:, :, :], inputs["pifc0_inv"])
        zifc0_inv = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a"
        )
        safe_assign_array(zifc0_inv.view[:, :, :], inputs["zifc0_inv"])
        pmid0_inv = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(pmid0_inv.view[:, :, :], inputs["pmid0_inv"])
        zmid0_inv = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(zmid0_inv.view[:, :, :], inputs["zmid0_inv"])
        kpbl_inv = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM], units="n/a"
        )
        safe_assign_array(kpbl_inv.view[:, :], inputs["kpbl_inv"])
        exnmid0_inv = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(exnmid0_inv.view[:, :, :], inputs["exnmid0_inv"])
        exnifc0_inv = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a"
        )
        safe_assign_array(exnifc0_inv.view[:, :, :], inputs["exnifc0_inv"])
        dp0_inv = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(dp0_inv.view[:, :, :], inputs["dp0_inv"])
        u0_inv = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(u0_inv.view[:, :, :], inputs["u0_inv"])
        v0_inv = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(v0_inv.view[:, :, :], inputs["v0_inv"])
        qv0_inv = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(qv0_inv.view[:, :, :], inputs["qv0_inv"])
        ql0_inv = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(ql0_inv.view[:, :, :], inputs["ql0_inv"])
        qi0_inv = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(qi0_inv.view[:, :, :], inputs["qi0_inv"])
        t0_inv = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(t0_inv.view[:, :, :], inputs["t0_inv"])
        frland = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM], units="n/a")
        safe_assign_array(frland.view[:, :], inputs["frland"])
        tke_inv = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a"
        )
        safe_assign_array(tke_inv.view[:, :, :], inputs["tke_inv"])
        rkfre = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM], units="n/a")
        safe_assign_array(rkfre.view[:, :], inputs["rkfre"])
        cush = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM], units="n/a")
        safe_assign_array(cush.view[:, :], inputs["cush"])
        shfx = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM], units="n/a")
        safe_assign_array(shfx.view[:, :], inputs["shfx"])
        evap = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM], units="n/a")
        safe_assign_array(evap.view[:, :], inputs["evap"])
        cnvtr = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM], units="n/a")
        safe_assign_array(cnvtr.view[:, :], inputs["cnvtr"])

        CNV_Tracers = self.make_ntracers_ijk_field(inputs["CNV_Tracers"])

        # Outputs
        # Z_interface fields
        umf_inv = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a"
        )
        qtflx_inv = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a"
        )
        slflx_inv = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a"
        )
        uflx_inv = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a"
        )
        vflx_inv = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a"
        )

        # FloatFields
        dcm_inv = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        qvten_inv = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        qlten_inv = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        qiten_inv = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        tten_inv = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        uten_inv = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        vten_inv = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        qrten_inv = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        qsten_inv = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        cufrc_inv = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        fer_inv = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        fdr_inv = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        ndrop_inv = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        nice_inv = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        qldet_inv = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        qlsub_inv = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        qidet_inv = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        qisub_inv = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")

        # FloatFieldIJs
        tpert_out = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM], units="n/a"
        )
        qpert_out = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM], units="n/a"
        )

        # Test variables
        testvar3D = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a"
        )

        testvar2D = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM], units="n/a"
        )

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
            cnvtr=cnvtr,
            CNV_Tracers=CNV_Tracers,
            # Float/Int inputs
            dotransport=dotransport,
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
            iter_cin=iter_cin,
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
            tpert_out=tpert_out,
            qpert_out=qpert_out,
            # Testvars
            # testvar3D=testvar3D,
            # testvar2D=testvar2D,
        )

        # with xr.open_dataset("/Users/kfandric/netcdf/ComputeUwshcu-Out2.nc") as ds:
        #     # Load in netcdf test var
        #     testvar = "thl0bot"
        #     var = testvar3D

        #     testvar_nan = ds.variables[testvar].data[0, 1, :, :-1, 0, 0]
        #     # Replace nans with zero
        #     testvar_zeros = np.nan_to_num(testvar_nan, nan=0)

        #     # Reshape and make testvar quantity
        #     i, j, k = self.grid.nic, self.grid.njc, self.grid.npz
        #     testvar_reshaped = np.reshape(testvar_zeros, (i, j, k)).astype(np.float32)
        #     testvar_transposed = testvar_reshaped.transpose(1, 0, 2)
        #     # testvar_out = self.make_zinterface_field(testvar_reshaped)
        #     testvar_out = self.make_ijk_field(testvar_transposed)
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

        # print(testvar3D.view[15, 16, :])
        # print(reference_variable.view[15, 16, :])

        # sys.exit()

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
            "tpert_out": tpert_out.view[:],
            "qpert_out": qpert_out.view[:],
            "dotransport": dotransport,
        }
