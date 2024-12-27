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
            "plfc_test": self.grid.compute_dict(),
            "cin_test": self.grid.compute_dict(),
            "thvubot_test": self.grid.compute_dict(),
            "thvlsrc": self.grid.compute_dict(),
            "thv0bot_test": self.grid.compute_dict(),
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
        ncnst = np.int64(inputs_reshaped["ncnst"])
        k0 = np.int64(inputs_reshaped["k0"])
        windsrcavg = np.int64(inputs_reshaped["windsrcavg"])
        qtsrchgt = inputs_reshaped["qtsrchgt"]
        qtsrc_fac = inputs_reshaped["qtsrc_fac"]
        thlsrc_fac = inputs_reshaped["thlsrc_fac"]
        frc_rasn = inputs_reshaped["frc_rasn"]
        rbuoy = inputs_reshaped["rbuoy"]
        epsvarw = inputs_reshaped["epsvarw"]
        use_CINcin = np.int64(inputs_reshaped["use_CINcin"])
        mumin1 = inputs_reshaped["mumin1"]
        rmaxfrac = inputs_reshaped["rmaxfrac"]
        PGFc = inputs_reshaped["PGFc"]
        dt = inputs_reshaped["dt"]

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
        thvl0top = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        tr0_o = self.make_ntracers_ijk_field(inputs_reshaped["tr0_inout"])
        sstr0_o = self.make_ntracers_ijk_field(inputs_reshaped["tr0_inout"])
        trflx = self.make_ntracers_zdim_field(np.zeros(shape=[24, 24, 73, 23]))
        trten = self.make_ntracers_ijk_field(inputs_reshaped["tr0_inout"])
        tru = self.make_ntracers_zdim_field(np.zeros(shape=[24, 24, 73, 23]))
        tru_emf = self.make_ntracers_zdim_field(np.zeros(shape=[24, 24, 73, 23]))
        kinv = self.make_ijk_field(inputs_reshaped["pmid0_in"], dtype=Int)
        umf_out = self.make_zinterface_field(inputs_reshaped["tke_in"])
        dcm_out = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        qvten_out = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        qlten_out = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        qiten_out = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        sten_out = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        uten_out = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        vten_out = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        qrten_out = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        qsten_out = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        cufrc_out = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        cush_inout = self.make_ij_field(inputs_reshaped["cush_inout"])
        qtflx_out = self.make_zinterface_field(inputs_reshaped["tke_in"])
        slflx_out = self.make_zinterface_field(inputs_reshaped["tke_in"])
        uflx_out = self.make_zinterface_field(inputs_reshaped["tke_in"])
        vflx_out = self.make_zinterface_field(inputs_reshaped["tke_in"])
        fer_out = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        fdr_out = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        thvlavg = self.make_ij_field(np.zeros(shape=[24, 24]))
        tkeavg = self.make_ij_field(np.zeros(shape=[24, 24]))
        uavg = self.make_ij_field(np.zeros(shape=[24, 24]))
        vavg = self.make_ij_field(np.zeros(shape=[24, 24]))
        thvlmin = self.make_ij_field(np.zeros(shape=[24, 24]))
        qtavg = self.make_ij_field(np.zeros(shape=[24, 24]))
        zmid0 = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        qt0 = self.make_ijk_field(np.zeros(shape=[24, 24, 72]))
        thvl0 = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        thvl0bot = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        dpi = self.make_ij_field(np.zeros(shape=[24, 24]))
        t0 = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        qv0 = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        pmid0 = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        thl0 = self.make_ijk_field(inputs_reshaped["pmid0_in"])
        thlsrc = self.make_ij_field(np.zeros(shape=[24, 24]))
        usrc = self.make_ij_field(np.zeros(shape=[24, 24]))
        vsrc = self.make_ij_field(np.zeros(shape=[24, 24]))
        trsrc = self.make_ntracers_ijk_field(np.zeros(shape=[24, 24, 72, 23]))
        plcl = self.make_ij_field(np.zeros(shape=[24, 24]))
        klcl = self.make_ijk_field(np.zeros(shape=[24, 24, 72]), dtype=Int)
        thl0lcl = self.make_ij_field(np.zeros(shape=[24, 24]))
        qt0lcl = self.make_ij_field(np.zeros(shape=[24, 24]))
        thv0lcl = self.make_ij_field(np.zeros(shape=[24, 24]))
        plfc = self.make_ij_field(np.zeros(shape=[24, 24]))
        cin = self.make_ij_field(np.zeros(shape=[24, 24]))
        thvubot = self.make_ij_field(np.zeros(shape=[24, 24]))
        thvlsrc = self.make_ij_field(np.zeros(shape=[24, 24]))

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
            trsrc=trsrc,
            plcl=plcl,
            klcl=klcl,
            thl0lcl=thl0lcl,
            qt0lcl=qt0lcl,
            thv0lcl=thv0lcl,
            thv0bot=thv0bot,
            plfc=plfc,
            cin=cin,
            thvubot=thvubot,
            thvlsrc=thvlsrc,
        )
        print("Performed compute_uwshcu on reshaped inputs")
        print(cin.view[:])
        # print("Reshaped outputs back to original shape")

        with xr.open_dataset("/Users/kfandric/netcdf/ComputeUwshcu-Out.nc") as ds:
            # Load in netcdf test var
            testvar = "cin"
            var = cin
            testvar_nan = ds.variables[testvar].data[0, 0, :, 0, 0, 0]
            # Replace nans with zero
            testvar_zeros = np.nan_to_num(testvar_nan, nan=0)

            # Reshape and make testvar quantity
            i, j, k = self.grid.nic, self.grid.njc, self.grid.npz
            testvar_reshaped = np.reshape(testvar_zeros, (i, j)).astype(np.float32)
            # testvar_out = self.make_zinterface_field(testvar_reshaped)
            testvar_out = self.make_ij_field(testvar_reshaped)
            # testvar_out = self.make_ntracers_ijk_field(testvar_reshaped)
            # testvar_out = self.make_ntracers_zdim_field(testvar_reshaped)

        # Run translate test by hand
        print("TESTING: " + testvar)
        tot_indicies = 0
        failed_indicies = 0
        testing_variable = var
        reference_variable = testvar_out
        # print(testing_variable.shape)
        # print(reference_variable.shape)

        for i in range(testing_variable.view[:].shape[0]):
            for j in range(testing_variable.view[:].shape[1]):
                # for k in range(testing_variable.view[:].shape[2]):
                diff = testing_variable.view[i, j] - reference_variable.view[i, j]

                tot_indicies = tot_indicies + 1
                if testing_variable.view[i, j] != reference_variable.view[i, j]:
                    failed_indicies = failed_indicies + 1
                    print(
                        "DIFF: ",
                        i,
                        j,
                        "computed: ",
                        testing_variable.view[i, j],
                        "reference: ",
                        reference_variable.view[i, j],
                        "difference: ",
                        testing_variable.view[i, j] - reference_variable.view[i, j],
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

        # return {
        #     "thvlavg": tkeavg_4D,
        # }
