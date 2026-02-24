from f90nml import Namelist

from ndsl import Quantity, QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import FloatField
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array
from pyMoist.UW.compute_uwshcu import ComputeUwshcuInv
from pyMoist.UW.config import UWConfiguration


# Dev NOTE: The data for this translate test comes from combining several ncfiles
# Run ComputeUwshcuInv_postprocess_netcdfs.py before running translate test.


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

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "PLE": {},
            "ZLE": {},
            "QLLS": {},
            "QILS": {},
            "QLCN": {},
            "QICN": {},
            "kpbl_inv": {},
            "u0_inv": {},
            "v0_inv": {},
            "qv0_inv": {},
            "t0_inv": {},
            "frland": {},
            "tke_inv": {},
            "cush": {},
            "shfx": {},
            "evap": {},
            "cnvtr": {},
            "CNV_Tracers": {},
        }

        # FloatField Outputs
        self.out_vars = {
            "CNV_Tracers": self.grid.compute_dict(),
            "CNPCRATE": self.grid.compute_dict(),
            "RKFRE": self.grid.compute_dict(),
            "DQADT_SC": self.grid.compute_dict(),
            "MFD_SC": self.grid.compute_dict(),
            "PTR3D": self.grid.compute_dict(),
            "Q": self.grid.compute_dict(),
            "QIDET_SC": self.grid.compute_dict(),
            "QIENT_SC": self.grid.compute_dict(),
            "QLDET_SC": self.grid.compute_dict(),
            "QLENT_SC": self.grid.compute_dict(),
            "T": self.grid.compute_dict(),
            "U": self.grid.compute_dict(),
            "V": self.grid.compute_dict(),
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
            "cush": self.grid.compute_dict(),
        }

    def make_ntracers_ijk_field(self, data) -> Quantity:
        qty = self.quantity_factory.empty(
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

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("ComputeUwshcuInv-constants")

    def compute(self, inputs):
        config = UWConfiguration(**self.constants)

        compute_uwshcu = ComputeUwshcuInv(
            self.stencil_factory,
            self.quantity_factory,
            config,
        )

        # Field inputs
        PLE = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(PLE.view[:, :, :], inputs["PLE"])
        ZLE = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(ZLE.view[:, :, :], inputs["ZLE"])
        QLLS = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(QLLS.view[:, :, :], inputs["QLLS"])
        QILS = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(QILS.view[:, :, :], inputs["QILS"])
        QLCN = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(QLCN.view[:, :, :], inputs["QLCN"])
        QICN = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(QICN.view[:, :, :], inputs["QICN"])
        kpbl_inv = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM], units="n/a")
        safe_assign_array(kpbl_inv.view[:, :], inputs["kpbl_inv"])
        u0_inv = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(u0_inv.view[:, :, :], inputs["u0_inv"])
        v0_inv = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(v0_inv.view[:, :, :], inputs["v0_inv"])
        qv0_inv = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(qv0_inv.view[:, :, :], inputs["qv0_inv"])
        t0_inv = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(t0_inv.view[:, :, :], inputs["t0_inv"])
        frland = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM], units="n/a")
        safe_assign_array(frland.view[:, :], inputs["frland"])
        tke_inv = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a"
        )
        safe_assign_array(tke_inv.view[:, :, :], inputs["tke_inv"])
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
        DQADT_SC = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        MFD_SC = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        PTR3D = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        QLENT_SC = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        QIENT_SC = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")

        # FloatFieldIJs
        RKFRE = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM], units="n/a")
        tpert_out = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM], units="n/a")
        qpert_out = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM], units="n/a")

        compute_uwshcu(
            # Field inputs
            PLE=PLE,
            ZLE=ZLE,
            QLLS=QLLS,
            QILS=QILS,
            QLCN=QLCN,
            QICN=QICN,
            kpbl_inv=kpbl_inv,
            u0_inv=u0_inv,
            v0_inv=v0_inv,
            qv0_inv=qv0_inv,
            t0_inv=t0_inv,
            frland=frland,
            tke_inv=tke_inv,
            cush=cush,
            shfx=shfx,
            evap=evap,
            cnvtr=cnvtr,
            CNV_Tracers=CNV_Tracers,
            # Outputs
            RKFRE=RKFRE,
            MFD_SC=MFD_SC,
            DQADT_SC=DQADT_SC,
            PTR3D=PTR3D,
            QLENT_SC=QLENT_SC,
            QIENT_SC=QIENT_SC,
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
        )

        return {
            "CNV_Tracers": CNV_Tracers.view[:],
            "CNPCRATE": cnvtr.view[:],
            "RKFRE": RKFRE.view[:],
            "Q": qv0_inv.view[:],
            "T": t0_inv.view[:],
            "U": u0_inv.view[:],
            "V": v0_inv.view[:],
            "QIDET_SC": qidet_inv.view[:],
            "QLDET_SC": qldet_inv.view[:],
            "MFD_SC": MFD_SC.view[:],
            "DQADT_SC": DQADT_SC.view[:],
            "PTR3D": PTR3D.view[:],
            "QLENT_SC": QLENT_SC.view[:],
            "QIENT_SC": QIENT_SC.view[:],
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
            "dotransport": config.dotransport,
            "cush": cush.view[:],
        }
