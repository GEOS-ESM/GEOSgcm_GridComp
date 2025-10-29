from ndsl import Namelist, QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float, Int
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array
from pyMoist.convection.GF_2020.GF_functions.cloud_dissipation import (
    CloudDissipation,
)


class TranslateCloudDissipation(TranslateFortranData2Py):
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
            "COUPL_MPHYSICS": {},
            "dtime": {},
            "heso_cup": {},
            "ierr": {},
            "kbcon": {},
            "ktop": {},
            "outq": {},
            "outqc": {},
            "outt": {},
            "qeso_cup": {},
            "qo_cup": {},
            "qrco": {},
            "rho_hydr": {},
            "sig": {},
            "tempco": {},
            "tn_cup": {},
            "use_cloud_dissipation": {},
            "vvel2d": {},
            "xmb": {},
            "zo": {},
            "zuo": {},
        }

        # Float/Int Inputs
        self.in_vars["parameters"] = []

        # FloatField Outputs
        self.out_vars = {
            "outq": self.grid.compute_dict(),
            "outqc": self.grid.compute_dict(),
            "outt": self.grid.compute_dict(),
            "qrco": self.grid.compute_dict(),
        }

    def compute(self, inputs):

        cloud_dissipation = CloudDissipation(
            self.stencil_factory,
            self.grid.quantity_factory,
        )

        # # Field inputs

        COUPL_MPHYSICS = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=Int,
        )
        safe_assign_array(COUPL_MPHYSICS.view[:, :, :], inputs["COUPL_MPHYSICS"])

        dtime = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a"
        )
        safe_assign_array(dtime.view[:, :, :], inputs["dtime"])

        heso_cup = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a"
        )
        safe_assign_array(heso_cup.view[:, :, :], inputs["heso_cup"])

        ierr = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=Int,
        )
        safe_assign_array(ierr.view[:, :, :], inputs["ierr"])

        kbcon = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int
        )
        safe_assign_array(kbcon.view[:, :, :], inputs["kbcon"])

        ktop = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int
        )
        safe_assign_array(ktop.view[:, :, :], inputs["ktop"])

        outq = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(outq.view[:, :, :], inputs["outq"])

        outqc = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(outqc.view[:, :, :], inputs["outqc"])

        outt = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(outt.view[:, :, :], inputs["outt"])

        qeso_cup = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(qeso_cup.view[:, :, :], inputs["qeso_cup"])

        qo_cup = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(qo_cup.view[:, :, :], inputs["qo_cup"])

        qrco = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(qrco.view[:, :, :], inputs["qrco"])

        rho_hydr = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(rho_hydr.view[:, :, :], inputs["rho_hydr"])

        sig = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(sig.view[:, :, :], inputs["sig"])

        tempco = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(tempco.view[:, :, :], inputs["tempco"])

        tn_cup = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(tn_cup.view[:, :, :], inputs["tn_cup"])

        use_cloud_dissipation = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(
            use_cloud_dissipation.view[:, :, :], inputs["use_cloud_dissipation"]
        )

        vvel2d = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(vvel2d.view[:, :, :], inputs["vvel2d"])

        xmb = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(xmb.view[:, :, :], inputs["xmb"])

        zo = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(zo.view[:, :, :], inputs["zo"])

        zuo = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(zuo.view[:, :, :], inputs["zuo"])

        cloud_dissipation(
            # In
            COUPL_MPHYSICS=COUPL_MPHYSICS,
            dtime=dtime,
            heso_cup=heso_cup,
            ierr=ierr,
            kbcon=kbcon,
            ktop=ktop,
            qeso_cup=qeso_cup,
            qo_cup=qo_cup,
            rho_hydr=rho_hydr,
            sig=sig,
            tempco=tempco,
            tn_cup=tn_cup,
            use_cloud_dissipation=use_cloud_dissipation,
            vvel2d=vvel2d,
            zo=zo,
            zuo=zuo,
            xmb=xmb,
            # Out
            qrco=qrco,
            outqc=outqc,
            outq=outq,
            outt=outt,
        )

        return {
            "qrco": qrco.view[:],
            "outqc": outqc.view[:],
            "outq": outq.view[:],
            "outt": outt.view[:],
        }
