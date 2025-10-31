from f90nml import Namelist
from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float, Int
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array
from pyMoist.convection.GF_2020.GF_functions.rain_evap_below_cloudbase import (
    RainEvapBelowCloudbase,
)


class TranslateRainEvapBelowCloudbase(TranslateFortranData2Py):
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
            "cumulus": {},
            "edto": {},
            "evap_flx": {},
            "ierr": {},
            "kbcon": {},
            "ktop": {},
            "outbuoy": {},
            "outq": {},
            "outt": {},
            "po_cup": {},
            "pre": {},
            "prec_flx": {},
            "psur": {},
            "pwavo": {},
            "pwdo": {},
            "pwevo": {},
            "pwo": {},
            "qes_cup": {},
            "qo_cup": {},
            "t_cup": {},
            "xland": {},
            "xmb": {},
        }

        # Float/Int Inputs
        self.in_vars["parameters"] = []

        # FloatField Outputs
        self.out_vars = {
            "evap_bcb": self.grid.compute_dict(),
            "evap_flx": self.grid.compute_dict(),
            "outbuoy": self.grid.compute_dict(),
            "outq": self.grid.compute_dict(),
            "outt": self.grid.compute_dict(),
            "pre": self.grid.compute_dict(),
            "prec_flx": self.grid.compute_dict(),
        }

    def compute(self, inputs):

        rain_evap_below_cloudbase = RainEvapBelowCloudbase(
            self.stencil_factory,
            self.grid.quantity_factory,
        )

        # # Field inputs

        cumulus = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=Int,
        )
        safe_assign_array(cumulus.view[:, :, :], inputs["cumulus"])

        edto = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(edto.view[:, :, :], inputs["edto"])

        evap_flx = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(evap_flx.view[:, :, :], inputs["evap_flx"])

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

        outbuoy = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(outbuoy.view[:, :, :], inputs["outbuoy"])

        outq = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(outq.view[:, :, :], inputs["outq"])

        outt = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(outt.view[:, :, :], inputs["outt"])

        po_cup = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(po_cup.view[:, :, :], inputs["po_cup"])

        pre = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(pre.view[:, :, :], inputs["pre"])

        prec_flx = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(prec_flx.view[:, :, :], inputs["prec_flx"])

        psur = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(psur.view[:, :, :], inputs["psur"])

        pwavo = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(pwavo.view[:, :, :], inputs["pwavo"])

        pwdo = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(pwdo.view[:, :, :], inputs["pwdo"])

        pwevo = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(pwevo.view[:, :, :], inputs["pwevo"])

        pwo = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(pwo.view[:, :, :], inputs["pwo"])

        qes_cup = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(qes_cup.view[:, :, :], inputs["qes_cup"])

        qo_cup = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(qo_cup.view[:, :, :], inputs["qo_cup"])

        t_cup = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(t_cup.view[:, :, :], inputs["t_cup"])

        xland = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(xland.view[:, :, :], inputs["xland"])

        xmb = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(xmb.view[:, :, :], inputs["xmb"])

        evap_bcb = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )

        rain_evap_below_cloudbase(
            # In
            cumulus=cumulus,
            edto=edto,
            ierr=ierr,
            kbcon=kbcon,
            ktop=ktop,
            po_cup=po_cup,
            psur=psur,
            pwavo=pwavo,
            pwdo=pwdo,
            pwevo=pwevo,
            pwo=pwo,
            qes_cup=qes_cup,
            qo_cup=qo_cup,
            t_cup=t_cup,
            xland=xland,
            xmb=xmb,
            # Out
            evap_bcb=evap_bcb,
            evap_flx=evap_flx,
            outbuoy=outbuoy,
            outq=outq,
            outt=outt,
            pre=pre,
            prec_flx=prec_flx,
        )

        return {
            "evap_bcb": evap_bcb.view[:],
            "evap_flx": evap_flx.view[:],
            "outbuoy": outbuoy.view[:],
            "outq": outq.view[:],
            "outt": outt.view[:],
            "pre": pre.view[:],
            "prec_flx": prec_flx.view[:],
        }
