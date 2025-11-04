from f90nml import Namelist
from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float, Int
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array
from pyMoist.convection.GF_2020.GF_functions.get_precip_fluxes import GetPrecipFluxes


class TranslateGetPrecipFluxes(TranslateFortranData2Py):
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
            "edto": {},
            "ierr": {},
            "ktop": {},
            "pwdo": {},
            "pwo": {},
            "xmb": {},
        }

        # Float/Int Inputs
        self.in_vars["parameters"] = []

        # FloatField Outputs
        self.out_vars = {
            "prec_flx": self.grid.compute_dict(),
            "evap_flx": self.grid.compute_dict(),
        }

    def compute(self, inputs):

        get_precip_fluxes = GetPrecipFluxes(
            self.stencil_factory,
            self.grid.quantity_factory,
        )

        # # Field inputs
        edto = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(edto.view[:, :, :], inputs["edto"])

        ierr = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=Int,
        )
        safe_assign_array(ierr.view[:, :, :], inputs["ierr"])

        ktop = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int
        )
        safe_assign_array(ktop.view[:, :, :], inputs["ktop"])

        pwdo = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(pwdo.view[:, :, :], inputs["pwdo"])

        pwo = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(pwo.view[:, :, :], inputs["pwo"])

        xmb = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(xmb.view[:, :, :], inputs["xmb"])

        prec_flx = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")

        evap_flx = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")

        get_precip_fluxes(
            # In
            edto=edto,
            ierr=ierr,
            ktop=ktop,
            pwdo=pwdo,
            pwo=pwo,
            xmb=xmb,
            # Out
            prec_flx=prec_flx,
            evap_flx=evap_flx,
        )

        return {
            "prec_flx": prec_flx.view[:],
            "evap_flx": evap_flx.view[:],
        }
