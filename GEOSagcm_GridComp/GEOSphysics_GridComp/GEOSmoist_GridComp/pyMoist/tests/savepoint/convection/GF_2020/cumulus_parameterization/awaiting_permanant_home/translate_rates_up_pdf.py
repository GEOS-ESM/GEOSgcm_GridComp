from f90nml import Namelist
from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float, Int
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array
from pyMoist.convection.GF_2020.GF_functions.rates_up_pdf import RatesUpPdf


class TranslateRatesUpPdf(TranslateFortranData2Py):
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
            "OVERSHOOT": {},
            "entr_rate": {},
            "heo": {},
            "heso_cup": {},
            "hkbo": {},
            "ierr": {},
            "kbcon": {},
            "ktop": {},
            "name": {},
            "z_cup": {},
        }

        # Float/Int Inputs
        self.in_vars["parameters"] = []

        # FloatField Outputs
        self.out_vars = {
            "ierr": self.grid.compute_dict(),
            "ktop": self.grid.compute_dict(),
        }

    def compute(self, inputs):

        rates_up_pdf = RatesUpPdf(
            self.stencil_factory,
            self.grid.quantity_factory,
        )

        # # Field inputs
        OVERSHOOT = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(OVERSHOOT.view[:, :, :], inputs["OVERSHOOT"])

        entr_rate = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(entr_rate.view[:, :, :], inputs["entr_rate"])

        heo = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(heo.view[:, :, :], inputs["heo"])

        heso_cup = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(heso_cup.view[:, :, :], inputs["heso_cup"])

        hkbo = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(hkbo.view[:, :, :], inputs["hkbo"])

        ierr = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int
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

        name = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int
        )
        safe_assign_array(name.view[:, :, :], inputs["name"])

        z_cup = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(z_cup.view[:, :, :], inputs["z_cup"])

        ktopIJ = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM], units="n/a", dtype=Int)

        rates_up_pdf(
            # In
            OVERSHOOT=OVERSHOOT,
            entr_rate=entr_rate,
            heo=heo,
            heso_cup=heso_cup,
            hkbo=hkbo,
            kbcon=kbcon,
            name=name,
            z_cup=z_cup,
            # Out
            ierr=ierr,
            ktopIJ=ktopIJ,
            ktop=ktop,
        )

        return {
            "ierr": ierr.view[:],
            "ktop": ktop.view[:] + 1,  # Add 1 to K-level array for testing
        }
