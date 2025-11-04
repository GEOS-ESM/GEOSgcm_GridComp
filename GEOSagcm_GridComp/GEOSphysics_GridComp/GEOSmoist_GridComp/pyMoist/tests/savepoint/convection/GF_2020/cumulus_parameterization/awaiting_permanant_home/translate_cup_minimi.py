from f90nml import Namelist
from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Float, Int
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array
from pyMoist.convection.GF_2020.GF_functions.cup_minimi import CupMinimi


class TranslateCupMinimi(TranslateFortranData2Py):
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
            "HEso_cup": {},
            "Kbcon": {},
            "ierr": {},
            "kstabm": {},
        }

        # Float/Int Inputs
        self.in_vars["parameters"] = []

        # FloatField Outputs
        self.out_vars = {
            "kstabi": self.grid.compute_dict(),
        }

    def compute(self, inputs):

        cup_minimi = CupMinimi(
            self.stencil_factory,
            self.grid.quantity_factory,
        )

        # Field inputs
        HEso_cup = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(HEso_cup.view[:, :, :], inputs["HEso_cup"])

        Kbcon = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int
        )
        safe_assign_array(Kbcon.view[:, :, :], inputs["Kbcon"])

        ierr = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int
        )
        safe_assign_array(ierr.view[:, :, :], inputs["ierr"])

        kstabm = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int
        )
        safe_assign_array(kstabm.view[:, :, :], inputs["kstabm"])

        kstabi = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int
        )

        cup_minimi(
            # In
            array=HEso_cup,
            ks=Kbcon,
            kend=kstabm,
            ierr=ierr,
            # Out
            kt=kstabi,
        )

        return {
            "kstabi": kstabi.view[:] + 1,  # kstabi is an array of klevels so need to adjust by 1
        }
