from f90nml import Namelist
from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float, Int
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array
from pyMoist.convection.GF_2020.GF_functions.cup_up_aa0 import CupUpAa0


class TranslateCupUpAa0(TranslateFortranData2Py):
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
            "dby": {},
            "gamma_cup": {},
            "k22": {},
            "kbcon": {},
            "klcl": {},
            "kbcon": {},
            "ktop": {},
            "t_cup": {},
            "z_cup": {},
            "zu": {},
            "ierr": {},
            "integ": {},
            "integ_interval": {},
        }

        # Float/Int Inputs
        self.in_vars["parameters"] = []

        # FloatField Outputs
        self.out_vars = {
            "ierr": self.grid.compute_dict(),
            "aa0": self.grid.compute_dict(),
        }

    def compute(self, inputs):

        cup_up_aa0 = CupUpAa0(
            self.stencil_factory,
            self.grid.quantity_factory,
        )

        # # Field inputs
        aa0 = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )

        dby = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(dby.view[:, :, :], inputs["dby"])

        gamma_cup = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(gamma_cup.view[:, :, :], inputs["gamma_cup"])

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

        k22 = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int
        )
        safe_assign_array(k22.view[:, :, :], inputs["k22"])

        t_cup = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(t_cup.view[:, :, :], inputs["t_cup"])

        z_cup = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(z_cup.view[:, :, :], inputs["z_cup"])

        zu = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(zu.view[:, :, :], inputs["zu"])

        integ = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=Int,
        )
        safe_assign_array(integ.view[:, :, :], inputs["integ"])

        integ_interval = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=Int,
        )
        safe_assign_array(integ_interval.view[:, :, :], inputs["integ_interval"])

        cup_up_aa0(
            # In
            dby=dby,
            gamma_cup=gamma_cup,
            ktop=ktop,
            kbcon=kbcon,
            k22=k22,
            z_cup=z_cup,
            t_cup=t_cup,
            zu=zu,
            integ=integ,
            integ_interval=integ_interval,
            # Out
            ierr=ierr,
            aa0=aa0,
        )

        return {
            "ierr": ierr.view[:],
            "aa0": aa0.view[:],
        }
