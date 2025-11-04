from f90nml import Namelist
from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Float, Int
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array
from pyMoist.convection.GF_2020.GF_functions.ke_to_heating import KeToHeating


class TranslateKeToHeating(TranslateFortranData2Py):
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
            "dellu": {},
            "dellv": {},
            "ierr": {},
            "ktop": {},
            "po_cup": {},
            "us": {},
            "vs": {},
            "dellat": {},
        }

        # Float/Int Inputs
        self.in_vars["parameters"] = []

        # FloatField Outputs
        self.out_vars = {
            "dellat": self.grid.compute_dict(),
        }

    def compute(self, inputs):

        ke_to_heating = KeToHeating(
            self.stencil_factory,
            self.grid.quantity_factory,
        )

        # # Field inputs
        dellu = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(dellu.view[:, :, :], inputs["dellu"])

        dellv = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(dellv.view[:, :, :], inputs["dellv"])

        po_cup = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(po_cup.view[:, :, :], inputs["po_cup"])

        us = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(us.view[:, :, :], inputs["us"])

        vs = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(vs.view[:, :, :], inputs["vs"])

        ierr = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int
        )
        safe_assign_array(ierr.view[:, :, :], inputs["ierr"])

        ktop = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int
        )
        safe_assign_array(ktop.view[:, :, :], inputs["ktop"])

        dellat = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(dellat.view[:, :, :], inputs["dellat"])

        ke_to_heating(
            # In
            dellu=dellu,
            dellv=dellv,
            ierr=ierr,
            ktop=ktop,
            po_cup=po_cup,
            us=us,
            vs=vs,
            # Out
            dellat=dellat,
        )

        return {
            "dellat": dellat.view[:],
        }
