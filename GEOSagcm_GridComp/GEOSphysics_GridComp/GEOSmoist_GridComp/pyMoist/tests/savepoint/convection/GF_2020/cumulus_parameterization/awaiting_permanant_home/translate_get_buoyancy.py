from ndsl import Namelist, QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Float, Int
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array
from pyMoist.convection.GF_2020.GF_functions.get_buoyancy import GetBuoyancy


class TranslateGetBuoyancy(TranslateFortranData2Py):
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
            "hc": {},
            "he_cup": {},
            "hes_cup": {},
            "ierr": {},
            "kbcon": {},
            "klcl": {},
            "ktop": {},

        }

        # Float/Int Inputs
        self.in_vars["parameters"] = [
        ]

        # FloatField Outputs
        self.out_vars = {
            "dby": self.grid.compute_dict(),
  
        }

    def compute(self, inputs):

        get_buoyancy = GetBuoyancy(
            self.stencil_factory,
            self.grid.quantity_factory,
        )


        # Field inputs
        hc = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM,Z_DIM], units="n/a", 
        )
        safe_assign_array(hc.view[:, :, :], inputs["hc"])

        he_cup = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM,Z_DIM], units="n/a", 
        )
        safe_assign_array(he_cup.view[:, :, :], inputs["he_cup"])

        hes_cup = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM,Z_DIM], units="n/a",
        )
        safe_assign_array(hes_cup.view[:, :, :], inputs["hes_cup"])

        ierr = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM,Z_DIM], units="n/a",dtype=Int
        )
        safe_assign_array(ierr.view[:, :, :], inputs["ierr"])

        kbcon = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM,Z_DIM], units="n/a", dtype=Int
        )
        safe_assign_array(kbcon.view[:, :, :], inputs["kbcon"])

        klcl = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM,Z_DIM], units="n/a",dtype=Int
        )
        safe_assign_array(klcl.view[:, :, :], inputs["klcl"])

        ktop = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM,Z_DIM], units="n/a",dtype=Int
        )
        safe_assign_array(ktop.view[:, :, :], inputs["ktop"])

       
        dby = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM,Z_DIM], units="n/a"
        )
      

        get_buoyancy(
            # In
            hc=hc,
            he_cup=he_cup,
            hes_cup=hes_cup,
            ierr=ierr,
            kbcon=kbcon,
            klcl=klcl,
            ktop=ktop,
            # Out
            dby=dby,
        )

        return {
            "dby": dby.view[:], 
        }
