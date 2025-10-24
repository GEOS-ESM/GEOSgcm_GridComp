from ndsl import Namelist, QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float, Int
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array
from pyMoist.convection.GF_2020.GF_functions.cup_env import CupEnv


class TranslateCupEnv(TranslateFortranData2Py):
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
            "SATUR_CALC": {},
            "ierr": {},
            "itest": {},
            "p": {},
            "psur": {},
            "q": {},
            "t": {},
            "z": {},
            "z1": {},
        }

        # Float/Int Inputs
        self.in_vars["parameters"] = []

        # FloatField Outputs
        self.out_vars = {
            "he": self.grid.compute_dict(),
            "hes": self.grid.compute_dict(),
            "ierr": self.grid.compute_dict(),
            "qes": self.grid.compute_dict(),
            "z": self.grid.compute_dict(),
        }

    def compute(self, inputs):

        cup_env = CupEnv(
            self.stencil_factory,
            self.grid.quantity_factory,
        )

        # # Field inputs
        SATUR_CALC = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int
        )
        safe_assign_array(SATUR_CALC.view[:, :, :], inputs["SATUR_CALC"])

        ierr = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=Int,
        )
        safe_assign_array(ierr.view[:, :, :], inputs["ierr"])

        itest = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int
        )
        safe_assign_array(itest.view[:, :, :], inputs["itest"])

        p = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a"
        )
        safe_assign_array(p.view[:, :, :], inputs["p"])

        psur = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(psur.view[:, :, :], inputs["psur"])

        q = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(q.view[:, :, :], inputs["q"])

        t = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(t.view[:, :, :], inputs["t"])

        z = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(z.view[:, :, :], inputs["z"])

        z1 = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(z1.view[:, :, :], inputs["z1"])

        he = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )

        hes = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )

        qes = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )

        cup_env(
            # In
            SATUR_CALC=SATUR_CALC,
            itest=itest,
            p=p,
            psur=psur,
            q=q,
            t=t,
            z1=z1,
            # Out
            he=he,
            hes=hes,
            ierr=ierr,
            qes=qes,
            z=z,
        )

        return {
            "he": he.view[:],
            "hes": hes.view[:],
            "ierr": ierr.view[:],
            "qes": qes.view[:],
            "z": z.view[:],
        }
