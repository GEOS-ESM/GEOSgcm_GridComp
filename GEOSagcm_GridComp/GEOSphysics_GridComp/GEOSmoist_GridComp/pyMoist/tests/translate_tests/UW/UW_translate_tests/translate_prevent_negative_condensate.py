from f90nml import Namelist
from gt4py.cartesian.gtscript import int32

import pyMoist.constants as constants
from ndsl import StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float, Int
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array
from pyMoist.UW.compute_uwshcu import prevent_negative_condensate
from pyMoist.UW.config import UWConfiguration


class TranslatePreventNegativeCondensate(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "condensation": {},
            "dp0": {},
            "qi0": {},
            "qiten": {},
            "ql0": {},
            "qlten": {},
            "qtten": {},
            "qv0": {},
            "qvten": {},
            "s0": {},
            "slten": {},
            "sten": {},
        }

        # FloatField Outputs
        self.out_vars = {
            "qiten": self.grid.compute_dict(),
            "qlten": self.grid.compute_dict(),
            "qvten": self.grid.compute_dict(),
        }

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("ComputeUwshcuInv-constants")

    def compute(self, inputs):
        config = UWConfiguration(**self.constants)

        self.quantity_factory.add_data_dimensions(
            {
                "ntracers": constants.NCNST,
            }
        )

        self._prevent_negative_condensate = self.stencil_factory.from_dims_halo(
            func=prevent_negative_condensate,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "k0": config.k0,
                "dt": config.dt,
            },
        )

        # Inputs
        condensation = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a", dtype=bool)

        for i in range(0, 24):
            for j in range(0, 24):
                if inputs["condensation"][i, j] == 1:
                    condensation.view[i, j] = False
                else:
                    condensation.view[i, j] = True

        dp0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(dp0.view[:], inputs["dp0"])
        qi0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(qi0.view[:], inputs["qi0"])
        qiten = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(qiten.view[:], inputs["qiten"])
        ql0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(ql0.view[:], inputs["ql0"])
        qlten = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(qlten.view[:], inputs["qlten"])
        qtten = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(qtten.view[:], inputs["qtten"])
        qv0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(qv0.view[:], inputs["qv0"])
        qvten = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(qvten.view[:], inputs["qvten"])
        s0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(s0.view[:], inputs["s0"])
        slten = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(slten.view[:], inputs["slten"])
        sten = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(sten.view[:], inputs["sten"])

        # Outputs
        qmin = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        # The iteration you want to test
        iter_test = int32(0)

        # # Call stencils
        self._prevent_negative_condensate(
            condensation=condensation,
            qv0=qv0,
            qvten=qvten,
            ql0=ql0,
            qlten=qlten,
            qi0=qi0,
            s0=s0,
            sten=sten,
            dp0=dp0,
            qiten=qiten,
            qmin=qmin,
            iteration=iter_test,
        )

        return {
            "qvten": qvten.view[:],
            "qiten": qiten.view[:],
            "qlten": qlten.view[:],
        }
