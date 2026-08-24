from f90nml import Namelist
from gt4py.cartesian.gtscript import int32
from ndsl import StencilFactory
from ndsl.constants import I_DIM, J_DIM, K_DIM, K_INTERFACE_DIM
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array

import pyMoist.constants as constants
from pyMoist.convection.UW.compute_uwshcu import calc_tracer_tendencies
from pyMoist.convection.UW.config import UWConfiguration


class TranslateTracerTendencies(TranslateFortranData2Py):
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
            "tr0_TracerTendencies": {},
            "trflx": {},
        }

        # FloatField Outputs
        self.out_vars = {
            "trten": self.grid.compute_dict(),
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

        self._calc_tracer_tendencies = self.stencil_factory.from_dims_halo(
            func=calc_tracer_tendencies,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "ncnst": config.NCNST,
                "dt": config.dt,
                "dotransport": config.dotransport,
            },
        )

        # Inputs
        condensation = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a", dtype=bool)

        for i in range(0, 24):
            for j in range(0, 24):
                if inputs["condensation"][i, j] == 1:
                    condensation.view[i, j] = False
                else:
                    condensation.view[i, j] = True

        dp0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(dp0.view[:], inputs["dp0"])
        tr0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM, "ntracers"], units="n/a")
        safe_assign_array(tr0.view[:], inputs["tr0_TracerTendencies"])
        trflx = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM, "ntracers"], units="n/a")
        safe_assign_array(trflx.view[:], inputs["trflx"])

        # Outputs
        trten = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM, "ntracers"], units="n/a")
        trflx_d = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM, "ntracers"], units="n/a")
        trflx_u = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM, "ntracers"], units="n/a")
        trmin = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, "ntracers"], units="n/a")

        # The iteration you want to test
        iter_test = int32(0)

        # # Call stencils
        self._calc_tracer_tendencies(
            condensation=condensation,
            dp0=dp0,
            trflx_d=trflx_d,
            trflx_u=trflx_u,
            trmin=trmin,
            tr0=tr0,
            trflx=trflx,
            trten=trten,
            iteration=iter_test,
        )

        return {
            "trten": trten.view[:],
        }
