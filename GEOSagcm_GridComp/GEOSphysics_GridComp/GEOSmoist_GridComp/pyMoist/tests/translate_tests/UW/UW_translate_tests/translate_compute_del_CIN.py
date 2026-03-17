from f90nml import Namelist
from gt4py.cartesian.gtscript import int32
from ndsl import StencilFactory
from ndsl.constants import I_DIM, J_DIM, K_DIM
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array

import pyMoist.constants as constants
from pyMoist.UW.compute_uwshcu import compute_del_CIN
from pyMoist.UW.config import UWConfiguration


class TranslateComputeDelCIN(TranslateFortranData2Py):
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
            "cin_i": {},
            "cinlcl_i": {},
            "cin": {},
            "cinlcl": {},
        }

        # FloatField Outputs
        self.out_vars = {
            "del_CIN": self.grid.compute_dict(),
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

        self._compute_del_CIN = self.stencil_factory.from_dims_halo(
            func=compute_del_CIN,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={"use_CINcin": config.use_CINcin},
        )

        # Inputs
        condensation = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a", dtype=bool)

        for i in range(0, 24):
            for j in range(0, 24):
                if inputs["condensation"][i, j] == 1:
                    condensation.view[i, j] = False
                else:
                    condensation.view[i, j] = True

        cin = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(cin.view[:], inputs["cin"])
        cinlcl = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(cinlcl.view[:], inputs["cinlcl"])
        cin_i = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(cin_i.view[:], inputs["cin_i"])
        cinlcl_i = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(cinlcl_i.view[:], inputs["cinlcl_i"])

        # Outputs
        del_CIN = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")

        # The iteration you want to test
        iter_test = int32(1)

        # # Call stencils
        self._compute_del_CIN(
            condensation=condensation,
            cin_IJ=cin,
            cinlcl_IJ=cinlcl,
            cin_i=cin_i,
            cinlcl_i=cinlcl_i,
            del_CIN=del_CIN,
        )

        return {
            "del_CIN": del_CIN.view[:],
        }
