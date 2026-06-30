from f90nml import Namelist
from gt4py.cartesian.gtscript import int32
from ndsl import StencilFactory
from ndsl.constants import I_DIM, J_DIM, K_DIM, K_INTERFACE_DIM
from ndsl.dsl.typing import Int
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array

import pyMoist.constants as constants
from pyMoist.convection.UW.compute_uwshcu import calc_ppen
from pyMoist.convection.UW.config import UWConfiguration


class TranslateCalcPpen(TranslateFortranData2Py):
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
            "bogbot": {},
            "bogtop": {},
            "condensation": {},
            "dp0": {},
            "drage": {},
            "kpen": {},
            "pifc0": {},
            "rhomid0j": {},
            "wtwb": {},
            "wu": {},
        }

        # FloatField Outputs
        self.out_vars = {
            "ppen": self.grid.compute_dict(),
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

        self._calc_ppen = self.stencil_factory.from_dims_halo(
            func=calc_ppen,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )

        # Field inputs
        condensation = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a", dtype=bool)

        for i in range(0, 24):
            for j in range(0, 24):
                if inputs["condensation"][i, j] == 1:
                    condensation.view[i, j] = False
                else:
                    condensation.view[i, j] = True

        kpen_IJ = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a", dtype=Int)
        safe_assign_array(kpen_IJ.view[:], inputs["kpen"] - 1)
        bogbot = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(bogbot.view[:], inputs["bogbot"])
        bogtop = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(bogtop.view[:], inputs["bogtop"])
        dp0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(dp0.view[:], inputs["dp0"])
        drage = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(drage.view[:], inputs["drage"])
        pifc0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(pifc0.view[:], inputs["pifc0"])
        rhomid0j = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(rhomid0j.view[:], inputs["rhomid0j"])
        wtwb = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(wtwb.view[:], inputs["wtwb"])
        wu = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(wu.view[:], inputs["wu"])

        # Outputs
        ppen = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        kpen = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a", dtype=Int)

        # The iteration you want to test
        iter_test = int32(0)

        # # Call stencils
        self._calc_ppen(
            condensation=condensation,
            drage=drage,
            bogbot=bogbot,
            bogtop=bogtop,
            pifc0=pifc0,
            kpen_IJ=kpen_IJ,
            kpen=kpen,
            wu=wu,
            rhomid0j=rhomid0j,
            dp0=dp0,
            wtwb=wtwb,
            ppen=ppen,
            iteration=iter_test,
        )

        return {
            "ppen": ppen.view[:],
        }
