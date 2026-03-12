from f90nml import Namelist
from gt4py.cartesian.gtscript import int32

import pyMoist.constants as constants
from ndsl import StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Float, Int
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array
from pyMoist.saturation_tables import get_saturation_vapor_pressure_table
from pyMoist.UW.compute_uwshcu import compute_diagnostic_outputs
from pyMoist.UW.config import UWConfiguration


class TranslateComputeDiagnosticOutputs(TranslateFortranData2Py):
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
            "prel": {},
            "qtu": {},
            "thlu": {},
            "krel": {},
        }

        # FloatField Outputs
        self.out_vars = {
            "qcubelow": self.grid.compute_dict(),
            "qiubelow": self.grid.compute_dict(),
            "qlubelow": self.grid.compute_dict(),
            "rcwp": self.grid.compute_dict(),
            "riwp": self.grid.compute_dict(),
            "rlwp": self.grid.compute_dict(),
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

        self._compute_diagnostic_outputs = self.stencil_factory.from_dims_halo(
            func=compute_diagnostic_outputs,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        # Inputs
        condensation = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a", dtype=bool)

        for i in range(0, 24):
            for j in range(0, 24):
                if inputs["condensation"][i, j] == 1:
                    condensation.view[i, j] = False
                else:
                    condensation.view[i, j] = True

        prel = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(prel.view[:], inputs["prel"])
        qtu = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(qtu.view[:], inputs["qtu"])
        thlu = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(thlu.view[:], inputs["thlu"])
        krel = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int)
        safe_assign_array(krel.view[:], inputs["krel"] - 1)

        # Outputs
        qcubelow = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        qiubelow = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        qlubelow = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        rcwp = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        riwp = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        rlwp = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")

        saturation_vapor_pressure_table = get_saturation_vapor_pressure_table(self.stencil_factory.backend)
        self.ese = saturation_vapor_pressure_table.ese
        self.esx = saturation_vapor_pressure_table.esx

        umf_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        qtflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        slflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        uflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        vflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        cush_inout = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")

        # The iteration you want to test
        iter_test = int32(0)

        # # Call stencils
        self._compute_diagnostic_outputs(
            condensation=condensation,
            prel=prel,
            thlu=thlu,
            qtu=qtu,
            krel=krel,
            ese=self.ese,
            esx=self.esx,
            umf_out=umf_out,
            qtflx_out=qtflx_out,
            slflx_out=slflx_out,
            uflx_out=uflx_out,
            vflx_out=vflx_out,
            qcubelow=qcubelow,
            qlubelow=qlubelow,
            qiubelow=qiubelow,
            rcwp=rcwp,
            rlwp=rlwp,
            riwp=riwp,
            cush_inout=cush_inout,
        )

        return {
            "qcubelow": qcubelow.view[:],
            "qiubelow": qiubelow.view[:],
            "qlubelow": qlubelow.view[:],
            "rcwp": rcwp.view[:],
            "riwp": riwp.view[:],
            "rlwp": rlwp.view[:],
        }
