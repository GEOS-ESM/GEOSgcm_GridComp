from f90nml import Namelist

from ndsl import StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.GFDL_1M.locals import GFDL1MLocals
from pyMoist.GFDL_1M.PhaseChange.rh_calculations import rh_calculations
from pyMoist.GFDL_1M.state import GFDL1MState


class TranslateGFDL_1M_RHCalculations(TranslateFortranData2Py):
    def __init__(self, grid: Grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "estimated_inversion_strength": {},
            "local_lcl_level": {},
            "area": {},
            "top_of_local_p_interface_mb": {},
            "local_p_mb": {},
            "local_alpha": {},
        }

        self.out_vars = self.in_vars["data_vars"].copy()

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GFDL_1M-constants")

    def compute(self, inputs):
        # initalize constants
        config = GFDL1MConfig(**self.constants)

        # initalize dataclasses
        state = GFDL1MState.zeros(self.quantity_factory)
        locals_ = GFDL1MLocals.make_as_state(self.quantity_factory)

        # Internal from wrapper class needed for this test
        alpha = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        # fill relavent parts of dataclasses
        state.estimated_inversion_strength.field[:] = inputs["estimated_inversion_strength"]
        locals_.lcl_level.field[:] = inputs["local_lcl_level"]
        state.area.field[:] = inputs["area"]
        locals_.p_interface_mb.field[:, :, -1] = inputs["top_of_local_p_interface_mb"]
        locals_.p_mb.field[:] = inputs["local_p_mb"]
        alpha.field[:] = inputs["local_alpha"]

        # construct test stencil
        code = self.stencil_factory.from_dims_halo(
            func=rh_calculations,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "DW_LAND": config.DW_LAND,
                "DW_OCEAN": config.DW_OCEAN,
                "TURNRHCRIT_PARAM": config.TURNRHCRIT_PARAM,
            },
        )
        code(
            estimated_inversion_strength=state.estimated_inversion_strength,
            p_mb=locals_.p_mb,
            p_interface_mb=locals_.p_interface_mb,
            area=state.area,
            lcl_level=locals_.lcl_level,
            alpha=alpha,
        )

        return {
            "estimated_inversion_strength": state.estimated_inversion_strength.field[:],
            "local_lcl_level": locals_.lcl_level.field[:],
            "area": state.area.field[:],
            "top_of_local_p_interface_mb": locals_.p_interface_mb.field[:, :, -1],
            "local_p_mb": locals_.p_mb.field[:],
            "local_alpha": alpha.field[:],
        }
