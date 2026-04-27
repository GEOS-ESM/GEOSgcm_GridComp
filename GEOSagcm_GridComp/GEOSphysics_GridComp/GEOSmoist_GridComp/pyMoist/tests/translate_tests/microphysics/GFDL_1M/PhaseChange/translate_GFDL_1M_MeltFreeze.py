from f90nml import Namelist
from ndsl import StencilFactory
from ndsl.constants import I_DIM, J_DIM, K_DIM
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py

from pyMoist.microphysics.GFDL_1M.config import GFDL1MConfig
from pyMoist.microphysics.GFDL_1M.PhaseChange.melt_freeze import melt_freeze
from pyMoist.microphysics.GFDL_1M.state import GFDL1MState
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable


class TranslateGFDL_1M_MeltFreeze(TranslateFortranData2Py):
    def __init__(self, grid: Grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.quantity_factory = grid.quantity_factory

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "t": {},
            "mixing_ratio_convective_liquid": {},
            "mixing_ratio_convective_ice": {},
            "mixing_ratio_large_scale_liquid": {},
            "mixing_ratio_large_scale_ice": {},
            "convection_fraction": {},
            "surface_type": {},
        }

        self.out_vars = self.in_vars["data_vars"].copy()

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GFDL_1M-constants")

    def compute(self, inputs):
        # initialize constants
        config = GFDL1MConfig(**self.constants)

        # initialize dataclasses
        state = GFDL1MState.zeros(self.quantity_factory)

        # Initialize saturation tables
        self.saturation_tables = SaturationVaporPressureTable(self.stencil_factory.backend)

        state.t.field[:] = inputs["t"]
        state.mixing_ratio.convective_liquid.field[:] = inputs["mixing_ratio_convective_liquid"]
        state.mixing_ratio.convective_ice.field[:] = inputs["mixing_ratio_convective_ice"]
        state.mixing_ratio.large_scale_liquid.field[:] = inputs["mixing_ratio_large_scale_liquid"]
        state.mixing_ratio.large_scale_ice.field[:] = inputs["mixing_ratio_large_scale_ice"]
        state.convection_fraction.field[:] = inputs["convection_fraction"]
        state.surface_type.field[:] = inputs["surface_type"]

        # construct test stencil
        code = self.stencil_factory.from_dims_halo(
            func=melt_freeze,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={"DT_MOIST": config.DT_MOIST},
        )
        code(
            convection_fraction=state.convection_fraction,
            surface_type=state.surface_type,
            t=state.t,
            liquid=state.mixing_ratio.convective_liquid,
            ice=state.mixing_ratio.convective_ice,
        )
        code(
            convection_fraction=state.convection_fraction,
            surface_type=state.surface_type,
            t=state.t,
            liquid=state.mixing_ratio.large_scale_liquid,
            ice=state.mixing_ratio.large_scale_ice,
        )

        return {
            "t": state.t.field[:],
            "mixing_ratio_convective_liquid": state.mixing_ratio.convective_liquid.field[:],
            "mixing_ratio_convective_ice": state.mixing_ratio.convective_ice.field[:],
            "mixing_ratio_large_scale_liquid": state.mixing_ratio.large_scale_liquid.field[:],
            "mixing_ratio_large_scale_ice": state.mixing_ratio.large_scale_ice.field[:],
            "convection_fraction": state.convection_fraction.field[:],
            "surface_type": state.surface_type.field[:],
        }
