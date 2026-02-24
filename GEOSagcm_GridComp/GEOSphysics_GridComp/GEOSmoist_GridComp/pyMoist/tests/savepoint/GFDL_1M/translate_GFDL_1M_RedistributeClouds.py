from f90nml import Namelist

from ndsl import StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.GFDL_1M.state import GFDL1MState
from pyMoist.redistribute_clouds import redistribute_clouds


class TranslateGFDL_1M_RedistributeClouds(TranslateFortranData2Py):
    def __init__(self, grid: Grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "radiation_cloud_fraction": {},
            "radiation_liquid": {},
            "radiation_ice": {},
            "cloud_fraction_convective": {},
            "cloud_fraction_large_scale": {},
            "mixing_ratio_convective_liquid": {},
            "mixing_ratio_large_scale_liquid": {},
            "mixing_ratio_convective_ice": {},
            "mixing_ratio_large_scale_ice": {},
            "radiation_vapor": {},
            "t": {},
        }

        self.out_vars = self.in_vars["data_vars"].copy()

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GFDL_1M-constants")

    def compute(self, inputs):
        # initialize dataclasses
        state = GFDL1MState.zeros(self.quantity_factory)

        state.radiation_field.cloud_fraction.field[:] = inputs["radiation_cloud_fraction"]
        state.radiation_field.liquid.field[:] = inputs["radiation_liquid"]
        state.radiation_field.ice.field[:] = inputs["radiation_ice"]
        state.cloud_fraction.convective.field[:] = inputs["cloud_fraction_convective"]
        state.cloud_fraction.large_scale.field[:] = inputs["cloud_fraction_large_scale"]
        state.mixing_ratio.convective_liquid.field[:] = inputs["mixing_ratio_convective_liquid"]
        state.mixing_ratio.large_scale_liquid.field[:] = inputs["mixing_ratio_large_scale_liquid"]
        state.mixing_ratio.convective_ice.field[:] = inputs["mixing_ratio_convective_ice"]
        state.mixing_ratio.large_scale_ice.field[:] = inputs["mixing_ratio_large_scale_ice"]
        state.radiation_field.vapor.field[:] = inputs["radiation_vapor"]
        state.t.field[:] = inputs["t"]

        # construct test stencil
        code = self.stencil_factory.from_dims_halo(
            func=redistribute_clouds,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )
        code(
            cloud_fraction=state.radiation_field.cloud_fraction,
            convective_cloud_fraction=state.cloud_fraction.convective,
            large_scale_cloud_fraction=state.cloud_fraction.large_scale,
            liquid=state.radiation_field.liquid,
            convective_liquid=state.mixing_ratio.convective_liquid,
            large_scale_liquid=state.mixing_ratio.large_scale_liquid,
            ice=state.radiation_field.ice,
            convective_ice=state.mixing_ratio.convective_ice,
            large_scale_ice=state.mixing_ratio.large_scale_ice,
            vapor=state.radiation_field.vapor,
            temperature=state.t,
        )

        return {
            "radiation_cloud_fraction": state.radiation_field.cloud_fraction.field[:],
            "radiation_liquid": state.radiation_field.liquid.field[:],
            "radiation_ice": state.radiation_field.ice.field[:],
            "cloud_fraction_convective": state.cloud_fraction.convective.field[:],
            "cloud_fraction_large_scale": state.cloud_fraction.large_scale.field[:],
            "mixing_ratio_convective_liquid": state.mixing_ratio.convective_liquid.field[:],
            "mixing_ratio_large_scale_liquid": state.mixing_ratio.large_scale_liquid.field[:],
            "mixing_ratio_convective_ice": state.mixing_ratio.convective_ice.field[:],
            "mixing_ratio_large_scale_ice": state.mixing_ratio.large_scale_ice.field[:],
            "radiation_vapor": state.radiation_field.vapor.field[:],
            "t": state.t.field[:],
        }
