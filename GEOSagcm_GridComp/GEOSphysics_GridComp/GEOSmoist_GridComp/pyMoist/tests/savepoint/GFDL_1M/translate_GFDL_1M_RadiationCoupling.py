from f90nml import Namelist

from ndsl import StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.GFDL_1M.state import GFDL1MState
from pyMoist.GFDL_1M.locals import GFDL1MLocals
from pyMoist.GFDL_1M.config import GFDL1MConfig
from ndsl.stencils.testing.savepoint import DataLoader
from pyMoist.radiation_coupling import GFDL1MRadiationCoupling


class TranslateGFDL_1M_RadiationCoupling(TranslateFortranData2Py):
    def __init__(self, grid: Grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "mixing_ratio_vapor": {},
            "t": {},
            "mixing_ratio_large_scale_liquid": {},
            "mixing_ratio_large_scale_ice": {},
            "cloud_fraction_large_scale": {},
            "mixing_ratio_convective_liquid": {},
            "mixing_ratio_convective_ice": {},
            "cloud_fraction_convective": {},
            "local_p_mb": {},
            "mixing_ratio_rain": {},
            "mixing_ratio_snow": {},
            "mixing_ratio_graupel": {},
            "concentration_liquid": {},
            "concentration_ice": {},
            "radiation_vapor": {},
            "radiation_liquid": {},
            "radiation_ice": {},
            "radiation_rain": {},
            "radiation_snow": {},
            "radiation_graupel": {},
            "radiation_cloud_fraction": {},
            "cloud_particle_effective_radius_liquid": {},
            "cloud_particle_effective_radius_ice": {},
            "relative_humidity_after_pdf": {},
        }

        self.out_vars = self.in_vars["data_vars"].copy()

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GFDL_1M-constants")

    def compute(self, inputs):

        # initalize constants
        config = GFDL1MConfig(**self.constants)

        # initalize dataclasses
        state = GFDL1MState.zeros(self.quantity_factory)
        locals = GFDL1MLocals.zeros(self.quantity_factory)

        # Initalize saturation tables
        saturation_tables = SaturationVaporPressureTable(self.stencil_factory.backend)

        state.mixing_ratio.vapor.field[:] = inputs["mixing_ratio_vapor"]
        state.t.field[:] = inputs["t"]
        state.mixing_ratio.large_scale_liquid.field[:] = inputs["mixing_ratio_large_scale_liquid"]
        state.mixing_ratio.large_scale_ice.field[:] = inputs["mixing_ratio_large_scale_ice"]
        state.cloud_fraction.large_scale.field[:] = inputs["cloud_fraction_large_scale"]
        state.mixing_ratio.convective_liquid.field[:] = inputs["mixing_ratio_convective_liquid"]
        state.mixing_ratio.convective_ice.field[:] = inputs["mixing_ratio_convective_ice"]
        state.cloud_fraction.convective.field[:] = inputs["cloud_fraction_convective"]
        locals.p_mb.field[:] = inputs["local_p_mb"]
        state.mixing_ratio.rain.field[:] = inputs["mixing_ratio_rain"]
        state.mixing_ratio.snow.field[:] = inputs["mixing_ratio_snow"]
        state.mixing_ratio.graupel.field[:] = inputs["mixing_ratio_graupel"]
        state.concentration.liquid.field[:] = inputs["concentration_liquid"]
        state.concentration.ice.field[:] = inputs["concentration_ice"]
        state.radiation_field.vapor.field[:] = inputs["radiation_vapor"]
        state.radiation_field.liquid.field[:] = inputs["radiation_liquid"]
        state.radiation_field.ice.field[:] = inputs["radiation_ice"]
        state.radiation_field.rain.field[:] = inputs["radiation_rain"]
        state.radiation_field.snow.field[:] = inputs["radiation_snow"]
        state.radiation_field.graupel.field[:] = inputs["radiation_graupel"]
        state.radiation_field.cloud_fraction.field[:] = inputs["radiation_cloud_fraction"]
        state.cloud_particle_effective_radius.liquid.field[:] = inputs[
            "cloud_particle_effective_radius_liquid"
        ]
        state.cloud_particle_effective_radius.ice.field[:] = inputs["cloud_particle_effective_radius_ice"]
        state.relative_humidity_after_pdf.field[:] = inputs["relative_humidity_after_pdf"]

        # construct test stencil
        code = GFDL1MRadiationCoupling(
            stencil_factory=self.stencil_factory,
            config=config,
            saturation_tables=saturation_tables,
        )
        code(
            state=state,
            locals=locals,
        )

        return {
            "mixing_ratio_vapor": state.mixing_ratio.vapor.field[:],
            "t": state.t.field[:],
            "mixing_ratio_large_scale_liquid": state.mixing_ratio.large_scale_liquid.field[:],
            "mixing_ratio_large_scale_ice": state.mixing_ratio.large_scale_ice.field[:],
            "cloud_fraction_large_scale": state.cloud_fraction.large_scale.field[:],
            "mixing_ratio_convective_liquid": state.mixing_ratio.convective_liquid.field[:],
            "mixing_ratio_convective_ice": state.mixing_ratio.convective_ice.field[:],
            "cloud_fraction_convective": state.cloud_fraction.convective.field[:],
            "local_p_mb": locals.p_mb.field[:],
            "mixing_ratio_rain": state.mixing_ratio.rain.field[:],
            "mixing_ratio_snow": state.mixing_ratio.snow.field[:],
            "mixing_ratio_graupel": state.mixing_ratio.graupel.field[:],
            "concentration_liquid": state.concentration.liquid.field[:],
            "concentration_ice": state.concentration.ice.field[:],
            "radiation_vapor": state.radiation_field.vapor.field[:],
            "radiation_liquid": state.radiation_field.liquid.field[:],
            "radiation_ice": state.radiation_field.ice.field[:],
            "radiation_rain": state.radiation_field.rain.field[:],
            "radiation_snow": state.radiation_field.snow.field[:],
            "radiation_graupel": state.radiation_field.graupel.field[:],
            "radiation_cloud_fraction": state.radiation_field.cloud_fraction.field[:],
            "cloud_particle_effective_radius_liquid": state.cloud_particle_effective_radius.liquid.field[:],
            "cloud_particle_effective_radius_ice": state.cloud_particle_effective_radius.ice.field[:],
            "relative_humidity_after_pdf": state.relative_humidity_after_pdf.field[:],
        }
