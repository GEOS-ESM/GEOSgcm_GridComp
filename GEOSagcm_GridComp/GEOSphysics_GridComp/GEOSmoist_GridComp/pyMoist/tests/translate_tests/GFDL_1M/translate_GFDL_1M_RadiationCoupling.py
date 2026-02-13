from f90nml import Namelist
from ndsl import StencilFactory
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array

from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.GFDL_1M.locals import GFDL1MLocals
from pyMoist.GFDL_1M.radiation_coupling import GFDL1MRadiationCoupling
from pyMoist.GFDL_1M.state import GFDL1MState
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable


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
        # initialize constants
        config = GFDL1MConfig(**self.constants)

        # initialize dataclasses
        state = GFDL1MState.zeros(self.quantity_factory)
        locals_ = GFDL1MLocals.make_as_state(self.quantity_factory)

        # Initialize saturation tables
        saturation_tables = SaturationVaporPressureTable(self.stencil_factory.backend)

        safe_assign_array(state.mixing_ratio.vapor.field[:], inputs["mixing_ratio_vapor"])
        safe_assign_array(state.t.field[:], inputs["t"])
        safe_assign_array(
            state.mixing_ratio.large_scale_liquid.field[:], inputs["mixing_ratio_large_scale_liquid"]
        )
        safe_assign_array(state.mixing_ratio.large_scale_ice.field[:], inputs["mixing_ratio_large_scale_ice"])
        safe_assign_array(state.cloud_fraction.large_scale.field[:], inputs["cloud_fraction_large_scale"])
        safe_assign_array(
            state.mixing_ratio.convective_liquid.field[:], inputs["mixing_ratio_convective_liquid"]
        )
        safe_assign_array(state.mixing_ratio.convective_ice.field[:], inputs["mixing_ratio_convective_ice"])
        safe_assign_array(state.cloud_fraction.convective.field[:], inputs["cloud_fraction_convective"])
        safe_assign_array(locals_.p_mb.field[:], inputs["local_p_mb"])
        safe_assign_array(state.mixing_ratio.rain.field[:], inputs["mixing_ratio_rain"])
        safe_assign_array(state.mixing_ratio.snow.field[:], inputs["mixing_ratio_snow"])
        safe_assign_array(state.mixing_ratio.graupel.field[:], inputs["mixing_ratio_graupel"])
        safe_assign_array(state.concentration.liquid.field[:], inputs["concentration_liquid"])
        safe_assign_array(state.concentration.ice.field[:], inputs["concentration_ice"])
        safe_assign_array(state.radiation_field.vapor.field[:], inputs["radiation_vapor"])
        safe_assign_array(state.radiation_field.liquid.field[:], inputs["radiation_liquid"])
        safe_assign_array(state.radiation_field.ice.field[:], inputs["radiation_ice"])
        safe_assign_array(state.radiation_field.rain.field[:], inputs["radiation_rain"])
        safe_assign_array(state.radiation_field.snow.field[:], inputs["radiation_snow"])
        safe_assign_array(state.radiation_field.graupel.field[:], inputs["radiation_graupel"])
        safe_assign_array(state.radiation_field.cloud_fraction.field[:], inputs["radiation_cloud_fraction"])
        safe_assign_array(
            state.cloud_particle_effective_radius.liquid.field[:],
            inputs["cloud_particle_effective_radius_liquid"],
        )
        safe_assign_array(
            state.cloud_particle_effective_radius.ice.field[:], inputs["cloud_particle_effective_radius_ice"]
        )
        safe_assign_array(state.relative_humidity_after_pdf.field[:], inputs["relative_humidity_after_pdf"])

        # construct test stencil
        code = GFDL1MRadiationCoupling(
            stencil_factory=self.stencil_factory,
            config=config,
            saturation_tables=saturation_tables,
        )
        code(
            t=state.t,
            mixing_ratio_vapor=state.mixing_ratio.vapor,
            mixing_ratio_large_scale_liquid=state.mixing_ratio.large_scale_liquid,
            mixing_ratio_large_scale_ice=state.mixing_ratio.large_scale_ice,
            mixing_ratio_convective_liquid=state.mixing_ratio.convective_liquid,
            mixing_ratio_rain=state.mixing_ratio.rain,
            mixing_ratio_snow=state.mixing_ratio.snow,
            mixing_ratio_graupel=state.mixing_ratio.graupel,
            mixing_ratio_convective_ice=state.mixing_ratio.convective_ice,
            cloud_fraction_large_scale=state.cloud_fraction.large_scale,
            cloud_fraction_convective=state.cloud_fraction.convective,
            concentration_liquid=state.concentration.liquid,
            concentration_ice=state.concentration.ice,
            liquid_radius=state.cloud_particle_effective_radius.liquid,
            ice_radius=state.cloud_particle_effective_radius.ice,
            relative_humidity_after_pdf=state.relative_humidity_after_pdf,
            radiation_vapor=state.radiation_field.vapor,
            radiation_liquid=state.radiation_field.liquid,
            radiation_ice=state.radiation_field.ice,
            radiation_rain=state.radiation_field.rain,
            radiation_snow=state.radiation_field.snow,
            radiation_graupel=state.radiation_field.graupel,
            radiation_cloud_fraction=state.radiation_field.cloud_fraction,
            local_p_mb=locals_.p_mb,
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
            "local_p_mb": locals_.p_mb.field[:],
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
