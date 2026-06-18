from typing import Any

from f90nml import Namelist
from ndsl import StencilFactory
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py

from pyMoist.microphysics.GFDL_1M.config import GFDL1MConfig
from pyMoist.microphysics.GFDL_1M.locals import GFDL1MLocals
from pyMoist.microphysics.GFDL_1M.setup import GFDL1MSetup
from pyMoist.microphysics.GFDL_1M.state import GFDL1MState
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable


class TranslateGFDL_1M_Setup(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.quantity_factory = grid.quantity_factory

        # grid.compute_dict is workaround to remove grid halo, which is hardcoded to 3
        self.in_vars["data_vars"] = {
            "p_interface": {},
            "z_interface": {},
            "t": {},
            "u": {},
            "v": {},
            "mixing_ratio_vapor": {},
            "mixing_ratio_rain": {},
            "mixing_ratio_snow": {},
            "mixing_ratio_graupel": {},
            "mixing_ratio_large_scale_liquid": {},
            "mixing_ratio_convective_liquid": {},
            "cloud_fraction_convective": {},
            "cloud_fraction_large_scale": {},
            "mixing_ratio_large_scale_ice": {},
            "mixing_ratio_convective_ice": {},
        }

        self.out_vars: dict[str, Any] = {
            "local_p_interface_mb": {},
            "local_p_mb": {},
            "local_edge_height_above_surface": {},
            "local_layer_height_above_surface": {},
            "local_layer_thickness": {},
            "local_dp": {},
            "local_mass": {},
            "local_mass_inverse": {},
            "local_u_unmodified": {},
            "local_v_unmodified": {},
            "local_saturation_specific_humidity": {},
            "local_dsaturation_specific_humidity": {},
            "local_lcl_level": {},
            "lower_tropospheric_stability": {},
            "estimated_inversion_strength": {},
            "mixing_ratio_rain": {},
            "mixing_ratio_snow": {},
            "dudt_macro": {},
            "dvdt_macro": {},
            "dtdt_macro": {},
            "dvapordt_macro": {},
            "dliquiddt_macro": {},
            "dicedt_macro": {},
            "dcloud_fractiondt_macro": {},
            "draindt_macro": {},
            "dsnowdt_macro": {},
            "dgraupeldt_macro": {},
        }

        # Initialize saturation tables
        self.saturation_tables = SaturationVaporPressureTable(self.stencil_factory.backend)

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GFDL_1M-constants")

    def compute(self, inputs):
        # initialize constants
        config = GFDL1MConfig(**self.constants)

        # initialize dataclasses
        state = GFDL1MState.zeros(self.quantity_factory)
        locals_ = GFDL1MLocals.make_as_state(self.quantity_factory)

        # fill relevant parts of dataclasses
        state.p_interface.field[:] = inputs["p_interface"]
        state.z_interface.field[:] = inputs["z_interface"]
        state.t.field[:] = inputs["t"]
        state.u.field[:] = inputs["u"]
        state.v.field[:] = inputs["v"]
        state.mixing_ratio.vapor.field[:] = inputs["mixing_ratio_vapor"]
        state.mixing_ratio.rain.field[:] = inputs["mixing_ratio_rain"]
        state.mixing_ratio.snow.field[:] = inputs["mixing_ratio_snow"]
        state.mixing_ratio.graupel.field[:] = inputs["mixing_ratio_graupel"]
        state.mixing_ratio.large_scale_liquid.field[:] = inputs["mixing_ratio_large_scale_liquid"]
        state.mixing_ratio.convective_liquid.field[:] = inputs["mixing_ratio_convective_liquid"]
        state.cloud_fraction.convective.field[:] = inputs["cloud_fraction_convective"]
        state.cloud_fraction.large_scale.field[:] = inputs["cloud_fraction_large_scale"]
        state.mixing_ratio.large_scale_ice.field[:] = inputs["mixing_ratio_large_scale_ice"]
        state.mixing_ratio.convective_ice.field[:] = inputs["mixing_ratio_convective_ice"]

        # initialize test class
        code = GFDL1MSetup(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            saturation_tables=self.saturation_tables,
        )

        # execute test code
        code(
            p_interface=state.p_interface,
            z_interface=state.z_interface,
            u=state.u,
            v=state.v,
            t=state.t,
            lcl_height=state.lcl_height,
            lower_tropospheric_stability=state.lower_tropospheric_stability,
            estimated_inversion_strength=state.estimated_inversion_strength,
            mixing_ratio_vapor=state.mixing_ratio.vapor,
            mixing_ratio_rain=state.mixing_ratio.rain,
            mixing_ratio_snow=state.mixing_ratio.snow,
            mixing_ratio_graupel=state.mixing_ratio.graupel,
            mixing_ratio_convective_liquid=state.mixing_ratio.convective_liquid,
            mixing_ratio_convective_ice=state.mixing_ratio.convective_ice,
            mixing_ratio_large_scale_liquid=state.mixing_ratio.large_scale_liquid,
            mixing_ratio_large_scale_ice=state.mixing_ratio.large_scale_ice,
            cloud_fraction_convective=state.cloud_fraction.convective,
            cloud_fraction_large_scale=state.cloud_fraction.large_scale,
            shallow_convection_rain=state.shallow_convection_rain,
            shallow_convection_snow=state.shallow_convection_snow,
            dudt_macro=state.tendencies.dudt_macro,
            dvdt_macro=state.tendencies.dvdt_macro,
            dtdt_macro=state.tendencies.dtdt_macro,
            dvapordt_macro=state.tendencies.dvapordt_macro,
            dliquiddt_macro=state.tendencies.dliquiddt_macro,
            dicedt_macro=state.tendencies.dicedt_macro,
            dcloud_fractiondt_macro=state.tendencies.dcloud_fractiondt_macro,
            draindt_macro=state.tendencies.draindt_macro,
            dsnowdt_macro=state.tendencies.dsnowdt_macro,
            dgraupeldt_macro=state.tendencies.dgraupeldt_macro,
            shallow_convective_precipitation=state.precipitation_at_surface.shallow_convective_precipitation,
            deep_convective_precipitation=state.precipitation_at_surface.deep_convective_precipitation,
            anvil_precipitation=state.precipitation_at_surface.anvil_precipitation,
            shallow_convective_snow=state.precipitation_at_surface.shallow_convective_snow,
            deep_convective_snow=state.precipitation_at_surface.deep_convective_snow,
            anvil_snow=state.precipitation_at_surface.anvil_snow,
            local_p_mb=locals_.p_mb,
            local_p_interface_mb=locals_.p_interface_mb,
            local_edge_height_above_surface=locals_.edge_height_above_surface,
            local_layer_height_above_surface=locals_.layer_height_above_surface,
            local_layer_thickness=locals_.layer_thickness,
            local_layer_thickness_negative=locals_.layer_thickness_negative,
            local_dp=locals_.dp,
            local_mass=locals_.mass,
            local_mass_inverse=locals_.mass_inverse,
            local_saturation_specific_humidity=locals_.saturation_specific_humidity,
            local_dsaturation_specific_humidity=locals_.dsaturation_specific_humidity,
            local_u_unmodified=locals_.u_unmodified,
            local_v_unmodified=locals_.v_unmodified,
            local_lcl_level=locals_.lcl_level,
        )

        return {
            "local_p_interface_mb": locals_.p_interface_mb.field[:],
            "local_p_mb": locals_.p_mb.field[:],
            "local_edge_height_above_surface": locals_.edge_height_above_surface.field[:],
            "local_layer_height_above_surface": locals_.layer_height_above_surface.field[:],
            "local_layer_thickness": locals_.layer_thickness.field[:],
            "local_dp": locals_.dp.field[:],
            "local_mass": locals_.mass.field[:],
            "local_mass_inverse": locals_.mass_inverse.field[:],
            "local_u_unmodified": locals_.u_unmodified.field[:],
            "local_v_unmodified": locals_.v_unmodified.field[:],
            "local_saturation_specific_humidity": locals_.saturation_specific_humidity.field[:],
            "local_dsaturation_specific_humidity": locals_.dsaturation_specific_humidity.field[:],
            "local_lcl_level": locals_.lcl_level.field[:] + 1,  # add one to shift back to fortran indexing
            "lower_tropospheric_stability": state.lower_tropospheric_stability.field[:],
            "estimated_inversion_strength": state.estimated_inversion_strength.field[:],
            "mixing_ratio_rain": state.mixing_ratio.rain.field[:],
            "mixing_ratio_snow": state.mixing_ratio.snow.field[:],
            "dudt_macro": state.tendencies.dudt_macro.field[:],
            "dvdt_macro": state.tendencies.dvdt_macro.field[:],
            "dtdt_macro": state.tendencies.dtdt_macro.field[:],
            "dvapordt_macro": state.tendencies.dvapordt_macro.field[:],
            "dliquiddt_macro": state.tendencies.dliquiddt_macro.field[:],
            "dicedt_macro": state.tendencies.dicedt_macro.field[:],
            "dcloud_fractiondt_macro": state.tendencies.dcloud_fractiondt_macro.field[:],
            "draindt_macro": state.tendencies.draindt_macro.field[:],
            "dsnowdt_macro": state.tendencies.dsnowdt_macro.field[:],
            "dgraupeldt_macro": state.tendencies.dgraupeldt_macro.field[:],
        }
