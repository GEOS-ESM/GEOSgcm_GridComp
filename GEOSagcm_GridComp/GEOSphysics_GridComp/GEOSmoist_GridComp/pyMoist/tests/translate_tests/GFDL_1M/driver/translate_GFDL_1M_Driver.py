import copy
import time

from f90nml import Namelist
from ndsl import StencilFactory, ndsl_log
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array

from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.GFDL_1M.driver.driver import GFDL1MDriver
from pyMoist.GFDL_1M.locals import GFDL1MLocals
from pyMoist.GFDL_1M.state import GFDL1MState


class TranslateGFDL_1M_Driver(TranslateFortranData2Py):
    def __init__(self, grid: Grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "radiation_vapor": {},
            "radiation_liquid": {},
            "radiation_rain": {},
            "radiation_ice": {},
            "radiation_snow": {},
            "radiation_graupel": {},
            "radiation_cloud_fraction": {},
            "local_total_concentration": {},
            "local_dvapordt_driver": {},
            "local_dliquiddt_driver": {},
            "local_draindt_driver": {},
            "local_dicedt_driver": {},
            "local_dsnowdt_driver": {},
            "local_dgraupeldt_driver": {},
            "local_dcloudfractiondt_driver": {},
            "local_dtdt_driver": {},
            "t": {},
            "w": {},
            "u": {},
            "v": {},
            "local_dudt_driver": {},
            "local_dvdt_driver": {},
            "local_dz": {},
            "local_dp": {},
            "area": {},
            "land_fraction": {},
            "convection_fraction": {},
            "surface_type": {},
            "estimated_inversion_strength": {},
            "critical_relative_humidity_for_pdf": {},
            "non_anvil_large_scale_evaporation": {},
            "non_anvil_large_scale_sublimation": {},
            "surface_precip_rain": {},
            "surface_precip_snow": {},
            "surface_precip_ice": {},
            "surface_precip_graupel": {},
            "non_anvil_large_scale_liquid_precip_flux": {},
            "non_anvil_large_scale_ice_precip_flux": {},
        }

        self.out_vars = self.in_vars["data_vars"].copy()
        self.out_vars.update(
            {
                "DEBUG_vapor_unmodified_setup": {},
                "DEBUG_vapor_modified_setup": {},
                "DEBUG_vapor_unmodified_fallspeed": {},
                "DEBUG_vapor_modified_fallspeed": {},
                "DEBUG_vapor_unmodified_terminalfall": {},
                "DEBUG_vapor_modified_terminalfall": {},
                "DEBUG_vapor_unmodified_warmrain": {},
                "DEBUG_vapor_modified_warmrain": {},
                "DEBUG_vapor_unmodified_icecloud": {},
                "DEBUG_vapor_modified_icecloud": {},
                #
                "DEBUG_TERMINALFALL_IN_driver_local_t_terminalfall": {},
                "DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_vapor_terminalfall": {},
                "DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_liquid_terminalfall": {},
                "DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_rain_terminalfall": {},
                "DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_graupel_terminalfall": {},
                "DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_snow_terminalfall": {},
                "DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_ice_terminalfall": {},
                "DEBUG_TERMINALFALL_IN_driver_local_ice_precip_flux_terminalfall": {},
                "DEBUG_TERMINALFALL_IN_driver_local_w_terminalfall": {},
                "DEBUG_TERMINALFALL_IN_driver_local_dz_terminalfall": {},
                "DEBUG_TERMINALFALL_IN_driver_local_dp_terminalfall": {},
                "DEBUG_TERMINALFALL_IN_driver_local_terminal_speed_ice_terminalfall": {},
                "DEBUG_TERMINALFALL_IN_driver_local_terminal_speed_snow_terminalfall": {},
                "DEBUG_TERMINALFALL_IN_driver_local_terminal_speed_graupel_terminalfall": {},
                "DEBUG_TERMINALFALL_IN_surface_precip_rain_terminalfall": {},
                "DEBUG_TERMINALFALL_IN_surface_precip_snow_terminalfall": {},
                "DEBUG_TERMINALFALL_IN_surface_precip_graupel_terminalfall": {},
                "DEBUG_TERMINALFALL_IN_surface_precip_ice_terminalfall": {},
                #
                "DEBUG_TERMINALFALL_OUT_driver_local_t_terminalfall": {},
                "DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_vapor_terminalfall": {},
                "DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_liquid_terminalfall": {},
                "DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_rain_terminalfall": {},
                "DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_graupel_terminalfall": {},
                "DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_snow_terminalfall": {},
                "DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_ice_terminalfall": {},
                "DEBUG_TERMINALFALL_OUT_driver_local_ice_precip_flux_terminalfall": {},
                "DEBUG_TERMINALFALL_OUT_driver_local_w_terminalfall": {},
                "DEBUG_TERMINALFALL_OUT_driver_local_dz_terminalfall": {},
                "DEBUG_TERMINALFALL_OUT_driver_local_dp_terminalfall": {},
                "DEBUG_TERMINALFALL_OUT_driver_local_terminal_speed_ice_terminalfall": {},
                "DEBUG_TERMINALFALL_OUT_driver_local_terminal_speed_snow_terminalfall": {},
                "DEBUG_TERMINALFALL_OUT_driver_local_terminal_speed_graupel_terminalfall": {},
                "DEBUG_TERMINALFALL_OUT_surface_precip_rain_terminalfall": {},
                "DEBUG_TERMINALFALL_OUT_surface_precip_snow_terminalfall": {},
                "DEBUG_TERMINALFALL_OUT_surface_precip_graupel_terminalfall": {},
                "DEBUG_TERMINALFALL_OUT_surface_precip_ice_terminalfall": {},
                #
                "DEBUG_WARMRAIN_IN_driver_local_dp_warmrain": {},
                "DEBUG_WARMRAIN_IN_driver_local_dz_warmrain": {},
                "DEBUG_WARMRAIN_IN_driver_local_t_warmrain": {},
                "DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_vapor_warmrain": {},
                "DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_liquid_warmrain": {},
                "DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_rain_warmrain": {},
                "DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_ice_warmrain": {},
                "DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_snow_warmrain": {},
                "DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_graupel_warmrain": {},
                "DEBUG_WARMRAIN_IN_driver_local_cloud_fraction_warmrain": {},
                "DEBUG_WARMRAIN_IN_driver_local_ccn_warmrain": {},
                "DEBUG_WARMRAIN_IN_driver_local_density_warmrain": {},
                "DEBUG_WARMRAIN_IN_driver_local_density_factor_warmrain": {},
                "DEBUG_WARMRAIN_IN_driver_local_c_praut_warmrain": {},
                "DEBUG_WARMRAIN_IN_driver_local_terminal_speed_rain_warmrain": {},
                "DEBUG_WARMRAIN_IN_driver_local_evaporation_warmrain": {},
                "DEBUG_WARMRAIN_IN_driver_local_liquid_precip_flux_warmrain": {},
                "DEBUG_WARMRAIN_IN_driver_local_w_warmrain": {},
                "DEBUG_WARMRAIN_IN_driver_local_rh_limited_warmrain": {},
                "DEBUG_WARMRAIN_IN_non_anvil_large_scale_evaporation_warmrain": {},
                "DEBUG_WARMRAIN_IN_non_anvil_large_scale_liquid_precip_flux_warmrain": {},
                "DEBUG_WARMRAIN_IN_non_anvil_large_scale_ice_precip_flux_warmrain": {},
                "DEBUG_WARMRAIN_IN_driver_local_mass_warmrain": {},
                "DEBUG_WARMRAIN_IN_driver_local_ice_precip_flux_warmrain": {},
                "DEBUG_WARMRAIN_IN_driver_local_rain_warmrain": {},
                "DEBUG_WARMRAIN_IN_surface_precip_rain_warmrain": {},
                "DEBUG_WARMRAIN_IN_estimated_inversion_strength_warmrain": {},
                "DEBUG_WARMRAIN_IN_driver_local_one_minus_sigma_warmrain": {},
                #
                "DEBUG_WARMRAIN_OUT_driver_local_dp_warmrain": {},
                "DEBUG_WARMRAIN_OUT_driver_local_dz_warmrain": {},
                "DEBUG_WARMRAIN_OUT_driver_local_t_warmrain": {},
                "DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_vapor_warmrain": {},
                "DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_liquid_warmrain": {},
                "DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_rain_warmrain": {},
                "DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_ice_warmrain": {},
                "DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_snow_warmrain": {},
                "DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_graupel_warmrain": {},
                "DEBUG_WARMRAIN_OUT_driver_local_cloud_fraction_warmrain": {},
                "DEBUG_WARMRAIN_OUT_driver_local_ccn_warmrain": {},
                "DEBUG_WARMRAIN_OUT_driver_local_density_warmrain": {},
                "DEBUG_WARMRAIN_OUT_driver_local_density_factor_warmrain": {},
                "DEBUG_WARMRAIN_OUT_driver_local_c_praut_warmrain": {},
                "DEBUG_WARMRAIN_OUT_driver_local_terminal_speed_rain_warmrain": {},
                "DEBUG_WARMRAIN_OUT_driver_local_evaporation_warmrain": {},
                "DEBUG_WARMRAIN_OUT_driver_local_liquid_precip_flux_warmrain": {},
                "DEBUG_WARMRAIN_OUT_driver_local_w_warmrain": {},
                "DEBUG_WARMRAIN_OUT_driver_local_rh_limited_warmrain": {},
                "DEBUG_WARMRAIN_OUT_non_anvil_large_scale_evaporation_warmrain": {},
                "DEBUG_WARMRAIN_OUT_non_anvil_large_scale_liquid_precip_flux_warmrain": {},
                "DEBUG_WARMRAIN_OUT_non_anvil_large_scale_ice_precip_flux_warmrain": {},
                "DEBUG_WARMRAIN_OUT_driver_local_mass_warmrain": {},
                "DEBUG_WARMRAIN_OUT_driver_local_ice_precip_flux_warmrain": {},
                "DEBUG_WARMRAIN_OUT_driver_local_rain_warmrain": {},
                "DEBUG_WARMRAIN_OUT_surface_precip_rain_warmrain": {},
                "DEBUG_WARMRAIN_OUT_estimated_inversion_strength_warmrain": {},
                "DEBUG_WARMRAIN_OUT_driver_local_one_minus_sigma_warmrain": {},
                #
                "DEBUG_ICECLOUD_IN_driver_local_t_icecloud": {},
                "DEBUG_ICECLOUD_IN_driver_local_p_dry_icecloud": {},
                "DEBUG_ICECLOUD_IN_driver_local_dp_icecloud": {},
                "DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_vapor_icecloud": {},
                "DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_liquid_icecloud": {},
                "DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_rain_icecloud": {},
                "DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_ice_icecloud": {},
                "DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_snow_icecloud": {},
                "DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_graupel_icecloud": {},
                "DEBUG_ICECLOUD_IN_driver_local_cloud_fraction_icecloud": {},
                "DEBUG_ICECLOUD_IN_driver_local_terminal_speed_snow_icecloud": {},
                "DEBUG_ICECLOUD_IN_driver_local_terminal_speed_graupel_icecloud": {},
                "DEBUG_ICECLOUD_IN_driver_local_terminal_speed_rain_icecloud": {},
                "DEBUG_ICECLOUD_IN_driver_local_density_icecloud": {},
                "DEBUG_ICECLOUD_IN_driver_local_density_factor_icecloud": {},
                "DEBUG_ICECLOUD_IN_driver_local_rh_limited_icecloud": {},
                "DEBUG_ICECLOUD_IN_non_anvil_large_scale_sublimation_icecloud": {},
                "DEBUG_ICECLOUD_IN_driver_local_ccn_icecloud": {},
                "DEBUG_ICECLOUD_IN_convection_fraction_icecloud": {},
                "DEBUG_ICECLOUD_IN_surface_type_icecloud": {},
                #
                "DEBUG_ICECLOUD_OUT_driver_local_t_icecloud": {},
                "DEBUG_ICECLOUD_OUT_driver_local_p_dry_icecloud": {},
                "DEBUG_ICECLOUD_OUT_driver_local_dp_icecloud": {},
                "DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_vapor_icecloud": {},
                "DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_liquid_icecloud": {},
                "DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_rain_icecloud": {},
                "DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_ice_icecloud": {},
                "DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_snow_icecloud": {},
                "DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_graupel_icecloud": {},
                "DEBUG_ICECLOUD_OUT_driver_local_cloud_fraction_icecloud": {},
                "DEBUG_ICECLOUD_OUT_driver_local_terminal_speed_snow_icecloud": {},
                "DEBUG_ICECLOUD_OUT_driver_local_terminal_speed_graupel_icecloud": {},
                "DEBUG_ICECLOUD_OUT_driver_local_terminal_speed_rain_icecloud": {},
                "DEBUG_ICECLOUD_OUT_driver_local_density_icecloud": {},
                "DEBUG_ICECLOUD_OUT_driver_local_density_factor_icecloud": {},
                "DEBUG_ICECLOUD_OUT_driver_local_rh_limited_icecloud": {},
                "DEBUG_ICECLOUD_OUT_non_anvil_large_scale_sublimation_icecloud": {},
                "DEBUG_ICECLOUD_OUT_driver_local_ccn_icecloud": {},
                "DEBUG_ICECLOUD_OUT_convection_fraction_icecloud": {},
                "DEBUG_ICECLOUD_OUT_surface_type_icecloud": {},
            }
        )

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GFDL_1M-constants")

    def compute(self, inputs):
        # initialize dataclasses
        state = GFDL1MState.zeros(self.quantity_factory)
        locals_ = GFDL1MLocals.make_as_state(self.quantity_factory)

        # initialize constants
        config = GFDL1MConfig(**self.constants)

        # fill relevant parts of dataclasses with input data
        safe_assign_array(state.radiation_field.vapor.field[:], inputs["radiation_vapor"])
        safe_assign_array(state.radiation_field.liquid.field[:], inputs["radiation_liquid"])
        safe_assign_array(state.radiation_field.rain.field[:], inputs["radiation_rain"])
        safe_assign_array(state.radiation_field.ice.field[:], inputs["radiation_ice"])
        safe_assign_array(state.radiation_field.snow.field[:], inputs["radiation_snow"])
        safe_assign_array(state.radiation_field.graupel.field[:], inputs["radiation_graupel"])
        safe_assign_array(state.radiation_field.cloud_fraction.field[:], inputs["radiation_cloud_fraction"])
        safe_assign_array(locals_.total_concentration.field[:], inputs["local_total_concentration"])
        safe_assign_array(locals_.driver_tendencies.dvapordt.field[:], inputs["local_dvapordt_driver"])
        safe_assign_array(locals_.driver_tendencies.dliquiddt.field[:], inputs["local_dliquiddt_driver"])
        safe_assign_array(locals_.driver_tendencies.draindt.field[:], inputs["local_draindt_driver"])
        safe_assign_array(locals_.driver_tendencies.dicedt.field[:], inputs["local_dicedt_driver"])
        safe_assign_array(locals_.driver_tendencies.dsnowdt.field[:], inputs["local_dsnowdt_driver"])
        safe_assign_array(locals_.driver_tendencies.dgraupeldt.field[:], inputs["local_dgraupeldt_driver"])
        safe_assign_array(
            locals_.driver_tendencies.dcloudfractiondt.field[:], inputs["local_dcloudfractiondt_driver"]
        )
        safe_assign_array(locals_.driver_tendencies.dtdt.field[:], inputs["local_dtdt_driver"])
        safe_assign_array(locals_.driver_tendencies.dudt.field[:], inputs["local_dudt_driver"])
        safe_assign_array(locals_.driver_tendencies.dvdt.field[:], inputs["local_dvdt_driver"])
        safe_assign_array(state.t.field[:], inputs["t"])
        safe_assign_array(state.u.field[:], inputs["u"])
        safe_assign_array(state.v.field[:], inputs["v"])
        safe_assign_array(state.vertical_motion.velocity.field[:], inputs["w"])
        safe_assign_array(locals_.layer_thickness_negative.field[:], inputs["local_dz"])
        safe_assign_array(locals_.dp.field[:], inputs["local_dp"])
        safe_assign_array(state.area.field[:], inputs["area"])
        safe_assign_array(state.land_fraction.field[:], inputs["land_fraction"])
        safe_assign_array(state.convection_fraction.field[:], inputs["convection_fraction"])
        safe_assign_array(state.surface_type.field[:], inputs["surface_type"])
        safe_assign_array(state.estimated_inversion_strength.field[:], inputs["estimated_inversion_strength"])
        safe_assign_array(
            state.critical_relative_humidity_for_pdf.field[:], inputs["critical_relative_humidity_for_pdf"]
        )
        safe_assign_array(
            state.non_anvil_large_scale.evaporation.field[:], inputs["non_anvil_large_scale_evaporation"]
        )
        safe_assign_array(
            state.non_anvil_large_scale.sublimation.field[:], inputs["non_anvil_large_scale_sublimation"]
        )
        safe_assign_array(
            state.non_anvil_large_scale.liquid_precip_flux.field[:],
            inputs["non_anvil_large_scale_liquid_precip_flux"],
        )
        safe_assign_array(
            state.non_anvil_large_scale.ice_precip_flux.field[:],
            inputs["non_anvil_large_scale_ice_precip_flux"],
        )
        safe_assign_array(state.precipitation_at_surface.rain.field[:], inputs["surface_precip_rain"])
        safe_assign_array(state.precipitation_at_surface.snow.field[:], inputs["surface_precip_snow"])
        safe_assign_array(state.precipitation_at_surface.ice.field[:], inputs["surface_precip_ice"])
        safe_assign_array(state.precipitation_at_surface.graupel.field[:], inputs["surface_precip_graupel"])

        # construct test stencil
        code = GFDL1MDriver(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
        )

        from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
        DEBUG_vapor_unmodified_setup = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_vapor_modified_setup = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_vapor_unmodified_fallspeed = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_vapor_modified_fallspeed = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_vapor_unmodified_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_vapor_modified_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_vapor_unmodified_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_vapor_modified_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_vapor_unmodified_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_vapor_modified_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        #
        DEBUG_TERMINALFALL_IN_driver_local_t_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_vapor_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_liquid_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_rain_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_graupel_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_snow_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_ice_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_TERMINALFALL_IN_driver_local_ice_precip_flux_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_TERMINALFALL_IN_driver_local_w_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_TERMINALFALL_IN_driver_local_dz_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_TERMINALFALL_IN_driver_local_dp_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_TERMINALFALL_IN_driver_local_terminal_speed_ice_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_TERMINALFALL_IN_driver_local_terminal_speed_snow_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_TERMINALFALL_IN_driver_local_terminal_speed_graupel_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_TERMINALFALL_IN_surface_precip_rain_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_TERMINALFALL_IN_surface_precip_snow_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_TERMINALFALL_IN_surface_precip_graupel_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_TERMINALFALL_IN_surface_precip_ice_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        #
        DEBUG_TERMINALFALL_OUT_driver_local_t_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_vapor_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_liquid_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_rain_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_graupel_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_snow_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_ice_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_TERMINALFALL_OUT_driver_local_ice_precip_flux_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_TERMINALFALL_OUT_driver_local_w_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_TERMINALFALL_OUT_driver_local_dz_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_TERMINALFALL_OUT_driver_local_dp_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_TERMINALFALL_OUT_driver_local_terminal_speed_ice_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_TERMINALFALL_OUT_driver_local_terminal_speed_snow_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_TERMINALFALL_OUT_driver_local_terminal_speed_graupel_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_TERMINALFALL_OUT_surface_precip_rain_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_TERMINALFALL_OUT_surface_precip_snow_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_TERMINALFALL_OUT_surface_precip_graupel_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_TERMINALFALL_OUT_surface_precip_ice_terminalfall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        #
        DEBUG_WARMRAIN_IN_driver_local_dp_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_IN_driver_local_dz_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_IN_driver_local_t_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_vapor_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_liquid_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_rain_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_ice_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_snow_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_graupel_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_IN_driver_local_cloud_fraction_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_IN_driver_local_ccn_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_IN_driver_local_density_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_IN_driver_local_density_factor_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_IN_driver_local_c_praut_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_IN_driver_local_terminal_speed_rain_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_IN_driver_local_evaporation_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_IN_driver_local_liquid_precip_flux_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_IN_driver_local_w_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_IN_driver_local_rh_limited_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_IN_non_anvil_large_scale_evaporation_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_IN_non_anvil_large_scale_liquid_precip_flux_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a")
        DEBUG_WARMRAIN_IN_non_anvil_large_scale_ice_precip_flux_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a")
        DEBUG_WARMRAIN_IN_driver_local_mass_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_IN_driver_local_ice_precip_flux_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_IN_driver_local_rain_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_IN_surface_precip_rain_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_IN_estimated_inversion_strength_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_IN_driver_local_one_minus_sigma_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        #
        DEBUG_WARMRAIN_OUT_driver_local_dp_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_OUT_driver_local_dz_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_OUT_driver_local_t_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_vapor_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_liquid_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_rain_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_ice_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_snow_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_graupel_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_OUT_driver_local_cloud_fraction_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_OUT_driver_local_ccn_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_OUT_driver_local_density_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_OUT_driver_local_density_factor_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_OUT_driver_local_c_praut_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_OUT_driver_local_terminal_speed_rain_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_OUT_driver_local_evaporation_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_OUT_driver_local_liquid_precip_flux_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_OUT_driver_local_w_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_OUT_driver_local_rh_limited_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_OUT_non_anvil_large_scale_evaporation_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_OUT_non_anvil_large_scale_liquid_precip_flux_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a")
        DEBUG_WARMRAIN_OUT_non_anvil_large_scale_ice_precip_flux_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a")
        DEBUG_WARMRAIN_OUT_driver_local_mass_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_OUT_driver_local_ice_precip_flux_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_OUT_driver_local_rain_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_OUT_surface_precip_rain_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_OUT_estimated_inversion_strength_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_WARMRAIN_OUT_driver_local_one_minus_sigma_warmrain = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        #
        DEBUG_ICECLOUD_IN_driver_local_t_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_IN_driver_local_p_dry_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_IN_driver_local_dp_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_vapor_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_liquid_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_rain_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_ice_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_snow_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_graupel_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_IN_driver_local_cloud_fraction_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_IN_driver_local_terminal_speed_snow_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_IN_driver_local_terminal_speed_graupel_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_IN_driver_local_terminal_speed_rain_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_IN_driver_local_density_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_IN_driver_local_density_factor_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_IN_driver_local_rh_limited_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_IN_non_anvil_large_scale_sublimation_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_IN_driver_local_ccn_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_IN_convection_fraction_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_IN_surface_type_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        #
        DEBUG_ICECLOUD_OUT_driver_local_t_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_OUT_driver_local_p_dry_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_OUT_driver_local_dp_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_vapor_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_liquid_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_rain_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_ice_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_snow_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_graupel_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_OUT_driver_local_cloud_fraction_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_OUT_driver_local_terminal_speed_snow_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_OUT_driver_local_terminal_speed_graupel_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_OUT_driver_local_terminal_speed_rain_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_OUT_driver_local_density_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_OUT_driver_local_density_factor_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_OUT_driver_local_rh_limited_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_OUT_non_anvil_large_scale_sublimation_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_OUT_driver_local_ccn_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_OUT_convection_fraction_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DEBUG_ICECLOUD_OUT_surface_type_icecloud = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        code(
            t=state.t,
            u=state.u,
            v=state.v,
            w=state.vertical_motion.velocity,
            dz=locals_.layer_thickness_negative,
            dp=locals_.dp,
            area=state.area,
            land_fraction=state.land_fraction,
            convection_fraction=state.convection_fraction,
            surface_type=state.surface_type,
            estimated_inversion_strength=state.estimated_inversion_strength,
            critical_relative_humidity_for_pdf=state.critical_relative_humidity_for_pdf,
            vapor=state.radiation_field.vapor,
            liquid=state.radiation_field.liquid,
            rain=state.radiation_field.rain,
            ice=state.radiation_field.ice,
            snow=state.radiation_field.snow,
            graupel=state.radiation_field.graupel,
            cloud_fraction=state.radiation_field.cloud_fraction,
            total_concentration=locals_.total_concentration,
            dvapordt=locals_.driver_tendencies.dvapordt,
            dliquiddt=locals_.driver_tendencies.dliquiddt,
            draindt=locals_.driver_tendencies.draindt,
            dicedt=locals_.driver_tendencies.dicedt,
            dsnowdt=locals_.driver_tendencies.dsnowdt,
            dgraupeldt=locals_.driver_tendencies.dgraupeldt,
            dcloudfractiondt=locals_.driver_tendencies.dcloudfractiondt,
            dtdt=locals_.driver_tendencies.dtdt,
            dudt=locals_.driver_tendencies.dudt,
            dvdt=locals_.driver_tendencies.dvdt,
            liquid_precip_flux=state.non_anvil_large_scale.liquid_precip_flux,
            ice_precip_flux=state.non_anvil_large_scale.ice_precip_flux,
            evaporation=state.non_anvil_large_scale.evaporation,
            sublimation=state.non_anvil_large_scale.sublimation,
            surface_precip_rain=state.precipitation_at_surface.rain,
            surface_precip_snow=state.precipitation_at_surface.snow,
            surface_precip_ice=state.precipitation_at_surface.ice,
            surface_precip_graupel=state.precipitation_at_surface.graupel,
            #
            DEBUG_vapor_unmodified_setup=DEBUG_vapor_unmodified_setup,
            DEBUG_vapor_modified_setup=DEBUG_vapor_modified_setup,
            DEBUG_vapor_unmodified_fallspeed=DEBUG_vapor_unmodified_fallspeed,
            DEBUG_vapor_modified_fallspeed=DEBUG_vapor_modified_fallspeed,
            DEBUG_vapor_unmodified_terminalfall=DEBUG_vapor_unmodified_terminalfall,
            DEBUG_vapor_modified_terminalfall=DEBUG_vapor_modified_terminalfall,
            DEBUG_vapor_unmodified_warmrain=DEBUG_vapor_unmodified_warmrain,
            DEBUG_vapor_modified_warmrain=DEBUG_vapor_modified_warmrain,
            DEBUG_vapor_unmodified_icecloud=DEBUG_vapor_unmodified_icecloud,
            DEBUG_vapor_modified_icecloud=DEBUG_vapor_modified_icecloud,
            #
            DEBUG_TERMINALFALL_IN_driver_local_t_terminalfall=DEBUG_TERMINALFALL_IN_driver_local_t_terminalfall,
            DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_vapor_terminalfall=DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_vapor_terminalfall,
            DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_liquid_terminalfall=DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_liquid_terminalfall,
            DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_rain_terminalfall=DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_rain_terminalfall,
            DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_graupel_terminalfall=DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_graupel_terminalfall,
            DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_snow_terminalfall=DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_snow_terminalfall,
            DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_ice_terminalfall=DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_ice_terminalfall,
            DEBUG_TERMINALFALL_IN_driver_local_ice_precip_flux_terminalfall=DEBUG_TERMINALFALL_IN_driver_local_ice_precip_flux_terminalfall,
            DEBUG_TERMINALFALL_IN_driver_local_w_terminalfall=DEBUG_TERMINALFALL_IN_driver_local_w_terminalfall,
            DEBUG_TERMINALFALL_IN_driver_local_dz_terminalfall=DEBUG_TERMINALFALL_IN_driver_local_dz_terminalfall,
            DEBUG_TERMINALFALL_IN_driver_local_dp_terminalfall=DEBUG_TERMINALFALL_IN_driver_local_dp_terminalfall,
            DEBUG_TERMINALFALL_IN_driver_local_terminal_speed_ice_terminalfall=DEBUG_TERMINALFALL_IN_driver_local_terminal_speed_ice_terminalfall,
            DEBUG_TERMINALFALL_IN_driver_local_terminal_speed_snow_terminalfall=DEBUG_TERMINALFALL_IN_driver_local_terminal_speed_snow_terminalfall,
            DEBUG_TERMINALFALL_IN_driver_local_terminal_speed_graupel_terminalfall=DEBUG_TERMINALFALL_IN_driver_local_terminal_speed_graupel_terminalfall,
            DEBUG_TERMINALFALL_IN_surface_precip_rain_terminalfall=DEBUG_TERMINALFALL_IN_surface_precip_rain_terminalfall,
            DEBUG_TERMINALFALL_IN_surface_precip_snow_terminalfall=DEBUG_TERMINALFALL_IN_surface_precip_snow_terminalfall,
            DEBUG_TERMINALFALL_IN_surface_precip_graupel_terminalfall=DEBUG_TERMINALFALL_IN_surface_precip_graupel_terminalfall,
            DEBUG_TERMINALFALL_IN_surface_precip_ice_terminalfall=DEBUG_TERMINALFALL_IN_surface_precip_ice_terminalfall,
            #
            DEBUG_TERMINALFALL_OUT_driver_local_t_terminalfall=DEBUG_TERMINALFALL_OUT_driver_local_t_terminalfall,
            DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_vapor_terminalfall=DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_vapor_terminalfall,
            DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_liquid_terminalfall=DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_liquid_terminalfall,
            DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_rain_terminalfall=DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_rain_terminalfall,
            DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_graupel_terminalfall=DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_graupel_terminalfall,
            DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_snow_terminalfall=DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_snow_terminalfall,
            DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_ice_terminalfall=DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_ice_terminalfall,
            DEBUG_TERMINALFALL_OUT_driver_local_ice_precip_flux_terminalfall=DEBUG_TERMINALFALL_OUT_driver_local_ice_precip_flux_terminalfall,
            DEBUG_TERMINALFALL_OUT_driver_local_w_terminalfall=DEBUG_TERMINALFALL_OUT_driver_local_w_terminalfall,
            DEBUG_TERMINALFALL_OUT_driver_local_dz_terminalfall=DEBUG_TERMINALFALL_OUT_driver_local_dz_terminalfall,
            DEBUG_TERMINALFALL_OUT_driver_local_dp_terminalfall=DEBUG_TERMINALFALL_OUT_driver_local_dp_terminalfall,
            DEBUG_TERMINALFALL_OUT_driver_local_terminal_speed_ice_terminalfall=DEBUG_TERMINALFALL_OUT_driver_local_terminal_speed_ice_terminalfall,
            DEBUG_TERMINALFALL_OUT_driver_local_terminal_speed_snow_terminalfall=DEBUG_TERMINALFALL_OUT_driver_local_terminal_speed_snow_terminalfall,
            DEBUG_TERMINALFALL_OUT_driver_local_terminal_speed_graupel_terminalfall=DEBUG_TERMINALFALL_OUT_driver_local_terminal_speed_graupel_terminalfall,
            DEBUG_TERMINALFALL_OUT_surface_precip_rain_terminalfall=DEBUG_TERMINALFALL_OUT_surface_precip_rain_terminalfall,
            DEBUG_TERMINALFALL_OUT_surface_precip_snow_terminalfall=DEBUG_TERMINALFALL_OUT_surface_precip_snow_terminalfall,
            DEBUG_TERMINALFALL_OUT_surface_precip_graupel_terminalfall=DEBUG_TERMINALFALL_OUT_surface_precip_graupel_terminalfall,
            DEBUG_TERMINALFALL_OUT_surface_precip_ice_terminalfall=DEBUG_TERMINALFALL_OUT_surface_precip_ice_terminalfall,
            #
            DEBUG_WARMRAIN_IN_driver_local_dp_warmrain=DEBUG_WARMRAIN_IN_driver_local_dp_warmrain,
            DEBUG_WARMRAIN_IN_driver_local_dz_warmrain=DEBUG_WARMRAIN_IN_driver_local_dz_warmrain,
            DEBUG_WARMRAIN_IN_driver_local_t_warmrain=DEBUG_WARMRAIN_IN_driver_local_t_warmrain,
            DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_vapor_warmrain=DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_vapor_warmrain,
            DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_liquid_warmrain=DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_liquid_warmrain,
            DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_rain_warmrain=DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_rain_warmrain,
            DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_ice_warmrain=DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_ice_warmrain,
            DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_snow_warmrain=DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_snow_warmrain,
            DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_graupel_warmrain=DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_graupel_warmrain,
            DEBUG_WARMRAIN_IN_driver_local_cloud_fraction_warmrain=DEBUG_WARMRAIN_IN_driver_local_cloud_fraction_warmrain,
            DEBUG_WARMRAIN_IN_driver_local_ccn_warmrain=DEBUG_WARMRAIN_IN_driver_local_ccn_warmrain,
            DEBUG_WARMRAIN_IN_driver_local_density_warmrain=DEBUG_WARMRAIN_IN_driver_local_density_warmrain,
            DEBUG_WARMRAIN_IN_driver_local_density_factor_warmrain=DEBUG_WARMRAIN_IN_driver_local_density_factor_warmrain,
            DEBUG_WARMRAIN_IN_driver_local_c_praut_warmrain=DEBUG_WARMRAIN_IN_driver_local_c_praut_warmrain,
            DEBUG_WARMRAIN_IN_driver_local_terminal_speed_rain_warmrain=DEBUG_WARMRAIN_IN_driver_local_terminal_speed_rain_warmrain,
            DEBUG_WARMRAIN_IN_driver_local_evaporation_warmrain=DEBUG_WARMRAIN_IN_driver_local_evaporation_warmrain,
            DEBUG_WARMRAIN_IN_driver_local_liquid_precip_flux_warmrain=DEBUG_WARMRAIN_IN_driver_local_liquid_precip_flux_warmrain,
            DEBUG_WARMRAIN_IN_driver_local_w_warmrain=DEBUG_WARMRAIN_IN_driver_local_w_warmrain,
            DEBUG_WARMRAIN_IN_driver_local_rh_limited_warmrain=DEBUG_WARMRAIN_IN_driver_local_rh_limited_warmrain,
            DEBUG_WARMRAIN_IN_non_anvil_large_scale_evaporation_warmrain=DEBUG_WARMRAIN_IN_non_anvil_large_scale_evaporation_warmrain,
            DEBUG_WARMRAIN_IN_non_anvil_large_scale_liquid_precip_flux_warmrain=DEBUG_WARMRAIN_IN_non_anvil_large_scale_liquid_precip_flux_warmrain,
            DEBUG_WARMRAIN_IN_non_anvil_large_scale_ice_precip_flux_warmrain=DEBUG_WARMRAIN_IN_non_anvil_large_scale_ice_precip_flux_warmrain,
            DEBUG_WARMRAIN_IN_driver_local_mass_warmrain=DEBUG_WARMRAIN_IN_driver_local_mass_warmrain,
            DEBUG_WARMRAIN_IN_driver_local_ice_precip_flux_warmrain=DEBUG_WARMRAIN_IN_driver_local_ice_precip_flux_warmrain,
            DEBUG_WARMRAIN_IN_driver_local_rain_warmrain=DEBUG_WARMRAIN_IN_driver_local_rain_warmrain,
            DEBUG_WARMRAIN_IN_surface_precip_rain_warmrain=DEBUG_WARMRAIN_IN_surface_precip_rain_warmrain,
            DEBUG_WARMRAIN_IN_estimated_inversion_strength_warmrain=DEBUG_WARMRAIN_IN_estimated_inversion_strength_warmrain,
            DEBUG_WARMRAIN_IN_driver_local_one_minus_sigma_warmrain=DEBUG_WARMRAIN_IN_driver_local_one_minus_sigma_warmrain,
            #
            DEBUG_WARMRAIN_OUT_driver_local_dp_warmrain=DEBUG_WARMRAIN_OUT_driver_local_dp_warmrain,
            DEBUG_WARMRAIN_OUT_driver_local_dz_warmrain=DEBUG_WARMRAIN_OUT_driver_local_dz_warmrain,
            DEBUG_WARMRAIN_OUT_driver_local_t_warmrain=DEBUG_WARMRAIN_OUT_driver_local_t_warmrain,
            DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_vapor_warmrain=DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_vapor_warmrain,
            DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_liquid_warmrain=DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_liquid_warmrain,
            DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_rain_warmrain=DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_rain_warmrain,
            DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_ice_warmrain=DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_ice_warmrain,
            DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_snow_warmrain=DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_snow_warmrain,
            DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_graupel_warmrain=DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_graupel_warmrain,
            DEBUG_WARMRAIN_OUT_driver_local_cloud_fraction_warmrain=DEBUG_WARMRAIN_OUT_driver_local_cloud_fraction_warmrain,
            DEBUG_WARMRAIN_OUT_driver_local_ccn_warmrain=DEBUG_WARMRAIN_OUT_driver_local_ccn_warmrain,
            DEBUG_WARMRAIN_OUT_driver_local_density_warmrain=DEBUG_WARMRAIN_OUT_driver_local_density_warmrain,
            DEBUG_WARMRAIN_OUT_driver_local_density_factor_warmrain=DEBUG_WARMRAIN_OUT_driver_local_density_factor_warmrain,
            DEBUG_WARMRAIN_OUT_driver_local_c_praut_warmrain=DEBUG_WARMRAIN_OUT_driver_local_c_praut_warmrain,
            DEBUG_WARMRAIN_OUT_driver_local_terminal_speed_rain_warmrain=DEBUG_WARMRAIN_OUT_driver_local_terminal_speed_rain_warmrain,
            DEBUG_WARMRAIN_OUT_driver_local_evaporation_warmrain=DEBUG_WARMRAIN_OUT_driver_local_evaporation_warmrain,
            DEBUG_WARMRAIN_OUT_driver_local_liquid_precip_flux_warmrain=DEBUG_WARMRAIN_OUT_driver_local_liquid_precip_flux_warmrain,
            DEBUG_WARMRAIN_OUT_driver_local_w_warmrain=DEBUG_WARMRAIN_OUT_driver_local_w_warmrain,
            DEBUG_WARMRAIN_OUT_driver_local_rh_limited_warmrain=DEBUG_WARMRAIN_OUT_driver_local_rh_limited_warmrain,
            DEBUG_WARMRAIN_OUT_non_anvil_large_scale_evaporation_warmrain=DEBUG_WARMRAIN_OUT_non_anvil_large_scale_evaporation_warmrain,
            DEBUG_WARMRAIN_OUT_non_anvil_large_scale_liquid_precip_flux_warmrain=DEBUG_WARMRAIN_OUT_non_anvil_large_scale_liquid_precip_flux_warmrain,
            DEBUG_WARMRAIN_OUT_non_anvil_large_scale_ice_precip_flux_warmrain=DEBUG_WARMRAIN_OUT_non_anvil_large_scale_ice_precip_flux_warmrain,
            DEBUG_WARMRAIN_OUT_driver_local_mass_warmrain=DEBUG_WARMRAIN_OUT_driver_local_mass_warmrain,
            DEBUG_WARMRAIN_OUT_driver_local_ice_precip_flux_warmrain=DEBUG_WARMRAIN_OUT_driver_local_ice_precip_flux_warmrain,
            DEBUG_WARMRAIN_OUT_driver_local_rain_warmrain=DEBUG_WARMRAIN_OUT_driver_local_rain_warmrain,
            DEBUG_WARMRAIN_OUT_surface_precip_rain_warmrain=DEBUG_WARMRAIN_OUT_surface_precip_rain_warmrain,
            DEBUG_WARMRAIN_OUT_estimated_inversion_strength_warmrain=DEBUG_WARMRAIN_OUT_estimated_inversion_strength_warmrain,
            DEBUG_WARMRAIN_OUT_driver_local_one_minus_sigma_warmrain=DEBUG_WARMRAIN_OUT_driver_local_one_minus_sigma_warmrain,
            #
            DEBUG_ICECLOUD_IN_driver_local_t_icecloud=DEBUG_ICECLOUD_IN_driver_local_t_icecloud,
            DEBUG_ICECLOUD_IN_driver_local_p_dry_icecloud=DEBUG_ICECLOUD_IN_driver_local_p_dry_icecloud,
            DEBUG_ICECLOUD_IN_driver_local_dp_icecloud=DEBUG_ICECLOUD_IN_driver_local_dp_icecloud,
            DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_vapor_icecloud=DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_vapor_icecloud,
            DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_liquid_icecloud=DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_liquid_icecloud,
            DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_rain_icecloud=DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_rain_icecloud,
            DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_ice_icecloud=DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_ice_icecloud,
            DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_snow_icecloud=DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_snow_icecloud,
            DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_graupel_icecloud=DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_graupel_icecloud,
            DEBUG_ICECLOUD_IN_driver_local_cloud_fraction_icecloud=DEBUG_ICECLOUD_IN_driver_local_cloud_fraction_icecloud,
            DEBUG_ICECLOUD_IN_driver_local_terminal_speed_snow_icecloud=DEBUG_ICECLOUD_IN_driver_local_terminal_speed_snow_icecloud,
            DEBUG_ICECLOUD_IN_driver_local_terminal_speed_graupel_icecloud=DEBUG_ICECLOUD_IN_driver_local_terminal_speed_graupel_icecloud,
            DEBUG_ICECLOUD_IN_driver_local_terminal_speed_rain_icecloud=DEBUG_ICECLOUD_IN_driver_local_terminal_speed_rain_icecloud,
            DEBUG_ICECLOUD_IN_driver_local_density_icecloud=DEBUG_ICECLOUD_IN_driver_local_density_icecloud,
            DEBUG_ICECLOUD_IN_driver_local_density_factor_icecloud=DEBUG_ICECLOUD_IN_driver_local_density_factor_icecloud,
            DEBUG_ICECLOUD_IN_driver_local_rh_limited_icecloud=DEBUG_ICECLOUD_IN_driver_local_rh_limited_icecloud,
            DEBUG_ICECLOUD_IN_non_anvil_large_scale_sublimation_icecloud=DEBUG_ICECLOUD_IN_non_anvil_large_scale_sublimation_icecloud,
            DEBUG_ICECLOUD_IN_driver_local_ccn_icecloud=DEBUG_ICECLOUD_IN_driver_local_ccn_icecloud,
            DEBUG_ICECLOUD_IN_convection_fraction_icecloud=DEBUG_ICECLOUD_IN_convection_fraction_icecloud,
            DEBUG_ICECLOUD_IN_surface_type_icecloud=DEBUG_ICECLOUD_IN_surface_type_icecloud,
            #
            DEBUG_ICECLOUD_OUT_driver_local_t_icecloud=DEBUG_ICECLOUD_OUT_driver_local_t_icecloud,
            DEBUG_ICECLOUD_OUT_driver_local_p_dry_icecloud=DEBUG_ICECLOUD_OUT_driver_local_p_dry_icecloud,
            DEBUG_ICECLOUD_OUT_driver_local_dp_icecloud=DEBUG_ICECLOUD_OUT_driver_local_dp_icecloud,
            DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_vapor_icecloud=DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_vapor_icecloud,
            DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_liquid_icecloud=DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_liquid_icecloud,
            DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_rain_icecloud=DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_rain_icecloud,
            DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_ice_icecloud=DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_ice_icecloud,
            DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_snow_icecloud=DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_snow_icecloud,
            DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_graupel_icecloud=DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_graupel_icecloud,
            DEBUG_ICECLOUD_OUT_driver_local_cloud_fraction_icecloud=DEBUG_ICECLOUD_OUT_driver_local_cloud_fraction_icecloud,
            DEBUG_ICECLOUD_OUT_driver_local_terminal_speed_snow_icecloud=DEBUG_ICECLOUD_OUT_driver_local_terminal_speed_snow_icecloud,
            DEBUG_ICECLOUD_OUT_driver_local_terminal_speed_graupel_icecloud=DEBUG_ICECLOUD_OUT_driver_local_terminal_speed_graupel_icecloud,
            DEBUG_ICECLOUD_OUT_driver_local_terminal_speed_rain_icecloud=DEBUG_ICECLOUD_OUT_driver_local_terminal_speed_rain_icecloud,
            DEBUG_ICECLOUD_OUT_driver_local_density_icecloud=DEBUG_ICECLOUD_OUT_driver_local_density_icecloud,
            DEBUG_ICECLOUD_OUT_driver_local_density_factor_icecloud=DEBUG_ICECLOUD_OUT_driver_local_density_factor_icecloud,
            DEBUG_ICECLOUD_OUT_driver_local_rh_limited_icecloud=DEBUG_ICECLOUD_OUT_driver_local_rh_limited_icecloud,
            DEBUG_ICECLOUD_OUT_non_anvil_large_scale_sublimation_icecloud=DEBUG_ICECLOUD_OUT_non_anvil_large_scale_sublimation_icecloud,
            DEBUG_ICECLOUD_OUT_driver_local_ccn_icecloud=DEBUG_ICECLOUD_OUT_driver_local_ccn_icecloud,
            DEBUG_ICECLOUD_OUT_convection_fraction_icecloud=DEBUG_ICECLOUD_OUT_convection_fraction_icecloud,
            DEBUG_ICECLOUD_OUT_surface_type_icecloud=DEBUG_ICECLOUD_OUT_surface_type_icecloud,
        )

        validation_outputs = {
            "radiation_vapor": copy.copy(state.radiation_field.vapor.field[:]),
            "radiation_liquid": copy.copy(state.radiation_field.liquid.field[:]),
            "radiation_rain": copy.copy(state.radiation_field.rain.field[:]),
            "radiation_ice": copy.copy(state.radiation_field.ice.field[:]),
            "radiation_snow": copy.copy(state.radiation_field.snow.field[:]),
            "radiation_graupel": copy.copy(state.radiation_field.graupel.field[:]),
            "radiation_cloud_fraction": copy.copy(state.radiation_field.cloud_fraction.field[:]),
            "local_total_concentration": copy.copy(locals_.total_concentration.field[:]),
            "local_dvapordt_driver": copy.copy(locals_.driver_tendencies.dvapordt.field[:]),
            "local_dliquiddt_driver": copy.copy(locals_.driver_tendencies.dliquiddt.field[:]),
            "local_draindt_driver": copy.copy(locals_.driver_tendencies.draindt.field[:]),
            "local_dicedt_driver": copy.copy(locals_.driver_tendencies.dicedt.field[:]),
            "local_dsnowdt_driver": copy.copy(locals_.driver_tendencies.dsnowdt.field[:]),
            "local_dgraupeldt_driver": copy.copy(locals_.driver_tendencies.dgraupeldt.field[:]),
            "local_dcloudfractiondt_driver": copy.copy(locals_.driver_tendencies.dcloudfractiondt.field[:]),
            "local_dtdt_driver": copy.copy(locals_.driver_tendencies.dtdt.field[:]),
            "local_dudt_driver": copy.copy(locals_.driver_tendencies.dudt.field[:]),
            "local_dvdt_driver": copy.copy(locals_.driver_tendencies.dvdt.field[:]),
            "t": copy.copy(state.t.field[:]),
            "u": copy.copy(state.u.field[:]),
            "v": copy.copy(state.v.field[:]),
            "w": copy.copy(state.vertical_motion.velocity.field[:]),
            "local_dz": copy.copy(locals_.layer_thickness_negative.field[:]),
            "local_dp": copy.copy(locals_.dp.field[:]),
            "area": copy.copy(state.area.field[:]),
            "land_fraction": copy.copy(state.land_fraction.field[:]),
            "convection_fraction": copy.copy(state.convection_fraction.field[:]),
            "surface_type": copy.copy(state.surface_type.field[:]),
            "estimated_inversion_strength": copy.copy(state.estimated_inversion_strength.field[:]),
            "critical_relative_humidity_for_pdf": copy.copy(
                state.critical_relative_humidity_for_pdf.field[:]
            ),
            "non_anvil_large_scale_evaporation": copy.copy(state.non_anvil_large_scale.evaporation.field[:]),
            "non_anvil_large_scale_sublimation": copy.copy(state.non_anvil_large_scale.sublimation.field[:]),
            "non_anvil_large_scale_liquid_precip_flux": copy.copy(
                state.non_anvil_large_scale.liquid_precip_flux.field
            )[:],
            "non_anvil_large_scale_ice_precip_flux": copy.copy(
                state.non_anvil_large_scale.ice_precip_flux.field[:]
            ),
            "surface_precip_rain": copy.copy(state.precipitation_at_surface.rain.field[:]),
            "surface_precip_snow": copy.copy(state.precipitation_at_surface.snow.field[:]),
            "surface_precip_ice": copy.copy(state.precipitation_at_surface.ice.field[:]),
            "surface_precip_graupel": copy.copy(state.precipitation_at_surface.graupel.field[:]),
            #
            "DEBUG_vapor_unmodified_setup": DEBUG_vapor_unmodified_setup.field[:],
            "DEBUG_vapor_modified_setup": DEBUG_vapor_modified_setup.field[:],
            "DEBUG_vapor_unmodified_fallspeed": DEBUG_vapor_unmodified_fallspeed.field[:],
            "DEBUG_vapor_modified_fallspeed": DEBUG_vapor_modified_fallspeed.field[:],
            "DEBUG_vapor_unmodified_terminalfall": DEBUG_vapor_unmodified_terminalfall.field[:],
            "DEBUG_vapor_modified_terminalfall": DEBUG_vapor_modified_terminalfall.field[:],
            "DEBUG_vapor_unmodified_warmrain": DEBUG_vapor_unmodified_warmrain.field[:],
            "DEBUG_vapor_modified_warmrain": DEBUG_vapor_modified_warmrain.field[:],
            "DEBUG_vapor_unmodified_icecloud": DEBUG_vapor_unmodified_icecloud.field[:],
            "DEBUG_vapor_modified_icecloud": DEBUG_vapor_modified_icecloud.field[:],
            #
            "DEBUG_TERMINALFALL_IN_driver_local_t_terminalfall": DEBUG_TERMINALFALL_IN_driver_local_t_terminalfall.field[:],
            "DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_vapor_terminalfall": DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_vapor_terminalfall.field[:],
            "DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_liquid_terminalfall": DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_liquid_terminalfall.field[:],
            "DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_rain_terminalfall": DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_rain_terminalfall.field[:],
            "DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_graupel_terminalfall": DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_graupel_terminalfall.field[:],
            "DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_snow_terminalfall": DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_snow_terminalfall.field[:],
            "DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_ice_terminalfall": DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_ice_terminalfall.field[:],
            "DEBUG_TERMINALFALL_IN_driver_local_ice_precip_flux_terminalfall": DEBUG_TERMINALFALL_IN_driver_local_ice_precip_flux_terminalfall.field[:],
            "DEBUG_TERMINALFALL_IN_driver_local_w_terminalfall": DEBUG_TERMINALFALL_IN_driver_local_w_terminalfall.field[:],
            "DEBUG_TERMINALFALL_IN_driver_local_dz_terminalfall": DEBUG_TERMINALFALL_IN_driver_local_dz_terminalfall.field[:],
            "DEBUG_TERMINALFALL_IN_driver_local_dp_terminalfall": DEBUG_TERMINALFALL_IN_driver_local_dp_terminalfall.field[:],
            "DEBUG_TERMINALFALL_IN_driver_local_terminal_speed_ice_terminalfall": DEBUG_TERMINALFALL_IN_driver_local_terminal_speed_ice_terminalfall.field[:],
            "DEBUG_TERMINALFALL_IN_driver_local_terminal_speed_snow_terminalfall": DEBUG_TERMINALFALL_IN_driver_local_terminal_speed_snow_terminalfall.field[:],
            "DEBUG_TERMINALFALL_IN_driver_local_terminal_speed_graupel_terminalfall": DEBUG_TERMINALFALL_IN_driver_local_terminal_speed_graupel_terminalfall.field[:],
            "DEBUG_TERMINALFALL_IN_surface_precip_rain_terminalfall": DEBUG_TERMINALFALL_IN_surface_precip_rain_terminalfall.field[:],
            "DEBUG_TERMINALFALL_IN_surface_precip_snow_terminalfall": DEBUG_TERMINALFALL_IN_surface_precip_snow_terminalfall.field[:],
            "DEBUG_TERMINALFALL_IN_surface_precip_graupel_terminalfall": DEBUG_TERMINALFALL_IN_surface_precip_graupel_terminalfall.field[:],
            "DEBUG_TERMINALFALL_IN_surface_precip_ice_terminalfall": DEBUG_TERMINALFALL_IN_surface_precip_ice_terminalfall.field[:],
            #
            "DEBUG_TERMINALFALL_OUT_driver_local_t_terminalfall": DEBUG_TERMINALFALL_OUT_driver_local_t_terminalfall.field[:],
            "DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_vapor_terminalfall": DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_vapor_terminalfall.field[:],
            "DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_liquid_terminalfall": DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_liquid_terminalfall.field[:],
            "DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_rain_terminalfall": DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_rain_terminalfall.field[:],
            "DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_graupel_terminalfall": DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_graupel_terminalfall.field[:],
            "DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_snow_terminalfall": DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_snow_terminalfall.field[:],
            "DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_ice_terminalfall": DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_ice_terminalfall.field[:],
            "DEBUG_TERMINALFALL_OUT_driver_local_ice_precip_flux_terminalfall": DEBUG_TERMINALFALL_OUT_driver_local_ice_precip_flux_terminalfall.field[:],
            "DEBUG_TERMINALFALL_OUT_driver_local_w_terminalfall": DEBUG_TERMINALFALL_OUT_driver_local_w_terminalfall.field[:],
            "DEBUG_TERMINALFALL_OUT_driver_local_dz_terminalfall": DEBUG_TERMINALFALL_OUT_driver_local_dz_terminalfall.field[:],
            "DEBUG_TERMINALFALL_OUT_driver_local_dp_terminalfall": DEBUG_TERMINALFALL_OUT_driver_local_dp_terminalfall.field[:],
            "DEBUG_TERMINALFALL_OUT_driver_local_terminal_speed_ice_terminalfall": DEBUG_TERMINALFALL_OUT_driver_local_terminal_speed_ice_terminalfall.field[:],
            "DEBUG_TERMINALFALL_OUT_driver_local_terminal_speed_snow_terminalfall": DEBUG_TERMINALFALL_OUT_driver_local_terminal_speed_snow_terminalfall.field[:],
            "DEBUG_TERMINALFALL_OUT_driver_local_terminal_speed_graupel_terminalfall": DEBUG_TERMINALFALL_OUT_driver_local_terminal_speed_graupel_terminalfall.field[:],
            "DEBUG_TERMINALFALL_OUT_surface_precip_rain_terminalfall": DEBUG_TERMINALFALL_OUT_surface_precip_rain_terminalfall.field[:],
            "DEBUG_TERMINALFALL_OUT_surface_precip_snow_terminalfall": DEBUG_TERMINALFALL_OUT_surface_precip_snow_terminalfall.field[:],
            "DEBUG_TERMINALFALL_OUT_surface_precip_graupel_terminalfall": DEBUG_TERMINALFALL_OUT_surface_precip_graupel_terminalfall.field[:],
            "DEBUG_TERMINALFALL_OUT_surface_precip_ice_terminalfall": DEBUG_TERMINALFALL_OUT_surface_precip_ice_terminalfall.field[:],
            #
            "DEBUG_WARMRAIN_IN_driver_local_dp_warmrain": DEBUG_WARMRAIN_IN_driver_local_dp_warmrain.field[:],
            "DEBUG_WARMRAIN_IN_driver_local_dz_warmrain": DEBUG_WARMRAIN_IN_driver_local_dz_warmrain.field[:],
            "DEBUG_WARMRAIN_IN_driver_local_t_warmrain": DEBUG_WARMRAIN_IN_driver_local_t_warmrain.field[:],
            "DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_vapor_warmrain": DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_vapor_warmrain.field[:],
            "DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_liquid_warmrain": DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_liquid_warmrain.field[:],
            "DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_rain_warmrain": DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_rain_warmrain.field[:],
            "DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_ice_warmrain": DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_ice_warmrain.field[:],
            "DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_snow_warmrain": DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_snow_warmrain.field[:],
            "DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_graupel_warmrain": DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_graupel_warmrain.field[:],
            "DEBUG_WARMRAIN_IN_driver_local_cloud_fraction_warmrain": DEBUG_WARMRAIN_IN_driver_local_cloud_fraction_warmrain.field[:],
            "DEBUG_WARMRAIN_IN_driver_local_ccn_warmrain": DEBUG_WARMRAIN_IN_driver_local_ccn_warmrain.field[:],
            "DEBUG_WARMRAIN_IN_driver_local_density_warmrain": DEBUG_WARMRAIN_IN_driver_local_density_warmrain.field[:],
            "DEBUG_WARMRAIN_IN_driver_local_density_factor_warmrain": DEBUG_WARMRAIN_IN_driver_local_density_factor_warmrain.field[:],
            "DEBUG_WARMRAIN_IN_driver_local_c_praut_warmrain": DEBUG_WARMRAIN_IN_driver_local_c_praut_warmrain.field[:],
            "DEBUG_WARMRAIN_IN_driver_local_terminal_speed_rain_warmrain": DEBUG_WARMRAIN_IN_driver_local_terminal_speed_rain_warmrain.field[:],
            "DEBUG_WARMRAIN_IN_driver_local_evaporation_warmrain": DEBUG_WARMRAIN_IN_driver_local_evaporation_warmrain.field[:],
            "DEBUG_WARMRAIN_IN_driver_local_liquid_precip_flux_warmrain": DEBUG_WARMRAIN_IN_driver_local_liquid_precip_flux_warmrain.field[:],
            "DEBUG_WARMRAIN_IN_driver_local_w_warmrain": DEBUG_WARMRAIN_IN_driver_local_w_warmrain.field[:],
            "DEBUG_WARMRAIN_IN_driver_local_rh_limited_warmrain": DEBUG_WARMRAIN_IN_driver_local_rh_limited_warmrain.field[:],
            "DEBUG_WARMRAIN_IN_non_anvil_large_scale_evaporation_warmrain": DEBUG_WARMRAIN_IN_non_anvil_large_scale_evaporation_warmrain.field[:],
            "DEBUG_WARMRAIN_IN_non_anvil_large_scale_liquid_precip_flux_warmrain": DEBUG_WARMRAIN_IN_non_anvil_large_scale_liquid_precip_flux_warmrain.field[:,:,1:],
            "DEBUG_WARMRAIN_IN_non_anvil_large_scale_ice_precip_flux_warmrain": DEBUG_WARMRAIN_IN_non_anvil_large_scale_ice_precip_flux_warmrain.field[:,:,1:],
            "DEBUG_WARMRAIN_IN_driver_local_mass_warmrain": DEBUG_WARMRAIN_IN_driver_local_mass_warmrain.field[:],
            "DEBUG_WARMRAIN_IN_driver_local_ice_precip_flux_warmrain": DEBUG_WARMRAIN_IN_driver_local_ice_precip_flux_warmrain.field[:],
            "DEBUG_WARMRAIN_IN_driver_local_rain_warmrain": DEBUG_WARMRAIN_IN_driver_local_rain_warmrain.field[:],
            "DEBUG_WARMRAIN_IN_surface_precip_rain_warmrain": DEBUG_WARMRAIN_IN_surface_precip_rain_warmrain.field[:],
            "DEBUG_WARMRAIN_IN_estimated_inversion_strength_warmrain": DEBUG_WARMRAIN_IN_estimated_inversion_strength_warmrain.field[:],
            "DEBUG_WARMRAIN_IN_driver_local_one_minus_sigma_warmrain": DEBUG_WARMRAIN_IN_driver_local_one_minus_sigma_warmrain.field[:],
            #
            "DEBUG_WARMRAIN_OUT_driver_local_dp_warmrain": DEBUG_WARMRAIN_OUT_driver_local_dp_warmrain.field[:],
            "DEBUG_WARMRAIN_OUT_driver_local_dz_warmrain": DEBUG_WARMRAIN_OUT_driver_local_dz_warmrain.field[:],
            "DEBUG_WARMRAIN_OUT_driver_local_t_warmrain": DEBUG_WARMRAIN_OUT_driver_local_t_warmrain.field[:],
            "DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_vapor_warmrain": DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_vapor_warmrain.field[:],
            "DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_liquid_warmrain": DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_liquid_warmrain.field[:],
            "DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_rain_warmrain": DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_rain_warmrain.field[:],
            "DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_ice_warmrain": DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_ice_warmrain.field[:],
            "DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_snow_warmrain": DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_snow_warmrain.field[:],
            "DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_graupel_warmrain": DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_graupel_warmrain.field[:],
            "DEBUG_WARMRAIN_OUT_driver_local_cloud_fraction_warmrain": DEBUG_WARMRAIN_OUT_driver_local_cloud_fraction_warmrain.field[:],
            "DEBUG_WARMRAIN_OUT_driver_local_ccn_warmrain": DEBUG_WARMRAIN_OUT_driver_local_ccn_warmrain.field[:],
            "DEBUG_WARMRAIN_OUT_driver_local_density_warmrain": DEBUG_WARMRAIN_OUT_driver_local_density_warmrain.field[:],
            "DEBUG_WARMRAIN_OUT_driver_local_density_factor_warmrain": DEBUG_WARMRAIN_OUT_driver_local_density_factor_warmrain.field[:],
            "DEBUG_WARMRAIN_OUT_driver_local_c_praut_warmrain": DEBUG_WARMRAIN_OUT_driver_local_c_praut_warmrain.field[:],
            "DEBUG_WARMRAIN_OUT_driver_local_terminal_speed_rain_warmrain": DEBUG_WARMRAIN_OUT_driver_local_terminal_speed_rain_warmrain.field[:],
            "DEBUG_WARMRAIN_OUT_driver_local_evaporation_warmrain": DEBUG_WARMRAIN_OUT_driver_local_evaporation_warmrain.field[:],
            "DEBUG_WARMRAIN_OUT_driver_local_liquid_precip_flux_warmrain": DEBUG_WARMRAIN_OUT_driver_local_liquid_precip_flux_warmrain.field[:],
            "DEBUG_WARMRAIN_OUT_driver_local_w_warmrain": DEBUG_WARMRAIN_OUT_driver_local_w_warmrain.field[:],
            "DEBUG_WARMRAIN_OUT_driver_local_rh_limited_warmrain": DEBUG_WARMRAIN_OUT_driver_local_rh_limited_warmrain.field[:],
            "DEBUG_WARMRAIN_OUT_non_anvil_large_scale_evaporation_warmrain": DEBUG_WARMRAIN_OUT_non_anvil_large_scale_evaporation_warmrain.field[:],
            "DEBUG_WARMRAIN_OUT_non_anvil_large_scale_liquid_precip_flux_warmrain": DEBUG_WARMRAIN_OUT_non_anvil_large_scale_liquid_precip_flux_warmrain.field[:,:,1:],
            "DEBUG_WARMRAIN_OUT_non_anvil_large_scale_ice_precip_flux_warmrain": DEBUG_WARMRAIN_OUT_non_anvil_large_scale_ice_precip_flux_warmrain.field[:,:,1:],
            "DEBUG_WARMRAIN_OUT_driver_local_mass_warmrain": DEBUG_WARMRAIN_OUT_driver_local_mass_warmrain.field[:],
            "DEBUG_WARMRAIN_OUT_driver_local_ice_precip_flux_warmrain": DEBUG_WARMRAIN_OUT_driver_local_ice_precip_flux_warmrain.field[:],
            "DEBUG_WARMRAIN_OUT_driver_local_rain_warmrain": DEBUG_WARMRAIN_OUT_driver_local_rain_warmrain.field[:],
            "DEBUG_WARMRAIN_OUT_surface_precip_rain_warmrain": DEBUG_WARMRAIN_OUT_surface_precip_rain_warmrain.field[:],
            "DEBUG_WARMRAIN_OUT_estimated_inversion_strength_warmrain": DEBUG_WARMRAIN_OUT_estimated_inversion_strength_warmrain.field[:],
            "DEBUG_WARMRAIN_OUT_driver_local_one_minus_sigma_warmrain": DEBUG_WARMRAIN_OUT_driver_local_one_minus_sigma_warmrain.field[:],
            #
            "DEBUG_ICECLOUD_IN_driver_local_t_icecloud": DEBUG_ICECLOUD_IN_driver_local_t_icecloud.field[:],
            "DEBUG_ICECLOUD_IN_driver_local_p_dry_icecloud": DEBUG_ICECLOUD_IN_driver_local_p_dry_icecloud.field[:],
            "DEBUG_ICECLOUD_IN_driver_local_dp_icecloud": DEBUG_ICECLOUD_IN_driver_local_dp_icecloud.field[:],
            "DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_vapor_icecloud": DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_vapor_icecloud.field[:],
            "DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_liquid_icecloud": DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_liquid_icecloud.field[:],
            "DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_rain_icecloud": DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_rain_icecloud.field[:],
            "DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_ice_icecloud": DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_ice_icecloud.field[:],
            "DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_snow_icecloud": DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_snow_icecloud.field[:],
            "DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_graupel_icecloud": DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_graupel_icecloud.field[:],
            "DEBUG_ICECLOUD_IN_driver_local_cloud_fraction_icecloud": DEBUG_ICECLOUD_IN_driver_local_cloud_fraction_icecloud.field[:],
            "DEBUG_ICECLOUD_IN_driver_local_terminal_speed_snow_icecloud": DEBUG_ICECLOUD_IN_driver_local_terminal_speed_snow_icecloud.field[:],
            "DEBUG_ICECLOUD_IN_driver_local_terminal_speed_graupel_icecloud": DEBUG_ICECLOUD_IN_driver_local_terminal_speed_graupel_icecloud.field[:],
            "DEBUG_ICECLOUD_IN_driver_local_terminal_speed_rain_icecloud": DEBUG_ICECLOUD_IN_driver_local_terminal_speed_rain_icecloud.field[:],
            "DEBUG_ICECLOUD_IN_driver_local_density_icecloud": DEBUG_ICECLOUD_IN_driver_local_density_icecloud.field[:],
            "DEBUG_ICECLOUD_IN_driver_local_density_factor_icecloud": DEBUG_ICECLOUD_IN_driver_local_density_factor_icecloud.field[:],
            "DEBUG_ICECLOUD_IN_driver_local_rh_limited_icecloud": DEBUG_ICECLOUD_IN_driver_local_rh_limited_icecloud.field[:],
            "DEBUG_ICECLOUD_IN_non_anvil_large_scale_sublimation_icecloud": DEBUG_ICECLOUD_IN_non_anvil_large_scale_sublimation_icecloud.field[:],
            "DEBUG_ICECLOUD_IN_driver_local_ccn_icecloud": DEBUG_ICECLOUD_IN_driver_local_ccn_icecloud.field[:],
            "DEBUG_ICECLOUD_IN_convection_fraction_icecloud": DEBUG_ICECLOUD_IN_convection_fraction_icecloud.field[:],
            "DEBUG_ICECLOUD_IN_surface_type_icecloud": DEBUG_ICECLOUD_IN_surface_type_icecloud.field[:],
            #
            "DEBUG_ICECLOUD_OUT_driver_local_t_icecloud": DEBUG_ICECLOUD_OUT_driver_local_t_icecloud.field[:],
            "DEBUG_ICECLOUD_OUT_driver_local_p_dry_icecloud": DEBUG_ICECLOUD_OUT_driver_local_p_dry_icecloud.field[:],
            "DEBUG_ICECLOUD_OUT_driver_local_dp_icecloud": DEBUG_ICECLOUD_OUT_driver_local_dp_icecloud.field[:],
            "DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_vapor_icecloud": DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_vapor_icecloud.field[:],
            "DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_liquid_icecloud": DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_liquid_icecloud.field[:],
            "DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_rain_icecloud": DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_rain_icecloud.field[:],
            "DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_ice_icecloud": DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_ice_icecloud.field[:],
            "DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_snow_icecloud": DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_snow_icecloud.field[:],
            "DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_graupel_icecloud": DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_graupel_icecloud.field[:],
            "DEBUG_ICECLOUD_OUT_driver_local_cloud_fraction_icecloud": DEBUG_ICECLOUD_OUT_driver_local_cloud_fraction_icecloud.field[:],
            "DEBUG_ICECLOUD_OUT_driver_local_terminal_speed_snow_icecloud": DEBUG_ICECLOUD_OUT_driver_local_terminal_speed_snow_icecloud.field[:],
            "DEBUG_ICECLOUD_OUT_driver_local_terminal_speed_graupel_icecloud": DEBUG_ICECLOUD_OUT_driver_local_terminal_speed_graupel_icecloud.field[:],
            "DEBUG_ICECLOUD_OUT_driver_local_terminal_speed_rain_icecloud": DEBUG_ICECLOUD_OUT_driver_local_terminal_speed_rain_icecloud.field[:],
            "DEBUG_ICECLOUD_OUT_driver_local_density_icecloud": DEBUG_ICECLOUD_OUT_driver_local_density_icecloud.field[:],
            "DEBUG_ICECLOUD_OUT_driver_local_density_factor_icecloud": DEBUG_ICECLOUD_OUT_driver_local_density_factor_icecloud.field[:],
            "DEBUG_ICECLOUD_OUT_driver_local_rh_limited_icecloud": DEBUG_ICECLOUD_OUT_driver_local_rh_limited_icecloud.field[:],
            "DEBUG_ICECLOUD_OUT_non_anvil_large_scale_sublimation_icecloud": DEBUG_ICECLOUD_OUT_non_anvil_large_scale_sublimation_icecloud.field[:],
            "DEBUG_ICECLOUD_OUT_driver_local_ccn_icecloud": DEBUG_ICECLOUD_OUT_driver_local_ccn_icecloud.field[:],
            "DEBUG_ICECLOUD_OUT_convection_fraction_icecloud": DEBUG_ICECLOUD_OUT_convection_fraction_icecloud.field[:],
            "DEBUG_ICECLOUD_OUT_surface_type_icecloud": DEBUG_ICECLOUD_OUT_surface_type_icecloud.field[:],
        }

        # Micro-bench
        ts = 0
        if ts > 0:
            s = time.perf_counter()
            for _ in range(ts):
                code(
                    t=state.t,
                    u=state.u,
                    v=state.v,
                    w=state.vertical_motion.velocity,
                    dz=locals_.layer_thickness_negative,
                    dp=locals_.dp,
                    area=state.area,
                    land_fraction=state.land_fraction,
                    convection_fraction=state.convection_fraction,
                    surface_type=state.surface_type,
                    estimated_inversion_strength=state.estimated_inversion_strength,
                    critical_relative_humidity_for_pdf=state.critical_relative_humidity_for_pdf,
                    vapor=state.radiation_field.vapor,
                    liquid=state.radiation_field.liquid,
                    rain=state.radiation_field.rain,
                    ice=state.radiation_field.ice,
                    snow=state.radiation_field.snow,
                    graupel=state.radiation_field.graupel,
                    cloud_fraction=state.radiation_field.cloud_fraction,
                    total_concentration=locals_.total_concentration,
                    dvapordt=locals_.driver_tendencies.dvapordt,
                    dliquiddt=locals_.driver_tendencies.dliquiddt,
                    draindt=locals_.driver_tendencies.draindt,
                    dicedt=locals_.driver_tendencies.dicedt,
                    dsnowdt=locals_.driver_tendencies.dsnowdt,
                    dgraupeldt=locals_.driver_tendencies.dgraupeldt,
                    dcloudfractiondt=locals_.driver_tendencies.dcloudfractiondt,
                    dtdt=locals_.driver_tendencies.dtdt,
                    dudt=locals_.driver_tendencies.dudt,
                    dvdt=locals_.driver_tendencies.dvdt,
                    liquid_precip_flux=state.non_anvil_large_scale.liquid_precip_flux,
                    ice_precip_flux=state.non_anvil_large_scale.ice_precip_flux,
                    evaporation=state.non_anvil_large_scale.evaporation,
                    sublimation=state.non_anvil_large_scale.sublimation,
                    surface_precip_rain=state.precipitation_at_surface.rain,
                    surface_precip_snow=state.precipitation_at_surface.snow,
                    surface_precip_ice=state.precipitation_at_surface.ice,
                    surface_precip_graupel=state.precipitation_at_surface.graupel,
                )
            e = time.perf_counter()
            ndsl_log.info(f"GFDL1M Driver micro bench: {(e - s) / ts:.4f}s")

        return validation_outputs
