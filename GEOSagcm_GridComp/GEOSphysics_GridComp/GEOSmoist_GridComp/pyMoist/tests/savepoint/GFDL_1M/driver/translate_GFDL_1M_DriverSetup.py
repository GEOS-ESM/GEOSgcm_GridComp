from f90nml import Namelist

from ndsl import StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.GFDL_1M.state import GFDL1MState
from pyMoist.GFDL_1M.locals import GFDL1MLocals
from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.GFDL_1M.driver.config_constants import GFDL1MDriverConfigDependentConstants
from pyMoist.GFDL_1M.driver.locals import GFDL1MDriverLocals
from pyMoist.GFDL_1M.driver.setup import GFDL1MDriverSetup
from ndsl.stencils.testing.savepoint import DataLoader
import numpy as np


class TranslateGFDL_1M_DriverSetup(TranslateFortranData2Py):
    def __init__(self, grid: Grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "t_driversetup": {},
            "local_dp_driversetup": {},
            "critical_relative_humidity_for_pdf_driversetup": {},
            "radiation_vapor_driversetup": {},
            "radiation_liquid_driversetup": {},
            "radiation_ice_driversetup": {},
            "radiation_rain_driversetup": {},
            "radiation_snow_driversetup": {},
            "radiation_graupel_driversetup": {},
            "radiation_cloud_fraction_driversetup": {},
            "total_concentration_driversetup": {},
            "driver_local_dry_mixing_ratio_vapor_unmodified_driversetup": {},
            "driver_local_dry_mixing_ratio_liquid_unmodified_driversetup": {},
            "driver_local_dry_mixing_ratio_rain_unmodified_driversetup": {},
            "driver_local_dry_mixing_ratio_ice_unmodified_driversetup": {},
            "driver_local_dry_mixing_ratio_snow_unmodified_driversetup": {},
            "driver_local_dry_mixing_ratio_graupel_unmodified_driversetup": {},
            "driver_local_dry_mixing_ratio_vapor_driversetup": {},
            "driver_local_dry_mixing_ratio_liquid_driversetup": {},
            "driver_local_dry_mixing_ratio_rain_driversetup": {},
            "driver_local_dry_mixing_ratio_ice_driversetup": {},
            "driver_local_dry_mixing_ratio_snow_driversetup": {},
            "driver_local_dry_mixing_ratio_graupel_driversetup": {},
            "driver_local_cloud_fraction_driversetup": {},
            "local_dz_drivesetup": {},
            "u_driversetup": {},
            "v_driversetup": {},
            "w_driversetup": {},
            "driver_local_t_driversetup": {},
            "driver_local_dp_driversetup": {},
            "driver_local_density_unmodified_driversetup": {},
            "driver_local_p_dry_driversetup": {},
            "driver_local_mass_driversetup": {},
            "driver_local_u_driversetup": {},
            "driver_local_v_driversetup": {},
            "driver_local_w_driversetup": {},
            "driver_local_ccn_driversetup": {},
            "driver_local_c_praut_driversetup": {},
            "driver_local_rh_limited_driversetup": {},
            "non_anvil_large_scale_liquid_precip_flux_driversetup": {},
            "non_anvil_large_scale_ice_precip_flux_driversetup": {},
            "non_anvil_large_scale_evaporation_driversetup": {},
            "non_anvil_large_scale_sublimation_driversetup": {},
            "driver_local_one_minus_sigma_driversetup": {},
            "area_driversetup": {},
            "surface_precip_rain_driversetup": {},
            "surface_precip_snow_driversetup": {},
            "surface_precip_graupel_driversetup": {},
            "surface_precip_ice_driversetup": {},
        }

        self.out_vars = self.in_vars["data_vars"].copy()

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GFDL_1M-constants")

    def compute(self, inputs):
        # initalize dataclasses
        state = GFDL1MState.zeros(self.quantity_factory)
        locals = GFDL1MLocals.zeros(self.quantity_factory)
        driver_locals = GFDL1MDriverLocals.zeros(self.quantity_factory)

        # initalize constants
        config = GFDL1MConfig(**self.constants)
        config_dependent_constants = GFDL1MDriverConfigDependentConstants.make(config)

        state.t.field[:] = inputs["t_driversetup"][:, :, :, 0]
        locals.dp.field[:] = inputs["local_dp_driversetup"][:, :, :, 0]
        state.critical_relative_humidity_for_pdf.field[:] = inputs[
            "critical_relative_humidity_for_pdf_driversetup"
        ][:, :, :, 0]
        state.radiation_field.vapor.field[:] = inputs["radiation_vapor_driversetup"][:, :, :, 0]
        state.radiation_field.liquid.field[:] = inputs["radiation_liquid_driversetup"][:, :, :, 0]
        state.radiation_field.ice.field[:] = inputs["radiation_ice_driversetup"][:, :, :, 0]
        state.radiation_field.rain.field[:] = inputs["radiation_rain_driversetup"][:, :, :, 0]
        state.radiation_field.snow.field[:] = inputs["radiation_snow_driversetup"][:, :, :, 0]
        state.radiation_field.graupel.field[:] = inputs["radiation_graupel_driversetup"][:, :, :, 0]
        state.radiation_field.cloud_fraction.field[:] = inputs["radiation_cloud_fraction_driversetup"][
            :, :, :, 0
        ]
        locals.total_concentration.field[:] = inputs["total_concentration_driversetup"][:, :, :, 0]
        driver_locals.unmodified.mixing_ratio.vapor.field[:] = inputs[
            "driver_local_dry_mixing_ratio_vapor_unmodified_driversetup"
        ][:, :, :, 0]
        driver_locals.unmodified.mixing_ratio.liquid.field[:] = inputs[
            "driver_local_dry_mixing_ratio_liquid_unmodified_driversetup"
        ][:, :, :, 0]
        driver_locals.unmodified.mixing_ratio.rain.field[:] = inputs[
            "driver_local_dry_mixing_ratio_rain_unmodified_driversetup"
        ][:, :, :, 0]
        driver_locals.unmodified.mixing_ratio.ice.field[:] = inputs[
            "driver_local_dry_mixing_ratio_ice_unmodified_driversetup"
        ][:, :, :, 0]
        driver_locals.unmodified.mixing_ratio.snow.field[:] = inputs[
            "driver_local_dry_mixing_ratio_snow_unmodified_driversetup"
        ][:, :, :, 0]
        driver_locals.unmodified.mixing_ratio.graupel.field[:] = inputs[
            "driver_local_dry_mixing_ratio_graupel_unmodified_driversetup"
        ][:, :, :, 0]
        driver_locals.dry_air_mixing_ratio.vapor.field[:] = inputs[
            "driver_local_dry_mixing_ratio_vapor_driversetup"
        ][:, :, :, 0]
        driver_locals.dry_air_mixing_ratio.liquid.field[:] = inputs[
            "driver_local_dry_mixing_ratio_liquid_driversetup"
        ][:, :, :, 0]
        driver_locals.dry_air_mixing_ratio.rain.field[:] = inputs[
            "driver_local_dry_mixing_ratio_rain_driversetup"
        ][:, :, :, 0]
        driver_locals.dry_air_mixing_ratio.ice.field[:] = inputs[
            "driver_local_dry_mixing_ratio_ice_driversetup"
        ][:, :, :, 0]
        driver_locals.dry_air_mixing_ratio.snow.field[:] = inputs[
            "driver_local_dry_mixing_ratio_snow_driversetup"
        ][:, :, :, 0]
        driver_locals.dry_air_mixing_ratio.graupel.field[:] = inputs[
            "driver_local_dry_mixing_ratio_graupel_driversetup"
        ][:, :, :, 0]
        driver_locals.cloud_fraction.field[:] = inputs["driver_local_cloud_fraction_driversetup"][:, :, :, 0]
        locals.dz.field[:] = inputs["local_dz_drivesetup"][:, :, :, 0]
        state.u.field[:] = inputs["u_driversetup"][:, :, :, 0]
        state.v.field[:] = inputs["v_driversetup"][:, :, :, 0]
        state.vertical_motion.velocity.field[:] = inputs["w_driversetup"][:, :, :, 0]
        driver_locals.t.field[:] = inputs["driver_local_t_driversetup"][:, :, :, 0]
        driver_locals.dp.field[:] = inputs["driver_local_dp_driversetup"][:, :, :, 0]
        driver_locals.density_unmodified.field[:] = inputs["driver_local_density_unmodified_driversetup"][
            :, :, :, 0
        ]
        driver_locals.p_dry.field[:] = inputs["driver_local_p_dry_driversetup"][:, :, :, 0]
        driver_locals.mass.field[:] = inputs["driver_local_mass_driversetup"][:, :, :, 0]
        driver_locals.u.field[:] = inputs["driver_local_u_driversetup"][:, :, :, 0]
        driver_locals.v.field[:] = inputs["driver_local_v_driversetup"][:, :, :, 0]
        driver_locals.w.field[:] = inputs["driver_local_w_driversetup"][:, :, :, 0]
        driver_locals.ccn.field[:] = inputs["driver_local_ccn_driversetup"][:, :, :, 0]
        driver_locals.c_praut.field[:] = inputs["driver_local_c_praut_driversetup"][:, :, :, 0]
        driver_locals.rh_limited.field[:] = inputs["driver_local_rh_limited_driversetup"][:, :, :, 0]
        state.non_anvil_large_scale.liquid_precip_flux.field[:, :, 0:-1] = inputs[
            "non_anvil_large_scale_liquid_precip_flux_driversetup"
        ][:, :, :, 0]
        state.non_anvil_large_scale.ice_precip_flux.field[:, :, 0:-1] = inputs[
            "non_anvil_large_scale_ice_precip_flux_driversetup"
        ][:, :, :, 0]
        state.non_anvil_large_scale.evaporation.field[:] = inputs[
            "non_anvil_large_scale_evaporation_driversetup"
        ][:, :, :, 0]
        state.non_anvil_large_scale.sublimation.field[:] = inputs[
            "non_anvil_large_scale_sublimation_driversetup"
        ][:, :, :, 0]
        driver_locals.one_minus_sigma.field[:] = inputs["driver_local_one_minus_sigma_driversetup"][
            :, :, 0, 0
        ]
        state.area.field[:] = inputs["area_driversetup"][:, :, 0, 0]
        state.precipitation_at_surface.rain.field[:] = inputs["surface_precip_rain_driversetup"][:, :, 0, 0]
        state.precipitation_at_surface.snow.field[:] = inputs["surface_precip_snow_driversetup"][:, :, 0, 0]
        state.precipitation_at_surface.graupel.field[:] = inputs["surface_precip_graupel_driversetup"][
            :, :, 0, 0
        ]
        state.precipitation_at_surface.ice.field[:] = inputs["surface_precip_ice_driversetup"][:, :, 0, 0]

        # construct test stencil
        code = GFDL1MDriverSetup(
            stencil_factory=self.stencil_factory,
            config=config,
            config_dependent_constants=config_dependent_constants,
        )

        code(
            unmodified_t=state.t,
            t=driver_locals.t,
            unmodified_dp=locals.dp,
            dp=driver_locals.dp,
            critical_relative_humidity_for_pdf=state.critical_relative_humidity_for_pdf,
            radiation_field_vapor=state.radiation_field.vapor,
            radiation_field_liquid=state.radiation_field.liquid,
            radiation_field_ice=state.radiation_field.ice,
            radiation_field_rain=state.radiation_field.rain,
            radiation_field_snow=state.radiation_field.snow,
            radiation_field_graupel=state.radiation_field.graupel,
            radiation_field_cloud_fraction=state.radiation_field.cloud_fraction,
            total_concentration=locals.total_concentration,
            unmodified_mixing_ratio_vapor=driver_locals.unmodified.mixing_ratio.vapor,
            unmodified_mixing_ratio_liquid=driver_locals.unmodified.mixing_ratio.liquid,
            unmodified_mixing_ratio_rain=driver_locals.unmodified.mixing_ratio.rain,
            unmodified_mixing_ratio_ice=driver_locals.unmodified.mixing_ratio.ice,
            unmodified_mixing_ratio_snow=driver_locals.unmodified.mixing_ratio.snow,
            unmodified_mixing_ratio_graupel=driver_locals.unmodified.mixing_ratio.graupel,
            dry_air_mixing_ratio_vapor=driver_locals.dry_air_mixing_ratio.vapor,
            dry_air_mixing_ratio_liquid=driver_locals.dry_air_mixing_ratio.liquid,
            dry_air_mixing_ratio_rain=driver_locals.dry_air_mixing_ratio.rain,
            dry_air_mixing_ratio_ice=driver_locals.dry_air_mixing_ratio.ice,
            dry_air_mixing_ratio_snow=driver_locals.dry_air_mixing_ratio.snow,
            dry_air_mixing_ratio_graupel=driver_locals.dry_air_mixing_ratio.graupel,
            cloud_fraction=driver_locals.cloud_fraction,
            dz=locals.dz,
            u_unmodified=state.u,
            u=driver_locals.u,
            v_unmodified=state.v,
            v=driver_locals.v,
            w_unmodified=state.vertical_motion.velocity,
            w=driver_locals.w,
            area=state.area,
            density_unmodified=driver_locals.density_unmodified,
            p_dry=driver_locals.p_dry,
            mass=driver_locals.mass,
            one_minus_sigma=driver_locals.one_minus_sigma,
            ccn=driver_locals.ccn,
            c_praut=driver_locals.c_praut,
            rh_limited=driver_locals.rh_limited,
            rain=state.precipitation_at_surface.rain,
            snow=state.precipitation_at_surface.snow,
            graupel=state.precipitation_at_surface.graupel,
            ice=state.precipitation_at_surface.ice,
            liquid_precip_flux=state.non_anvil_large_scale.liquid_precip_flux,
            ice_precip_flux=state.non_anvil_large_scale.ice_precip_flux,
            evaporation=state.non_anvil_large_scale.evaporation,
            sublimation=state.non_anvil_large_scale.sublimation,
        )

        # get the shape of the field
        nx, ny, nz, ntimes = inputs["t_driversetup"].shape

        # prefill output array with nans
        outputs = {}
        for key in self.out_vars:
            outputs[key] = np.full((nx, ny, nz, ntimes), np.nan)

        outputs["t_driversetup"][:, :, :, 0] = state.t.field[:]
        outputs["local_dp_driversetup"][:, :, :, 0] = locals.dp.field[:]
        outputs["critical_relative_humidity_for_pdf_driversetup"][:, :, :, 0] = (
            state.critical_relative_humidity_for_pdf.field[:]
        )
        outputs["radiation_vapor_driversetup"][:, :, :, 0] = state.radiation_field.vapor.field[:]
        outputs["radiation_liquid_driversetup"][:, :, :, 0] = state.radiation_field.liquid.field[:]
        outputs["radiation_ice_driversetup"][:, :, :, 0] = state.radiation_field.ice.field[:]
        outputs["radiation_rain_driversetup"][:, :, :, 0] = state.radiation_field.rain.field[:]
        outputs["radiation_snow_driversetup"][:, :, :, 0] = state.radiation_field.snow.field[:]
        outputs["radiation_graupel_driversetup"][:, :, :, 0] = state.radiation_field.graupel.field[:]
        outputs["radiation_cloud_fraction_driversetup"][:, :, :, 0] = (
            state.radiation_field.cloud_fraction.field[:]
        )
        outputs["total_concentration_driversetup"][:, :, :, 0] = locals.total_concentration.field[:]
        outputs["driver_local_dry_mixing_ratio_vapor_unmodified_driversetup"][:, :, :, 0] = (
            driver_locals.unmodified.mixing_ratio.vapor.field[:]
        )
        outputs["driver_local_dry_mixing_ratio_liquid_unmodified_driversetup"][:, :, :, 0] = (
            driver_locals.unmodified.mixing_ratio.liquid.field[:]
        )
        outputs["driver_local_dry_mixing_ratio_rain_unmodified_driversetup"][:, :, :, 0] = (
            driver_locals.unmodified.mixing_ratio.rain.field[:]
        )
        outputs["driver_local_dry_mixing_ratio_ice_unmodified_driversetup"][:, :, :, 0] = (
            driver_locals.unmodified.mixing_ratio.ice.field[:]
        )
        outputs["driver_local_dry_mixing_ratio_snow_unmodified_driversetup"][:, :, :, 0] = (
            driver_locals.unmodified.mixing_ratio.snow.field[:]
        )
        outputs["driver_local_dry_mixing_ratio_graupel_unmodified_driversetup"][:, :, :, 0] = (
            driver_locals.unmodified.mixing_ratio.graupel.field[:]
        )
        outputs["driver_local_dry_mixing_ratio_vapor_driversetup"][:, :, :, 0] = (
            driver_locals.dry_air_mixing_ratio.vapor.field[:]
        )
        outputs["driver_local_dry_mixing_ratio_liquid_driversetup"][:, :, :, 0] = (
            driver_locals.dry_air_mixing_ratio.liquid.field[:]
        )
        outputs["driver_local_dry_mixing_ratio_rain_driversetup"][:, :, :, 0] = (
            driver_locals.dry_air_mixing_ratio.rain.field[:]
        )
        outputs["driver_local_dry_mixing_ratio_ice_driversetup"][:, :, :, 0] = (
            driver_locals.dry_air_mixing_ratio.ice.field[:]
        )
        outputs["driver_local_dry_mixing_ratio_snow_driversetup"][:, :, :, 0] = (
            driver_locals.dry_air_mixing_ratio.snow.field[:]
        )
        outputs["driver_local_dry_mixing_ratio_graupel_driversetup"][:, :, :, 0] = (
            driver_locals.dry_air_mixing_ratio.graupel.field[:]
        )
        outputs["driver_local_cloud_fraction_driversetup"][:, :, :, 0] = driver_locals.cloud_fraction.field[:]
        outputs["local_dz_drivesetup"][:, :, :, 0] = locals.dz.field[:]
        outputs["u_driversetup"][:, :, :, 0] = state.u.field[:]
        outputs["v_driversetup"][:, :, :, 0] = state.v.field[:]
        outputs["w_driversetup"][:, :, :, 0] = state.vertical_motion.velocity.field[:]
        outputs["driver_local_t_driversetup"][:, :, :, 0] = driver_locals.t.field[:]
        outputs["driver_local_dp_driversetup"][:, :, :, 0] = driver_locals.dp.field[:]
        outputs["driver_local_density_unmodified_driversetup"][:, :, :, 0] = (
            driver_locals.density_unmodified.field[:]
        )
        outputs["driver_local_p_dry_driversetup"][:, :, :, 0] = driver_locals.p_dry.field[:]
        outputs["driver_local_mass_driversetup"][:, :, :, 0] = driver_locals.mass.field[:]
        outputs["driver_local_u_driversetup"][:, :, :, 0] = driver_locals.u.field[:]
        outputs["driver_local_v_driversetup"][:, :, :, 0] = driver_locals.v.field[:]
        outputs["driver_local_w_driversetup"][:, :, :, 0] = driver_locals.w.field[:]
        outputs["driver_local_ccn_driversetup"][:, :, :, 0] = driver_locals.ccn.field[:]
        outputs["driver_local_c_praut_driversetup"][:, :, :, 0] = driver_locals.c_praut.field[:]
        outputs["driver_local_rh_limited_driversetup"][:, :, :, 0] = driver_locals.rh_limited.field[:]
        outputs["non_anvil_large_scale_liquid_precip_flux_driversetup"][:, :, :, 0] = (
            state.non_anvil_large_scale.liquid_precip_flux.field[:, :, 0:-1]
        )
        outputs["non_anvil_large_scale_ice_precip_flux_driversetup"][:, :, :, 0] = (
            state.non_anvil_large_scale.ice_precip_flux.field[:, :, 0:-1]
        )
        outputs["non_anvil_large_scale_evaporation_driversetup"][:, :, :, 0] = (
            state.non_anvil_large_scale.evaporation.field[:]
        )
        outputs["non_anvil_large_scale_sublimation_driversetup"][:, :, :, 0] = (
            state.non_anvil_large_scale.sublimation.field[:]
        )

        for k in range(nz):
            outputs["driver_local_one_minus_sigma_driversetup"][:, :, k, 0] = (
                driver_locals.one_minus_sigma.field[:]
            )
            outputs["area_driversetup"][:, :, k, 0] = state.area.field[:]
            outputs["surface_precip_rain_driversetup"][:, :, k, 0] = (
                state.precipitation_at_surface.rain.field[:]
            )
            outputs["surface_precip_snow_driversetup"][:, :, k, 0] = (
                state.precipitation_at_surface.snow.field[:]
            )
            outputs["surface_precip_graupel_driversetup"][:, :, k, 0] = (
                state.precipitation_at_surface.graupel.field[:]
            )
            outputs["surface_precip_ice_driversetup"][:, :, k, 0] = state.precipitation_at_surface.ice.field[
                :
            ]

        return outputs
