import numpy as np
from f90nml import Namelist

from ndsl import StencilFactory
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array
from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.GFDL_1M.driver.config_constants import GFDL1MDriverConfigDependentConstants
from pyMoist.GFDL_1M.driver.locals import GFDL1MDriverLocals
from pyMoist.GFDL_1M.driver.sat_tables import get_tables
from pyMoist.GFDL_1M.driver.warm_rain import GFDL1MWarmRain
from pyMoist.GFDL_1M.state import GFDL1MState


class TranslateGFDL_1M_WarmRain(TranslateFortranData2Py):
    def __init__(self, grid: Grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "driver_local_dp_warmrain": {},
            "driver_local_dz_warmrain": {},
            "driver_local_t_warmrain": {},
            "driver_local_dry_mixing_ratio_vapor_warmrain": {},
            "driver_local_dry_mixing_ratio_liquid_warmrain": {},
            "driver_local_dry_mixing_ratio_rain_warmrain": {},
            "driver_local_dry_mixing_ratio_ice_warmrain": {},
            "driver_local_dry_mixing_ratio_snow_warmrain": {},
            "driver_local_dry_mixing_ratio_graupel_warmrain": {},
            "driver_local_cloud_fraction_warmrain": {},
            "driver_local_ccn_warmrain": {},
            "driver_local_density_warmrain": {},
            "driver_local_density_factor_warmrain": {},
            "driver_local_c_praut_warmrain": {},
            "driver_local_terminal_speed_rain_warmrain": {},
            "driver_local_evaporation_warmrain": {},
            "driver_local_liquid_precip_flux_warmrain": {},
            "driver_local_w_warmrain": {},
            "driver_local_rh_limited_warmrain": {},
            "non_anvil_large_scale_evaporation_warmrain": {},
            "non_anvil_large_scale_liquid_precip_flux_warmrain": {},
            "non_anvil_large_scale_ice_precip_flux_warmrain": {},
            "driver_local_mass_warmrain": {},
            "driver_local_ice_precip_flux_warmrain": {},
            "driver_local_rain_warmrain": {},
            "surface_precip_rain_warmrain": {},
            "estimated_inversion_strength_warmrain": {},
            "driver_local_one_minus_sigma_warmrain": {},
        }

        self.out_vars = self.in_vars["data_vars"].copy()

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GFDL_1M-constants")

    def compute(self, inputs):
        # initialize dataclasses
        state = GFDL1MState.zeros(self.quantity_factory)
        driver_locals = GFDL1MDriverLocals.make_as_state(self.quantity_factory)

        # initialize constants
        config = GFDL1MConfig(**self.constants)
        config_dependent_constants = GFDL1MDriverConfigDependentConstants.make(config)

        # initialize saturation tables
        saturation_tables = get_tables(
            backend=self.stencil_factory.backend,
            dace_config=self.stencil_factory.config.dace_config,
        )

        # get the shape of the field
        nx, ny, nz, ntimes = inputs["driver_local_dp_warmrain"].shape

        # preset output dictionary to be filled inside the for loop
        outputs = {}
        for key in self.out_vars:
            outputs[key] = np.full((nx, ny, nz, ntimes), np.nan)

        # construct test stencil
        code = GFDL1MWarmRain(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            config_dependent_constants=config_dependent_constants,
            saturation_tables=saturation_tables,
        )

        for n in range(1, 2):
            safe_assign_array(
                driver_locals.dp.field[:],
                inputs["driver_local_dp_warmrain"][:, :, :, n],
            )
            safe_assign_array(
                driver_locals.dz.field[:],
                inputs["driver_local_dz_warmrain"][:, :, :, n],
            )
            safe_assign_array(driver_locals.t.field[:], inputs["driver_local_t_warmrain"][:, :, :, n])
            safe_assign_array(
                driver_locals.dry_air_mixing_ratio.vapor.field[:],
                inputs["driver_local_dry_mixing_ratio_vapor_warmrain"][:, :, :, n],
            )
            safe_assign_array(
                driver_locals.dry_air_mixing_ratio.liquid.field[:],
                inputs["driver_local_dry_mixing_ratio_liquid_warmrain"][:, :, :, n],
            )
            safe_assign_array(
                driver_locals.dry_air_mixing_ratio.rain.field[:],
                inputs["driver_local_dry_mixing_ratio_rain_warmrain"][:, :, :, n],
            )
            safe_assign_array(
                driver_locals.dry_air_mixing_ratio.ice.field[:],
                inputs["driver_local_dry_mixing_ratio_ice_warmrain"][:, :, :, n],
            )
            safe_assign_array(
                driver_locals.dry_air_mixing_ratio.snow.field[:],
                inputs["driver_local_dry_mixing_ratio_snow_warmrain"][:, :, :, n],
            )
            safe_assign_array(
                driver_locals.dry_air_mixing_ratio.graupel.field[:],
                inputs["driver_local_dry_mixing_ratio_graupel_warmrain"][:, :, :, n],
            )
            safe_assign_array(
                driver_locals.cloud_fraction.field[:],
                inputs["driver_local_cloud_fraction_warmrain"][:, :, :, n],
            )
            safe_assign_array(
                driver_locals.ccn.field[:],
                inputs["driver_local_ccn_warmrain"][:, :, :, n],
            )
            safe_assign_array(
                driver_locals.density.field[:],
                inputs["driver_local_density_warmrain"][:, :, :, n],
            )
            safe_assign_array(
                driver_locals.density_factor.field[:],
                inputs["driver_local_density_factor_warmrain"][:, :, :, n],
            )
            safe_assign_array(
                driver_locals.c_praut.field[:],
                inputs["driver_local_c_praut_warmrain"][:, :, :, n],
            )
            safe_assign_array(
                driver_locals.terminal_speed.rain.field[:],
                inputs["driver_local_terminal_speed_rain_warmrain"][:, :, :, n],
            )
            safe_assign_array(
                driver_locals.evaporation.field[:],
                inputs["driver_local_evaporation_warmrain"][:, :, :, n],
            )
            safe_assign_array(
                driver_locals.liquid_precip_flux.field[:],
                inputs["driver_local_liquid_precip_flux_warmrain"][:, :, :, n],
            )
            safe_assign_array(driver_locals.w.field[:], inputs["driver_local_w_warmrain"][:, :, :, n])
            safe_assign_array(
                driver_locals.rh_limited.field[:],
                inputs["driver_local_rh_limited_warmrain"][:, :, :, n],
            )
            safe_assign_array(
                state.non_anvil_large_scale.evaporation.field[:],
                inputs["non_anvil_large_scale_evaporation_warmrain"][:, :, :, n],
            )
            safe_assign_array(
                state.non_anvil_large_scale.liquid_precip_flux.field[:, :, 0:-1],
                inputs["non_anvil_large_scale_liquid_precip_flux_warmrain"][:, :, :, n],
            )
            safe_assign_array(
                state.non_anvil_large_scale.ice_precip_flux.field[:, :, 0:-1],
                inputs["non_anvil_large_scale_ice_precip_flux_warmrain"][:, :, :, n],
            )
            safe_assign_array(
                driver_locals.mass.field[:],
                inputs["driver_local_mass_warmrain"][:, :, :, n],
            )
            safe_assign_array(
                driver_locals.ice_precip_flux.field[:],
                inputs["driver_local_ice_precip_flux_warmrain"][:, :, :, n],
            )
            safe_assign_array(
                driver_locals.rain.field[:],
                inputs["driver_local_rain_warmrain"][:, :, 0, n],
            )
            safe_assign_array(
                state.precipitation_at_surface.rain.field[:],
                inputs["surface_precip_rain_warmrain"][:, :, 0, n],
            )
            safe_assign_array(
                state.estimated_inversion_strength.field[:],
                inputs["estimated_inversion_strength_warmrain"][:, :, 0, n],
            )
            safe_assign_array(
                driver_locals.one_minus_sigma.field[:],
                inputs["driver_local_one_minus_sigma_warmrain"][:, :, 0, n],
            )

            # run the test code
            code(
                t=driver_locals.t,
                dp=driver_locals.dp,
                dz=driver_locals.dz,
                w=driver_locals.w,
                mixing_ratio_vapor=driver_locals.dry_air_mixing_ratio.vapor,
                mixing_ratio_liquid=driver_locals.dry_air_mixing_ratio.liquid,
                mixing_ratio_rain=driver_locals.dry_air_mixing_ratio.rain,
                mixing_ratio_ice=driver_locals.dry_air_mixing_ratio.ice,
                mixing_ratio_snow=driver_locals.dry_air_mixing_ratio.snow,
                mixing_ratio_graupel=driver_locals.dry_air_mixing_ratio.graupel,
                cloud_fraction=driver_locals.cloud_fraction,
                ccn=driver_locals.ccn,
                density=driver_locals.density,
                density_factor=driver_locals.density_factor,
                c_praut=driver_locals.c_praut,
                terminal_speed_rain=driver_locals.terminal_speed.rain,
                rh_limited=driver_locals.rh_limited,
                estimated_inversion_strength=state.estimated_inversion_strength,
                one_minus_sigma=driver_locals.one_minus_sigma,
                mass=driver_locals.mass,
                rain=state.precipitation_at_surface.rain,
                driver_rain=driver_locals.rain,
                ice_precip_flux=state.non_anvil_large_scale.ice_precip_flux,
                driver_ice_precip_flux=driver_locals.ice_precip_flux,
                liquid_precip_flux=state.non_anvil_large_scale.liquid_precip_flux,
                driver_liquid_precip_flux=driver_locals.liquid_precip_flux,
                evaporation=state.non_anvil_large_scale.evaporation,
                driver_evaporation=driver_locals.evaporation,
            )

            # fill the output arrays so that all calls are tested
            safe_assign_array(
                outputs["driver_local_dp_warmrain"][:, :, :, n],
                driver_locals.dp.field[:],
            )
            safe_assign_array(
                outputs["driver_local_dz_warmrain"][:, :, :, n],
                driver_locals.dz.field[:],
            )
            safe_assign_array(outputs["driver_local_t_warmrain"][:, :, :, n], driver_locals.t.field[:])
            safe_assign_array(
                outputs["driver_local_dry_mixing_ratio_vapor_warmrain"][:, :, :, n],
                driver_locals.dry_air_mixing_ratio.vapor.field[:],
            )
            safe_assign_array(
                outputs["driver_local_dry_mixing_ratio_liquid_warmrain"][:, :, :, n],
                driver_locals.dry_air_mixing_ratio.liquid.field[:],
            )
            safe_assign_array(
                outputs["driver_local_dry_mixing_ratio_rain_warmrain"][:, :, :, n],
                driver_locals.dry_air_mixing_ratio.rain.field[:],
            )
            safe_assign_array(
                outputs["driver_local_dry_mixing_ratio_ice_warmrain"][:, :, :, n],
                driver_locals.dry_air_mixing_ratio.ice.field[:],
            )
            safe_assign_array(
                outputs["driver_local_dry_mixing_ratio_snow_warmrain"][:, :, :, n],
                driver_locals.dry_air_mixing_ratio.snow.field[:],
            )
            safe_assign_array(
                outputs["driver_local_dry_mixing_ratio_graupel_warmrain"][:, :, :, n],
                driver_locals.dry_air_mixing_ratio.graupel.field[:],
            )
            safe_assign_array(
                outputs["driver_local_cloud_fraction_warmrain"][:, :, :, n],
                driver_locals.cloud_fraction.field[:],
            )
            safe_assign_array(
                outputs["driver_local_ccn_warmrain"][:, :, :, n],
                driver_locals.ccn.field[:],
            )
            safe_assign_array(
                outputs["driver_local_density_warmrain"][:, :, :, n],
                driver_locals.density.field[:],
            )
            safe_assign_array(
                outputs["driver_local_density_factor_warmrain"][:, :, :, n],
                driver_locals.density_factor.field[:],
            )
            safe_assign_array(
                outputs["driver_local_c_praut_warmrain"][:, :, :, n],
                driver_locals.c_praut.field[:],
            )
            safe_assign_array(
                outputs["driver_local_terminal_speed_rain_warmrain"][:, :, :, n],
                driver_locals.terminal_speed.rain.field[:],
            )
            safe_assign_array(
                outputs["driver_local_evaporation_warmrain"][:, :, :, n],
                driver_locals.evaporation.field[:],
            )
            safe_assign_array(
                outputs["driver_local_liquid_precip_flux_warmrain"][:, :, :, n],
                driver_locals.liquid_precip_flux.field[:],
            )
            safe_assign_array(outputs["driver_local_w_warmrain"][:, :, :, n], driver_locals.w.field[:])
            safe_assign_array(
                outputs["driver_local_rh_limited_warmrain"][:, :, :, n],
                driver_locals.rh_limited.field[:],
            )
            safe_assign_array(
                outputs["non_anvil_large_scale_evaporation_warmrain"][:, :, :, n],
                state.non_anvil_large_scale.evaporation.field[:],
            )
            safe_assign_array(
                outputs["non_anvil_large_scale_liquid_precip_flux_warmrain"][:, :, :, n],
                state.non_anvil_large_scale.liquid_precip_flux.field[:, :, 0:-1],
            )
            safe_assign_array(
                outputs["non_anvil_large_scale_ice_precip_flux_warmrain"][:, :, :, n],
                state.non_anvil_large_scale.ice_precip_flux.field[:, :, 0:-1],
            )
            safe_assign_array(
                outputs["driver_local_mass_warmrain"][:, :, :, n],
                driver_locals.mass.field[:],
            )
            safe_assign_array(
                outputs["driver_local_ice_precip_flux_warmrain"][:, :, :, n],
                driver_locals.ice_precip_flux.field[:],
            )

            for k in range(nz):
                safe_assign_array(
                    outputs["driver_local_rain_warmrain"][:, :, k, n],
                    driver_locals.rain.field[:],
                )
                safe_assign_array(
                    outputs["surface_precip_rain_warmrain"][:, :, k, n],
                    state.precipitation_at_surface.rain.field[:],
                )
                safe_assign_array(
                    outputs["estimated_inversion_strength_warmrain"][:, :, k, n],
                    state.estimated_inversion_strength.field[:],
                )
                safe_assign_array(
                    outputs["driver_local_one_minus_sigma_warmrain"][:, :, k, n],
                    driver_locals.one_minus_sigma.field[:],
                )

        return outputs
