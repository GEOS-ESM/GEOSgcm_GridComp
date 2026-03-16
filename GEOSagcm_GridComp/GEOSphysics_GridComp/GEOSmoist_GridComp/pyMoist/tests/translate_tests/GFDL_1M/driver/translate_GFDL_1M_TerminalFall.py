import ndsl.xumpy as xp
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
from pyMoist.GFDL_1M.driver.terminal_fall import GFDL1MTerminalFall
from pyMoist.GFDL_1M.state import GFDL1MState


class TranslateGFDL_1M_TerminalFall(TranslateFortranData2Py):
    def __init__(self, grid: Grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.backend = self.stencil_factory.backend
        self.quantity_factory = grid.quantity_factory

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "driver_local_t_terminalfall": {},
            "driver_local_dry_mixing_ratio_vapor_terminalfall": {},
            "driver_local_dry_mixing_ratio_liquid_terminalfall": {},
            "driver_local_dry_mixing_ratio_rain_terminalfall": {},
            "driver_local_dry_mixing_ratio_graupel_terminalfall": {},
            "driver_local_dry_mixing_ratio_snow_terminalfall": {},
            "driver_local_dry_mixing_ratio_ice_terminalfall": {},
            "driver_local_ice_precip_flux_terminalfall": {},
            "driver_local_w_terminalfall": {},
            "driver_local_dz_terminalfall": {},
            "driver_local_dp_terminalfall": {},
            "driver_local_terminal_speed_ice_terminalfall": {},
            "driver_local_terminal_speed_snow_terminalfall": {},
            "driver_local_terminal_speed_graupel_terminalfall": {},
            "surface_precip_rain_terminalfall": {},
            "surface_precip_snow_terminalfall": {},
            "surface_precip_graupel_terminalfall": {},
            "surface_precip_ice_terminalfall": {},
        }

        self.out_vars = self.in_vars["data_vars"].copy()

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GFDL_1M-constants")

    def compute(self, inputs):
        # initialize dataclasses
        driver_locals = GFDL1MDriverLocals.make_as_state(self.quantity_factory)
        # initialize constants
        config = GFDL1MConfig(**self.constants)
        config_dependent_constants = GFDL1MDriverConfigDependentConstants.make(config)

        # get the shape of the field
        state = GFDL1MState.zeros(self.quantity_factory)
        nx, ny, nz = inputs["driver_local_t_terminalfall"].shape

        # preset output dictionary to be filled inside the for loop
        outputs = {}
        for key in self.out_vars:
            outputs[key] = xp.full((nx, ny, nz), self.backend, np.nan)

        # construct test stencil
        # NOTE: required `extra_data_load` means we can't allocate code in __init__
        code = GFDL1MTerminalFall(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            config_dependent_constants=config_dependent_constants,
        )

        safe_assign_array(
            driver_locals.t.field[:],
            inputs["driver_local_t_terminalfall"][:, :, :],
        )
        safe_assign_array(
            driver_locals.dry_air_mixing_ratio.vapor.field[:],
            inputs["driver_local_dry_mixing_ratio_vapor_terminalfall"][:, :, :],
        )
        safe_assign_array(
            driver_locals.dry_air_mixing_ratio.liquid.field[:],
            inputs["driver_local_dry_mixing_ratio_liquid_terminalfall"][:, :, :],
        )
        safe_assign_array(
            driver_locals.dry_air_mixing_ratio.rain.field[:],
            inputs["driver_local_dry_mixing_ratio_rain_terminalfall"][:, :, :],
        )
        safe_assign_array(
            driver_locals.dry_air_mixing_ratio.graupel.field[:],
            inputs["driver_local_dry_mixing_ratio_graupel_terminalfall"][:, :, :],
        )
        safe_assign_array(
            driver_locals.dry_air_mixing_ratio.snow.field[:],
            inputs["driver_local_dry_mixing_ratio_snow_terminalfall"][:, :, :],
        )
        safe_assign_array(
            driver_locals.dry_air_mixing_ratio.ice.field[:],
            inputs["driver_local_dry_mixing_ratio_ice_terminalfall"][:, :, :],
        )
        safe_assign_array(
            driver_locals.ice_precip_flux.field[:],
            inputs["driver_local_ice_precip_flux_terminalfall"][:, :, :],
        )
        safe_assign_array(
            driver_locals.w.field[:],
            inputs["driver_local_w_terminalfall"][:, :, :],
        )
        safe_assign_array(
            driver_locals.dz.field[:],
            inputs["driver_local_dz_terminalfall"][:, :, :],
        )
        safe_assign_array(
            driver_locals.dp.field[:],
            inputs["driver_local_dp_terminalfall"][:, :, :],
        )
        safe_assign_array(
            driver_locals.terminal_speed.ice.field[:],
            inputs["driver_local_terminal_speed_ice_terminalfall"][:, :, :],
        )
        safe_assign_array(
            driver_locals.terminal_speed.snow.field[:],
            inputs["driver_local_terminal_speed_snow_terminalfall"][:, :, :],
        )
        safe_assign_array(
            driver_locals.terminal_speed.graupel.field[:],
            inputs["driver_local_terminal_speed_graupel_terminalfall"][:, :, :],
        )
        safe_assign_array(
            state.precipitation_at_surface.rain.field[:],
            inputs["surface_precip_rain_terminalfall"][:, :, 0],
        )
        safe_assign_array(
            state.precipitation_at_surface.snow.field[:],
            inputs["surface_precip_snow_terminalfall"][:, :, 0],
        )
        safe_assign_array(
            state.precipitation_at_surface.graupel.field[:],
            inputs["surface_precip_graupel_terminalfall"][:, :, 0],
        )
        safe_assign_array(
            state.precipitation_at_surface.ice.field[:],
            inputs["surface_precip_ice_terminalfall"][:, :, 0],
        )

        # run the test code
        code(
            t=driver_locals.t,
            w=driver_locals.w,
            mixing_ratio_vapor=driver_locals.dry_air_mixing_ratio.vapor,
            mixing_ratio_liquid=driver_locals.dry_air_mixing_ratio.liquid,
            mixing_ratio_rain=driver_locals.dry_air_mixing_ratio.rain,
            mixing_ratio_graupel=driver_locals.dry_air_mixing_ratio.graupel,
            mixing_ratio_snow=driver_locals.dry_air_mixing_ratio.snow,
            mixing_ratio_ice=driver_locals.dry_air_mixing_ratio.ice,
            dz=driver_locals.dz,
            dp=driver_locals.dp,
            terminal_velocity_graupel=driver_locals.terminal_speed.graupel,
            terminal_velocity_snow=driver_locals.terminal_speed.snow,
            terminal_velocity_ice=driver_locals.terminal_speed.ice,
            rain=state.precipitation_at_surface.rain,
            graupel=state.precipitation_at_surface.graupel,
            snow=state.precipitation_at_surface.snow,
            ice=state.precipitation_at_surface.ice,
            ice_precip_flux=driver_locals.ice_precip_flux,
        )

        # fill the output arrays so that all calls are tested
        safe_assign_array(
            outputs["driver_local_t_terminalfall"][:, :, :],
            driver_locals.t.field[:],
        )
        safe_assign_array(
            outputs["driver_local_dry_mixing_ratio_vapor_terminalfall"][:, :, :],
            driver_locals.dry_air_mixing_ratio.vapor.field[:],
        )
        safe_assign_array(
            outputs["driver_local_dry_mixing_ratio_liquid_terminalfall"][:, :, :],
            driver_locals.dry_air_mixing_ratio.liquid.field[:],
        )
        safe_assign_array(
            outputs["driver_local_dry_mixing_ratio_rain_terminalfall"][:, :, :],
            driver_locals.dry_air_mixing_ratio.rain.field[:],
        )
        safe_assign_array(
            outputs["driver_local_dry_mixing_ratio_graupel_terminalfall"][:, :, :],
            driver_locals.dry_air_mixing_ratio.graupel.field[:],
        )
        safe_assign_array(
            outputs["driver_local_dry_mixing_ratio_snow_terminalfall"][:, :, :],
            driver_locals.dry_air_mixing_ratio.snow.field[:],
        )
        safe_assign_array(
            outputs["driver_local_dry_mixing_ratio_ice_terminalfall"][:, :, :],
            driver_locals.dry_air_mixing_ratio.ice.field[:],
        )
        safe_assign_array(
            outputs["driver_local_ice_precip_flux_terminalfall"][:, :, :],
            driver_locals.ice_precip_flux.field[:],
        )
        safe_assign_array(
            outputs["driver_local_w_terminalfall"][:, :, :],
            driver_locals.w.field[:],
        )
        safe_assign_array(
            outputs["driver_local_dz_terminalfall"][:, :, :],
            driver_locals.dz.field[:],
        )
        safe_assign_array(
            outputs["driver_local_dp_terminalfall"][:, :, :],
            driver_locals.dp.field[:],
        )
        safe_assign_array(
            outputs["driver_local_terminal_speed_ice_terminalfall"][:, :, :],
            driver_locals.terminal_speed.ice.field[:],
        )
        safe_assign_array(
            outputs["driver_local_terminal_speed_snow_terminalfall"][:, :, :],
            driver_locals.terminal_speed.snow.field[:],
        )
        safe_assign_array(
            outputs["driver_local_terminal_speed_graupel_terminalfall"][:, :, :],
            driver_locals.terminal_speed.graupel.field[:],
        )

        for k in range(nz):
            safe_assign_array(
                outputs["surface_precip_rain_terminalfall"][:, :, k],
                state.precipitation_at_surface.rain.field[:],
            )
            safe_assign_array(
                outputs["surface_precip_snow_terminalfall"][:, :, k],
                state.precipitation_at_surface.snow.field[:],
            )
            safe_assign_array(
                outputs["surface_precip_graupel_terminalfall"][:, :, k],
                state.precipitation_at_surface.graupel.field[:],
            )
            safe_assign_array(
                outputs["surface_precip_ice_terminalfall"][:, :, k],
                state.precipitation_at_surface.ice.field[:],
            )

        return outputs
