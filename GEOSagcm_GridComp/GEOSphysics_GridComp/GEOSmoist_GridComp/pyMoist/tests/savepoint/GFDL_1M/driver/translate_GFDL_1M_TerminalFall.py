from f90nml import Namelist

from ndsl import StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.GFDL_1M.state import GFDL1MState
from pyMoist.GFDL_1M.locals import GFDL1MLocals
from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.GFDL_1M.driver.config_constants import GFDL1MDriverConfigDependentConstants
from pyMoist.GFDL_1M.driver.locals import GFDL1MDriverLocals
from pyMoist.GFDL_1M.driver.terminal_fall import TerminalFall
from ndsl.stencils.testing.savepoint import DataLoader
import numpy as np


class TranslateGFDL_1M_TerminalFall(TranslateFortranData2Py):
    def __init__(self, grid: Grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
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
            "driver_local_rain_terminalfall": {},
            "driver_local_ice_terminalfall": {},
            "driver_local_snow_terminalfall": {},
            "driver_local_graupel_terminalfall": {},
            "surface_precip_rain_terminalfall": {},
            "surface_precip_snow_terminalfall": {},
            "surface_precip_graupel_terminalfall": {},
            "surface_precip_ice_terminalfall": {},
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

        # get the shape of the field
        nx, ny, nz, ntimes = inputs["driver_local_t_terminalfall"].shape

        # preset output dictionary to be filled inside the for loop
        outputs = {}
        for key in self.out_vars:
            outputs[key] = np.full((nx, ny, nz, ntimes), np.nan)

        # construct test stencil
        code = TerminalFall(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            config_dependent_constants=config_dependent_constants,
        )

        for n in range(1):  # range(ntimes):
            driver_locals.t.field[:] = inputs["driver_local_t_terminalfall"][:, :, :, n]
            driver_locals.dry_air_mixing_ratio.vapor.field[:] = inputs[
                "driver_local_dry_mixing_ratio_vapor_terminalfall"
            ][:, :, :, n]
            driver_locals.dry_air_mixing_ratio.liquid.field[:] = inputs[
                "driver_local_dry_mixing_ratio_liquid_terminalfall"
            ][:, :, :, n]
            driver_locals.dry_air_mixing_ratio.rain.field[:] = inputs[
                "driver_local_dry_mixing_ratio_rain_terminalfall"
            ][:, :, :, n]
            driver_locals.dry_air_mixing_ratio.graupel.field[:] = inputs[
                "driver_local_dry_mixing_ratio_graupel_terminalfall"
            ][:, :, :, n]
            driver_locals.dry_air_mixing_ratio.snow.field[:] = inputs[
                "driver_local_dry_mixing_ratio_snow_terminalfall"
            ][:, :, :, n]
            driver_locals.dry_air_mixing_ratio.ice.field[:] = inputs[
                "driver_local_dry_mixing_ratio_ice_terminalfall"
            ][:, :, :, n]
            driver_locals.ice_precip_flux.field[:] = inputs["driver_local_ice_precip_flux_terminalfall"][
                :, :, :, n
            ]
            driver_locals.w.field[:] = inputs["driver_local_w_terminalfall"][:, :, :, n]
            driver_locals.dz.field[:] = inputs["driver_local_dz_terminalfall"][:, :, :, n]
            driver_locals.dp.field[:] = inputs["driver_local_dp_terminalfall"][:, :, :, n]
            driver_locals.terminal_speed.ice.field[:] = inputs[
                "driver_local_terminal_speed_ice_terminalfall"
            ][:, :, :, n]
            driver_locals.terminal_speed.snow.field[:] = inputs[
                "driver_local_terminal_speed_snow_terminalfall"
            ][:, :, :, n]
            driver_locals.terminal_speed.graupel.field[:] = inputs[
                "driver_local_terminal_speed_graupel_terminalfall"
            ][:, :, :, n]
            state.precipitation_at_surface.rain.field[:] = inputs["surface_precip_rain_terminalfall"][
                :, :, 0, n
            ]
            state.precipitation_at_surface.snow.field[:] = inputs["surface_precip_snow_terminalfall"][
                :, :, 0, n
            ]
            state.precipitation_at_surface.graupel.field[:] = inputs["surface_precip_graupel_terminalfall"][
                :, :, 0, n
            ]
            state.precipitation_at_surface.ice.field[:] = inputs["surface_precip_ice_terminalfall"][
                :, :, 0, n
            ]

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
                graupel=state.precipitation_at_surface.snow,
                snow=state.precipitation_at_surface.graupel,
                ice=state.precipitation_at_surface.ice,
                ice_precip_flux=driver_locals.ice_precip_flux,
            )

            # fill the output arrays so that all calls are tested
            outputs["driver_local_t_terminalfall"][:, :, :, n] = driver_locals.t.field[:]
            outputs["driver_local_dry_mixing_ratio_vapor_terminalfall"][:, :, :, n] = (
                driver_locals.dry_air_mixing_ratio.vapor.field[:]
            )
            outputs["driver_local_dry_mixing_ratio_liquid_terminalfall"][:, :, :, n] = (
                driver_locals.dry_air_mixing_ratio.liquid.field[:]
            )
            outputs["driver_local_dry_mixing_ratio_rain_terminalfall"][:, :, :, n] = (
                driver_locals.dry_air_mixing_ratio.rain.field[:]
            )
            outputs["driver_local_dry_mixing_ratio_graupel_terminalfall"][:, :, :, n] = (
                driver_locals.dry_air_mixing_ratio.graupel.field[:]
            )
            outputs["driver_local_dry_mixing_ratio_snow_terminalfall"][:, :, :, n] = (
                driver_locals.dry_air_mixing_ratio.snow.field[:]
            )
            outputs["driver_local_dry_mixing_ratio_ice_terminalfall"][:, :, :, n] = (
                driver_locals.dry_air_mixing_ratio.ice.field[:]
            )
            outputs["driver_local_ice_precip_flux_terminalfall"][:, :, :, n] = (
                driver_locals.ice_precip_flux.field[:]
            )
            outputs["driver_local_w_terminalfall"][:, :, :, n] = driver_locals.w.field[:]
            outputs["driver_local_dz_terminalfall"][:, :, :, n] = driver_locals.dz.field[:]
            outputs["driver_local_dp_terminalfall"][:, :, :, n] = driver_locals.dp.field[:]
            outputs["driver_local_terminal_speed_ice_terminalfall"][:, :, :, n] = (
                driver_locals.terminal_speed.ice.field[:]
            )
            outputs["driver_local_terminal_speed_snow_terminalfall"][:, :, :, n] = (
                driver_locals.terminal_speed.snow.field[:]
            )
            outputs["driver_local_terminal_speed_graupel_terminalfall"][:, :, :, n] = (
                driver_locals.terminal_speed.graupel.field[:]
            )

            for k in range(nz):
                outputs["driver_local_rain_terminalfall"][:, :, k, n] = code._locals.rain.field[:]
                outputs["driver_local_ice_terminalfall"][:, :, k, n] = code._locals.ice.field[:]
                outputs["driver_local_snow_terminalfall"][:, :, k, n] = code._locals.snow.field[:]
                outputs["driver_local_graupel_terminalfall"][:, :, k, n] = code._locals.graupel.field[:]
                outputs["surface_precip_rain_terminalfall"][:, :, k, n] = (
                    state.precipitation_at_surface.rain.field[:]
                )
                outputs["surface_precip_snow_terminalfall"][:, :, k, n] = (
                    state.precipitation_at_surface.snow.field[:]
                )
                outputs["surface_precip_graupel_terminalfall"][:, :, k, n] = (
                    state.precipitation_at_surface.graupel.field[:]
                )
                outputs["surface_precip_ice_terminalfall"][:, :, k, n] = (
                    state.precipitation_at_surface.ice.field[:]
                )

        return outputs
