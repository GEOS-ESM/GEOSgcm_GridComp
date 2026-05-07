import numpy as np
from f90nml import Namelist
from ndsl import StencilFactory
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py

from pyMoist.microphysics.GFDL_1M.config import GFDL1MConfig
from pyMoist.microphysics.GFDL_1M.driver.config_constants import GFDL1MDriverConfigDependentConstants
from pyMoist.microphysics.GFDL_1M.driver.ice_cloud import GFDL1MIceCloud
from pyMoist.microphysics.GFDL_1M.driver.locals import GFDL1MDriverLocals
from pyMoist.microphysics.GFDL_1M.driver.sat_tables import get_tables
from pyMoist.microphysics.GFDL_1M.state import GFDL1MState


class TranslateGFDL_1M_IceCloud(TranslateFortranData2Py):
    def __init__(self, grid: Grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.quantity_factory = grid.quantity_factory

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "driver_local_t_icecloud": {},
            "driver_local_p_dry_icecloud": {},
            "driver_local_dp_icecloud": {},
            "driver_local_dry_mixing_ratio_vapor_icecloud": {},
            "driver_local_dry_mixing_ratio_liquid_icecloud": {},
            "driver_local_dry_mixing_ratio_rain_icecloud": {},
            "driver_local_dry_mixing_ratio_ice_icecloud": {},
            "driver_local_dry_mixing_ratio_snow_icecloud": {},
            "driver_local_dry_mixing_ratio_graupel_icecloud": {},
            "driver_local_cloud_fraction_icecloud": {},
            "driver_local_terminal_speed_snow_icecloud": {},
            "driver_local_terminal_speed_graupel_icecloud": {},
            "driver_local_terminal_speed_rain_icecloud": {},
            "driver_local_density_icecloud": {},
            "driver_local_density_factor_icecloud": {},
            "driver_local_rh_limited_icecloud": {},
            "non_anvil_large_scale_sublimation_icecloud": {},
            "driver_local_ccn_icecloud": {},
            "convection_fraction_icecloud": {},
            "surface_type_icecloud": {},
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
        nx, ny, nz, _ntimes = inputs["driver_local_t_icecloud"].shape

        # preset output dictionary to be filled inside the for loop
        outputs = {}
        for key in self.out_vars:
            outputs[key] = np.full((nx, ny, nz), np.nan)

        # construct test stencil
        code = GFDL1MIceCloud(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            config_dependent_constants=config_dependent_constants,
            saturation_tables=saturation_tables,
        )

        driver_locals.t.field[:] = inputs["driver_local_t_icecloud"][:, :, :]
        driver_locals.p_dry.field[:] = inputs["driver_local_p_dry_icecloud"][:, :, :]
        driver_locals.dp.field[:] = inputs["driver_local_dp_icecloud"][:, :, :]
        driver_locals.dry_air_mixing_ratio.vapor.field[:] = inputs["driver_local_dry_mixing_ratio_vapor_icecloud"][:, :, :]
        driver_locals.dry_air_mixing_ratio.liquid.field[:] = inputs["driver_local_dry_mixing_ratio_liquid_icecloud"][:, :, :]
        driver_locals.dry_air_mixing_ratio.rain.field[:] = inputs["driver_local_dry_mixing_ratio_rain_icecloud"][:, :, :]
        driver_locals.dry_air_mixing_ratio.ice.field[:] = inputs["driver_local_dry_mixing_ratio_ice_icecloud"][:, :, :]
        driver_locals.dry_air_mixing_ratio.snow.field[:] = inputs["driver_local_dry_mixing_ratio_snow_icecloud"][:, :, :]
        driver_locals.dry_air_mixing_ratio.graupel.field[:] = inputs["driver_local_dry_mixing_ratio_graupel_icecloud"][:, :, :]
        driver_locals.cloud_fraction.field[:] = inputs["driver_local_cloud_fraction_icecloud"][:, :, :]
        driver_locals.terminal_speed.snow.field[:] = inputs["driver_local_terminal_speed_snow_icecloud"][:, :, :]
        driver_locals.terminal_speed.graupel.field[:] = inputs["driver_local_terminal_speed_graupel_icecloud"][:, :, :]
        driver_locals.terminal_speed.rain.field[:] = inputs["driver_local_terminal_speed_rain_icecloud"][:, :, :]
        driver_locals.density.field[:] = inputs["driver_local_density_icecloud"][:, :, :]
        driver_locals.density_factor.field[:] = inputs["driver_local_density_factor_icecloud"][:, :, :]
        driver_locals.rh_limited.field[:] = inputs["driver_local_rh_limited_icecloud"][:, :, :]
        state.non_anvil_large_scale.sublimation.field[:] = inputs["non_anvil_large_scale_sublimation_icecloud"][:, :, :]
        driver_locals.ccn.field[:] = inputs["driver_local_ccn_icecloud"][:, :, :]
        state.convection_fraction.field[:] = inputs["convection_fraction_icecloud"][:, :, 0]
        state.surface_type.field[:] = inputs["surface_type_icecloud"][:, :, 0]

        # run the test code
        code(
            t=driver_locals.t,
            p_dry=driver_locals.p_dry,
            dp=driver_locals.dp,
            mixing_ratio_vapor=driver_locals.dry_air_mixing_ratio.vapor,
            mixing_ratio_liquid=driver_locals.dry_air_mixing_ratio.liquid,
            mixing_ratio_rain=driver_locals.dry_air_mixing_ratio.rain,
            mixing_ratio_ice=driver_locals.dry_air_mixing_ratio.ice,
            mixing_ratio_snow=driver_locals.dry_air_mixing_ratio.snow,
            mixing_ratio_graupel=driver_locals.dry_air_mixing_ratio.graupel,
            cloud_fraction=driver_locals.cloud_fraction,
            density=driver_locals.density,
            density_factor=driver_locals.density_factor,
            terminal_fall_snow=driver_locals.terminal_speed.snow,
            terminal_fall_graupel=driver_locals.terminal_speed.graupel,
            terminal_fall_rain=driver_locals.terminal_speed.rain,
            sublimation=state.non_anvil_large_scale.sublimation,
            rh_limited=driver_locals.rh_limited,
            ccn=driver_locals.ccn,
            convection_fraction=state.convection_fraction,
            surface_type=state.surface_type,
        )

        # fill the output arrays so that all calls are tested
        outputs["driver_local_t_icecloud"][:, :, :] = driver_locals.t.field[:]
        outputs["driver_local_p_dry_icecloud"][:, :, :] = driver_locals.p_dry.field[:]
        outputs["driver_local_dp_icecloud"][:, :, :] = driver_locals.dp.field[:]
        outputs["driver_local_dry_mixing_ratio_vapor_icecloud"][:, :, :] = driver_locals.dry_air_mixing_ratio.vapor.field[:]
        outputs["driver_local_dry_mixing_ratio_liquid_icecloud"][:, :, :] = driver_locals.dry_air_mixing_ratio.liquid.field[:]
        outputs["driver_local_dry_mixing_ratio_rain_icecloud"][:, :, :] = driver_locals.dry_air_mixing_ratio.rain.field[:]
        outputs["driver_local_dry_mixing_ratio_ice_icecloud"][:, :, :] = driver_locals.dry_air_mixing_ratio.ice.field[:]
        outputs["driver_local_dry_mixing_ratio_snow_icecloud"][:, :, :] = driver_locals.dry_air_mixing_ratio.snow.field[:]
        outputs["driver_local_dry_mixing_ratio_graupel_icecloud"][:, :, :] = driver_locals.dry_air_mixing_ratio.graupel.field[:]
        outputs["driver_local_cloud_fraction_icecloud"][:, :, :] = driver_locals.cloud_fraction.field[:]
        outputs["driver_local_terminal_speed_snow_icecloud"][:, :, :] = driver_locals.terminal_speed.snow.field[:]
        outputs["driver_local_terminal_speed_graupel_icecloud"][:, :, :] = driver_locals.terminal_speed.graupel.field[:]
        outputs["driver_local_terminal_speed_rain_icecloud"][:, :, :] = driver_locals.terminal_speed.rain.field[:]
        outputs["driver_local_density_icecloud"][:, :, :] = driver_locals.density.field[:]
        outputs["driver_local_density_factor_icecloud"][:, :, :] = driver_locals.density_factor.field[:]
        outputs["driver_local_rh_limited_icecloud"][:, :, :] = driver_locals.rh_limited.field[:]
        outputs["non_anvil_large_scale_sublimation_icecloud"][:, :, :] = state.non_anvil_large_scale.sublimation.field[:]
        outputs["driver_local_ccn_icecloud"][:, :, :] = driver_locals.ccn.field[:]

        for k in range(nz):
            outputs["convection_fraction_icecloud"][:, :, k] = state.convection_fraction.field[:]
            outputs["surface_type_icecloud"][:, :, k] = state.surface_type.field[:]

        return outputs
