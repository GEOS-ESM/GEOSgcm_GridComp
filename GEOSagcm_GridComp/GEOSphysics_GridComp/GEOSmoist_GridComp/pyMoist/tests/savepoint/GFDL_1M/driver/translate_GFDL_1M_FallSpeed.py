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
from pyMoist.GFDL_1M.driver.fall_speed import fall_speed
from ndsl.stencils.testing.savepoint import DataLoader
from pyMoist.redistribute_clouds import redistribute_clouds
import numpy as np


class TranslateGFDL_1M_FallSpeed(TranslateFortranData2Py):
    def __init__(self, grid: Grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "driver_local_p_dry_fallspeed": {},
            "driver_local_density_fallspeed": {},
            "driver_local_dry_mixing_ratio_snow_fallspeed": {},
            "driver_local_dry_mixing_ratio_ice_fallspeed": {},
            "driver_local_dry_mixing_ratio_graupel_fallspeed": {},
            "driver_local_dry_mixing_ratio_liquid_fallspeed": {},
            "driver_local_terminal_speed_ice_fallspeed": {},
            "driver_local_terminal_speed_snow_fallspeed": {},
            "driver_local_terminal_speed_graupel_fallspeed": {},
            "driver_local_t_fallspeed": {},
            "driver_local_dz_unmodified_fallspeed": {},
            "driver_local_dz_fallspeed": {},
            "driver_local_density_unmodified_fallspeed": {},
            "driver_local_density_factor_fallspeed": {},
            "driver_local_t_unmodified_fallspeed": {},
            "convection_fraction_fallspeed": {},
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
        nx, ny, nz, ntimes = inputs["driver_local_p_dry_fallspeed"].shape

        # preset output dictionary to be filled inside the for loop
        outputs = {}
        for key in self.out_vars:
            outputs[key] = np.zeros((nx, ny, nz, ntimes), dtype=inputs[key].dtype)

        # construct test stencil
        code = self.stencil_factory.from_dims_halo(
            func=fall_speed,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "p_nonhydro": config_dependent_constants.P_NONHYDRO,
                "const_vi": config.CONST_VI,
                "const_vs": config.CONST_VS,
                "const_vg": config.CONST_VG,
                "vi_fac": config.VI_FAC,
                "vi_max": config.VI_MAX,
                "vs_fac": config.VS_FAC,
                "vs_max": config.VS_MAX,
                "vg_fac": config.VG_FAC,
                "vg_max": config.VG_MAX,
                "anv_icefall": config.ANV_ICEFALL,
                "ls_icefall": config.LS_ICEFALL,
            },
        )

        for n in range(ntimes):

            driver_locals.p_dry.field[:] = inputs["driver_local_p_dry_fallspeed"][:, :, :, n]
            driver_locals.density.field[:] = inputs["driver_local_density_fallspeed"][:, :, :, n]
            driver_locals.dry_air_mixing_ratio.snow.field[:] = inputs[
                "driver_local_dry_mixing_ratio_snow_fallspeed"
            ][:, :, :, n]
            driver_locals.dry_air_mixing_ratio.ice.field[:] = inputs[
                "driver_local_dry_mixing_ratio_ice_fallspeed"
            ][:, :, :, n]
            driver_locals.dry_air_mixing_ratio.graupel.field[:] = inputs[
                "driver_local_dry_mixing_ratio_graupel_fallspeed"
            ][:, :, :, n]
            driver_locals.dry_air_mixing_ratio.liquid.field[:] = inputs[
                "driver_local_dry_mixing_ratio_liquid_fallspeed"
            ][:, :, :, n]
            driver_locals.terminal_speed.ice.field[:] = inputs["driver_local_terminal_speed_ice_fallspeed"][
                :, :, :, n
            ]
            driver_locals.terminal_speed.snow.field[:] = inputs["driver_local_terminal_speed_snow_fallspeed"][
                :, :, :, n
            ]
            driver_locals.terminal_speed.graupel.field[:] = inputs[
                "driver_local_terminal_speed_graupel_fallspeed"
            ][:, :, :, n]
            driver_locals.t.field[:] = inputs["driver_local_t_fallspeed"][:, :, :, n]
            driver_locals.unmodified.dz.field[:] = inputs["driver_local_dz_unmodified_fallspeed"][:, :, :, n]
            driver_locals.dz.field[:] = inputs["driver_local_dz_fallspeed"][:, :, :, n]
            driver_locals.unmodified.density.field[:] = inputs["driver_local_density_unmodified_fallspeed"][
                :, :, :, n
            ]
            driver_locals.density_factor.field[:] = inputs["driver_local_density_factor_fallspeed"][
                :, :, :, n
            ]
            driver_locals.unmodified.t.field[:] = inputs["driver_local_t_unmodified_fallspeed"][:, :, :, n]
            state.convection_fraction.field[:] = inputs["convection_fraction_fallspeed"][:, :, 0, n]

            code(
                liquid=driver_locals.dry_air_mixing_ratio.liquid,
                ice=driver_locals.dry_air_mixing_ratio.ice,
                snow=driver_locals.dry_air_mixing_ratio.snow,
                graupel=driver_locals.dry_air_mixing_ratio.graupel,
                t_unmodified=driver_locals.unmodified.t,
                t=driver_locals.t,
                dz_unmodified=driver_locals.unmodified.dz,
                dz=driver_locals.dz,
                density_unmodified=driver_locals.unmodified.density,
                density=driver_locals.density,
                density_factor=driver_locals.density_factor,
                ice_terminal_velocity=driver_locals.terminal_speed.ice,
                snow_terminal_velocity=driver_locals.terminal_speed.snow,
                graupel_terminal_velosity=driver_locals.terminal_speed.graupel,
                convection_fraction=state.convection_fraction,
            )

            # fill the output arrays so that all calls are tested
            outputs["driver_local_p_dry_fallspeed"][:, :, :, n] = driver_locals.p_dry.field[:]
            outputs["driver_local_density_fallspeed"][:, :, :, n] = driver_locals.density.field[:]
            outputs["driver_local_dry_mixing_ratio_snow_fallspeed"][:, :, :, n] = (
                driver_locals.dry_air_mixing_ratio.snow.field[:]
            )
            outputs["driver_local_dry_mixing_ratio_ice_fallspeed"][:, :, :, n] = (
                driver_locals.dry_air_mixing_ratio.ice.field[:]
            )
            outputs["driver_local_dry_mixing_ratio_graupel_fallspeed"][:, :, :, n] = (
                driver_locals.dry_air_mixing_ratio.graupel.field[:]
            )
            outputs["driver_local_dry_mixing_ratio_liquid_fallspeed"][:, :, :, n] = (
                driver_locals.dry_air_mixing_ratio.liquid.field[:]
            )
            outputs["driver_local_terminal_speed_ice_fallspeed"][:, :, :, n] = (
                driver_locals.terminal_speed.ice.field[:]
            )
            outputs["driver_local_terminal_speed_snow_fallspeed"][:, :, :, n] = (
                driver_locals.terminal_speed.snow.field[:]
            )
            outputs["driver_local_terminal_speed_graupel_fallspeed"][:, :, :, n] = (
                driver_locals.terminal_speed.graupel.field[:]
            )
            outputs["driver_local_t_fallspeed"][:, :, :, n] = driver_locals.t.field[:]
            outputs["driver_local_dz_unmodified_fallspeed"][:, :, :, n] = driver_locals.unmodified.dz.field[:]
            outputs["driver_local_dz_fallspeed"][:, :, :, n] = driver_locals.dz.field[:]
            outputs["driver_local_density_unmodified_fallspeed"][:, :, :, n] = (
                driver_locals.unmodified.density.field[:]
            )
            outputs["driver_local_density_factor_fallspeed"][:, :, :, n] = driver_locals.density_factor.field[
                :
            ]
            outputs["driver_local_t_unmodified_fallspeed"][:, :, :, n] = driver_locals.unmodified.t.field[:]

            for k in range(nz):
                outputs["convection_fraction_fallspeed"][:, :, k, n] = state.convection_fraction.field[:]

        return outputs
