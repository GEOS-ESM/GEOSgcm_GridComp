import numpy as np
from f90nml import Namelist
from ndsl import StencilFactory
from ndsl.constants import I_DIM, J_DIM, K_DIM
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py

from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.GFDL_1M.driver.config_constants import GFDL1MDriverConfigDependentConstants
from pyMoist.GFDL_1M.driver.finish import update_tendencies
from pyMoist.GFDL_1M.driver.locals import GFDL1MDriverLocals
from pyMoist.GFDL_1M.locals import GFDL1MLocals
from pyMoist.GFDL_1M.state import GFDL1MState


class TranslateGFDL_1M_DriverFinish(TranslateFortranData2Py):
    def __init__(self, grid: Grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "driver_local_dry_mixing_ratio_vapor_unmodified_driverfinish": {},
            "driver_local_dry_mixing_ratio_liquid_unmodified_driverfinish": {},
            "driver_local_dry_mixing_ratio_rain_unmodified_driverfinish": {},
            "driver_local_dry_mixing_ratio_ice_unmodified_driverfinish": {},
            "driver_local_dry_mixing_ratio_graupel_unmodified_driverfinish": {},
            "driver_local_cloud_fraciton_unmodified_driverfinish": {},
            "driver_local_dry_mixing_ratio_vapor_driverfinish": {},
            "driver_local_dry_mixing_ratio_liquid_driverfinish": {},
            "driver_local_dry_mixing_ratio_rain_driverfinish": {},
            "driver_local_dry_mixing_ratio_ice_driverfinish": {},
            "driver_local_dry_mixing_ratio_snow_driverfinish": {},
            "driver_local_dry_mixing_ratio_graupel_driverfinish": {},
            "local_dvapordt_driver_driverfinish": {},
            "local_dliquiddt_driver_driverfinish": {},
            "local_draindt_driver_driverfinish": {},
            "local_dicedt_driver_driverfinish": {},
            "local_dsnowdt_driver_driverfinish": {},
            "local_dgraupeldt_driver_driverfinish": {},
            "local_dcloudfractiondt_driver_driverfinish": {},
            "driver_local_t_unmodified_driverfinish": {},
            "driver_local_t_driverfinish": {},
            "local_dtdt_driver_driverfinish": {},
            "w_driverfinish": {},
            "driver_local_w_driverfinish": {},
            "u_driverfinish": {},
            "driver_local_u_driverfinish": {},
            "local_dudt_driver_driverfinish": {},
            "v_driverfinish": {},
            "driver_local_v_driverfinish": {},
            "local_dvdt_driver_driverfinish": {},
            "local_dp_driverfinish": {},
            "driver_local_dp_driverfinish": {},
            "driver_local_mass_driverfinish": {},
            "surface_precip_rain_driverfinish": {},
            "surface_precip_snow_driverfinish": {},
            "surface_precip_ice_driverfinish": {},
            "surface_precip_graupel_driverfinish": {},
        }

        self.out_vars = self.in_vars["data_vars"].copy()

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GFDL_1M-constants")

    def compute(self, inputs):
        # initialize dataclasses
        state = GFDL1MState.zeros(self.quantity_factory)
        locals_ = GFDL1MLocals.make_as_state(self.quantity_factory)
        driver_locals = GFDL1MDriverLocals.make_as_state(self.quantity_factory)

        # initialize constants
        config = GFDL1MConfig(**self.constants)
        config_dependent_constants = GFDL1MDriverConfigDependentConstants.make(config)

        driver_locals.unmodified.mixing_ratio.vapor.field[:] = inputs[
            "driver_local_dry_mixing_ratio_vapor_unmodified_driverfinish"
        ][:, :, :, 0]
        driver_locals.unmodified.mixing_ratio.liquid.field[:] = inputs[
            "driver_local_dry_mixing_ratio_liquid_unmodified_driverfinish"
        ][:, :, :, 0]
        driver_locals.unmodified.mixing_ratio.rain.field[:] = inputs[
            "driver_local_dry_mixing_ratio_rain_unmodified_driverfinish"
        ][:, :, :, 0]
        driver_locals.unmodified.mixing_ratio.ice.field[:] = inputs[
            "driver_local_dry_mixing_ratio_ice_unmodified_driverfinish"
        ][:, :, :, 0]
        driver_locals.unmodified.mixing_ratio.rain.field[:] = inputs[
            "driver_local_dry_mixing_ratio_rain_unmodified_driverfinish"
        ][:, :, :, 0]
        driver_locals.unmodified.mixing_ratio.graupel.field[:] = inputs[
            "driver_local_dry_mixing_ratio_graupel_unmodified_driverfinish"
        ][:, :, :, 0]
        state.radiation_field.cloud_fraction.field[:] = inputs[
            "driver_local_cloud_fraciton_unmodified_driverfinish"
        ][:, :, :, 0]
        driver_locals.dry_air_mixing_ratio.vapor.field[:] = inputs[
            "driver_local_dry_mixing_ratio_vapor_driverfinish"
        ][:, :, :, 0]
        driver_locals.dry_air_mixing_ratio.liquid.field[:] = inputs[
            "driver_local_dry_mixing_ratio_liquid_driverfinish"
        ][:, :, :, 0]
        driver_locals.dry_air_mixing_ratio.rain.field[:] = inputs[
            "driver_local_dry_mixing_ratio_rain_driverfinish"
        ][:, :, :, 0]
        driver_locals.dry_air_mixing_ratio.ice.field[:] = inputs[
            "driver_local_dry_mixing_ratio_ice_driverfinish"
        ][:, :, :, 0]
        driver_locals.dry_air_mixing_ratio.snow.field[:] = inputs[
            "driver_local_dry_mixing_ratio_snow_driverfinish"
        ][:, :, :, 0]
        driver_locals.dry_air_mixing_ratio.graupel.field[:] = inputs[
            "driver_local_dry_mixing_ratio_graupel_driverfinish"
        ][:, :, :, 0]
        locals_.driver_tendencies.dvapordt.field[:] = inputs["local_dvapordt_driver_driverfinish"][:, :, :, 0]
        locals_.driver_tendencies.dliquiddt.field[:] = inputs["local_dliquiddt_driver_driverfinish"][
            :, :, :, 0
        ]
        locals_.driver_tendencies.draindt.field[:] = inputs["local_draindt_driver_driverfinish"][:, :, :, 0]
        locals_.driver_tendencies.dicedt.field[:] = inputs["local_dicedt_driver_driverfinish"][:, :, :, 0]
        locals_.driver_tendencies.dsnowdt.field[:] = inputs["local_dsnowdt_driver_driverfinish"][:, :, :, 0]
        locals_.driver_tendencies.dgraupeldt.field[:] = inputs["local_dgraupeldt_driver_driverfinish"][
            :, :, :, 0
        ]
        locals_.driver_tendencies.dcloudfractiondt.field[:] = inputs[
            "local_dcloudfractiondt_driver_driverfinish"
        ][:, :, :, 0]
        state.t.field[:] = inputs["driver_local_t_unmodified_driverfinish"][:, :, :, 0]
        driver_locals.t.field[:] = inputs["driver_local_t_driverfinish"][:, :, :, 0]
        locals_.driver_tendencies.dtdt.field[:] = inputs["local_dtdt_driver_driverfinish"][:, :, :, 0]
        state.vertical_motion.velocity.field[:] = inputs["w_driverfinish"][:, :, :, 0]
        driver_locals.w.field[:] = inputs["driver_local_w_driverfinish"][:, :, :, 0]
        state.u.field[:] = inputs["u_driverfinish"][:, :, :, 0]
        driver_locals.u.field[:] = inputs["driver_local_u_driverfinish"][:, :, :, 0]
        locals_.driver_tendencies.dudt.field[:] = inputs["local_dudt_driver_driverfinish"][:, :, :, 0]
        state.v.field[:] = inputs["v_driverfinish"][:, :, :, 0]
        driver_locals.v.field[:] = inputs["driver_local_v_driverfinish"][:, :, :, 0]
        locals_.driver_tendencies.dvdt.field[:] = inputs["local_dvdt_driver_driverfinish"][:, :, :, 0]
        locals_.dp.field[:] = inputs["local_dp_driverfinish"][:, :, :, 0]
        driver_locals.dp.field[:] = inputs["driver_local_dp_driverfinish"][:, :, :, 0]
        driver_locals.mass.field[:] = inputs["driver_local_mass_driverfinish"][:, :, :, 0]
        state.precipitation_at_surface.rain.field[:] = inputs["surface_precip_rain_driverfinish"][:, :, 0, 0]
        state.precipitation_at_surface.snow.field[:] = inputs["surface_precip_snow_driverfinish"][:, :, 0, 0]
        state.precipitation_at_surface.ice.field[:] = inputs["surface_precip_ice_driverfinish"][:, :, 0, 0]
        state.precipitation_at_surface.graupel.field[:] = inputs["surface_precip_graupel_driverfinish"][
            :, :, 0, 0
        ]

        # construct test stencil
        code = self.stencil_factory.from_dims_halo(
            func=update_tendencies,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "c_air": config_dependent_constants.C_AIR,
                "c_vap": config_dependent_constants.C_VAP,
                "rdt": config_dependent_constants.RDT,
                "do_sedi_w": config.DO_SEDI_W,
                "sedi_transport": config.SEDI_TRANSPORT,
                "do_qa": config.DO_QA,
            },
        )

        code(
            mixing_ratio_vapor_unmodified=driver_locals.unmodified.mixing_ratio.vapor,
            mixing_ratio_liquid_unmodified=driver_locals.unmodified.mixing_ratio.liquid,
            mixing_ratio_rain_unmodified=driver_locals.unmodified.mixing_ratio.rain,
            mixing_ratio_ice_unmodified=driver_locals.unmodified.mixing_ratio.ice,
            mixing_ratio_snow_unmodified=driver_locals.unmodified.mixing_ratio.rain,
            mixing_ratio_graupel_unmodified=driver_locals.unmodified.mixing_ratio.graupel,
            cloud_fraction_unmodified=state.radiation_field.cloud_fraction,
            mixing_ratio_driver_vapor=driver_locals.dry_air_mixing_ratio.vapor,
            mixing_ratio_driver_liquid=driver_locals.dry_air_mixing_ratio.liquid,
            mixing_ratio_driver_rain=driver_locals.dry_air_mixing_ratio.rain,
            mixing_ratio_driver_ice=driver_locals.dry_air_mixing_ratio.ice,
            mixing_ratio_driver_snow=driver_locals.dry_air_mixing_ratio.snow,
            mixing_ratio_driver_graupel=driver_locals.dry_air_mixing_ratio.graupel,
            dvapordt=locals_.driver_tendencies.dvapordt,
            dliquiddt=locals_.driver_tendencies.dliquiddt,
            draindt=locals_.driver_tendencies.draindt,
            dicedt=locals_.driver_tendencies.dicedt,
            dsnowdt=locals_.driver_tendencies.dsnowdt,
            dgraupeldt=locals_.driver_tendencies.dgraupeldt,
            dcloudfractiondt=locals_.driver_tendencies.dcloudfractiondt,
            t_unmodified=state.t,
            driver_t=driver_locals.t,
            dtdt=locals_.driver_tendencies.dtdt,
            w_unmodified=state.vertical_motion.velocity,
            driver_w=driver_locals.w,
            u_unmodified=state.u,
            driver_u=driver_locals.u,
            dudt=locals_.driver_tendencies.dudt,
            v_unmodified=state.v,
            driver_v=driver_locals.v,
            dvdt=locals_.driver_tendencies.dvdt,
            dp_unmodified=locals_.dp,
            driver_dp=driver_locals.dp,
            driver_mass=driver_locals.mass,
            rain=state.precipitation_at_surface.rain,
            snow=state.precipitation_at_surface.snow,
            ice=state.precipitation_at_surface.ice,
            graupel=state.precipitation_at_surface.graupel,
        )

        # get the shape of the field
        nx, ny, nz, ntimes = inputs["driver_local_dry_mixing_ratio_vapor_unmodified_driverfinish"].shape

        # prefill output array with nans
        outputs = {}
        for key in self.out_vars:
            outputs[key] = np.full((nx, ny, nz, ntimes), np.nan)

        outputs["driver_local_dry_mixing_ratio_vapor_unmodified_driverfinish"][:, :, :, 0] = (
            driver_locals.unmodified.mixing_ratio.vapor.field[:]
        )
        outputs["driver_local_dry_mixing_ratio_liquid_unmodified_driverfinish"][:, :, :, 0] = (
            driver_locals.unmodified.mixing_ratio.liquid.field[:]
        )
        outputs["driver_local_dry_mixing_ratio_rain_unmodified_driverfinish"][:, :, :, 0] = (
            driver_locals.unmodified.mixing_ratio.rain.field[:]
        )
        outputs["driver_local_dry_mixing_ratio_ice_unmodified_driverfinish"][:, :, :, 0] = (
            driver_locals.unmodified.mixing_ratio.ice.field[:]
        )
        outputs["driver_local_dry_mixing_ratio_rain_unmodified_driverfinish"][:, :, :, 0] = (
            driver_locals.unmodified.mixing_ratio.rain.field[:]
        )
        outputs["driver_local_dry_mixing_ratio_graupel_unmodified_driverfinish"][:, :, :, 0] = (
            driver_locals.unmodified.mixing_ratio.graupel.field[:]
        )
        outputs["driver_local_cloud_fraciton_unmodified_driverfinish"][:, :, :, 0] = (
            state.radiation_field.cloud_fraction.field[:]
        )
        outputs["driver_local_dry_mixing_ratio_vapor_driverfinish"][:, :, :, 0] = (
            driver_locals.dry_air_mixing_ratio.vapor.field[:]
        )
        outputs["driver_local_dry_mixing_ratio_liquid_driverfinish"][:, :, :, 0] = (
            driver_locals.dry_air_mixing_ratio.liquid.field[:]
        )
        outputs["driver_local_dry_mixing_ratio_rain_driverfinish"][:, :, :, 0] = (
            driver_locals.dry_air_mixing_ratio.rain.field[:]
        )
        outputs["driver_local_dry_mixing_ratio_ice_driverfinish"][:, :, :, 0] = (
            driver_locals.dry_air_mixing_ratio.ice.field[:]
        )
        outputs["driver_local_dry_mixing_ratio_snow_driverfinish"][:, :, :, 0] = (
            driver_locals.dry_air_mixing_ratio.snow.field[:]
        )
        outputs["driver_local_dry_mixing_ratio_graupel_driverfinish"][:, :, :, 0] = (
            driver_locals.dry_air_mixing_ratio.graupel.field[:]
        )
        outputs["local_dvapordt_driver_driverfinish"][:, :, :, 0] = locals_.driver_tendencies.dvapordt.field[
            :
        ]
        outputs["local_dliquiddt_driver_driverfinish"][:, :, :, 0] = (
            locals_.driver_tendencies.dliquiddt.field[:]
        )
        outputs["local_draindt_driver_driverfinish"][:, :, :, 0] = locals_.driver_tendencies.draindt.field[:]
        outputs["local_dicedt_driver_driverfinish"][:, :, :, 0] = locals_.driver_tendencies.dicedt.field[:]
        outputs["local_dsnowdt_driver_driverfinish"][:, :, :, 0] = locals_.driver_tendencies.dsnowdt.field[:]
        outputs["local_dgraupeldt_driver_driverfinish"][:, :, :, 0] = (
            locals_.driver_tendencies.dgraupeldt.field[:]
        )
        outputs["local_dcloudfractiondt_driver_driverfinish"][:, :, :, 0] = (
            locals_.driver_tendencies.dcloudfractiondt.field[:]
        )
        outputs["driver_local_t_unmodified_driverfinish"][:, :, :, 0] = state.t.field[:]
        outputs["driver_local_t_driverfinish"][:, :, :, 0] = driver_locals.t.field[:]
        outputs["local_dtdt_driver_driverfinish"][:, :, :, 0] = locals_.driver_tendencies.dtdt.field[:]
        outputs["w_driverfinish"][:, :, :, 0] = state.vertical_motion.velocity.field[:]
        outputs["driver_local_w_driverfinish"][:, :, :, 0] = driver_locals.w.field[:]
        outputs["u_driverfinish"][:, :, :, 0] = state.u.field[:]
        outputs["driver_local_u_driverfinish"][:, :, :, 0] = driver_locals.u.field[:]
        outputs["local_dudt_driver_driverfinish"][:, :, :, 0] = locals_.driver_tendencies.dudt.field[:]
        outputs["v_driverfinish"][:, :, :, 0] = state.v.field[:]
        outputs["driver_local_v_driverfinish"][:, :, :, 0] = driver_locals.v.field[:]
        outputs["local_dvdt_driver_driverfinish"][:, :, :, 0] = locals_.driver_tendencies.dvdt.field[:]
        outputs["local_dp_driverfinish"][:, :, :, 0] = locals_.dp.field[:]
        outputs["driver_local_dp_driverfinish"][:, :, :, 0] = driver_locals.dp.field[:]
        outputs["driver_local_mass_driverfinish"][:, :, :, 0] = driver_locals.mass.field[:]
        outputs["surface_precip_rain_driverfinish"][:, :, 0, 0] = state.precipitation_at_surface.rain.field[:]
        outputs["surface_precip_snow_driverfinish"][:, :, 0, 0] = state.precipitation_at_surface.snow.field[:]
        outputs["surface_precip_ice_driverfinish"][:, :, 0, 0] = state.precipitation_at_surface.ice.field[:]
        outputs["surface_precip_graupel_driverfinish"][:, :, 0, 0] = (
            state.precipitation_at_surface.graupel.field[:]
        )

        return outputs
