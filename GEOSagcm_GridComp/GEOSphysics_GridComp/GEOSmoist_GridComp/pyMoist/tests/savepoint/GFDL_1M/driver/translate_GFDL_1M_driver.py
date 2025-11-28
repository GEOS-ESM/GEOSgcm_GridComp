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
from pyMoist.GFDL_1M.driver.driver import GFDL1MDriver
from ndsl.stencils.testing.savepoint import DataLoader
import numpy as np


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

        # fill relevant parts of dataclasses with input data
        state.radiation_field.vapor.field[:] = inputs["radiation_vapor"]
        state.radiation_field.liquid.field[:] = inputs["radiation_liquid"]
        state.radiation_field.rain.field[:] = inputs["radiation_rain"]
        state.radiation_field.ice.field[:] = inputs["radiation_ice"]
        state.radiation_field.snow.field[:] = inputs["radiation_snow"]
        state.radiation_field.graupel.field[:] = inputs["radiation_graupel"]
        state.radiation_field.cloud_fraction.field[:] = inputs["radiation_cloud_fraction"]
        locals.total_concentration.field[:] = inputs["local_total_concentration"]
        locals.driver_tencencies.dvapordt.field[:] = inputs["local_dvapordt_driver"]
        locals.driver_tencencies.dliquiddt.field[:] = inputs["local_dliquiddt_driver"]
        locals.driver_tencencies.draindt.field[:] = inputs["local_draindt_driver"]
        locals.driver_tencencies.dicedt.field[:] = inputs["local_dicedt_driver"]
        locals.driver_tencencies.dsnowdt.field[:] = inputs["local_dsnowdt_driver"]
        locals.driver_tencencies.dgraupeldt.field[:] = inputs["local_dgraupeldt_driver"]
        locals.driver_tencencies.dcloudfractiondt.field[:] = inputs["local_dcloudfractiondt_driver"]
        locals.driver_tencencies.dtdt.field[:] = inputs["local_dtdt_driver"]
        locals.driver_tencencies.dudt.field[:] = inputs["local_dudt_driver"]
        locals.driver_tencencies.dvdt.field[:] = inputs["local_dvdt_driver"]
        state.t.field[:] = inputs["t"]
        state.u.field[:] = inputs["u"]
        state.v.field[:] = inputs["v"]
        state.vertical_motion.velocity.field[:] = inputs["w"]
        locals.dz.field[:] = inputs["local_dz"]
        locals.dp.field[:] = inputs["local_dp"]
        state.area.field[:] = inputs["area"]
        state.land_fraction.field[:] = inputs["land_fraction"]
        state.convection_fraction.field[:] = inputs["convection_fraction"]
        state.surface_type.field[:] = inputs["surface_type"]
        state.estimated_inversion_strength.field[:] = inputs["estimated_inversion_strength"]
        state.critical_relative_humidity_for_pdf.field[:] = inputs["critical_relative_humidity_for_pdf"]
        state.non_anvil_large_scale.evaporation.field[:] = inputs["non_anvil_large_scale_evaporation"]
        state.non_anvil_large_scale.sublimation.field[:] = inputs["non_anvil_large_scale_sublimation"]
        state.non_anvil_large_scale.liquid_precip_flux.field[:] = inputs[
            "non_anvil_large_scale_liquid_precip_flux"
        ]
        state.non_anvil_large_scale.ice_precip_flux.field[:] = inputs["non_anvil_large_scale_ice_precip_flux"]
        state.precipitation_at_surface.rain.field[:] = inputs["surface_precip_rain"]
        state.precipitation_at_surface.snow.field[:] = inputs["surface_precip_snow"]
        state.precipitation_at_surface.ice.field[:] = inputs["surface_precip_ice"]
        state.precipitation_at_surface.graupel.field[:] = inputs["surface_precip_graupel"]

        # construct test stencil
        code = GFDL1MDriver(
            stencil_factory=self.stencil_factory, quantity_factory=self.quantity_factory, config=config
        )

        code(
            t=state.t,
            u=state.u,
            v=state.v,
            w=state.vertical_motion.velocity,
            dz=locals.dz,
            dp=locals.dp,
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
            total_concentration=locals.total_concentration,
            dvapordt=locals.driver_tencencies.dvapordt,
            dliquiddt=locals.driver_tencencies.dliquiddt,
            draindt=locals.driver_tencencies.draindt,
            dicedt=locals.driver_tencencies.dicedt,
            dsnowdt=locals.driver_tencencies.dsnowdt,
            dgraupeldt=locals.driver_tencencies.dgraupeldt,
            dcloudfractiondt=locals.driver_tencencies.dcloudfractiondt,
            dtdt=locals.driver_tencencies.dtdt,
            dudt=locals.driver_tencencies.dudt,
            dvdt=locals.driver_tencencies.dvdt,
            liquid_precip_flux=state.non_anvil_large_scale.liquid_precip_flux,
            ice_precip_flux=state.non_anvil_large_scale.ice_precip_flux,
            evaporation=state.non_anvil_large_scale.evaporation,
            sublimation=state.non_anvil_large_scale.sublimation,
            surface_precip_rain=state.precipitation_at_surface.rain,
            surface_precip_snow=state.precipitation_at_surface.snow,
            surface_precip_ice=state.precipitation_at_surface.ice,
            surface_precip_graupel=state.precipitation_at_surface.graupel,
        )

        return {
            "radiation_vapor": state.radiation_field.vapor.field[:],
            "radiation_liquid": state.radiation_field.liquid.field[:],
            "radiation_rain": state.radiation_field.rain.field[:],
            "radiation_ice": state.radiation_field.ice.field[:],
            "radiation_snow": state.radiation_field.snow.field[:],
            "radiation_graupel": state.radiation_field.graupel.field[:],
            "radiation_cloud_fraction": state.radiation_field.cloud_fraction.field[:],
            "local_total_concentration": locals.total_concentration.field[:],
            "local_dvapordt_driver": locals.driver_tencencies.dvapordt.field[:],
            "local_dliquiddt_driver": locals.driver_tencencies.dliquiddt.field[:],
            "local_draindt_driver": locals.driver_tencencies.draindt.field[:],
            "local_dicedt_driver": locals.driver_tencencies.dicedt.field[:],
            "local_dsnowdt_driver": locals.driver_tencencies.dsnowdt.field[:],
            "local_dgraupeldt_driver": locals.driver_tencencies.dgraupeldt.field[:],
            "local_dcloudfractiondt_driver": locals.driver_tencencies.dcloudfractiondt.field[:],
            "local_dtdt_driver": locals.driver_tencencies.dtdt.field[:],
            "local_dudt_driver": locals.driver_tencencies.dudt.field[:],
            "local_dvdt_driver": locals.driver_tencencies.dvdt.field[:],
            "t": state.t.field[:],
            "u": state.u.field[:],
            "v": state.v.field[:],
            "w": state.vertical_motion.velocity.field[:],
            "local_dz": locals.dz.field[:],
            "local_dp": locals.dp.field[:],
            "area": state.area.field[:],
            "land_fraction": state.land_fraction.field[:],
            "convection_fraction": state.convection_fraction.field[:],
            "surface_type": state.surface_type.field[:],
            "estimated_inversion_strength": state.estimated_inversion_strength.field[:],
            "critical_relative_humidity_for_pdf": state.critical_relative_humidity_for_pdf.field[:],
            "non_anvil_large_scale_evaporation": state.non_anvil_large_scale.evaporation.field[:],
            "non_anvil_large_scale_sublimation": state.non_anvil_large_scale.sublimation.field[:],
            "non_anvil_large_scale_liquid_precip_flux": state.non_anvil_large_scale.liquid_precip_flux.field[
                :
            ],
            "non_anvil_large_scale_ice_precip_flux": state.non_anvil_large_scale.ice_precip_flux.field[:],
            "surface_precip_rain": state.precipitation_at_surface.rain.field[:],
            "surface_precip_snow": state.precipitation_at_surface.snow.field[:],
            "surface_precip_ice": state.precipitation_at_surface.ice.field[:],
            "surface_precip_graupel": state.precipitation_at_surface.graupel.field[:],
        }
