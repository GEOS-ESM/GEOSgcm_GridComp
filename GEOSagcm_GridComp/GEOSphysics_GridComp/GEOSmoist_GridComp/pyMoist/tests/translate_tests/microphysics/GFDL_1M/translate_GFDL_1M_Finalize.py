from f90nml import Namelist
from ndsl import StencilFactory
from ndsl.constants import I_DIM, J_DIM, K_DIM
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py

from pyMoist.microphysics.GFDL_1M.config import GFDL1MConfig
from pyMoist.microphysics.GFDL_1M.finalize import GFDL1MFinalize
from pyMoist.microphysics.GFDL_1M.locals import GFDL1MLocals
from pyMoist.microphysics.GFDL_1M.shared_stencils import update_tendencies
from pyMoist.microphysics.GFDL_1M.state import GFDL1MState
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable


class TranslateGFDL_1M_Finalize(TranslateFortranData2Py):
    def __init__(self, grid: Grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "t": {},
            "u": {},
            "v": {},
            "local_p_mb": {},
            "radiation_cloud_fraction": {},
            "radiation_ice": {},
            "radiation_liquid": {},
            "local_u_unmodified": {},
            "local_v_unmodified": {},
            "radiation_vapor": {},
            "radiation_rain": {},
            "radiation_snow": {},
            "radiation_graupel": {},
            "local_mass": {},
            "cloud_fraction_convective": {},
            "cloud_fraction_large_scale": {},
            "mixing_ratio_convective_liquid": {},
            "mixing_ratio_convective_ice": {},
            "mixing_ratio_large_scale_ice": {},
            "mixing_ratio_large_scale_liquid": {},
            "concentration_ice": {},
            "concentration_liquid": {},
            "mixing_ratio_vapor": {},
            "mixing_ratio_rain": {},
            "mixing_ratio_snow": {},
            "mixing_ratio_graupel": {},
            "non_anvil_large_scale_evaporation": {},
            "non_anvil_large_scale_sublimation": {},
            "surface_precip_rain": {},
            "surface_precip_snow": {},
            "surface_precip_ice": {},
            "surface_precip_graupel": {},
            "non_anvil_large_scale_precip": {},
            "non_anvil_large_scale_snow": {},
            "icefall": {},
            "freezing_rainfall": {},
            "cloud_particle_effective_radius_liquid": {},
            "cloud_particle_effective_radius_ice": {},
            "non_anvil_large_scale_liquid_precip_flux": {},
            "non_anvil_large_scale_ice_precip_flux": {},
            "anvil_liquid_precip_flux": {},
            "anvil_ice_precip_flux": {},
            "relative_humidity_after_pdf": {},
            "dvapordt_micro": {},
            "dliquiddt_micro": {},
            "dicedt_micro": {},
            "dcloud_fractiondt_micro": {},
            "draindt_micro": {},
            "dsnowdt_micro": {},
            "dgraupeldt_micro": {},
            "dudt_micro": {},
            "dvdt_micro": {},
            "dtdt_micro": {},
        }

        self.out_vars = self.in_vars["data_vars"].copy()

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GFDL_1M-constants")

    def compute(self, inputs):
        # initialize constants
        config = GFDL1MConfig(**self.constants)

        # initialize dataclasses
        state = GFDL1MState.zeros(self.quantity_factory)
        locals_ = GFDL1MLocals.make_as_state(self.quantity_factory)

        # Initialize saturation tables
        saturation_tables = SaturationVaporPressureTable(self.stencil_factory.backend)

        state.t.field[:] = inputs["t"]
        state.u.field[:] = inputs["u"]
        state.v.field[:] = inputs["v"]
        locals_.p_mb.field[:] = inputs["local_p_mb"]
        state.radiation_field.cloud_fraction.field[:] = inputs["radiation_cloud_fraction"]
        state.radiation_field.ice.field[:] = inputs["radiation_ice"]
        state.radiation_field.liquid.field[:] = inputs["radiation_liquid"]
        locals_.u_unmodified.field[:] = inputs["local_u_unmodified"]
        locals_.v_unmodified.field[:] = inputs["local_v_unmodified"]
        state.radiation_field.vapor.field[:] = inputs["radiation_vapor"]
        state.radiation_field.rain.field[:] = inputs["radiation_rain"]
        state.radiation_field.snow.field[:] = inputs["radiation_snow"]
        state.radiation_field.graupel.field[:] = inputs["radiation_graupel"]
        locals_.mass.field[:] = inputs["local_mass"]
        state.cloud_fraction.convective.field[:] = inputs["cloud_fraction_convective"]
        state.cloud_fraction.large_scale.field[:] = inputs["cloud_fraction_large_scale"]
        state.mixing_ratio.convective_liquid.field[:] = inputs["mixing_ratio_convective_liquid"]
        state.mixing_ratio.convective_ice.field[:] = inputs["mixing_ratio_convective_ice"]
        state.mixing_ratio.large_scale_ice.field[:] = inputs["mixing_ratio_large_scale_ice"]
        state.mixing_ratio.large_scale_liquid.field[:] = inputs["mixing_ratio_large_scale_liquid"]
        state.concentration.ice.field[:] = inputs["concentration_ice"]
        state.concentration.liquid.field[:] = inputs["concentration_liquid"]
        state.mixing_ratio.vapor.field[:] = inputs["mixing_ratio_vapor"]
        state.mixing_ratio.rain.field[:] = inputs["mixing_ratio_rain"]
        state.mixing_ratio.snow.field[:] = inputs["mixing_ratio_snow"]
        state.mixing_ratio.graupel.field[:] = inputs["mixing_ratio_graupel"]
        state.non_anvil_large_scale.evaporation.field[:] = inputs["non_anvil_large_scale_evaporation"]
        state.non_anvil_large_scale.sublimation.field[:] = inputs["non_anvil_large_scale_sublimation"]
        state.precipitation_at_surface.rain.field[:] = inputs["surface_precip_rain"]
        state.precipitation_at_surface.snow.field[:] = inputs["surface_precip_snow"]
        state.precipitation_at_surface.ice.field[:] = inputs["surface_precip_ice"]
        state.precipitation_at_surface.graupel.field[:] = inputs["surface_precip_graupel"]
        state.non_anvil_large_scale.precip.field[:] = inputs["non_anvil_large_scale_precip"]
        state.non_anvil_large_scale.snow.field[:] = inputs["non_anvil_large_scale_snow"]
        state.icefall.field[:] = inputs["icefall"]
        state.freezing_rainfall.field[:] = inputs["freezing_rainfall"]
        state.cloud_particle_effective_radius.liquid.field[:] = inputs["cloud_particle_effective_radius_liquid"]
        state.cloud_particle_effective_radius.ice.field[:] = inputs["cloud_particle_effective_radius_ice"]
        state.non_anvil_large_scale.liquid_precip_flux.field[:] = inputs["non_anvil_large_scale_liquid_precip_flux"]
        state.non_anvil_large_scale.ice_precip_flux.field[:] = inputs["non_anvil_large_scale_ice_precip_flux"]
        state.anvil.liquid_precip_flux.field[:] = inputs["anvil_liquid_precip_flux"]
        state.anvil.ice_precip_flux.field[:] = inputs["anvil_ice_precip_flux"]
        state.relative_humidity_after_pdf.field[:] = inputs["relative_humidity_after_pdf"]
        state.tendencies.dvapordt_micro.field[:] = inputs["dvapordt_micro"]
        state.tendencies.dliquiddt_micro.field[:] = inputs["dliquiddt_micro"]
        state.tendencies.dicedt_micro.field[:] = inputs["dicedt_micro"]
        state.tendencies.dcloud_fractiondt_micro.field[:] = inputs["dcloud_fractiondt_micro"]
        state.tendencies.draindt_micro.field[:] = inputs["draindt_micro"]
        state.tendencies.dsnowdt_micro.field[:] = inputs["dsnowdt_micro"]
        state.tendencies.dgraupeldt_micro.field[:] = inputs["dgraupeldt_micro"]
        state.tendencies.dudt_micro.field[:] = inputs["dudt_micro"]
        state.tendencies.dvdt_micro.field[:] = inputs["dvdt_micro"]
        state.tendencies.dtdt_micro.field[:] = inputs["dtdt_micro"]

        # construct test stencil
        _update_tendencies = self.stencil_factory.from_dims_halo(
            func=update_tendencies,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={"DT_MOIST": config.DT_MOIST},
        )

        code = GFDL1MFinalize(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            saturation_tables=saturation_tables,
            update_tendencies=_update_tendencies,
        )
        code(
            t=state.t,
            u=state.u,
            v=state.v,
            mixing_ratio_vapor=state.mixing_ratio.vapor,
            mixing_ratio_convective_liquid=state.mixing_ratio.convective_liquid,
            mixing_ratio_large_scale_liquid=state.mixing_ratio.large_scale_liquid,
            mixing_ratio_convective_ice=state.mixing_ratio.convective_ice,
            mixing_ratio_large_scale_ice=state.mixing_ratio.large_scale_ice,
            mixing_ratio_rain=state.mixing_ratio.rain,
            mixing_ratio_snow=state.mixing_ratio.snow,
            mixing_ratio_graupel=state.mixing_ratio.graupel,
            cloud_fraction_convective=state.cloud_fraction.convective,
            cloud_fraction_large_scale=state.cloud_fraction.large_scale,
            non_anvil_large_scale_precip=state.non_anvil_large_scale.precip,
            non_anvil_large_scale_snow=state.non_anvil_large_scale.snow,
            non_anvil_large_scale_ice_precip_flux=state.non_anvil_large_scale.ice_precip_flux,
            non_anvil_large_scale_liquid_precip_flux=state.non_anvil_large_scale.liquid_precip_flux,
            anvil_liquid_precip_flux=state.anvil.liquid_precip_flux,
            anvil_ice_precip_flux=state.anvil.ice_precip_flux,
            surface_rain=state.precipitation_at_surface.rain,
            surface_snow=state.precipitation_at_surface.snow,
            surface_ice=state.precipitation_at_surface.ice,
            surface_graupel=state.precipitation_at_surface.graupel,
            icefall=state.icefall,
            freezing_rainfall=state.freezing_rainfall,
            concentration_liquid=state.concentration.liquid,
            concentration_ice=state.concentration.ice,
            cloud_particle_effective_radius_liquid=state.cloud_particle_effective_radius.liquid,
            cloud_particle_effective_radius_ice=state.cloud_particle_effective_radius.ice,
            relative_humidity_after_pdf=state.relative_humidity_after_pdf,
            large_scale_rainwater_source=state.large_scale_rainwater_source,
            radiation_vapor=state.radiation_field.vapor,
            radiation_liquid=state.radiation_field.liquid,
            radiation_rain=state.radiation_field.rain,
            radiation_snow=state.radiation_field.snow,
            radiation_graupel=state.radiation_field.graupel,
            radiation_ice=state.radiation_field.ice,
            radiation_cloud_fraction=state.radiation_field.cloud_fraction,
            dudt_micro=state.tendencies.dudt_micro,
            dvdt_micro=state.tendencies.dvdt_micro,
            dtdt_micro=state.tendencies.dtdt_micro,
            dvapordt_micro=state.tendencies.dvapordt_micro,
            dliquiddt_micro=state.tendencies.dliquiddt_micro,
            dicedt_micro=state.tendencies.dicedt_micro,
            dcloud_fractiondt_micro=state.tendencies.dcloud_fractiondt_micro,
            draindt_micro=state.tendencies.draindt_micro,
            dsnowdt_micro=state.tendencies.dsnowdt_micro,
            dgraupeldt_micro=state.tendencies.dgraupeldt_micro,
            dudt_macro=state.tendencies.dudt_macro,
            dvdt_macro=state.tendencies.dvdt_macro,
            draindt_macro=state.tendencies.draindt_macro,
            dtdt_friction_pressure_weighted=state.tendencies.dtdt_friction_pressure_weighted,
            local_p_mb=locals_.p_mb,
            local_mass=locals_.mass,
            local_u_unmodified=locals_.u_unmodified,
            local_v_unmodified=locals_.v_unmodified,
            simulated_reflectivity=None,
            maximum_composite_reflectivity=None,
            base_1km_agl_reflectivity=None,
            echo_top_reflectivity=None,
            minus_10c_reflectivity=None,
            mass_fraction_suspended_rain=None,
            mass_fraction_suspended_snow=None,
            mass_fraction_suspended_graupel=None,
        )

        return {
            "t": state.t.field[:],
            "u": state.u.field[:],
            "v": state.v.field[:],
            "local_p_mb": locals_.p_mb.field[:],
            "radiation_cloud_fraction": state.radiation_field.cloud_fraction.field[:],
            "radiation_ice": state.radiation_field.ice.field[:],
            "radiation_liquid": state.radiation_field.liquid.field[:],
            "local_u_unmodified": locals_.u_unmodified.field[:],
            "local_v_unmodified": locals_.v_unmodified.field[:],
            "radiation_vapor": state.radiation_field.vapor.field[:],
            "radiation_rain": state.radiation_field.rain.field[:],
            "radiation_snow": state.radiation_field.snow.field[:],
            "radiation_graupel": state.radiation_field.graupel.field[:],
            "local_mass": locals_.mass.field[:],
            "cloud_fraction_convective": state.cloud_fraction.convective.field[:],
            "cloud_fraction_large_scale": state.cloud_fraction.large_scale.field[:],
            "mixing_ratio_convective_liquid": state.mixing_ratio.convective_liquid.field[:],
            "mixing_ratio_convective_ice": state.mixing_ratio.convective_ice.field[:],
            "mixing_ratio_large_scale_ice": state.mixing_ratio.large_scale_ice.field[:],
            "mixing_ratio_large_scale_liquid": state.mixing_ratio.large_scale_liquid.field[:],
            "concentration_ice": state.concentration.ice.field[:],
            "concentration_liquid": state.concentration.liquid.field[:],
            "mixing_ratio_vapor": state.mixing_ratio.vapor.field[:],
            "mixing_ratio_rain": state.mixing_ratio.rain.field[:],
            "mixing_ratio_snow": state.mixing_ratio.snow.field[:],
            "mixing_ratio_graupel": state.mixing_ratio.graupel.field[:],
            "non_anvil_large_scale_evaporation": state.non_anvil_large_scale.evaporation.field[:],
            "non_anvil_large_scale_sublimation": state.non_anvil_large_scale.sublimation.field[:],
            "surface_precip_rain": state.precipitation_at_surface.rain.field[:],
            "surface_precip_snow": state.precipitation_at_surface.snow.field[:],
            "surface_precip_ice": state.precipitation_at_surface.ice.field[:],
            "surface_precip_graupel": state.precipitation_at_surface.graupel.field[:],
            "non_anvil_large_scale_precip": state.non_anvil_large_scale.precip.field[:],
            "non_anvil_large_scale_snow": state.non_anvil_large_scale.snow.field[:],
            "icefall": state.icefall.field[:],
            "freezing_rainfall": state.freezing_rainfall.field[:],
            "cloud_particle_effective_radius_liquid": state.cloud_particle_effective_radius.liquid.field[:],
            "cloud_particle_effective_radius_ice": state.cloud_particle_effective_radius.ice.field[:],
            "non_anvil_large_scale_liquid_precip_flux": state.non_anvil_large_scale.liquid_precip_flux.field[:],
            "non_anvil_large_scale_ice_precip_flux": state.non_anvil_large_scale.ice_precip_flux.field[:],
            "anvil_liquid_precip_flux": state.anvil.liquid_precip_flux.field[:],
            "anvil_ice_precip_flux": state.anvil.ice_precip_flux.field[:],
            "relative_humidity_after_pdf": state.relative_humidity_after_pdf.field[:],
            "dvapordt_micro": state.tendencies.dvapordt_micro.field[:],
            "dliquiddt_micro": state.tendencies.dliquiddt_micro.field[:],
            "dicedt_micro": state.tendencies.dicedt_micro.field[:],
            "dcloud_fractiondt_micro": state.tendencies.dcloud_fractiondt_micro.field[:],
            "draindt_micro": state.tendencies.draindt_micro.field[:],
            "dsnowdt_micro": state.tendencies.dsnowdt_micro.field[:],
            "dgraupeldt_micro": state.tendencies.dgraupeldt_micro.field[:],
            "dudt_micro": state.tendencies.dudt_micro.field[:],
            "dvdt_micro": state.tendencies.dvdt_micro.field[:],
            "dtdt_micro": state.tendencies.dtdt_micro.field[:],
        }
