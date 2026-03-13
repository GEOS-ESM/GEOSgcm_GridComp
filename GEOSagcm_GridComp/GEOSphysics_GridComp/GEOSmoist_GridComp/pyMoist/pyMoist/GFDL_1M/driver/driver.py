"""GFDL_1M driver"""

import dace
from ndsl import NDSLRuntime, QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField, FloatFieldIJ

from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.GFDL_1M.driver.check_flags import check_flags
from pyMoist.GFDL_1M.driver.config_constants import GFDL1MDriverConfigDependentConstants
from pyMoist.GFDL_1M.driver.fall_speed import fall_speed
from pyMoist.GFDL_1M.driver.finish import update_tendencies
from pyMoist.GFDL_1M.driver.ice_cloud import GFDL1MIceCloud
from pyMoist.GFDL_1M.driver.locals import GFDL1MDriverLocals
from pyMoist.GFDL_1M.driver.sat_tables import get_tables
from pyMoist.GFDL_1M.driver.setup import GFDL1MDriverSetup
from pyMoist.GFDL_1M.driver.terminal_fall import GFDL1MTerminalFall
from pyMoist.GFDL_1M.driver.warm_rain import GFDL1MWarmRain


class GFDL1MDriver(NDSLRuntime):
    """
    Computes precipitates and microphysics tendencies using the following functions:
    __init__:
        - checks validity of constants and trigger parameters for unimplemented options
        - initializes internal fields
        - constructs stencils
        Arguments: StencilFactory, QuantityFactory, GFDL1MConfig

    __call__:
        Evaluate the microphysics driver. The driver call is broken into six parts:
        - Setup: fill temporaries, compute required intermediary fields from inputs
        - FallSpeed: compute real fall speed of precipitates
        - TerminalFall: compute terminal fall speed of precipitates
        - WarmRain: warm rain cloud microphysics
        - IceCloud: ice cloud microphysical processes
        - Finish: compute output tendencies
        Arguments: various state fields (pressure, temperature, wind, mixing ratios, etc)
    """

    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GFDL1MConfig,
    ):
        """
        Perform setup for the microphysics driver. Check flags for unimplemented options,
        initialize internal fields, and compile stencils.

        Arguments:
            stencil_factory: StencilFactory with model domain information
            quantity_factory: QuantityFactory with model domain information
            config: driver configuration
        """
        # init NDSLRuntime
        super().__init__(stencil_factory)

        self.config_dependent_constants = GFDL1MDriverConfigDependentConstants.make(config)

        constant_print_dict = {
            "dts": self.config_dependent_constants.DTS,
            "tau_imlt": config.TAU_IMLT,
            "ql_mlt": config.QL_MLT,
            "c_air": self.config_dependent_constants.C_AIR,
            "c_vap": self.config_dependent_constants.C_VAP,
            "d0_vap": self.config_dependent_constants.D0_VAP,
            "lv00": self.config_dependent_constants.LV00,
            "do_sedi_w": config.DO_SEDI_W,
        }
        print(f"RELEVANT CONSTANTS {constant_print_dict}")

        # Check constants for unimplemented and untested code paths
        check_flags(
            config,
            self.config_dependent_constants.DTS,
        )

        # initialize locals
        self._locals = GFDL1MDriverLocals.make_locals(quantity_factory)

        # pull saturation specific humidity tables, generate if first call
        self.driver_saturation_tables = get_tables(
            stencil_factory.backend,
            stencil_factory.config.dace_config,
        )

        # construct stencils
        self._setup = GFDL1MDriverSetup(
            stencil_factory=stencil_factory,
            config=config,
            config_dependent_constants=self.config_dependent_constants,
        )

        self._fall_speed = stencil_factory.from_dims_halo(
            func=fall_speed,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "p_nonhydro": self.config_dependent_constants.P_NONHYDRO,
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

        self._terminal_fall = GFDL1MTerminalFall(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            config_dependent_constants=self.config_dependent_constants,
        )

        self._warm_rain = GFDL1MWarmRain(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            config_dependent_constants=self.config_dependent_constants,
            saturation_tables=self.driver_saturation_tables,
        )

        self._ice_cloud = GFDL1MIceCloud(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            config_dependent_constants=self.config_dependent_constants,
            saturation_tables=self.driver_saturation_tables,
        )

        self._finish = stencil_factory.from_dims_halo(
            func=update_tendencies,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "c_air": self.config_dependent_constants.C_AIR,
                "c_vap": self.config_dependent_constants.C_VAP,
                "rdt": self.config_dependent_constants.RDT,
                "do_sedi_w": config.DO_SEDI_W,
                "sedi_transport": config.SEDI_TRANSPORT,
                "do_qa": config.DO_QA,
            },
        )

    def __call__(
        self,
        t: FloatField,
        u: FloatField,
        v: FloatField,
        w: FloatField,
        dz: FloatField,
        dp: FloatField,
        area: FloatFieldIJ,
        land_fraction: FloatFieldIJ,
        convection_fraction: FloatFieldIJ,
        surface_type: FloatFieldIJ,
        estimated_inversion_strength: FloatFieldIJ,
        critical_relative_humidity_for_pdf: FloatField,
        vapor: FloatField,
        liquid: FloatField,
        rain: FloatField,
        ice: FloatField,
        snow: FloatField,
        graupel: FloatField,
        cloud_fraction: FloatField,
        total_concentration: FloatField,
        dvapordt: FloatField,
        dliquiddt: FloatField,
        draindt: FloatField,
        dicedt: FloatField,
        dsnowdt: FloatField,
        dgraupeldt: FloatField,
        dcloudfractiondt: FloatField,
        dtdt: FloatField,
        dudt: FloatField,
        dvdt: FloatField,
        liquid_precip_flux: FloatField,
        ice_precip_flux: FloatField,
        evaporation: FloatField,
        sublimation: FloatField,
        surface_precip_rain: FloatFieldIJ,
        surface_precip_snow: FloatFieldIJ,
        surface_precip_ice: FloatFieldIJ,
        surface_precip_graupel: FloatFieldIJ,
        #
        DEBUG_vapor_unmodified_setup,
        DEBUG_vapor_modified_setup,
        DEBUG_vapor_unmodified_fallspeed,
        DEBUG_vapor_modified_fallspeed,
        DEBUG_vapor_unmodified_terminalfall,
        DEBUG_vapor_modified_terminalfall,
        DEBUG_vapor_unmodified_warmrain,
        DEBUG_vapor_modified_warmrain,
        DEBUG_vapor_unmodified_icecloud,
        DEBUG_vapor_modified_icecloud,
        #
        DEBUG_TERMINALFALL_IN_driver_local_t_terminalfall,
        DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_vapor_terminalfall,
        DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_liquid_terminalfall,
        DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_rain_terminalfall,
        DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_graupel_terminalfall,
        DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_snow_terminalfall,
        DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_ice_terminalfall,
        DEBUG_TERMINALFALL_IN_driver_local_ice_precip_flux_terminalfall,
        DEBUG_TERMINALFALL_IN_driver_local_w_terminalfall,
        DEBUG_TERMINALFALL_IN_driver_local_dz_terminalfall,
        DEBUG_TERMINALFALL_IN_driver_local_dp_terminalfall,
        DEBUG_TERMINALFALL_IN_driver_local_terminal_speed_ice_terminalfall,
        DEBUG_TERMINALFALL_IN_driver_local_terminal_speed_snow_terminalfall,
        DEBUG_TERMINALFALL_IN_driver_local_terminal_speed_graupel_terminalfall,
        DEBUG_TERMINALFALL_IN_surface_precip_rain_terminalfall,
        DEBUG_TERMINALFALL_IN_surface_precip_snow_terminalfall,
        DEBUG_TERMINALFALL_IN_surface_precip_graupel_terminalfall,
        DEBUG_TERMINALFALL_IN_surface_precip_ice_terminalfall,
        #
        DEBUG_TERMINALFALL_OUT_driver_local_t_terminalfall,
        DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_vapor_terminalfall,
        DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_liquid_terminalfall,
        DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_rain_terminalfall,
        DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_graupel_terminalfall,
        DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_snow_terminalfall,
        DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_ice_terminalfall,
        DEBUG_TERMINALFALL_OUT_driver_local_ice_precip_flux_terminalfall,
        DEBUG_TERMINALFALL_OUT_driver_local_w_terminalfall,
        DEBUG_TERMINALFALL_OUT_driver_local_dz_terminalfall,
        DEBUG_TERMINALFALL_OUT_driver_local_dp_terminalfall,
        DEBUG_TERMINALFALL_OUT_driver_local_terminal_speed_ice_terminalfall,
        DEBUG_TERMINALFALL_OUT_driver_local_terminal_speed_snow_terminalfall,
        DEBUG_TERMINALFALL_OUT_driver_local_terminal_speed_graupel_terminalfall,
        DEBUG_TERMINALFALL_OUT_surface_precip_rain_terminalfall,
        DEBUG_TERMINALFALL_OUT_surface_precip_snow_terminalfall,
        DEBUG_TERMINALFALL_OUT_surface_precip_graupel_terminalfall,
        DEBUG_TERMINALFALL_OUT_surface_precip_ice_terminalfall,
        #
        DEBUG_WARMRAIN_IN_driver_local_dp_warmrain,
        DEBUG_WARMRAIN_IN_driver_local_dz_warmrain,
        DEBUG_WARMRAIN_IN_driver_local_t_warmrain,
        DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_vapor_warmrain,
        DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_liquid_warmrain,
        DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_rain_warmrain,
        DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_ice_warmrain,
        DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_snow_warmrain,
        DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_graupel_warmrain,
        DEBUG_WARMRAIN_IN_driver_local_cloud_fraction_warmrain,
        DEBUG_WARMRAIN_IN_driver_local_ccn_warmrain,
        DEBUG_WARMRAIN_IN_driver_local_density_warmrain,
        DEBUG_WARMRAIN_IN_driver_local_density_factor_warmrain,
        DEBUG_WARMRAIN_IN_driver_local_c_praut_warmrain,
        DEBUG_WARMRAIN_IN_driver_local_terminal_speed_rain_warmrain,
        DEBUG_WARMRAIN_IN_driver_local_evaporation_warmrain,
        DEBUG_WARMRAIN_IN_driver_local_liquid_precip_flux_warmrain,
        DEBUG_WARMRAIN_IN_driver_local_w_warmrain,
        DEBUG_WARMRAIN_IN_driver_local_rh_limited_warmrain,
        DEBUG_WARMRAIN_IN_non_anvil_large_scale_evaporation_warmrain,
        DEBUG_WARMRAIN_IN_non_anvil_large_scale_liquid_precip_flux_warmrain,
        DEBUG_WARMRAIN_IN_non_anvil_large_scale_ice_precip_flux_warmrain,
        DEBUG_WARMRAIN_IN_driver_local_mass_warmrain,
        DEBUG_WARMRAIN_IN_driver_local_ice_precip_flux_warmrain,
        DEBUG_WARMRAIN_IN_driver_local_rain_warmrain,
        DEBUG_WARMRAIN_IN_surface_precip_rain_warmrain,
        DEBUG_WARMRAIN_IN_estimated_inversion_strength_warmrain,
        DEBUG_WARMRAIN_IN_driver_local_one_minus_sigma_warmrain,
        #
        DEBUG_WARMRAIN_OUT_driver_local_dp_warmrain,
        DEBUG_WARMRAIN_OUT_driver_local_dz_warmrain,
        DEBUG_WARMRAIN_OUT_driver_local_t_warmrain,
        DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_vapor_warmrain,
        DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_liquid_warmrain,
        DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_rain_warmrain,
        DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_ice_warmrain,
        DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_snow_warmrain,
        DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_graupel_warmrain,
        DEBUG_WARMRAIN_OUT_driver_local_cloud_fraction_warmrain,
        DEBUG_WARMRAIN_OUT_driver_local_ccn_warmrain,
        DEBUG_WARMRAIN_OUT_driver_local_density_warmrain,
        DEBUG_WARMRAIN_OUT_driver_local_density_factor_warmrain,
        DEBUG_WARMRAIN_OUT_driver_local_c_praut_warmrain,
        DEBUG_WARMRAIN_OUT_driver_local_terminal_speed_rain_warmrain,
        DEBUG_WARMRAIN_OUT_driver_local_evaporation_warmrain,
        DEBUG_WARMRAIN_OUT_driver_local_liquid_precip_flux_warmrain,
        DEBUG_WARMRAIN_OUT_driver_local_w_warmrain,
        DEBUG_WARMRAIN_OUT_driver_local_rh_limited_warmrain,
        DEBUG_WARMRAIN_OUT_non_anvil_large_scale_evaporation_warmrain,
        DEBUG_WARMRAIN_OUT_non_anvil_large_scale_liquid_precip_flux_warmrain,
        DEBUG_WARMRAIN_OUT_non_anvil_large_scale_ice_precip_flux_warmrain,
        DEBUG_WARMRAIN_OUT_driver_local_mass_warmrain,
        DEBUG_WARMRAIN_OUT_driver_local_ice_precip_flux_warmrain,
        DEBUG_WARMRAIN_OUT_driver_local_rain_warmrain,
        DEBUG_WARMRAIN_OUT_surface_precip_rain_warmrain,
        DEBUG_WARMRAIN_OUT_estimated_inversion_strength_warmrain,
        DEBUG_WARMRAIN_OUT_driver_local_one_minus_sigma_warmrain,
        #
        DEBUG_ICECLOUD_IN_driver_local_t_icecloud,
        DEBUG_ICECLOUD_IN_driver_local_p_dry_icecloud,
        DEBUG_ICECLOUD_IN_driver_local_dp_icecloud,
        DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_vapor_icecloud,
        DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_liquid_icecloud,
        DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_rain_icecloud,
        DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_ice_icecloud,
        DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_snow_icecloud,
        DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_graupel_icecloud,
        DEBUG_ICECLOUD_IN_driver_local_cloud_fraction_icecloud,
        DEBUG_ICECLOUD_IN_driver_local_terminal_speed_snow_icecloud,
        DEBUG_ICECLOUD_IN_driver_local_terminal_speed_graupel_icecloud,
        DEBUG_ICECLOUD_IN_driver_local_terminal_speed_rain_icecloud,
        DEBUG_ICECLOUD_IN_driver_local_density_icecloud,
        DEBUG_ICECLOUD_IN_driver_local_density_factor_icecloud,
        DEBUG_ICECLOUD_IN_driver_local_rh_limited_icecloud,
        DEBUG_ICECLOUD_IN_non_anvil_large_scale_sublimation_icecloud,
        DEBUG_ICECLOUD_IN_driver_local_ccn_icecloud,
        DEBUG_ICECLOUD_IN_convection_fraction_icecloud,
        DEBUG_ICECLOUD_IN_surface_type_icecloud,
        #
        DEBUG_ICECLOUD_OUT_driver_local_t_icecloud,
        DEBUG_ICECLOUD_OUT_driver_local_p_dry_icecloud,
        DEBUG_ICECLOUD_OUT_driver_local_dp_icecloud,
        DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_vapor_icecloud,
        DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_liquid_icecloud,
        DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_rain_icecloud,
        DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_ice_icecloud,
        DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_snow_icecloud,
        DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_graupel_icecloud,
        DEBUG_ICECLOUD_OUT_driver_local_cloud_fraction_icecloud,
        DEBUG_ICECLOUD_OUT_driver_local_terminal_speed_snow_icecloud,
        DEBUG_ICECLOUD_OUT_driver_local_terminal_speed_graupel_icecloud,
        DEBUG_ICECLOUD_OUT_driver_local_terminal_speed_rain_icecloud,
        DEBUG_ICECLOUD_OUT_driver_local_density_icecloud,
        DEBUG_ICECLOUD_OUT_driver_local_density_factor_icecloud,
        DEBUG_ICECLOUD_OUT_driver_local_rh_limited_icecloud,
        DEBUG_ICECLOUD_OUT_non_anvil_large_scale_sublimation_icecloud,
        DEBUG_ICECLOUD_OUT_driver_local_ccn_icecloud,
        DEBUG_ICECLOUD_OUT_convection_fraction_icecloud,
        DEBUG_ICECLOUD_OUT_surface_type_icecloud,
    ):
        """

        Arguments:
            t (in): atmospheric temperature (K)
            u (in): eastward winds (m/s)
            v (in): northward winds (m/s)
            w (in): vertical velocity (m/s)
            dz (in): layer thickness (m)
            dp (in): change in pressure between model levels (Pa)
            area (in): grid cell area
            land_fraction (in): land fraction
            convection_fraction (in): convection fraction
            surface_type (in): surface type
            estimated_inversion_strength (in): estimated inversion strength (K)
            critical_relative_humidity_for_pdf (in): critical relative humidity for pdf
            vapor: (in): water vapor mixing ratio (kg/kg)
            liquid (in): in cloud liquid mixing radio (kg/kg)
            rain (in): falling rain (kg/kg)
            ice (in): in cloud ice mixing radio (kg/kg)
            snow (in): in cloud snow mixing radio (kg/kg)
            graupel (in): in cloud graupel mixing radio (kg/kg)
            cloud_fraction (in): cloud fraction (convective + large scale)
            total_concentration (in): total ice + liquid concentration (m^-3)
            dvapordt (out): water vapor tendency
            dliquiddt (out): in cloud liquid water tendency
            draindt (out): falling rain tendency
            dicedt (out): in cloud frozen water tendency
            dsnowdt (out): in cloud snow tendency
            dgraupeldt (out): in cloud graupel tendency
            dcloudfractiondt (out): cloud fraction (convective + large scale) tendency
            dtdt (out): atmospheric temperature tendency
            dudt (out): eastward wind tendency
            dvdt (out): northward wind tendency
            liquid_precip_flux (out): non-anvil large scale liquid precip flux
            ice_precip_flux (out): non-anvil large scale ice precip flux
            evaporation (out): non-anvil large scale evaporation
            sublimation (out): non-anvil large scale sublimation
            surface_precip_rain (out): rain precip at surface (kg/m^2/s)
            surface_precip_snow (out): snow precip at surface (kg/m^2/s)
            surface_precip_ice (out): ice precip at surface (kg/m^2/s)
            surface_precip_graupel (out): graupel precip at surface (kg/m^2/s)
        """
        self._setup(
            unmodified_t=t,
            t=self._locals.t,  # driver_locals.t,
            unmodified_dp=dp,  # locals_.dp,
            dp=self._locals.dp,  # driver_locals.dp,
            critical_relative_humidity_for_pdf=critical_relative_humidity_for_pdf,  # state.critical_relative_humidity_for_pdf,
            radiation_field_vapor=vapor,  # state.radiation_field.vapor,
            radiation_field_liquid=liquid,  # state.radiation_field.liquid,
            radiation_field_ice=ice,  # state.radiation_field.ice,
            radiation_field_rain=rain,  # state.radiation_field.rain,
            radiation_field_snow=snow,  # state.radiation_field.snow,
            radiation_field_graupel=graupel,  # state.radiation_field.graupel,
            radiation_field_cloud_fraction=cloud_fraction,  # state.radiation_field.cloud_fraction,
            total_concentration=total_concentration,  # locals_.total_concentration,
            unmodified_mixing_ratio_vapor=self._locals.unmodified.mixing_ratio.vapor,  # driver_locals.unmodified.mixing_ratio.vapor,
            unmodified_mixing_ratio_liquid=self._locals.unmodified.mixing_ratio.liquid,  # driver_locals.unmodified.mixing_ratio.liquid,
            unmodified_mixing_ratio_rain=self._locals.unmodified.mixing_ratio.rain,  # driver_locals.unmodified.mixing_ratio.rain,
            unmodified_mixing_ratio_ice=self._locals.unmodified.mixing_ratio.ice,  # driver_locals.unmodified.mixing_ratio.ice,
            unmodified_mixing_ratio_snow=self._locals.unmodified.mixing_ratio.snow,  # driver_locals.unmodified.mixing_ratio.snow,
            unmodified_mixing_ratio_graupel=self._locals.unmodified.mixing_ratio.graupel,  # driver_locals.unmodified.mixing_ratio.graupel,
            dry_air_mixing_ratio_vapor=self._locals.dry_air_mixing_ratio.vapor,  # driver_locals.dry_air_mixing_ratio.vapor,
            dry_air_mixing_ratio_liquid=self._locals.dry_air_mixing_ratio.liquid,  # driver_locals.dry_air_mixing_ratio.liquid,
            dry_air_mixing_ratio_rain=self._locals.dry_air_mixing_ratio.rain,  # driver_locals.dry_air_mixing_ratio.rain,
            dry_air_mixing_ratio_ice=self._locals.dry_air_mixing_ratio.ice,  # driver_locals.dry_air_mixing_ratio.ice,
            dry_air_mixing_ratio_snow=self._locals.dry_air_mixing_ratio.snow,  # driver_locals.dry_air_mixing_ratio.snow,
            dry_air_mixing_ratio_graupel=self._locals.dry_air_mixing_ratio.graupel,  # driver_locals.dry_air_mixing_ratio.graupel,
            cloud_fraction=self._locals.cloud_fraction,  # driver_locals.cloud_fraction,
            dz=dz,  # locals_.layer_height_above_surface,
            u_unmodified=u,  # state.u,
            u=self._locals.u,  # driver_locals.u,
            v_unmodified=v,  # state.v,
            v=self._locals.v,  # driver_locals.v,
            w_unmodified=w,  # state.vertical_motion.velocity,
            w=self._locals.w,  # driver_locals.w,
            area=area,  # state.area,
            density_unmodified=self._locals.density_unmodified,  # driver_locals.density_unmodified,
            p_dry=self._locals.p_dry,  # driver_locals.p_dry,
            mass=self._locals.mass,  # driver_locals.mass,
            one_minus_sigma=self._locals.one_minus_sigma,  # driver_locals.one_minus_sigma,
            ccn=self._locals.ccn,  # driver_locals.ccn,
            c_praut=self._locals.c_praut,  # driver_locals.c_praut,
            rh_limited=self._locals.rh_limited,  # driver_locals.rh_limited,
            rain=surface_precip_rain,  # state.precipitation_at_surface.rain,
            snow=surface_precip_snow,  # state.precipitation_at_surface.snow,
            graupel=surface_precip_graupel,  # state.precipitation_at_surface.graupel,
            ice=surface_precip_ice,  # state.precipitation_at_surface.ice,
            liquid_precip_flux=liquid_precip_flux,  # state.non_anvil_large_scale.liquid_precip_flux,
            ice_precip_flux=ice_precip_flux,  # state.non_anvil_large_scale.ice_precip_flux,
            evaporation=evaporation,  # state.non_anvil_large_scale.evaporation,
            sublimation=sublimation,  # state.non_anvil_large_scale.sublimation,
        )

        DEBUG_vapor_unmodified_setup.field[:] = self._locals.unmodified.mixing_ratio.vapor.field[:]
        DEBUG_vapor_modified_setup.field[:] = self._locals.dry_air_mixing_ratio.vapor.field[:]

        for _ in range(1):  # dace.nounroll(range(self.config_dependent_constants.NTIMES)):
            self._fall_speed(
                liquid=self._locals.dry_air_mixing_ratio.liquid,  # driver_locals.dry_air_mixing_ratio.liquid,
                ice=self._locals.dry_air_mixing_ratio.ice,  # driver_locals.dry_air_mixing_ratio.ice,
                snow=self._locals.dry_air_mixing_ratio.snow,  # driver_locals.dry_air_mixing_ratio.snow,
                graupel=self._locals.dry_air_mixing_ratio.graupel,  # driver_locals.dry_air_mixing_ratio.graupel,
                t_unmodified=t,  # state.t,
                t=self._locals.t,  # driver_locals.t,
                dz_unmodified=dz,  # gfdl1m_locals.layer_thickness_negative,
                dz=self._locals.dz,  # driver_locals.dz,
                density_unmodified=self._locals.density_unmodified,  # driver_locals.density_unmodified,
                density=self._locals.density,  # driver_locals.density,
                density_factor=self._locals.density_factor,  # driver_locals.density_factor,
                ice_terminal_velocity=self._locals.terminal_speed.ice,  # driver_locals.terminal_speed.ice,
                snow_terminal_velocity=self._locals.terminal_speed.snow,  # driver_locals.terminal_speed.snow,
                graupel_terminal_velocity=self._locals.terminal_speed.graupel,  # driver_locals.terminal_speed.graupel,
                convection_fraction=convection_fraction,  # state.convection_fraction,
            )

            DEBUG_vapor_unmodified_fallspeed.field[:] = self._locals.unmodified.mixing_ratio.vapor.field[:]
            DEBUG_vapor_modified_fallspeed.field[:] = self._locals.dry_air_mixing_ratio.vapor.field[:]

            DEBUG_TERMINALFALL_IN_driver_local_t_terminalfall.field[:] = self._locals.t.field[:]
            DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_vapor_terminalfall.field[:] = (
                self._locals.dry_air_mixing_ratio.vapor.field[:]
            )
            DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_liquid_terminalfall.field[:] = (
                self._locals.dry_air_mixing_ratio.liquid.field[:]
            )
            DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_rain_terminalfall.field[:] = (
                self._locals.dry_air_mixing_ratio.rain.field[:]
            )
            DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_graupel_terminalfall.field[:] = (
                self._locals.dry_air_mixing_ratio.graupel.field[:]
            )
            DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_snow_terminalfall.field[:] = (
                self._locals.dry_air_mixing_ratio.snow.field[:]
            )
            DEBUG_TERMINALFALL_IN_driver_local_dry_mixing_ratio_ice_terminalfall.field[:] = (
                self._locals.dry_air_mixing_ratio.ice.field[:]
            )
            DEBUG_TERMINALFALL_IN_driver_local_ice_precip_flux_terminalfall.field[:] = (
                self._locals.ice_precip_flux.field[:]
            )
            DEBUG_TERMINALFALL_IN_driver_local_w_terminalfall.field[:] = self._locals.w.field[:]
            DEBUG_TERMINALFALL_IN_driver_local_dz_terminalfall.field[:] = self._locals.dz.field[:]
            DEBUG_TERMINALFALL_IN_driver_local_dp_terminalfall.field[:] = self._locals.dp.field[:]
            DEBUG_TERMINALFALL_IN_driver_local_terminal_speed_ice_terminalfall.field[:] = (
                self._locals.terminal_speed.ice.field[:]
            )
            DEBUG_TERMINALFALL_IN_driver_local_terminal_speed_snow_terminalfall.field[:] = (
                self._locals.terminal_speed.snow.field[:]
            )
            DEBUG_TERMINALFALL_IN_driver_local_terminal_speed_graupel_terminalfall.field[:] = (
                self._locals.terminal_speed.graupel.field[:]
            )
            for k in range(72):
                DEBUG_TERMINALFALL_IN_surface_precip_rain_terminalfall.field[:, :, k] = (
                    surface_precip_rain.field[:]
                )
                DEBUG_TERMINALFALL_IN_surface_precip_snow_terminalfall.field[:, :, k] = (
                    surface_precip_snow.field[:]
                )
                DEBUG_TERMINALFALL_IN_surface_precip_graupel_terminalfall.field[:, :, k] = (
                    surface_precip_graupel.field[:]
                )
                DEBUG_TERMINALFALL_IN_surface_precip_ice_terminalfall.field[:, :, k] = (
                    surface_precip_ice.field[:]
                )

            self._terminal_fall(
                t=self._locals.t,  # driver_locals.t,
                w=self._locals.w,  # driver_locals.w,
                mixing_ratio_vapor=self._locals.dry_air_mixing_ratio.vapor,  # driver_locals.dry_air_mixing_ratio.vapor,
                mixing_ratio_liquid=self._locals.dry_air_mixing_ratio.liquid,  # driver_locals.dry_air_mixing_ratio.liquid,
                mixing_ratio_rain=self._locals.dry_air_mixing_ratio.rain,  # driver_locals.dry_air_mixing_ratio.rain,
                mixing_ratio_graupel=self._locals.dry_air_mixing_ratio.graupel,  # driver_locals.dry_air_mixing_ratio.graupel,
                mixing_ratio_snow=self._locals.dry_air_mixing_ratio.snow,  # driver_locals.dry_air_mixing_ratio.snow,
                mixing_ratio_ice=self._locals.dry_air_mixing_ratio.ice,  # driver_locals.dry_air_mixing_ratio.ice,
                dz=self._locals.dz,  # driver_locals.dz,
                dp=self._locals.dp,  # driver_locals.dp,
                terminal_velocity_graupel=self._locals.terminal_speed.ice,  # driver_locals.terminal_speed.graupel,
                terminal_velocity_snow=self._locals.terminal_speed.snow,  # driver_locals.terminal_speed.snow,
                terminal_velocity_ice=self._locals.terminal_speed.graupel,  # driver_locals.terminal_speed.ice,
                rain=surface_precip_rain,  # state.precipitation_at_surface.rain,
                graupel=surface_precip_graupel,  # state.precipitation_at_surface.graupel,
                snow=surface_precip_snow,  # state.precipitation_at_surface.snow,
                ice=surface_precip_ice,  # state.precipitation_at_surface.ice,
                ice_precip_flux=self._locals.ice_precip_flux,  # driver_locals.ice_precip_flux,
            )

            DEBUG_TERMINALFALL_OUT_driver_local_t_terminalfall.field[:] = self._locals.t.field[:]
            DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_vapor_terminalfall.field[:] = (
                self._locals.dry_air_mixing_ratio.vapor.field[:]
            )
            DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_liquid_terminalfall.field[:] = (
                self._locals.dry_air_mixing_ratio.liquid.field[:]
            )
            DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_rain_terminalfall.field[:] = (
                self._locals.dry_air_mixing_ratio.rain.field[:]
            )
            DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_graupel_terminalfall.field[:] = (
                self._locals.dry_air_mixing_ratio.graupel.field[:]
            )
            DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_snow_terminalfall.field[:] = (
                self._locals.dry_air_mixing_ratio.snow.field[:]
            )
            DEBUG_TERMINALFALL_OUT_driver_local_dry_mixing_ratio_ice_terminalfall.field[:] = (
                self._locals.dry_air_mixing_ratio.ice.field[:]
            )
            DEBUG_TERMINALFALL_OUT_driver_local_ice_precip_flux_terminalfall.field[:] = (
                self._locals.ice_precip_flux.field[:]
            )
            DEBUG_TERMINALFALL_OUT_driver_local_w_terminalfall.field[:] = self._locals.w.field[:]
            DEBUG_TERMINALFALL_OUT_driver_local_dz_terminalfall.field[:] = self._locals.dz.field[:]
            DEBUG_TERMINALFALL_OUT_driver_local_dp_terminalfall.field[:] = self._locals.dp.field[:]
            DEBUG_TERMINALFALL_OUT_driver_local_terminal_speed_ice_terminalfall.field[:] = (
                self._locals.terminal_speed.ice.field[:]
            )
            DEBUG_TERMINALFALL_OUT_driver_local_terminal_speed_snow_terminalfall.field[:] = (
                self._locals.terminal_speed.snow.field[:]
            )
            DEBUG_TERMINALFALL_OUT_driver_local_terminal_speed_graupel_terminalfall.field[:] = (
                self._locals.terminal_speed.graupel.field[:]
            )
            for k in range(72):
                DEBUG_TERMINALFALL_OUT_surface_precip_rain_terminalfall.field[:, :, k] = (
                    surface_precip_rain.field[:]
                )
                DEBUG_TERMINALFALL_OUT_surface_precip_snow_terminalfall.field[:, :, k] = (
                    surface_precip_snow.field[:]
                )
                DEBUG_TERMINALFALL_OUT_surface_precip_graupel_terminalfall.field[:, :, k] = (
                    surface_precip_graupel.field[:]
                )
                DEBUG_TERMINALFALL_OUT_surface_precip_ice_terminalfall.field[:, :, k] = (
                    surface_precip_ice.field[:]
                )

            DEBUG_vapor_unmodified_terminalfall.field[:] = self._locals.unmodified.mixing_ratio.vapor.field[:]
            DEBUG_vapor_modified_terminalfall.field[:] = self._locals.dry_air_mixing_ratio.vapor.field[:]

            DEBUG_WARMRAIN_IN_driver_local_dp_warmrain.field[:] = self._locals.dp.field[:]
            DEBUG_WARMRAIN_IN_driver_local_dz_warmrain.field[:] = self._locals.dz.field[:]
            DEBUG_WARMRAIN_IN_driver_local_t_warmrain.field[:] = self._locals.t.field[:]
            DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_vapor_warmrain.field[:] = (
                self._locals.dry_air_mixing_ratio.vapor.field[:]
            )
            DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_liquid_warmrain.field[:] = (
                self._locals.dry_air_mixing_ratio.liquid.field[:]
            )
            DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_rain_warmrain.field[:] = (
                self._locals.dry_air_mixing_ratio.rain.field[:]
            )
            DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_ice_warmrain.field[:] = (
                self._locals.dry_air_mixing_ratio.ice.field[:]
            )
            DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_snow_warmrain.field[:] = (
                self._locals.dry_air_mixing_ratio.snow.field[:]
            )
            DEBUG_WARMRAIN_IN_driver_local_dry_mixing_ratio_graupel_warmrain.field[:] = (
                self._locals.dry_air_mixing_ratio.graupel.field[:]
            )
            DEBUG_WARMRAIN_IN_driver_local_cloud_fraction_warmrain.field[:] = (
                self._locals.cloud_fraction.field[:]
            )
            DEBUG_WARMRAIN_IN_driver_local_ccn_warmrain.field[:] = self._locals.ccn.field[:]
            DEBUG_WARMRAIN_IN_driver_local_density_warmrain.field[:] = self._locals.density.field[:]
            DEBUG_WARMRAIN_IN_driver_local_density_factor_warmrain.field[:] = (
                self._locals.density_factor.field[:]
            )
            DEBUG_WARMRAIN_IN_driver_local_c_praut_warmrain.field[:] = self._locals.c_praut.field[:]
            DEBUG_WARMRAIN_IN_driver_local_terminal_speed_rain_warmrain.field[:] = (
                self._locals.terminal_speed.rain.field[:]
            )
            DEBUG_WARMRAIN_IN_driver_local_evaporation_warmrain.field[:] = self._locals.evaporation.field[:]
            DEBUG_WARMRAIN_IN_driver_local_liquid_precip_flux_warmrain.field[:] = (
                self._locals.liquid_precip_flux.field[:]
            )
            DEBUG_WARMRAIN_IN_driver_local_w_warmrain.field[:] = self._locals.w.field[:]
            DEBUG_WARMRAIN_IN_driver_local_rh_limited_warmrain.field[:] = self._locals.rh_limited.field[:]
            DEBUG_WARMRAIN_IN_non_anvil_large_scale_evaporation_warmrain.field[:] = evaporation.field[:]
            DEBUG_WARMRAIN_IN_non_anvil_large_scale_liquid_precip_flux_warmrain.field[:] = (
                liquid_precip_flux.field[:]
            )
            DEBUG_WARMRAIN_IN_non_anvil_large_scale_ice_precip_flux_warmrain.field[:] = ice_precip_flux.field[
                :
            ]
            DEBUG_WARMRAIN_IN_driver_local_mass_warmrain.field[:] = self._locals.mass.field[:]
            DEBUG_WARMRAIN_IN_driver_local_ice_precip_flux_warmrain.field[:] = (
                self._locals.ice_precip_flux.field[:]
            )
            for k in range(72):
                DEBUG_WARMRAIN_IN_driver_local_rain_warmrain.field[:, :, k] = self._locals.rain.field[:]
                DEBUG_WARMRAIN_IN_surface_precip_rain_warmrain.field[:, :, k] = surface_precip_rain.field[:]
                DEBUG_WARMRAIN_IN_estimated_inversion_strength_warmrain.field[:, :, k] = (
                    estimated_inversion_strength.field[:]
                )
                DEBUG_WARMRAIN_IN_driver_local_one_minus_sigma_warmrain.field[:, :, k] = (
                    self._locals.one_minus_sigma.field[:]
                )

            self._warm_rain(
                t=self._locals.t,  # driver_locals.t,
                dp=self._locals.dp,  # driver_locals.dp,
                dz=self._locals.dz,  # driver_locals.dz,
                w=self._locals.w,  # driver_locals.w,
                mixing_ratio_vapor=self._locals.dry_air_mixing_ratio.vapor,  # driver_locals.dry_air_mixing_ratio.vapor,
                mixing_ratio_liquid=self._locals.dry_air_mixing_ratio.liquid,  # driver_locals.dry_air_mixing_ratio.liquid,
                mixing_ratio_rain=self._locals.dry_air_mixing_ratio.rain,  # driver_locals.dry_air_mixing_ratio.rain,
                mixing_ratio_ice=self._locals.dry_air_mixing_ratio.ice,  # driver_locals.dry_air_mixing_ratio.ice,
                mixing_ratio_snow=self._locals.dry_air_mixing_ratio.snow,  # driver_locals.dry_air_mixing_ratio.snow,
                mixing_ratio_graupel=self._locals.dry_air_mixing_ratio.graupel,  # driver_locals.dry_air_mixing_ratio.graupel,
                cloud_fraction=self._locals.cloud_fraction,  # driver_locals.cloud_fraction,
                ccn=self._locals.ccn,  # driver_locals.ccn,
                density=self._locals.density,  # driver_locals.density,
                density_factor=self._locals.density_factor,  # driver_locals.density_factor,
                c_praut=self._locals.c_praut,  # driver_locals.c_praut,
                terminal_speed_rain=self._locals.terminal_speed.rain,  # driver_locals.terminal_speed.rain,
                rh_limited=self._locals.rh_limited,  # driver_locals.rh_limited,
                estimated_inversion_strength=estimated_inversion_strength,  # state.estimated_inversion_strength,
                one_minus_sigma=self._locals.one_minus_sigma,  # driver_locals.one_minus_sigma,
                mass=self._locals.mass,  # driver_locals.mass,
                rain=surface_precip_rain,  # state.precipitation_at_surface.rain,
                driver_rain=self._locals.rain,  # driver_locals.rain,
                ice_precip_flux=ice_precip_flux,  # state.non_anvil_large_scale.ice_precip_flux,
                driver_ice_precip_flux=self._locals.ice_precip_flux,  # driver_locals.ice_precip_flux,
                liquid_precip_flux=liquid_precip_flux,  # state.non_anvil_large_scale.liquid_precip_flux,
                driver_liquid_precip_flux=self._locals.liquid_precip_flux,  # driver_locals.liquid_precip_flux,
                evaporation=evaporation,  # state.non_anvil_large_scale.evaporation,
                driver_evaporation=self._locals.evaporation,  # driver_locals.evaporation,
            )

            DEBUG_WARMRAIN_OUT_driver_local_dp_warmrain.field[:] = self._locals.dp.field[:]
            DEBUG_WARMRAIN_OUT_driver_local_dz_warmrain.field[:] = self._locals.dz.field[:]
            DEBUG_WARMRAIN_OUT_driver_local_t_warmrain.field[:] = self._locals.t.field[:]
            DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_vapor_warmrain.field[:] = (
                self._locals.dry_air_mixing_ratio.vapor.field[:]
            )
            DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_liquid_warmrain.field[:] = (
                self._locals.dry_air_mixing_ratio.liquid.field[:]
            )
            DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_rain_warmrain.field[:] = (
                self._locals.dry_air_mixing_ratio.rain.field[:]
            )
            DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_ice_warmrain.field[:] = (
                self._locals.dry_air_mixing_ratio.ice.field[:]
            )
            DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_snow_warmrain.field[:] = (
                self._locals.dry_air_mixing_ratio.snow.field[:]
            )
            DEBUG_WARMRAIN_OUT_driver_local_dry_mixing_ratio_graupel_warmrain.field[:] = (
                self._locals.dry_air_mixing_ratio.graupel.field[:]
            )
            DEBUG_WARMRAIN_OUT_driver_local_cloud_fraction_warmrain.field[:] = (
                self._locals.cloud_fraction.field[:]
            )
            DEBUG_WARMRAIN_OUT_driver_local_ccn_warmrain.field[:] = self._locals.ccn.field[:]
            DEBUG_WARMRAIN_OUT_driver_local_density_warmrain.field[:] = self._locals.density.field[:]
            DEBUG_WARMRAIN_OUT_driver_local_density_factor_warmrain.field[:] = (
                self._locals.density_factor.field[:]
            )
            DEBUG_WARMRAIN_OUT_driver_local_c_praut_warmrain.field[:] = self._locals.c_praut.field[:]
            DEBUG_WARMRAIN_OUT_driver_local_terminal_speed_rain_warmrain.field[:] = (
                self._locals.terminal_speed.rain.field[:]
            )
            DEBUG_WARMRAIN_OUT_driver_local_evaporation_warmrain.field[:] = self._locals.evaporation.field[:]
            DEBUG_WARMRAIN_OUT_driver_local_liquid_precip_flux_warmrain.field[:] = (
                self._locals.liquid_precip_flux.field[:]
            )
            DEBUG_WARMRAIN_OUT_driver_local_w_warmrain.field[:] = self._locals.w.field[:]
            DEBUG_WARMRAIN_OUT_driver_local_rh_limited_warmrain.field[:] = self._locals.rh_limited.field[:]
            DEBUG_WARMRAIN_OUT_non_anvil_large_scale_evaporation_warmrain.field[:] = evaporation.field[:]
            DEBUG_WARMRAIN_OUT_non_anvil_large_scale_liquid_precip_flux_warmrain.field[:] = (
                liquid_precip_flux.field[:]
            )
            DEBUG_WARMRAIN_OUT_non_anvil_large_scale_ice_precip_flux_warmrain.field[:] = (
                ice_precip_flux.field[:]
            )
            DEBUG_WARMRAIN_OUT_driver_local_mass_warmrain.field[:] = self._locals.mass.field[:]
            DEBUG_WARMRAIN_OUT_driver_local_ice_precip_flux_warmrain.field[:] = (
                self._locals.ice_precip_flux.field[:]
            )
            for k in range(72):
                DEBUG_WARMRAIN_OUT_driver_local_rain_warmrain.field[:, :, k] = self._locals.rain.field[:]
                DEBUG_WARMRAIN_OUT_surface_precip_rain_warmrain.field[:, :, k] = surface_precip_rain.field[:]
                DEBUG_WARMRAIN_OUT_estimated_inversion_strength_warmrain.field[:, :, k] = (
                    estimated_inversion_strength.field[:]
                )
                DEBUG_WARMRAIN_OUT_driver_local_one_minus_sigma_warmrain.field[:, :, k] = (
                    self._locals.one_minus_sigma.field[:]
                )

            DEBUG_vapor_unmodified_warmrain.field[:] = self._locals.unmodified.mixing_ratio.vapor.field[:]
            DEBUG_vapor_modified_warmrain.field[:] = self._locals.dry_air_mixing_ratio.vapor.field[:]

            DEBUG_ICECLOUD_IN_driver_local_t_icecloud.field[:] = self._locals.t.field[:]
            DEBUG_ICECLOUD_IN_driver_local_p_dry_icecloud.field[:] = self._locals.p_dry.field[:]
            DEBUG_ICECLOUD_IN_driver_local_dp_icecloud.field[:] = self._locals.dp.field[:]
            DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_vapor_icecloud.field[:] = (
                self._locals.dry_air_mixing_ratio.vapor.field[:]
            )
            DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_liquid_icecloud.field[:] = (
                self._locals.dry_air_mixing_ratio.liquid.field[:]
            )
            DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_rain_icecloud.field[:] = (
                self._locals.dry_air_mixing_ratio.rain.field[:]
            )
            DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_ice_icecloud.field[:] = (
                self._locals.dry_air_mixing_ratio.ice.field[:]
            )
            DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_snow_icecloud.field[:] = (
                self._locals.dry_air_mixing_ratio.snow.field[:]
            )
            DEBUG_ICECLOUD_IN_driver_local_dry_mixing_ratio_graupel_icecloud.field[:] = (
                self._locals.dry_air_mixing_ratio.graupel.field[:]
            )
            DEBUG_ICECLOUD_IN_driver_local_cloud_fraction_icecloud.field[:] = (
                self._locals.cloud_fraction.field[:]
            )
            DEBUG_ICECLOUD_IN_driver_local_terminal_speed_snow_icecloud.field[:] = (
                self._locals.terminal_speed.snow.field[:]
            )
            DEBUG_ICECLOUD_IN_driver_local_terminal_speed_graupel_icecloud.field[:] = (
                self._locals.terminal_speed.graupel.field[:]
            )
            DEBUG_ICECLOUD_IN_driver_local_terminal_speed_rain_icecloud.field[:] = (
                self._locals.terminal_speed.rain.field[:]
            )
            DEBUG_ICECLOUD_IN_driver_local_density_icecloud.field[:] = self._locals.density.field[:]
            DEBUG_ICECLOUD_IN_driver_local_density_factor_icecloud.field[:] = (
                self._locals.density_factor.field[:]
            )
            DEBUG_ICECLOUD_IN_driver_local_rh_limited_icecloud.field[:] = self._locals.rh_limited.field[:]
            DEBUG_ICECLOUD_IN_non_anvil_large_scale_sublimation_icecloud.field[:] = sublimation.field[:]
            DEBUG_ICECLOUD_IN_driver_local_ccn_icecloud.field[:] = self._locals.ccn.field[:]
            for k in range(72):
                DEBUG_ICECLOUD_IN_convection_fraction_icecloud.field[:, :, k] = convection_fraction.field[:]
                DEBUG_ICECLOUD_IN_surface_type_icecloud.field[:, :, k] = surface_type.field[:]

            self._ice_cloud(
                t=self._locals.t,  # driver_locals.t,
                p_dry=self._locals.p_dry,  # driver_locals.p_dry,
                dp=self._locals.dp,  # driver_locals.dp,
                vapor=self._locals.dry_air_mixing_ratio.vapor,  # driver_locals.dry_air_mixing_ratio.vapor,
                liquid=self._locals.dry_air_mixing_ratio.liquid,  # driver_locals.dry_air_mixing_ratio.liquid,
                rain=self._locals.dry_air_mixing_ratio.rain,  # driver_locals.dry_air_mixing_ratio.rain,
                ice=self._locals.dry_air_mixing_ratio.ice,  # driver_locals.dry_air_mixing_ratio.ice,
                snow=self._locals.dry_air_mixing_ratio.snow,  # driver_locals.dry_air_mixing_ratio.snow,
                graupel=self._locals.dry_air_mixing_ratio.graupel,  # driver_locals.dry_air_mixing_ratio.graupel,
                cloud_fraction=self._locals.cloud_fraction,  # driver_locals.cloud_fraction,
                density=self._locals.density,  # driver_locals.density,
                density_factor=self._locals.density_factor,  # driver_locals.density_factor,
                terminal_fall_snow=self._locals.terminal_speed.snow,  # driver_locals.terminal_speed.snow,
                terminal_fall_graupel=self._locals.terminal_speed.graupel,  # driver_locals.terminal_speed.graupel,
                terminal_fall_rain=self._locals.terminal_speed.rain,  # driver_locals.terminal_speed.rain,
                sublimation=sublimation,  # state.non_anvil_large_scale.sublimation,
                rh_limited=self._locals.rh_limited,  # driver_locals.rh_limited,
                ccn=self._locals.ccn,  # driver_locals.ccn,
                convection_fraction=convection_fraction,  # state.convection_fraction,
                surface_type=surface_type,  # state.surface_type,
            )

            DEBUG_vapor_unmodified_icecloud.field[:] = self._locals.unmodified.mixing_ratio.vapor.field[:]
            DEBUG_vapor_modified_icecloud.field[:] = self._locals.dry_air_mixing_ratio.vapor.field[:]

            DEBUG_ICECLOUD_OUT_driver_local_t_icecloud.field[:] = self._locals.t.field[:]
            DEBUG_ICECLOUD_OUT_driver_local_p_dry_icecloud.field[:] = self._locals.p_dry.field[:]
            DEBUG_ICECLOUD_OUT_driver_local_dp_icecloud.field[:] = self._locals.dp.field[:]
            DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_vapor_icecloud.field[:] = (
                self._locals.dry_air_mixing_ratio.vapor.field[:]
            )
            DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_liquid_icecloud.field[:] = (
                self._locals.dry_air_mixing_ratio.liquid.field[:]
            )
            DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_rain_icecloud.field[:] = (
                self._locals.dry_air_mixing_ratio.rain.field[:]
            )
            DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_ice_icecloud.field[:] = (
                self._locals.dry_air_mixing_ratio.ice.field[:]
            )
            DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_snow_icecloud.field[:] = (
                self._locals.dry_air_mixing_ratio.snow.field[:]
            )
            DEBUG_ICECLOUD_OUT_driver_local_dry_mixing_ratio_graupel_icecloud.field[:] = (
                self._locals.dry_air_mixing_ratio.graupel.field[:]
            )
            DEBUG_ICECLOUD_OUT_driver_local_cloud_fraction_icecloud.field[:] = (
                self._locals.cloud_fraction.field[:]
            )
            DEBUG_ICECLOUD_OUT_driver_local_terminal_speed_snow_icecloud.field[:] = (
                self._locals.terminal_speed.snow.field[:]
            )
            DEBUG_ICECLOUD_OUT_driver_local_terminal_speed_graupel_icecloud.field[:] = (
                self._locals.terminal_speed.graupel.field[:]
            )
            DEBUG_ICECLOUD_OUT_driver_local_terminal_speed_rain_icecloud.field[:] = (
                self._locals.terminal_speed.rain.field[:]
            )
            DEBUG_ICECLOUD_OUT_driver_local_density_icecloud.field[:] = self._locals.density.field[:]
            DEBUG_ICECLOUD_OUT_driver_local_density_factor_icecloud.field[:] = (
                self._locals.density_factor.field[:]
            )
            DEBUG_ICECLOUD_OUT_driver_local_rh_limited_icecloud.field[:] = self._locals.rh_limited.field[:]
            DEBUG_ICECLOUD_OUT_non_anvil_large_scale_sublimation_icecloud.field[:] = sublimation.field[:]
            DEBUG_ICECLOUD_OUT_driver_local_ccn_icecloud.field[:] = self._locals.ccn.field[:]
            for k in range(72):
                DEBUG_ICECLOUD_OUT_convection_fraction_icecloud.field[:, :, k] = convection_fraction.field[:]
                DEBUG_ICECLOUD_OUT_surface_type_icecloud.field[:, :, k] = surface_type.field[:]

        self._finish(
            mixing_ratio_vapor_unmodified=self._locals.unmodified.mixing_ratio.vapor,  # driver_locals.unmodified.mixing_ratio.vapor,
            mixing_ratio_liquid_unmodified=self._locals.unmodified.mixing_ratio.liquid,  # driver_locals.unmodified.mixing_ratio.liquid,
            mixing_ratio_rain_unmodified=self._locals.unmodified.mixing_ratio.rain,  # driver_locals.unmodified.mixing_ratio.rain,
            mixing_ratio_ice_unmodified=self._locals.unmodified.mixing_ratio.ice,  # driver_locals.unmodified.mixing_ratio.ice,
            mixing_ratio_snow_unmodified=self._locals.unmodified.mixing_ratio.snow,  # driver_locals.unmodified.mixing_ratio.rain,
            mixing_ratio_graupel_unmodified=self._locals.unmodified.mixing_ratio.graupel,  # driver_locals.unmodified.mixing_ratio.graupel,
            cloud_fraction_unmodified=cloud_fraction,  # state.radiation_field.cloud_fraction,
            mixing_ratio_driver_vapor=self._locals.dry_air_mixing_ratio.vapor,  # driver_locals.dry_air_mixing_ratio.vapor,
            mixing_ratio_driver_liquid=self._locals.dry_air_mixing_ratio.liquid,  # driver_locals.dry_air_mixing_ratio.liquid,
            mixing_ratio_driver_rain=self._locals.dry_air_mixing_ratio.rain,  # driver_locals.dry_air_mixing_ratio.rain,
            mixing_ratio_driver_ice=self._locals.dry_air_mixing_ratio.ice,  # driver_locals.dry_air_mixing_ratio.ice,
            mixing_ratio_driver_snow=self._locals.dry_air_mixing_ratio.snow,  # driver_locals.dry_air_mixing_ratio.snow,
            mixing_ratio_driver_graupel=self._locals.dry_air_mixing_ratio.graupel,  # driver_locals.dry_air_mixing_ratio.graupel,
            dvapordt=dvapordt,  # locals_.driver_tendencies.dvapordt,
            dliquiddt=dliquiddt,  # locals_.driver_tendencies.dliquiddt,
            draindt=draindt,  # locals_.driver_tendencies.draindt,
            dicedt=dicedt,  # locals_.driver_tendencies.dicedt,
            dsnowdt=dsnowdt,  # locals_.driver_tendencies.dsnowdt,
            dgraupeldt=dgraupeldt,  # locals_.driver_tendencies.dgraupeldt,
            dcloudfractiondt=dcloudfractiondt,  # locals_.driver_tendencies.dcloudfractiondt,
            t_unmodified=t,  # state.t,
            driver_t=self._locals.t,  # driver_locals.t,
            dtdt=dtdt,  # locals_.driver_tendencies.dtdt,
            w_unmodified=w,  # state.vertical_motion.velocity,
            driver_w=self._locals.w,  # driver_locals.w,
            u_unmodified=u,  # state.u,
            driver_u=self._locals.u,  # driver_locals.u,
            dudt=dudt,  # locals_.driver_tendencies.dudt,
            v_unmodified=v,  # state.v,
            driver_v=self._locals.v,  # driver_locals.v,
            dvdt=dvdt,  # locals_.driver_tendencies.dvdt,
            dp_unmodified=dp,  # locals_.dp,
            driver_dp=self._locals.dp,  # driver_locals.dp,
            driver_mass=self._locals.mass,  # driver_locals.mass,
            rain=surface_precip_rain,  # state.precipitation_at_surface.rain,
            snow=surface_precip_snow,  # state.precipitation_at_surface.snow,
            ice=surface_precip_ice,  # state.precipitation_at_surface.ice,
            graupel=surface_precip_graupel,  # state.precipitation_at_surface.graupel,
        )
