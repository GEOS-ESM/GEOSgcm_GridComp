"""GFDL_1M driver"""

import dace
from ndsl import NDSLRuntime, Quantity, QuantityFactory, StencilFactory
from ndsl.constants import I_DIM, J_DIM, K_DIM

from pyMoist.microphysics.GFDL_1M.config import GFDL1MConfig
from pyMoist.microphysics.GFDL_1M.driver.check_flags import check_flags
from pyMoist.microphysics.GFDL_1M.driver.config_constants import GFDL1MDriverConfigDependentConstants
from pyMoist.microphysics.GFDL_1M.driver.fall_speed import fall_speed
from pyMoist.microphysics.GFDL_1M.driver.finish import update_tendencies
from pyMoist.microphysics.GFDL_1M.driver.ice_cloud import GFDL1MIceCloud
from pyMoist.microphysics.GFDL_1M.driver.locals import GFDL1MDriverLocals
from pyMoist.microphysics.GFDL_1M.driver.sat_tables import get_tables
from pyMoist.microphysics.GFDL_1M.driver.setup import GFDL1MDriverSetup
from pyMoist.microphysics.GFDL_1M.driver.terminal_fall import GFDL1MTerminalFall
from pyMoist.microphysics.GFDL_1M.driver.warm_rain import GFDL1MWarmRain


class GFDL1MDriver(NDSLRuntime):
    """
    Computes precipitates and microphysics using the Geophysical Fluid Dynamics Labratory Single Moment
    Microphysics package tendencies using the following functions:
    __init__:
        - checks constants for options which will require unimplemented options
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
        """Perform setup for the microphysics driver. Check flags for unimplemented options,
        initialize internal fields, and compile stencils.

        Args:
            stencil_factory (StencilFactory): StencilFactory with model domain information
            quantity_factory (QuantityFactory): QuantityFactory with model domain information
            config (GFDL1MConfig): driver configuration
        """
        # init NDSLRuntime
        super().__init__(stencil_factory)

        self.config_dependent_constants = GFDL1MDriverConfigDependentConstants.make(config)

        # Check constants for unimplemented and untested code paths
        check_flags(
            config,
            self.config_dependent_constants,
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
            compute_dims=[I_DIM, J_DIM, K_DIM],
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
            compute_dims=[I_DIM, J_DIM, K_DIM],
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
        t: Quantity,
        u: Quantity,
        v: Quantity,
        w: Quantity,
        dz: Quantity,
        dp: Quantity,
        area: Quantity,
        land_fraction: Quantity,
        convection_fraction: Quantity,
        surface_type: Quantity,
        estimated_inversion_strength: Quantity,
        critical_relative_humidity_for_pdf: Quantity,
        vapor: Quantity,
        liquid: Quantity,
        rain: Quantity,
        ice: Quantity,
        snow: Quantity,
        graupel: Quantity,
        cloud_fraction: Quantity,
        total_concentration: Quantity,
        dvapordt: Quantity,
        dliquiddt: Quantity,
        draindt: Quantity,
        dicedt: Quantity,
        dsnowdt: Quantity,
        dgraupeldt: Quantity,
        dcloudfractiondt: Quantity,
        dtdt: Quantity,
        dudt: Quantity,
        dvdt: Quantity,
        liquid_precip_flux: Quantity,
        ice_precip_flux: Quantity,
        evaporation: Quantity,
        sublimation: Quantity,
        surface_precip_rain: Quantity,
        surface_precip_snow: Quantity,
        surface_precip_ice: Quantity,
        surface_precip_graupel: Quantity,
    ):
        """

        Arguments:
            t (Quantity): (in) atmospheric temperature (K)
            u (Quantity): (in) eastward winds (m/s)
            v (Quantity): (in) northward winds (m/s)
            w (Quantity): (in) vertical velocity (m/s)
            dz (Quantity): (in) layer thickness (m)
            dp (Quantity): (in) change in pressure between model levels (Pa)
            area (Quantity): (in) grid cell area
            land_fraction (Quantity): (in) land fraction
            convection_fraction (Quantity): (in) convection fraction
            surface_type (Quantity): (in) surface type
            estimated_inversion_strength (Quantity): (in) estimated inversion strength (K)
            critical_relative_humidity_for_pdf (Quantity): (in) critical relative humidity for pdf
            vapor:(Quantity):  (in: water vapor mixing ratio (kg/kg)
            liquid (Quantity): (in) in cloud liquid mixing radio (kg/kg)
            rain (Quantity): (in) falling rain (kg/kg)
            ice (Quantity): (in) in cloud ice mixing radio (kg/kg)
            snow (Quantity): (in) in cloud snow mixing radio (kg/kg)
            graupel (Quantity): (in) in cloud graupel mixing radio (kg/kg)
            cloud_fraction (Quantity): (in) cloud fraction (convective + large scale)
            total_concentration (Quantity): (in) total ice + liquid concentration (m^-3)
            dvapordt (Quantity): (out: water vapor tendency
            dliquiddt (Quantity): (out: in cloud liquid water tendency
            draindt (Quantity): (out: falling rain tendency
            dicedt (Quantity): (out: in cloud frozen water tendency
            dsnowdt (Quantity): (out: in cloud snow tendency
            dgraupeldt (Quantity): (out: in cloud graupel tendency
            dcloudfractiondt (Quantity): (out: cloud fraction (convective + large scale) tendency
            dtdt (Quantity): (out: atmospheric temperature tendency
            dudt (Quantity): (out: eastward wind tendency
            dvdt (Quantity): (out: northward wind tendency
            liquid_precip_flux (Quantity): (out: non-anvil large scale liquid precip flux
            ice_precip_flux (Quantity): (out: non-anvil large scale ice precip flux
            evaporation (Quantity): (out: non-anvil large scale evaporation
            sublimation (Quantity): (out: non-anvil large scale sublimation
            surface_precip_rain (Quantity): (out: rain precip at surface (kg/m^2/s)
            surface_precip_snow (Quantity): (out: snow precip at surface (kg/m^2/s)
            surface_precip_ice (Quantity): (out: ice precip at surface (kg/m^2/s)
            surface_precip_graupel (Quantity): (out: graupel precip at surface (kg/m^2/s)
        """
        self._setup(
            unmodified_t=t,
            t=self._locals.t,
            unmodified_dp=dp,
            dp=self._locals.dp,
            critical_relative_humidity_for_pdf=critical_relative_humidity_for_pdf,
            radiation_field_vapor=vapor,
            radiation_field_liquid=liquid,
            radiation_field_ice=ice,
            radiation_field_rain=rain,
            radiation_field_snow=snow,
            radiation_field_graupel=graupel,
            radiation_field_cloud_fraction=cloud_fraction,
            total_concentration=total_concentration,
            unmodified_mixing_ratio_vapor=self._locals.unmodified.mixing_ratio.vapor,
            unmodified_mixing_ratio_liquid=self._locals.unmodified.mixing_ratio.liquid,
            unmodified_mixing_ratio_rain=self._locals.unmodified.mixing_ratio.rain,
            unmodified_mixing_ratio_ice=self._locals.unmodified.mixing_ratio.ice,
            unmodified_mixing_ratio_snow=self._locals.unmodified.mixing_ratio.snow,
            unmodified_mixing_ratio_graupel=self._locals.unmodified.mixing_ratio.graupel,
            dry_air_mixing_ratio_vapor=self._locals.dry_air_mixing_ratio.vapor,
            dry_air_mixing_ratio_liquid=self._locals.dry_air_mixing_ratio.liquid,
            dry_air_mixing_ratio_rain=self._locals.dry_air_mixing_ratio.rain,
            dry_air_mixing_ratio_ice=self._locals.dry_air_mixing_ratio.ice,
            dry_air_mixing_ratio_snow=self._locals.dry_air_mixing_ratio.snow,
            dry_air_mixing_ratio_graupel=self._locals.dry_air_mixing_ratio.graupel,
            cloud_fraction=self._locals.cloud_fraction,
            dz=dz,
            u_unmodified=u,
            u=self._locals.u,
            v_unmodified=v,
            v=self._locals.v,
            w_unmodified=w,
            w=self._locals.w,
            area=area,
            density_unmodified=self._locals.density_unmodified,
            p_dry=self._locals.p_dry,
            mass=self._locals.mass,
            one_minus_sigma=self._locals.one_minus_sigma,
            ccn=self._locals.ccn,
            c_praut=self._locals.c_praut,
            rh_limited=self._locals.rh_limited,
            rain=surface_precip_rain,
            snow=surface_precip_snow,
            graupel=surface_precip_graupel,
            ice=surface_precip_ice,
            liquid_precip_flux=liquid_precip_flux,
            ice_precip_flux=ice_precip_flux,
            evaporation=evaporation,
            sublimation=sublimation,
        )

        for _ in dace.nounroll(range(self.config_dependent_constants.NTIMES)):
            self._fall_speed(
                mixing_ratio_liquid=self._locals.dry_air_mixing_ratio.liquid,
                mixing_ratio_ice=self._locals.dry_air_mixing_ratio.ice,
                mixing_ratio_snow=self._locals.dry_air_mixing_ratio.snow,
                mixing_ratio_graupel=self._locals.dry_air_mixing_ratio.graupel,
                t_unmodified=t,
                t=self._locals.t,
                dz_unmodified=dz,
                dz=self._locals.dz,
                density_unmodified=self._locals.density_unmodified,
                density=self._locals.density,
                density_factor=self._locals.density_factor,
                ice_terminal_velocity=self._locals.terminal_speed.ice,
                snow_terminal_velocity=self._locals.terminal_speed.snow,
                graupel_terminal_velocity=self._locals.terminal_speed.graupel,
                convection_fraction=convection_fraction,
            )

            self._terminal_fall(
                t=self._locals.t,
                w=self._locals.w,
                mixing_ratio_vapor=self._locals.dry_air_mixing_ratio.vapor,
                mixing_ratio_liquid=self._locals.dry_air_mixing_ratio.liquid,
                mixing_ratio_rain=self._locals.dry_air_mixing_ratio.rain,
                mixing_ratio_graupel=self._locals.dry_air_mixing_ratio.graupel,
                mixing_ratio_snow=self._locals.dry_air_mixing_ratio.snow,
                mixing_ratio_ice=self._locals.dry_air_mixing_ratio.ice,
                dz=self._locals.dz,
                dp=self._locals.dp,
                terminal_velocity_graupel=self._locals.terminal_speed.graupel,
                terminal_velocity_snow=self._locals.terminal_speed.snow,
                terminal_velocity_ice=self._locals.terminal_speed.ice,
                rain=surface_precip_rain,
                graupel=surface_precip_graupel,
                snow=surface_precip_snow,
                ice=surface_precip_ice,
                ice_precip_flux=self._locals.ice_precip_flux,
            )

            self._warm_rain(
                t=self._locals.t,
                dp=self._locals.dp,
                dz=self._locals.dz,
                w=self._locals.w,
                mixing_ratio_vapor=self._locals.dry_air_mixing_ratio.vapor,
                mixing_ratio_liquid=self._locals.dry_air_mixing_ratio.liquid,
                mixing_ratio_rain=self._locals.dry_air_mixing_ratio.rain,
                mixing_ratio_ice=self._locals.dry_air_mixing_ratio.ice,
                mixing_ratio_snow=self._locals.dry_air_mixing_ratio.snow,
                mixing_ratio_graupel=self._locals.dry_air_mixing_ratio.graupel,
                cloud_fraction=self._locals.cloud_fraction,
                ccn=self._locals.ccn,
                density=self._locals.density,
                density_factor=self._locals.density_factor,
                c_praut=self._locals.c_praut,
                terminal_speed_rain=self._locals.terminal_speed.rain,
                rh_limited=self._locals.rh_limited,
                estimated_inversion_strength=estimated_inversion_strength,
                one_minus_sigma=self._locals.one_minus_sigma,
                mass=self._locals.mass,
                rain=surface_precip_rain,
                driver_rain=self._locals.rain,
                ice_precip_flux=ice_precip_flux,
                driver_ice_precip_flux=self._locals.ice_precip_flux,
                liquid_precip_flux=liquid_precip_flux,
                driver_liquid_precip_flux=self._locals.liquid_precip_flux,
                evaporation=evaporation,
                driver_evaporation=self._locals.evaporation,
            )

            self._ice_cloud(
                t=self._locals.t,
                p_dry=self._locals.p_dry,
                dp=self._locals.dp,
                mixing_ratio_vapor=self._locals.dry_air_mixing_ratio.vapor,
                mixing_ratio_liquid=self._locals.dry_air_mixing_ratio.liquid,
                mixing_ratio_rain=self._locals.dry_air_mixing_ratio.rain,
                mixing_ratio_ice=self._locals.dry_air_mixing_ratio.ice,
                mixing_ratio_snow=self._locals.dry_air_mixing_ratio.snow,
                mixing_ratio_graupel=self._locals.dry_air_mixing_ratio.graupel,
                cloud_fraction=self._locals.cloud_fraction,
                density=self._locals.density,
                density_factor=self._locals.density_factor,
                terminal_fall_snow=self._locals.terminal_speed.snow,
                terminal_fall_graupel=self._locals.terminal_speed.graupel,
                terminal_fall_rain=self._locals.terminal_speed.rain,
                sublimation=sublimation,
                rh_limited=self._locals.rh_limited,
                ccn=self._locals.ccn,
                convection_fraction=convection_fraction,
                surface_type=surface_type,
            )

        self._finish(
            mixing_ratio_vapor_unmodified=self._locals.unmodified.mixing_ratio.vapor,
            mixing_ratio_liquid_unmodified=self._locals.unmodified.mixing_ratio.liquid,
            mixing_ratio_rain_unmodified=self._locals.unmodified.mixing_ratio.rain,
            mixing_ratio_ice_unmodified=self._locals.unmodified.mixing_ratio.ice,
            mixing_ratio_snow_unmodified=self._locals.unmodified.mixing_ratio.snow,
            mixing_ratio_graupel_unmodified=self._locals.unmodified.mixing_ratio.graupel,
            cloud_fraction_unmodified=cloud_fraction,
            mixing_ratio_driver_vapor=self._locals.dry_air_mixing_ratio.vapor,
            mixing_ratio_driver_liquid=self._locals.dry_air_mixing_ratio.liquid,
            mixing_ratio_driver_rain=self._locals.dry_air_mixing_ratio.rain,
            mixing_ratio_driver_ice=self._locals.dry_air_mixing_ratio.ice,
            mixing_ratio_driver_snow=self._locals.dry_air_mixing_ratio.snow,
            mixing_ratio_driver_graupel=self._locals.dry_air_mixing_ratio.graupel,
            dvapordt=dvapordt,
            dliquiddt=dliquiddt,
            draindt=draindt,
            dicedt=dicedt,
            dsnowdt=dsnowdt,
            dgraupeldt=dgraupeldt,
            dcloudfractiondt=dcloudfractiondt,
            t_unmodified=t,
            driver_t=self._locals.t,
            dtdt=dtdt,
            w_unmodified=w,
            driver_w=self._locals.w,
            u_unmodified=u,
            driver_u=self._locals.u,
            dudt=dudt,
            v_unmodified=v,
            driver_v=self._locals.v,
            dvdt=dvdt,
            dp_unmodified=dp,
            driver_dp=self._locals.dp,
            driver_mass=self._locals.mass,
            rain=surface_precip_rain,
            snow=surface_precip_snow,
            ice=surface_precip_ice,
            graupel=surface_precip_graupel,
        )
