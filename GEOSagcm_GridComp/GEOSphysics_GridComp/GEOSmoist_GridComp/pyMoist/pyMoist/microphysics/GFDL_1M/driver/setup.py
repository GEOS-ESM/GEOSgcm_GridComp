from ndsl import NDSLRuntime, StencilFactory, Quantity
from ndsl.constants import I_DIM, J_DIM, K_DIM
from ndsl.dsl.gt4py import FORWARD, PARALLEL, computation, function, interval, sqrt
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ

from pyMoist.microphysics.GFDL_1M.config import GFDL1MConfig
from pyMoist.microphysics.GFDL_1M.driver.config_constants import GFDL1MDriverConfigDependentConstants
from pyMoist.microphysics.GFDL_1M.driver.constants import constants
from pyMoist.shared.atmos_recipes import sigma


def init_temporaries(
    unmodified_t: FloatField,
    unmodified_dp: FloatField,
    critical_relative_humidity_for_pdf: FloatField,
    radiation_field_vapor: FloatField,
    radiation_field_liquid: FloatField,
    radiation_field_ice: FloatField,
    radiation_field_rain: FloatField,
    radiation_field_snow: FloatField,
    radiation_field_graupel: FloatField,
    radiation_field_cloud_fraction: FloatField,
    total_concentration: FloatField,
    unmodified_mixing_ratio_vapor: FloatField,
    unmodified_mixing_ratio_liquid: FloatField,
    unmodified_mixing_ratio_rain: FloatField,
    unmodified_mixing_ratio_ice: FloatField,
    unmodified_mixing_ratio_snow: FloatField,
    unmodified_mixing_ratio_graupel: FloatField,
    dry_air_mixing_ratio_vapor: FloatField,
    dry_air_mixing_ratio_liquid: FloatField,
    dry_air_mixing_ratio_rain: FloatField,
    dry_air_mixing_ratio_ice: FloatField,
    dry_air_mixing_ratio_snow: FloatField,
    dry_air_mixing_ratio_graupel: FloatField,
    cloud_fraction: FloatField,
    dz: FloatField,
    u_unmodified: FloatField,
    v_unmodified: FloatField,
    w_unmodified: FloatField,
    area: FloatFieldIJ,
    t: FloatField,
    dp: FloatField,
    density_unmodified: FloatField,
    p_dry: FloatField,
    mass: FloatField,
    u: FloatField,
    v: FloatField,
    w: FloatField,
    one_minus_sigma: FloatFieldIJ,
    ccn: FloatField,
    c_praut: FloatField,
    rh_limited: FloatField,
    rain: FloatFieldIJ,
    snow: FloatFieldIJ,
    graupel: FloatFieldIJ,
    ice: FloatFieldIJ,
    liquid_precip_flux: FloatField,
    ice_precip_flux: FloatField,
    evaporation: FloatField,
    sublimation: FloatField,
):
    """Initialize temporary copies of many fields. Data from these fields are not directly returned to the
    larger model, but it may be used to update another field which is returned to the model

    Args:
        unmodified_t (FloatField)
        unmodified_dp (FloatField)
        critical_relative_humidity_for_pdf (FloatField)
        radiation_field_vapor (FloatField)
        radiation_field_liquid (FloatField)
        radiation_field_ice (FloatField)
        radiation_field_rain (FloatField)
        radiation_field_snow (FloatField)
        radiation_field_graupel (FloatField)
        radiation_field_cloud_fraction (FloatField)
        total_concentration (FloatField)
        unmodified_mixing_ratio_vapor (FloatField)
        unmodified_mixing_ratio_liquid (FloatField)
        unmodified_mixing_ratio_rain (FloatField)
        unmodified_mixing_ratio_ice (FloatField)
        unmodified_mixing_ratio_snow (FloatField)
        unmodified_mixing_ratio_graupel (FloatField)
        dry_air_mixing_ratio_vapor (FloatField)
        dry_air_mixing_ratio_liquid (FloatField)
        dry_air_mixing_ratio_rain (FloatField)
        dry_air_mixing_ratio_ice (FloatField)
        dry_air_mixing_ratio_snow (FloatField)
        dry_air_mixing_ratio_graupel (FloatField)
        cloud_fraction (FloatField)
        dz (FloatField)
        u_unmodified (FloatField)
        v_unmodified (FloatField)
        w_unmodified (FloatField)
        area (FloatFieldIJ)
        t (FloatField)
        dp (FloatField)
        density_unmodified (FloatField)
        p_dry (FloatField)
        mass (FloatField)
        u (FloatField)
        v (FloatField)
        w (FloatField)
        one_minus_sigma (FloatFieldIJ)
        ccn (FloatField)
        c_praut (FloatField)
        rh_limited (FloatField)
        rain (FloatFieldIJ)
        snow (FloatFieldIJ)
        graupel (FloatFieldIJ)
        ice (FloatFieldIJ)
        liquid_precip_flux (FloatField)
        ice_precip_flux (FloatField)
        evaporation (FloatField)
        sublimation (FloatField)
    """
    from __externals__ import DO_SEDI_W, cpaut

    with computation(PARALLEL), interval(...):
        t = unmodified_t
        dp = unmodified_dp  # moist air mass * grav

        # -----------------------------------------------------------------------
        # import horizontal subgrid variability with pressure dependence
        # total water subgrid deviation in horizontal direction
        # default area dependent form: use dx ~ 100 km as the base
        # -----------------------------------------------------------------------
        rh_limited = min(0.30, 1.0 - critical_relative_humidity_for_pdf)  # restricted to 70%

        # -----------------------------------------------------------------------
        # convert moist mixing ratios to dry mixing ratios
        # -----------------------------------------------------------------------

        dp = dp * (1.0 - radiation_field_vapor)  # gfs
        omq = unmodified_dp / dp

        unmodified_mixing_ratio_vapor = radiation_field_vapor * omq
        unmodified_mixing_ratio_liquid = radiation_field_liquid * omq
        unmodified_mixing_ratio_rain = radiation_field_rain * omq
        unmodified_mixing_ratio_ice = radiation_field_ice * omq
        unmodified_mixing_ratio_snow = radiation_field_snow * omq
        unmodified_mixing_ratio_graupel = radiation_field_graupel * omq

        dry_air_mixing_ratio_vapor = unmodified_mixing_ratio_vapor
        dry_air_mixing_ratio_liquid = unmodified_mixing_ratio_liquid
        dry_air_mixing_ratio_rain = unmodified_mixing_ratio_rain
        dry_air_mixing_ratio_ice = unmodified_mixing_ratio_ice
        dry_air_mixing_ratio_snow = unmodified_mixing_ratio_snow
        dry_air_mixing_ratio_graupel = unmodified_mixing_ratio_graupel

        cloud_fraction = radiation_field_cloud_fraction

        density_unmodified = -dp / (constants.GRAV * dz)  # density of dry air
        p_dry = density_unmodified * constants.RDGAS * unmodified_t  # dry air pressure

        # -----------------------------------------------------------------------
        # for sedi_momentum
        # -----------------------------------------------------------------------

        mass = 0.0
        u = u_unmodified
        v = v_unmodified

        if DO_SEDI_W:
            w = w_unmodified

        # ccn needs units #/m^3
        ccn = total_concentration
        c_praut = cpaut * (ccn * constants.RHOR) ** (-1.0 / 3.0)

        # Reset precipitation aggregates to zero
        liquid_precip_flux = 0
        ice_precip_flux = 0
        evaporation = 0
        sublimation = 0

    with computation(FORWARD), interval(0, 1):
        # 1 minus sigma used to control minimum cloud
        # fraction needed to autoconvert ql->qr
        one_minus_sigma = 1.0 - sigma(sqrt(area))

        # Reset precipitation aggregates to zero
        rain = 0
        snow = 0
        graupel = 0
        ice = 0


@function
def fix_negative_core(
    t: Float,
    mixing_ratio_vapor: Float,
    mixing_ratio_liquid: Float,
    mixing_ratio_rain: Float,
    mixing_ratio_ice: Float,
    mixing_ratio_snow: Float,
    mixing_ratio_graupel: Float,
    c_air: Float,
    c_vap: Float,
    lv00: Float,
    d0_vap: Float,
):
    """Adjusts/removes negative mixing ratios

    reference Fortran: gfdl_cloud_microphys.F90: subroutine neg_adj

    Args:
        t (Float): temperature (Kelvin)
        mixing_ratio_vapor (Float): water vapor mixing ratio (kg/kg)
        mixing_ratio_liquid (Float): liquid water mixing ratio (kg/kg)
        mixing_ratio_rain (Float): rain mixing ratio (kg/kg)
        mixing_ratio_ice (Float): ice mixing ratio (kg/kg)
        mixing_ratio_snow (Float): snow mixing ratio (kg/kg)
        mixing_ratio_graupel (Float): graupel mixing ratio (kg/kg)
        c_air (Float)
        c_vap (Float)
        lv00 (Float)
        d0_vap (Float)

    Returns:
        Float: t
        Float: mixing_ratio_vapor
        Float: mixing_ratio_liquid
        Float: mixing_ratio_rain
        Float: mixing_ratio_ice
        Float: mixing_ratio_snow
        Float: mixing_ratio_graupel
    """
    # -----------------------------------------------------------------------
    # define heat capacity and latent heat coefficient
    # -----------------------------------------------------------------------

    cvm = (
        c_air
        + mixing_ratio_vapor * c_vap
        + (mixing_ratio_rain + mixing_ratio_liquid) * constants.C_LIQ
        + (mixing_ratio_ice + mixing_ratio_snow + mixing_ratio_graupel) * constants.C_ICE
    )
    lcpk = (lv00 + d0_vap * t) / cvm
    icpk = (constants.LI00 + constants.DC_ICE * t) / cvm

    # -----------------------------------------------------------------------
    # ice phase:
    # -----------------------------------------------------------------------

    # if cloud ice < 0, borrow from snow
    if mixing_ratio_ice < 0.0:
        mixing_ratio_snow = mixing_ratio_snow + mixing_ratio_ice
        mixing_ratio_ice = 0.0
    # if snow < 0, borrow from graupel
    if mixing_ratio_snow < 0.0:
        mixing_ratio_graupel = mixing_ratio_graupel + mixing_ratio_snow
        mixing_ratio_snow = 0.0
    # if graupel < 0, borrow from rain
    if mixing_ratio_graupel < 0.0:
        mixing_ratio_rain = mixing_ratio_rain + mixing_ratio_graupel
        t = t - mixing_ratio_graupel * icpk  # heating
        mixing_ratio_graupel = 0.0

    # -----------------------------------------------------------------------
    # liquid phase:
    # -----------------------------------------------------------------------

    # if rain < 0, borrow from cloud water
    if mixing_ratio_rain < 0.0:
        mixing_ratio_liquid = mixing_ratio_liquid + mixing_ratio_rain
        mixing_ratio_rain = 0.0
    # if cloud water < 0, borrow from water vapor
    if mixing_ratio_liquid < 0.0:
        mixing_ratio_vapor = mixing_ratio_vapor + mixing_ratio_liquid
        t = t - mixing_ratio_liquid * lcpk  # heating
        mixing_ratio_liquid = 0.0

    return (
        t,
        mixing_ratio_vapor,
        mixing_ratio_liquid,
        mixing_ratio_rain,
        mixing_ratio_ice,
        mixing_ratio_snow,
        mixing_ratio_graupel,
    )


def fix_negative_values(
    t: FloatField,
    dry_air_mixing_ratio_vapor: FloatField,
    dry_air_mixing_ratio_liquid: FloatField,
    dry_air_mixing_ratio_rain: FloatField,
    dry_air_mixing_ratio_ice: FloatField,
    dry_air_mixing_ratio_snow: FloatField,
    dry_air_mixing_ratio_graupel: FloatField,
    dp: FloatField,
):
    """Adjusts/removes negative mixing ratios and updates vapor and temperature
    where necessary. Core math is done in fix_negative_core

    reference Fortran: gfdl_cloud_microphys.F90: subroutine mpdrv

    Args:
        t (FloatField): temperature (Kelvin)
        dry_air_mixing_ratio_vapor (FloatField): water vapor mixing ratio (kg/kg)
        dry_air_mixing_ratio_liquid (FloatField): liquid water mixing ratio (kg/kg)
        dry_air_mixing_ratio_rain (FloatField): rain mixing ratio (kg/kg)
        dry_air_mixing_ratio_ice (FloatField): ice mixing ratio (kg/kg)
        dry_air_mixing_ratio_snow (FloatField): snow mixing ratio (kg/kg)
        dry_air_mixing_ratio_graupel (FloatField): graupel mixing ratio (kg/kg)
        dp (FloatField): pressure change between model layers (Pa)
    """
    from __externals__ import c_air, c_vap, d0_vap, lv00

    # -----------------------------------------------------------------------
    # fix all negative water species
    # -----------------------------------------------------------------------

    with computation(FORWARD), interval(0, -1):
        (
            t,
            dry_air_mixing_ratio_vapor,
            dry_air_mixing_ratio_liquid,
            dry_air_mixing_ratio_rain,
            dry_air_mixing_ratio_ice,
            dry_air_mixing_ratio_snow,
            dry_air_mixing_ratio_graupel,
        ) = fix_negative_core(
            t,
            dry_air_mixing_ratio_vapor,
            dry_air_mixing_ratio_liquid,
            dry_air_mixing_ratio_rain,
            dry_air_mixing_ratio_ice,
            dry_air_mixing_ratio_snow,
            dry_air_mixing_ratio_graupel,
            c_air,
            c_vap,
            lv00,
            d0_vap,
        )
        if dry_air_mixing_ratio_vapor < 0.0:
            dry_air_mixing_ratio_vapor[0, 0, 1] = (
                dry_air_mixing_ratio_vapor[0, 0, 1] + dry_air_mixing_ratio_vapor * dp / dp[0, 0, 1]
            )
            dry_air_mixing_ratio_vapor = 0.0

    with computation(FORWARD), interval(-1, None):
        (
            t,
            dry_air_mixing_ratio_vapor,
            dry_air_mixing_ratio_liquid,
            dry_air_mixing_ratio_rain,
            dry_air_mixing_ratio_ice,
            dry_air_mixing_ratio_snow,
            dry_air_mixing_ratio_graupel,
        ) = fix_negative_core(
            t,
            dry_air_mixing_ratio_vapor,
            dry_air_mixing_ratio_liquid,
            dry_air_mixing_ratio_rain,
            dry_air_mixing_ratio_ice,
            dry_air_mixing_ratio_snow,
            dry_air_mixing_ratio_graupel,
            c_air,
            c_vap,
            lv00,
            d0_vap,
        )

        if dry_air_mixing_ratio_vapor < 0.0 and dry_air_mixing_ratio_vapor[0, 0, -1] > 0.0:
            dq = min(-dry_air_mixing_ratio_vapor * dp, dry_air_mixing_ratio_vapor[0, 0, -1] * dp[0, 0, -1])
            dry_air_mixing_ratio_vapor[0, 0, -1] = dry_air_mixing_ratio_vapor[0, 0, -1] - dq / dp[0, 0, -1]
            dry_air_mixing_ratio_vapor = dry_air_mixing_ratio_vapor + dq / dp


class GFDL1MDriverSetup(NDSLRuntime):
    def __init__(
        self,
        stencil_factory: StencilFactory,
        config: GFDL1MConfig,
        config_dependent_constants: GFDL1MDriverConfigDependentConstants,
    ):
        """The driver modifies a number of variables (t, p, qX) but does not pass
        the changes back to the rest of the model. To replicate this behavior,
        temporary copies of these variables are used throughout the driver.

        Args:
            stencil_factory (StencilFactory)
            config (GFDL1MConfig)
            config_dependent_constants (GFDL1MDriverConfigDependentConstants)
        """
        # init NDSLRuntime
        super().__init__(stencil_factory)

        # construct stencils
        self._init_temporaries = stencil_factory.from_dims_halo(
            func=init_temporaries,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "DO_SEDI_W": config.DO_SEDI_W,
                "cpaut": config_dependent_constants.CPAUT,
            },
        )

        self._fix_negative_values = stencil_factory.from_dims_halo(
            func=fix_negative_values,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "c_air": config_dependent_constants.C_AIR,
                "c_vap": config_dependent_constants.C_VAP,
                "d0_vap": config_dependent_constants.D0_VAP,
                "lv00": config_dependent_constants.LV00,
            },
        )

    def __call__(
        self,
        unmodified_t: Quantity,
        t: Quantity,
        unmodified_dp: Quantity,
        dp: Quantity,
        critical_relative_humidity_for_pdf: Quantity,
        radiation_field_vapor: Quantity,
        radiation_field_liquid: Quantity,
        radiation_field_ice: Quantity,
        radiation_field_rain: Quantity,
        radiation_field_snow: Quantity,
        radiation_field_graupel: Quantity,
        radiation_field_cloud_fraction: Quantity,
        total_concentration: Quantity,
        unmodified_mixing_ratio_vapor: Quantity,
        unmodified_mixing_ratio_liquid: Quantity,
        unmodified_mixing_ratio_rain: Quantity,
        unmodified_mixing_ratio_ice: Quantity,
        unmodified_mixing_ratio_snow: Quantity,
        unmodified_mixing_ratio_graupel: Quantity,
        dry_air_mixing_ratio_vapor: Quantity,
        dry_air_mixing_ratio_liquid: Quantity,
        dry_air_mixing_ratio_rain: Quantity,
        dry_air_mixing_ratio_ice: Quantity,
        dry_air_mixing_ratio_snow: Quantity,
        dry_air_mixing_ratio_graupel: Quantity,
        cloud_fraction: Quantity,
        dz: Quantity,
        u_unmodified: Quantity,
        u: Quantity,
        v_unmodified: Quantity,
        v: Quantity,
        w_unmodified: Quantity,
        w: Quantity,
        area: Quantity,
        density_unmodified: Quantity,
        p_dry: Quantity,
        mass: Quantity,
        one_minus_sigma: Quantity,
        ccn: Quantity,
        c_praut: Quantity,
        rh_limited: Quantity,
        rain: Quantity,
        snow: Quantity,
        graupel: Quantity,
        ice: Quantity,
        liquid_precip_flux: Quantity,
        ice_precip_flux: Quantity,
        evaporation: Quantity,
        sublimation: Quantity,
    ):
        """Setup the driver: initialize/prefill locals and ensure there are no negative mixing ratios.

        Args:
            unmodified_t (Quantity): (in) temperature from the model state (K), unmodified within driver
            t (Quantity): (out) copy of temperature (K) - can be modified within driver
            unmodified_dp (Quantity): (in) pressure between model layers from model state (Pa), unmodified
                within driver
            dp (Quantity): (out) copy of pressure between model layers (Pa) - can be modified within driver
            critical_relative_humidity_for_pdf (Quantity): (in) critical relative humidity
            radiation_field_vapor (Quantity): (in) water vapor mixing ratio - used in radiation scheme (kg/kg)
            radiation_field_liquid (Quantity): (in) liquid water mixing ratio - used in radiation scheme (kg/kg)
            radiation_field_ice (Quantity): (in) ice mixing ratio - used in radiation scheme (kg/kg)
            radiation_field_rain (Quantity): (in) rain mixing ratio - used in radiation scheme (kg/kg)
            radiation_field_snow (Quantity): (in) snow mixing ratio - used in radiation scheme (kg/kg)
            radiation_field_graupel (Quantity): (in) graupel mixing ratio - used in radiation scheme (kg/kg)
            radiation_field_cloud_fraction (Quantity): (in) cloud fraction - used in radiation scheme
            total_concentration (Quantity): (in) total liquid + ice concentration (m^-3)
            unmodified_mixing_ratio_vapor (Quantity): (out) unit-converted copy of radiation_field_vapor
                mixing ratio (kg/kg), unmodified within the remaining driver
            unmodified_mixing_ratio_liquid (Quantity): (out) unit-converted copy of radiation_field_liquid
                mixing ratio (kg/kg), unmodified within the remaining driver
            unmodified_mixing_ratio_rain (Quantity): (out) unit-converted copy of radiation_field_rain
                mixing ratio (kg/kg), unmodified within the remaining driver
            unmodified_mixing_ratio_ice (Quantity): (out) unit-converted copy of radiation_field_ice
                mixing ratio (kg/kg), unmodified within the remaining driver
            unmodified_mixing_ratio_snow (Quantity): (out) unit-converted copy of radiation_field_snow
                mixing ratio (kg/kg), unmodified within the remaining driver
            unmodified_mixing_ratio_graupel (Quantity): (out) unit-converted copy of radiation_field_graupel
                mixing ratio (kg/kg), unmodified within the remaining driver
            dry_air_mixing_ratio_vapor (Quantity): (out) copy of radiation_field_vapor (kg/kg)
            dry_air_mixing_ratio_liquid (Quantity): (out) copy of radiation_field_liquid (kg/kg)
            dry_air_mixing_ratio_rain (Quantity): (out) copy of radiation_field_rain (kg/kg)
            dry_air_mixing_ratio_ice (Quantity): (out) copy of radiation_field_ice (kg/kg)
            dry_air_mixing_ratio_snow (Quantity): (out) copy of radiation_field_snow (kg/kg)
            dry_air_mixing_ratio_graupel (Quantity): (out) copy of radiation_field_graupel (kg/kg)
            cloud_fraction (Quantity): (out) copy of radiation_field_cloud_fraction
            dz (Quantity): (in) height between model layer (m)
            u_unmodified (Quantity): (in) zonal wind (m/s), unmodified within driver
            u (Quantity): (out) zonal wind (m), may be modified within driver
            v_unmodified (Quantity): (in) meridional wind (m/s), unmodified within driver
            v (Quantity): (out) meridional wind (m), may be modified within driver
            w_unmodified (Quantity): (in) vertical motion (m/s), unmodified within driver
            w (Quantity): (out) vertical motion (m/s), may be modified within driver
            area (Quantity): (in) area of grid cell
            density_unmodified (Quantity): (out) density of grid cell, unmodified in the rest of the driver
            p_dry (Quantity): (out) dry air pressure
            mass (Quantity): (out) mass of grid call
            one_minus_sigma (Quantity): (out) details unknown
            ccn (Quantity): (out) cloud condensation nuclei of a grid cell (m^-3)
            c_praut (Quantity): (out) details unknown
            rh_limited (Quantity): (out) relative humidity, with limits applied
            rain (Quantity): (out) precipitated rain (kg m^-2 s^-1)
            snow (Quantity): (out) precipitated snow (kg m^-2 s^-1)
            graupel (Quantity): (out) precipitated graupel (kg m^-2 s^-1)
            ice (Quantity): (out) precipitated ice (kg m^-2 s^-1)
            liquid_precip_flux (Quantity): (out) non-anvil large scale liquid precip flux (kg m^-2 s^-1)
            ice_precip_flux (Quantity): (out) non-anvil large scale ice precip flux (kg m^-2 s^-1)
            evaporation (Quantity): (out) non-anvil large scale evaporation (kg kg^-1 s^-1)
            sublimation (Quantity): (out) non-anvil large scale sublimation (kg kg^-1 s^-1)
        """
        # The driver modifies a number of variables but does not pass the changes back to
        # the rest of the model. To replicate this behavior, temporary copies of these
        # variables are used throughout the driver. Prefill them here.
        self._init_temporaries(
            unmodified_t=unmodified_t,
            unmodified_dp=unmodified_dp,
            critical_relative_humidity_for_pdf=critical_relative_humidity_for_pdf,
            radiation_field_vapor=radiation_field_vapor,
            radiation_field_liquid=radiation_field_liquid,
            radiation_field_ice=radiation_field_ice,
            radiation_field_rain=radiation_field_rain,
            radiation_field_snow=radiation_field_snow,
            radiation_field_graupel=radiation_field_graupel,
            radiation_field_cloud_fraction=radiation_field_cloud_fraction,
            total_concentration=total_concentration,
            unmodified_mixing_ratio_vapor=unmodified_mixing_ratio_vapor,
            unmodified_mixing_ratio_liquid=unmodified_mixing_ratio_liquid,
            unmodified_mixing_ratio_rain=unmodified_mixing_ratio_rain,
            unmodified_mixing_ratio_ice=unmodified_mixing_ratio_ice,
            unmodified_mixing_ratio_snow=unmodified_mixing_ratio_snow,
            unmodified_mixing_ratio_graupel=unmodified_mixing_ratio_graupel,
            dry_air_mixing_ratio_vapor=dry_air_mixing_ratio_vapor,
            dry_air_mixing_ratio_liquid=dry_air_mixing_ratio_liquid,
            dry_air_mixing_ratio_rain=dry_air_mixing_ratio_rain,
            dry_air_mixing_ratio_ice=dry_air_mixing_ratio_ice,
            dry_air_mixing_ratio_snow=dry_air_mixing_ratio_snow,
            dry_air_mixing_ratio_graupel=dry_air_mixing_ratio_graupel,
            cloud_fraction=cloud_fraction,
            dz=dz,
            u_unmodified=u_unmodified,
            v_unmodified=v_unmodified,
            w_unmodified=w_unmodified,
            area=area,
            t=t,
            dp=dp,
            density_unmodified=density_unmodified,
            p_dry=p_dry,
            mass=mass,
            u=u,
            v=v,
            w=w,
            one_minus_sigma=one_minus_sigma,
            ccn=ccn,
            c_praut=c_praut,
            rh_limited=rh_limited,
            rain=rain,
            snow=snow,
            graupel=graupel,
            ice=ice,
            liquid_precip_flux=liquid_precip_flux,
            ice_precip_flux=ice_precip_flux,
            evaporation=evaporation,
            sublimation=sublimation,
        )

        self._fix_negative_values(
            t=t,
            dry_air_mixing_ratio_vapor=dry_air_mixing_ratio_vapor,
            dry_air_mixing_ratio_liquid=dry_air_mixing_ratio_liquid,
            dry_air_mixing_ratio_rain=dry_air_mixing_ratio_rain,
            dry_air_mixing_ratio_ice=dry_air_mixing_ratio_ice,
            dry_air_mixing_ratio_snow=dry_air_mixing_ratio_snow,
            dry_air_mixing_ratio_graupel=dry_air_mixing_ratio_graupel,
            dp=dp,
        )
