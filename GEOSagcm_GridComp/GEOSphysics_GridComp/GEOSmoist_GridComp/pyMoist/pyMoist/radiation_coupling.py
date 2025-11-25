from ndsl import StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.gt4py import PARALLEL, computation, interval
from ndsl.dsl.typing import Float, FloatField
from pyMoist.saturation_tables import (
    GlobalTable_saturation_tables,
    saturation_specific_humidity,
    SaturationVaporPressureTable,
)
from pyMoist.shared_incloud_processes import (
    cloud_effective_radius_ice,
    cloud_effective_radius_liquid,
    fix_up_clouds,
)
from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.GFDL_1M.state import GFDL1MState
from pyMoist.GFDL_1M.locals import GFDL1MLocals


def _radiation_coupling(
    temperature: FloatField,
    pressure: FloatField,
    large_scale_cloud_fraction: FloatField,
    convective_cloud_fraction: FloatField,
    vapor: FloatField,
    large_scale_liquid: FloatField,
    large_scale_ice: FloatField,
    convective_liquid: FloatField,
    convective_ice: FloatField,
    rain: FloatField,
    snow: FloatField,
    graupel: FloatField,
    liquid_concentration: FloatField,
    ice_concentration: FloatField,
    radiation_vapor: FloatField,
    radiation_liquid: FloatField,
    radiation_ice: FloatField,
    radiation_rain: FloatField,
    radiation_snow: FloatField,
    radiation_graupel: FloatField,
    radiation_cloud_fraction: FloatField,
    liquid_radius: FloatField,
    ice_radius: FloatField,
) -> None:
    """
    Couple radiation with cloud variables to ensure physical consistency.
    """

    from __externals__ import FAC_RI, FAC_RL, MAX_RI, MAX_RL, MIN_RI, MIN_RL

    with computation(PARALLEL), interval(...):
        # water vapor
        radiation_vapor = vapor

        # total cloud fraction
        radiation_cloud_fraction = max(min(large_scale_cloud_fraction + convective_cloud_fraction, 1.0), 0.0)
        if radiation_cloud_fraction >= 1.0e-5:
            radiation_liquid = (
                (large_scale_liquid + convective_liquid) / radiation_cloud_fraction
                if (large_scale_liquid + convective_liquid) >= 1.0e-8
                else 0.0
            )
            radiation_ice = (
                (large_scale_ice + convective_ice) / radiation_cloud_fraction
                if (large_scale_ice + convective_ice) >= 1.0e-8
                else 0.0
            )
            radiation_rain = rain / radiation_cloud_fraction if rain >= 1.0e-8 else 0.0
            radiation_snow = snow / radiation_cloud_fraction if snow >= 1.0e-8 else 0.0
            radiation_graupel = graupel / radiation_cloud_fraction if graupel >= 1.0e-8 else 0.0
        else:
            radiation_cloud_fraction = 0.0
            radiation_liquid = 0.0
            radiation_ice = 0.0
            radiation_rain = 0.0
            radiation_snow = 0.0
            radiation_graupel = 0.0

        # Cap the high end of condensates
        radiation_liquid = min(radiation_liquid, 0.01)
        radiation_ice = min(radiation_ice, 0.01)
        radiation_rain = min(radiation_rain, 0.01)
        radiation_snow = min(radiation_snow, 0.01)
        radiation_graupel = min(radiation_graupel, 0.01)

        # Liquid radii - Brams formulation with limits
        liquid_radius = max(
            MIN_RL,
            min(
                cloud_effective_radius_liquid(pressure, temperature, radiation_liquid, liquid_concentration)
                * FAC_RL,
                MAX_RL,
            ),
        )
        # Ice radii - Brams formulation with limits
        ice_radius = max(
            MIN_RI,
            min(
                cloud_effective_radius_ice(pressure, temperature, radiation_ice) * FAC_RI,
                MAX_RI,
            ),
        )


def _update_humidity(
    temperature: FloatField,
    pressure: FloatField,
    vapor: FloatField,
    humidity: FloatField,
    ese: GlobalTable_saturation_tables,
    esx: GlobalTable_saturation_tables,
):
    with computation(PARALLEL), interval(...):
        qsat, _ = saturation_specific_humidity(temperature, pressure * 100, ese, esx)
        humidity = vapor * qsat


class GFDL1MRadiationCoupling:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        config: GFDL1MConfig,
        saturation_tables: SaturationVaporPressureTable,
    ):
        """
        Initialize GFDL radiation coupling class

        Arguments:
            stencil_factory (StencilFactory): Factory to create stencils.
            config (GFDL1MConfig): contains all constants for GFDL Single Moment Microphysics
        """

        # make config and saturation tables visible at runtime
        self.config = config
        self.saturation_tables = saturation_tables

        # construct stencils
        self.fix_up_clouds = stencil_factory.from_dims_halo(
            func=fix_up_clouds,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self.radiation_coupling = stencil_factory.from_dims_halo(
            func=_radiation_coupling,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "FAC_RL": config.FAC_RL,
                "MIN_RL": config.MIN_RL,
                "MAX_RL": config.MAX_RL,
                "FAC_RI": config.FAC_RI,
                "MIN_RI": config.MIN_RI,
                "MAX_RI": config.MAX_RI,
            },
        )
        self.update_humidity = stencil_factory.from_dims_halo(
            func=_update_humidity,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        if config.DO_QA:
            raise NotImplementedError(
                "[Radiation Coupling] DO_QA option implemented, but untested. "
                "Disable message manually to process."
            )

    def __call__(
        self,
        state: GFDL1MState,
        locals: GFDL1MLocals,
    ):
        """
        Perform the radiation coupling process. This prefills fields for the proper radiation scheme.
        Fields are (generally) copied cleanly from non-radiation storage to the radiation driven counterpart
        with only minor modifications. Exceptions include extreme value checks and unit conversions.
        """
        self.fix_up_clouds(
            vapor=state.mixing_ratio.vapor,
            t=state.t,
            large_scale_liquid=state.mixing_ratio.large_scale_liquid,
            large_scale_ice=state.mixing_ratio.large_scale_ice,
            large_scale_cloud_fraction=state.cloud_fraction.large_scale,
            convective_liquid=state.mixing_ratio.convective_liquid,
            convective_ice=state.mixing_ratio.convective_ice,
            convective_cloud_fraction=state.cloud_fraction.convective,
        )
        
        self.radiation_coupling(
            temperature=state.t,
            pressure=locals.p_mb,
            large_scale_cloud_fraction=state.cloud_fraction.large_scale,
            convective_cloud_fraction=state.cloud_fraction.convective,
            vapor=state.mixing_ratio.vapor,
            large_scale_liquid=state.mixing_ratio.large_scale_liquid,
            large_scale_ice=state.mixing_ratio.large_scale_ice,
            convective_liquid=state.mixing_ratio.convective_liquid,
            convective_ice=state.mixing_ratio.convective_ice,
            rain=state.mixing_ratio.rain,
            snow=state.mixing_ratio.snow,
            graupel=state.mixing_ratio.graupel,
            liquid_concentration=state.concentration.liquid,
            ice_concentration=state.concentration.ice,
            radiation_vapor=state.radiation_field.vapor,
            radiation_liquid=state.radiation_field.liquid,
            radiation_ice=state.radiation_field.ice,
            radiation_rain=state.radiation_field.rain,
            radiation_snow=state.radiation_field.snow,
            radiation_graupel=state.radiation_field.graupel,
            radiation_cloud_fraction=state.radiation_field.cloud_fraction,
            liquid_radius=state.cloud_particle_effective_radius.liquid,
            ice_radius=state.cloud_particle_effective_radius.ice,
        )

        if self.config.DO_QA:
            self.update_humidity(
                temperature=state.t,
                pressure=locals.p_mb,
                vapor=state.mixing_ratio.vapor,
                humidity=state.relative_humidity_after_pdf,
                ese=self.saturation_tables.ese,
                esx=self.saturation_tables.esx,
            )
