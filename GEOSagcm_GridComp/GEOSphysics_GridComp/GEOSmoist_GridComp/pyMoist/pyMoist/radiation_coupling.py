from ndsl.dsl.gt4py import PARALLEL, computation, interval

from ndsl import StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float, FloatField
from pyMoist.field_types import GlobalTable_saturaion_tables
from pyMoist.saturation_tables.qsat_functions import saturation_specific_humidity
from pyMoist.shared_incloud_processes import (
    cloud_effective_radius_ice,
    cloud_effective_radius_liquid,
    fix_up_clouds,
)


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
    ese: GlobalTable_saturaion_tables,
    esx: GlobalTable_saturaion_tables,
):
    with computation(PARALLEL), interval(...):
        qsat, _ = saturation_specific_humidity(temperature, pressure * 100, ese, esx)
        humidity = vapor * qsat


class GFDL1MRadiationCoupling:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        DO_QA: bool,
        FAC_RL: Float,
        MIN_RL: Float,
        MAX_RL: Float,
        FAC_RI: Float,
        MIN_RI: Float,
        MAX_RI: Float,
    ) -> None:
        """
        Initialize GFDL radiation coupling class.

        Arguments:
            stencil_factory (StencilFactory): Factory to create stencils.
            quantity_factory (QuantityFactory): Factory to create quantities.
            DO_QA (bool): Flow control parameter. Details unknown.
            FAC_RL (Float): Factor for liquid effective radius.
            MIN_RL (Float): Minimum liquid effective radius.
            MAX_RL (Float): Maximum liquid effective radius.
            FAC_RI (Float): Factor for ice effective radius.
            MIN_RI (Float): Minimum ice effective radius.
            MAX_RI (Float): Maximum ice effective radius.
        """
        self.DO_QA = DO_QA
        self.fix_up_clouds = stencil_factory.from_dims_halo(
            func=fix_up_clouds,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )
        self.radiation_coupling = stencil_factory.from_dims_halo(
            func=_radiation_coupling,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "FAC_RL": FAC_RL,
                "MIN_RL": MIN_RL,
                "MAX_RL": MAX_RL,
                "FAC_RI": FAC_RI,
                "MIN_RI": MIN_RI,
                "MAX_RI": MAX_RI,
            },
        )
        self.update_humidity = stencil_factory.from_dims_halo(
            func=_update_humidity,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        if self.DO_QA:
            raise NotImplementedError(
                "[Radiation Coupling] DO_QA option implemented, but untested. "
                "Disable message manually to process."
            )

    def __call__(
        self,
        vapor: FloatField,
        temperature: FloatField,
        large_scale_liquid: FloatField,
        large_scale_ice: FloatField,
        large_scale_cloud_fraction: FloatField,
        convective_liquid: FloatField,
        convective_ice: FloatField,
        convective_cloud_fraction: FloatField,
        pressure: FloatField,
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
        relative_humidity_after_pdf: FloatField,
        ese: GlobalTable_saturaion_tables,
        esx: GlobalTable_saturaion_tables,
    ):
        """
        Perform the radiation coupling process.

        Arguments:
            vapor (inout): water vapor mixing ratio (kg/kg)
            temperature (inout): temperature (Kelvin)
            large_scale_liquid (inout): large scale liquid mixing ratio (kg/kg)
            large_scale_ice (inout): large scale ice mixing ratio (kg/kg)
            large_scale_cloud_fraction (inout): large scale cloud fraction (unitless)
            convective_liquid (inout): convective liquid mixing ratio (kg/kg)
            convective_ice (inout): convective ice mixing ratio (kg/kg)
            convective_cloud_fraction (inout): convective cloud fraction
            pressure (in): pressure (millibars)
            rain (inout): rain mixing ratio (unitless)
            snow (inout): snow mixing ratio (unitless)
            graupel (inout): graupel mixing ratio (unitless)
            liquid_concentration (in): liquid cloud droplet concentration (m^-3)
            ice_concentration (in): ice cloud droplet concnetration (m^-3)
            radiation_vapor (inout): radiation water vapor mixing ratio (kg/kg)
            radiation_liquid (inout): radiation liquid cloud mixing ratio (kg/kg)
            radiation_ice (inout): radiation ice cloud mixing ratio (kg/kg)
            radiation_rain (inout): radiation rain mixing ratio (unitless)
            radiation_snow (inout): radiation snow mixing ratio (unitless)
            radiation_graupel (inout): radiation graupel mixing ratio (unitless)
            radiation_cloud_fraction (inout): radiation cloud fraction (unitless)
            liquid_radius (inout): radiation liquid effective radius (m)
            ice_radius (inout): radiation ice effective radius (m)
        """
        self.fix_up_clouds(
            vapor,
            temperature,
            large_scale_liquid,
            large_scale_ice,
            large_scale_cloud_fraction,
            convective_liquid,
            convective_ice,
            convective_cloud_fraction,
        )
        self.radiation_coupling(
            temperature,
            pressure,
            large_scale_cloud_fraction,
            convective_cloud_fraction,
            vapor,
            large_scale_liquid,
            large_scale_ice,
            convective_liquid,
            convective_ice,
            rain,
            snow,
            graupel,
            liquid_concentration,
            ice_concentration,
            radiation_vapor,
            radiation_liquid,
            radiation_ice,
            radiation_rain,
            radiation_snow,
            radiation_graupel,
            radiation_cloud_fraction,
            liquid_radius,
            ice_radius,
        )

        if self.DO_QA:
            self.update_humidity(
                temperature,
                pressure,
                vapor,
                relative_humidity_after_pdf,
                ese,
                esx,
            )
