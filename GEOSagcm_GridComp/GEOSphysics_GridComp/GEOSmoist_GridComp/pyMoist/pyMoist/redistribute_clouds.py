from ndsl import StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.gt4py import PARALLEL, computation, interval
from ndsl.dsl.typing import FloatField
from pyMoist.constants import ALHLBCP, ALHSBCP


def redistribute_clouds(
    cloud_fraction: FloatField,
    convective_cloud_fraction: FloatField,
    large_scale_cloud_fraction: FloatField,
    liquid: FloatField,
    convective_liquid: FloatField,
    large_scale_liquid: FloatField,
    ice: FloatField,
    convective_ice: FloatField,
    large_scale_ice: FloatField,
    vapor: FloatField,
    temperature: FloatField,
):
    with computation(PARALLEL), interval(...):

        # Fix cloud quants if too small
        if liquid + ice < 1e-8:
            vapor = vapor + liquid + ice
            temperature = temperature - ALHLBCP * liquid - ALHSBCP * ice
            cloud_fraction = 0.0
            liquid = 0.0
            ice = 0.0

        if cloud_fraction < 1e-5:
            vapor = vapor + liquid + ice
            temperature = temperature - (ALHLBCP * liquid) - (ALHSBCP * ice)
            cloud_fraction = 0.0
            liquid = 0.0
            ice = 0.0

        # Redistribute liquid CN/LS portions based on prior fractions
        local_cloud_fraction = 0.0
        if convective_liquid + large_scale_liquid > 0.0:
            local_cloud_fraction = min(
                max(convective_liquid / (convective_liquid + large_scale_liquid), 0.0), 1.0
            )

        # Put all new condensate into LS
        dqc = liquid - (convective_liquid + large_scale_liquid)
        if dqc > 0.0:
            large_scale_liquid = large_scale_liquid + dqc
            dqc = 0.0

        # Any loss of condensate uses the FCN ratio
        convective_liquid = convective_liquid + dqc * local_cloud_fraction
        large_scale_liquid = large_scale_liquid + dqc * (1.0 - local_cloud_fraction)

        # Redistribute ice CN/LS portions based on prior fractions
        local_cloud_fraction = 0.0
        if convective_ice + large_scale_ice > 0.0:
            local_cloud_fraction = min(max(convective_ice / (convective_ice + large_scale_ice), 0.0), 1.0)

        # Put all new condensate into LS
        dqc = ice - (convective_ice + large_scale_ice)
        if dqc > 0.0:
            large_scale_ice = large_scale_ice + dqc
            dqc = 0.0

        # Any loss of condensate uses the FCN ratio
        convective_ice = convective_ice + dqc * local_cloud_fraction
        large_scale_ice = large_scale_ice + dqc * (1.0 - local_cloud_fraction)

        # Redistribute cloud-fraction CN/LS portions based on prior fractions
        local_cloud_fraction = 0.0
        if convective_cloud_fraction + large_scale_cloud_fraction > 0.0:
            local_cloud_fraction = min(
                max(
                    convective_cloud_fraction / (convective_cloud_fraction + large_scale_cloud_fraction), 0.0
                ),
                1.0,
            )

        # Put all new condensate into LS
        dqc = cloud_fraction - (convective_cloud_fraction + large_scale_cloud_fraction)
        if dqc > 0.0:
            large_scale_cloud_fraction = large_scale_cloud_fraction + dqc
            dqc = 0.0

        # Any loss of condensate uses the FCN ratio
        convective_cloud_fraction = convective_cloud_fraction + dqc * local_cloud_fraction
        large_scale_cloud_fraction = large_scale_cloud_fraction + dqc * (1.0 - local_cloud_fraction)


class RedistributeClouds:
    def __init__(
        self,
        stencil_factory: StencilFactory,
    ) -> None:
        self._redistribute_clouds = stencil_factory.from_dims_halo(
            func=redistribute_clouds,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        cloud_fraction: FloatField,
        convective_cloud_fraction: FloatField,
        large_scale_cloud_fraction: FloatField,
        liquid: FloatField,
        convective_liquid: FloatField,
        large_scale_liquid: FloatField,
        ice: FloatField,
        convective_ice: FloatField,
        large_scale_ice: FloatField,
        vapor: FloatField,
        temperature: FloatField,
    ):
        """
        Redistribute non-physical in-cloud quantities.

        Parameters:
        cloud_fraction (3D inout): total cloud fraction (unitless)
        convective_cloud_fraction (3D inout): convective cloud fraction (unitless)
        large_scale_cloud_fraction (3D inout): large scale cloud fraction (unitless)
        liquid (3D inout): in-cloud liquid mixing ratio (kg/kg)
        convective_liquid (3D inout): convective cloud liquid mixing ratio (unitless)
        large_scale_liquid (3D inout): large scale cloud liquid mixing ratio (unitless)
        ice (3D inout): in-cloud liquid mixing ratio (kg/kg)
        convective_ice (3D inout): convective cloud ice mixing ratio (unitless)
        large_scale_ice (3D inout): large scale cloud liquid mixing ratio (unitless)
        vapor (3D inout): water vapor mixing ratio (kg/kg)
        temperature (3D inout): temperature (Kelvin)
        """
        self._redistribute_clouds(
            cloud_fraction=cloud_fraction,
            convective_cloud_fraction=convective_cloud_fraction,
            large_scale_cloud_fraction=large_scale_cloud_fraction,
            liquid=liquid,
            convective_liquid=convective_liquid,
            large_scale_liquid=large_scale_liquid,
            ice=ice,
            convective_ice=convective_ice,
            large_scale_ice=large_scale_ice,
            vapor=vapor,
            temperature=temperature,
        )
