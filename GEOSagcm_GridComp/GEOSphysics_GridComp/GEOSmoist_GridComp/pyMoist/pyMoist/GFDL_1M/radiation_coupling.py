from ndsl import NDSLRuntime, StencilFactory, ndsl_log
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.gt4py import PARALLEL, computation, interval
from ndsl.dsl.typing import FloatField

from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.saturation_tables import (
    GlobalTable_saturation_tables,
    SaturationVaporPressureTable,
    saturation_specific_humidity,
)
from pyMoist.shared.incloud_processes import fix_up_clouds
from pyMoist.shared.radiation_coupling import radiation_coupling


def update_humidity(
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


class GFDL1MRadiationCoupling(NDSLRuntime):
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
        # init NDSLRuntime
        super().__init__(stencil_factory)

        # make config and saturation tables visible at runtime
        self.config = config
        self.saturation_tables = saturation_tables

        # construct stencils
        self._fix_up_clouds = stencil_factory.from_dims_halo(
            func=fix_up_clouds,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._radiation_coupling = stencil_factory.from_dims_halo(
            func=radiation_coupling,
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
        self._update_humidity = stencil_factory.from_dims_halo(
            func=update_humidity,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        if config.DO_QA:
            ndsl_log.log(
                "[Radiation Coupling] DO_QA option implemented, but untested. "
                "Running untested code... proceed with caution"
            )

    def __call__(
        self,
        t,
        mixing_ratio_vapor,
        mixing_ratio_large_scale_liquid,
        mixing_ratio_large_scale_ice,
        mixing_ratio_convective_liquid,
        mixing_ratio_rain,
        mixing_ratio_snow,
        mixing_ratio_graupel,
        mixing_ratio_convective_ice,
        cloud_fraction_large_scale,
        cloud_fraction_convective,
        concentration_liquid,
        concentration_ice,
        liquid_radius,
        ice_radius,
        relative_humidity_after_pdf,
        radiation_vapor,
        radiation_liquid,
        radiation_ice,
        radiation_rain,
        radiation_snow,
        radiation_graupel,
        radiation_cloud_fraction,
        local_p_mb,
    ):
        """
        Perform the radiation coupling process. This prefills fields for the proper radiation scheme.
        Fields are (generally) copied cleanly from non-radiation storage to the radiation driven counterpart
        with only minor modifications. Exceptions include extreme value checks and unit conversions.
        """
        self._fix_up_clouds(
            vapor=mixing_ratio_vapor,
            t=t,
            large_scale_liquid=mixing_ratio_large_scale_liquid,
            large_scale_ice=mixing_ratio_large_scale_ice,
            large_scale_cloud_fraction=cloud_fraction_large_scale,
            convective_liquid=mixing_ratio_convective_liquid,
            convective_ice=mixing_ratio_convective_ice,
            convective_cloud_fraction=cloud_fraction_convective,
        )

        self._radiation_coupling(
            temperature=t,
            pressure=local_p_mb,
            large_scale_cloud_fraction=cloud_fraction_large_scale,
            convective_cloud_fraction=cloud_fraction_convective,
            vapor=mixing_ratio_vapor,
            large_scale_liquid=mixing_ratio_large_scale_liquid,
            large_scale_ice=mixing_ratio_large_scale_ice,
            convective_liquid=mixing_ratio_convective_liquid,
            convective_ice=mixing_ratio_convective_ice,
            rain=mixing_ratio_rain,
            snow=mixing_ratio_snow,
            graupel=mixing_ratio_graupel,
            liquid_concentration=concentration_liquid,
            ice_concentration=concentration_ice,
            radiation_vapor=radiation_vapor,
            radiation_liquid=radiation_liquid,
            radiation_ice=radiation_ice,
            radiation_rain=radiation_rain,
            radiation_snow=radiation_snow,
            radiation_graupel=radiation_graupel,
            radiation_cloud_fraction=radiation_cloud_fraction,
            liquid_radius=liquid_radius,
            ice_radius=ice_radius,
        )

        if self.config.DO_QA:
            self._update_humidity(
                temperature=t,
                pressure=local_p_mb,
                vapor=mixing_ratio_vapor,
                humidity=relative_humidity_after_pdf,
                ese=self.saturation_tables.ese,
                esx=self.saturation_tables.esx,
            )
