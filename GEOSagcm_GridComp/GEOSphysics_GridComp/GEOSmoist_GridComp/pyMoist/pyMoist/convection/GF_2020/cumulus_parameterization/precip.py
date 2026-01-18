from ndsl import StencilFactory, QuantityFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import (
    GF2020CumulusParameterizationConfig,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.state import (
    GF2020CumulusParameterizationState,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.locals import (
    GF2020CumulusParameterizationLocals,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)
from ndsl.logging import ndsl_log
from ndsl.dsl.gt4py import (
    PARALLEL,
    computation,
    interval,
    FORWARD,
    function,
    BACKWARD,
    K,
    sqrt,
)
from ndsl.dsl.typing import FloatField, FloatFieldIJ, IntField, Int, IntFieldIJ
import pyMoist.constants as constants
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import (
    FloatField_Plume,
    IntFieldIJ_Plume,
    FloatFieldIJ_Plume,
    Float,
)
from pyMoist.shared_incloud_processes import ice_fraction


@function
def liquid_fraction(
    t,
    convection_fraction,
    surface_type,
    FRAC_MODIS,
):
    """
    Get the fraction of liquid condensates

    Args:
        t (in): temperature
        convection_fraction (in)
        surface_type (in)
        FRAC_MODIS (in): use fraction liq/ice content derived from MODIS/CALIPO sensors
    """
    if FRAC_MODIS == 1:
        liquid_fraction = 1.0 - ice_fraction(t, convection_fraction, surface_type)
    else:
        liquid_fraction = min(
            1.0,
            (
                max(0.0, (t - cumulus_parameterization_constants.T_ICE))
                / (cumulus_parameterization_constants.T_0 - cumulus_parameterization_constants.T_ICE)
            )
            ** 2,
        )

    return liquid_fraction


def partition_liquid_ice(
    t: FloatField,
    p: FloatField_Plume,
    geopotential_height: FloatField,
    topography_height_no_negative: FloatFieldIJ,
    surface_type: FloatFieldIJ,
    convection_fraction: FloatFieldIJ,
    error_code: IntFieldIJ_Plume,
    melting_layer: FloatField,
    part_liquid_ice: FloatField,
    plume: Int,
):
    """
    Partition total condensate into liquid and ice phases

    Args:

    """
    from __externals__ import MELT_GLAC, FRAC_MODIS, k_end

    with computation(PARALLEL), interval(...):
        # constants, set internally because they may differ from global constants
        # and need to only exist inside this stencil
        t1 = 276.16
        z_meltlayer1 = 4000.0
        z_meltlayer2 = 6000.0
        del_t = 3.0

        # prefill some fields
        part_liquid_ice = 1.0
        melting_layer = 0.0

    with computation(PARALLEL), interval(0, -1):
        if MELT_GLAC == True and plume == 2:
            if error_code[0, 0][plume] == 0:
                # get function of T for partition of total condensate into liq and ice phases
                part_liquid_ice = liquid_fraction(t, convection_fraction, surface_type, FRAC_MODIS)

    with computation(PARALLEL), interval(0, -1):
        if MELT_GLAC == True and plume == 2:
            if error_code[0, 0][plume] == 0:
                # define the melting layer (the layer will be between T_0+1 < TEMP < T_1
                if t <= (cumulus_parameterization_constants.T_0 - del_t):
                    melting_layer = 0.0

                elif t < (cumulus_parameterization_constants.T_0 + del_t) and t > (
                    cumulus_parameterization_constants.T_0 - del_t
                ):
                    melting_layer = (
                        (t - (cumulus_parameterization_constants.T_0 - del_t)) / (2.0 * del_t)
                    ) ** 2

                else:
                    melting_layer = 1.0

                melting_layer = melting_layer * (1.0 - melting_layer)

    with computation(FORWARD), interval(0, 1):
        if MELT_GLAC == True and plume == 2:
            # normalize vertical integral of melting_layer to 1
            norm: FloatFieldIJ = 0.0

    with computation(FORWARD), interval(0, -2):
        if MELT_GLAC == True and plume == 2:
            if error_code[0, 0][plume] == 0:
                # normalize vertical integral of melting_layer to 1
                dp = 100.0 * (p[0, 0, 0][plume] - p[0, 0, 1][plume])
                norm = norm + melting_layer * dp / constants.MAPL_GRAV

    with computation(PARALLEL), interval(...):
        if MELT_GLAC == True and plume == 2:
            if error_code[0, 0][plume] == 0:
                # normalize vertical integral of melting_layer to 1
                melting_layer = (
                    melting_layer
                    / (norm + 1.0e-6)
                    * (
                        100
                        * (p.at(K=0, ddim=[plume]) - p.at(K=k_end - 1, ddim=[plume]))
                        / constants.MAPL_GRAV
                    )
                )


def get_precip_fluxes(
    error_code: IntFieldIJ_Plume,
    cloud_top_level: IntFieldIJ_Plume,
    cloud_base_mass_flux_modified: FloatFieldIJ_Plume,
    epsilon_forced: FloatFieldIJ_Plume,
    condensate_to_fall_forced: FloatField_Plume,
    evaporate_in_downdraft_forced: FloatField_Plume,
    precipitation_flux: FloatField,
    evaporation_flux: FloatField,
    plume: Int,
):
    with computation(PARALLEL), interval(...):
        precipitation_flux = 0.0
        evaporation_flux = 0.0

    with computation(BACKWARD), interval(...):
        if error_code[0, 0][plume] == 0:
            if K <= cloud_top_level[0, 0][plume]:
                precipitation_flux = precipitation_flux[0, 0, 1] + cloud_base_mass_flux_modified[0, 0][
                    plume
                ] * (
                    condensate_to_fall_forced[0, 0, 0][plume]
                    + epsilon_forced[0, 0][plume] * evaporate_in_downdraft_forced[0, 0, 0][plume]
                )
                precipitation_flux = max(0.0, precipitation_flux)

                evaporation_flux = (
                    evaporation_flux[0, 0, 1]
                    - cloud_base_mass_flux_modified[0, 0][plume]
                    * epsilon_forced[0, 0][plume]
                    * evaporate_in_downdraft_forced[0, 0, 0][plume]
                )
                evaporation_flux = max(0.0, evaporation_flux)


def output_evaporation_flux(
    error_code: IntFieldIJ_Plume,
    plume: Int,
    cloud_top_level: IntFieldIJ_Plume,
    p_cloud_levels_forced: FloatField_Plume,
    evap_flux: FloatField,
    evap_subl_tendency: FloatField,
):
    with computation(FORWARD), interval(...):
        if error_code[0, 0][plume] == 0:
            if K <= cloud_top_level[0, 0][plume]:
                dp: FloatFieldIJ = 100.0 * (
                    p_cloud_levels_forced[0, 0, 0][plume] - p_cloud_levels_forced[0, 0, 1][plume]
                )
                evap_subl_tendency = evap_subl_tendency + evap_flux * constants.MAPL_GRAV / dp


def output_deep_precipitation(
    error_code: IntFieldIJ_Plume,
    plume: Int,
    cloud_top_level: IntFieldIJ_Plume,
    precipitation_flux: FloatField,
    convective_precip_flux: FloatField,
):
    with computation(PARALLEL), interval(...):
        if plume == cumulus_parameterization_constants.deep:
            if error_code[0, 0][plume] == 0:
                if K <= cloud_top_level[0, 0][plume] + 1:
                    convective_precip_flux = precipitation_flux


def rain_evaporation_below_cloud_base(
    error_code: IntFieldIJ_Plume,
    updraft_lfc_level: IntFieldIJ_Plume,
    cloud_top_level: IntFieldIJ_Plume,
    ocean_fraction: FloatFieldIJ,
    p_cloud_levels_forced: FloatField_Plume,
    p_surface: FloatFieldIJ,
    t_cloud_levels: FloatField,
    vapor_cloud_levels_forced: FloatField,
    environment_saturation_mixing_ratio_cloud_levels: FloatField,
    epsilon_forced: FloatFieldIJ_Plume,
    cloud_base_mass_flux_modified: FloatFieldIJ_Plume,
    condensate_to_fall_forced: FloatField_Plume,
    evaporate_in_downdraft_forced: FloatField_Plume,
    precip: FloatFieldIJ_Plume,
    precipitation_flux: FloatField,
    evaporation_flux: FloatField,
    evaporation_below_cloud_base: FloatField,
    dtdt: FloatField_Plume,
    dvapordt: FloatField_Plume,
    dbuoyancydt: FloatField_Plume,
    plume: Int,
):
    with computation(FORWARD), interval(0, 1):
        # setup internal constants
        alpha1: FloatFieldIJ = 5.44e-4
        alpha2: FloatFieldIJ = 5.09e-3
        alpha3: FloatFieldIJ = 0.5777
        c_conv: FloatFieldIJ = 0.05

    with computation(FORWARD), interval(0, 1):
        if plume == cumulus_parameterization_constants.shallow:
            critical_rh_ocean: FloatFieldIJ = 1.0
            critical_rh_land: FloatFieldIJ = 1.0
            eff_c_conv: FloatFieldIJ = min(0.2, max(cloud_base_mass_flux_modified[0, 0][plume], c_conv))
        else:
            critical_rh_ocean: FloatFieldIJ = 0.95
            critical_rh_land: FloatFieldIJ = 0.85
            eff_c_conv: FloatFieldIJ = c_conv

        total_evaporation_below_cloud_base: FloatFieldIJ = 0.0

    with computation(PARALLEL), interval(...):
        precipitation_flux = 0.0
        evaporation_flux = 0.0

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            critical_rh: FloatFieldIJ = critical_rh_ocean * ocean_fraction + critical_rh_land * (
                1.0 - ocean_fraction
            )

    with computation(BACKWARD), interval(...):
        if error_code[0, 0][plume] == 0:
            if K <= cloud_top_level[0, 0][plume]:
                dp: FloatFieldIJ = 100.0 * (
                    p_cloud_levels_forced[0, 0, 0][plume] - p_cloud_levels_forced[0, 0, 1][plume]
                )

                if K <= updraft_lfc_level[0, 0][plume]:
                    vapor_deficit: FloatFieldIJ = max(
                        0.0,
                        (
                            critical_rh * environment_saturation_mixing_ratio_cloud_levels
                            - vapor_cloud_levels_forced
                        ),
                    )

                    evaporation_below_cloud_base = (
                        eff_c_conv
                        * alpha1
                        * vapor_deficit
                        * (
                            sqrt(p_cloud_levels_forced[0, 0, 0][plume] / p_surface)
                            / alpha2
                            * precipitation_flux[0, 0, 1]
                            / eff_c_conv
                        )
                        ** alpha3
                    )

                    evaporation_below_cloud_base = evaporation_below_cloud_base * dp / constants.MAPL_GRAV

                else:

                    evaporation_below_cloud_base = 0.0

                precipitation_flux = (
                    precipitation_flux[0, 0, 1]
                    - evaporation_below_cloud_base
                    + cloud_base_mass_flux_modified[0, 0][plume]
                    * (
                        condensate_to_fall_forced[0, 0, 0][plume]
                        + epsilon_forced[0, 0][plume] * evaporate_in_downdraft_forced[0, 0, 0][plume]
                    )
                )
                precipitation_flux = max(0.0, precipitation_flux)

                evaporation_flux = (
                    evaporation_flux[0, 0, 1]
                    + evaporation_below_cloud_base
                    - cloud_base_mass_flux_modified[0, 0][plume]
                    * epsilon_forced[0, 0][plume]
                    * evaporate_in_downdraft_forced[0, 0, 0][plume]
                )
                evaporation_flux = max(0.0, evaporation_flux)

                total_evaporation_below_cloud_base = (
                    total_evaporation_below_cloud_base + evaporation_below_cloud_base
                )

                del_vapor = evaporation_below_cloud_base * constants.MAPL_GRAV / dp
                del_t = (
                    -evaporation_below_cloud_base
                    * constants.MAPL_GRAV
                    / dp
                    * (cumulus_parameterization_constants.XLV / cumulus_parameterization_constants.CP)
                )

                dvapordt[0, 0, 0][plume] = dvapordt[0, 0, 0][plume] + del_vapor
                dtdt[0, 0, 0][plume] = dtdt[0, 0, 0][plume] + del_t
                dbuoyancydt[0, 0, 0][plume] = (
                    dbuoyancydt[0, 0, 0][plume]
                    + cumulus_parameterization_constants.CP * del_t
                    + cumulus_parameterization_constants.XLV * del_vapor
                )

                precip[0, 0][plume] = precip[0, 0][plume] - evaporation_below_cloud_base


class PrecipFactor:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GF2020Config,
        cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
    ):
        pass
        # # make configuration visible at runtime
        # self.config = config
        # self.cumulus_parameterization_config = cumulus_parameterization_config

        # # construct stencils and functions
        # self._ = stencil_factory.from_dims_halo(
        #     func=,
        #     compute_dims=[X_DIM, Y_DIM, Z_DIM],
        #     externals={
        #         "MELT_GLAC": cumulus_parameterization_config.MELT_GLAC,
        #         "FRAC_MODIS": cumulus_parameterization_config.FRAC_MODIS,
        #     },
        # )

    def __call__(self, *args, **kwds):
        pass


class RainEvapBelowCloudBase:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GF2020Config,
        cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
    ):
        # make configuration visible at runtime
        self.config = config
        self.cumulus_parameterization_config = cumulus_parameterization_config

        # construct stencils and functions
        self._rain_evap_below_cloudbase = stencil_factory.from_dims_halo(
            func=rain_evaporation_below_cloud_base,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._rain_evap_below_cloudbase(
            error_code=state.output.error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
            epsilon_forced=state.output.epsilon_forced,
            updraft_lfc_level=state.output.updraft_lfc_level,
            cloud_top_level=state.output.cloud_top_level,
            p_cloud_levels_forced=state.output.p_cloud_levels_forced,
            p_surface=state.input_output.p_surface,
            evaporate_in_downdraft_forced=state.output.evaporate_in_downdraft_forced,
            condensate_to_fall_forced=state.output.condensate_to_fall_forced,
            local_env_saturation_mixing_ratio_cloud_levels=locals.environment_saturation_mixing_ratio_cloud_levels,
            local_vapor_cloud_levels_forced=locals.vapor_cloud_levels_forced,
            ocean_fraction=locals.ocean_fraction,
            cloud_base_mass_flux=state.output.cloud_base_mass_flux,
            local_evap_bcb=locals.evap_bcb,
            evap_flux=locals.evap_flux,
            dbuoyancydt=state.output.dbuoyancydt,
            dvapordt=state.output.dvapordt,
            t=state.output.dtdt,
            precip=state.output.precip,
            prec_flux=locals.prec_flux,
        )


class CloudDissipation:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass


class OutputEvaporationFlux:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GF2020Config,
        cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
    ):
        # make configuration visible at runtime
        self.config = config
        self.cumulus_parameterization_config = cumulus_parameterization_config

        # construct stencils and functions
        self._output_evaporation_flux = stencil_factory.from_dims_halo(
            func=output_evaporation_flux,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._output_evaporation_flux(
            error_code=state.output.error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
            cloud_top_level=state.output.cloud_top_level,
            p_cloud_levels_forced=state.output.p_cloud_levels_forced,
            evap_flux=locals.evap_flux,
            evap_subl_tendency=state.output.evap_subl_tendency,
        )


class LightningFlassDensity:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass


class OutputDeepPrecipitation:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GF2020Config,
        cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
    ):
        # make configuration visible at runtime
        self.config = config
        self.cumulus_parameterization_config = cumulus_parameterization_config

        # construct stencils and functions
        self._output_deep_precipitation = stencil_factory.from_dims_halo(
            func=output_deep_precipitation,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._output_deep_precipitation(
            error_code=state.output.error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
            cloud_top_level=state.output.cloud_top_level,
            precipitation_flux=locals.precipitation_flux,
            convective_precip_flux=state.output.convective_precip_flux,
        )


class UpdateWorkfunctionsAndCondensates:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass
