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
from pyMoist.convection.GF_2020.cumulus_parameterization.shared_stencils import updraft_vertical_velocity


from ndsl.dsl.typing import FloatField, FloatFieldIJ, Float, IntFieldIJ, Int, IntField
from ndsl.dsl.gt4py import computation, PARALLEL, interval, FORWARD, K
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
import pyMoist.constants as constants
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import (
    IntFieldIJ_Plume,
    FloatFieldIJ_Plume,
    FloatField_Plume,
    FloatFieldIJ_Ensemble,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.shared_functions import (
    get_cloud_boundary_conditions,
    liquid_fraction,
)


def updraft_moisture_light(
    start_level: IntFieldIJ,
    error_code: IntFieldIJ_Plume,
    geopotential_height_cloud_levels_forced: FloatField,
    cloud_vapor_mixing_ratio_forced: FloatField,
    cloud_liquid_after_rain_forced: FloatField_Plume,
    precipitable_water_updraft_forced: FloatField_Plume,
    total_normalized_integrated_condensate_forced: FloatFieldIJ_Plume,
    cloud_moist_static_energy_forced: FloatField,
    unspecifid_temperature: FloatField,
    ocean_fraction: FloatFieldIJ,
    convection_fraction: FloatFieldIJ,
    surface_type: FloatFieldIJ,
    p_forced: FloatField,
    cloud_top_level: IntFieldIJ_Plume,
    buoyancy_forced: FloatField,
    cloud_liquid_before_rain_forced: FloatField,
    t_cloud_levels: FloatField,
    vapor_forced: FloatField,
    gamma_cloud_levels_forced: FloatField,
    normalized_massflux_updraft_forced: FloatField_Plume,
    environment_saturation_mixing_ratio_cloud_levels_forced: FloatField,
    updraft_origin_level: FloatFieldIJ_Plume,
    vapor_cloud_levels_forced: FloatField,
    vapor_excess: FloatFieldIJ,
    mass_entrainment_updraft: FloatField,
    mass_detrainment_updraft: FloatField,
    psum: FloatFieldIJ,
    psumh: FloatFieldIJ,
    c1d: FloatField,
    add_buoyancy: FloatFieldIJ,
    AVERAGE_LAYER_DEPTH: Float,
    C0: Float,
    plume: Int,
):
    from __externals__ import (
        BOUNDARY_CONDITION_METHOD,
        k_end,
        MELT_GLAC,
        FRAC_MODIS,
        QRC_CRIT_OCN,
        QRC_CRIT_LND,
    )

    with computation(FORWARD), interval(0, 1):
        total_normalized_integrated_condensate_forced[0, 0][plume] = 0.0
        psum = 0.0
        psumh = 0.0

    with computation(PARALLEL), interval(...):
        # make garbage field so the get_cloud_boundary_conditions call does not break (this is never touched)
        dummy_field_no_read = 0.0 + BOUNDARY_CONDITION_METHOD

    with computation(PARALLEL), interval(...):
        precipitable_water_updraft_forced[0, 0, 0][plume] = 0.0
        cloud_liquid_after_rain_forced[0, 0, 0][plume] = 0.0
        cloud_liquid_before_rain_forced = 0.0
        unspecifid_temperature = t_cloud_levels
        cloud_vapor_mixing_ratio_forced = vapor_cloud_levels_forced

    # with computation(FORWARD), interval(0, 1):
    #     if error_code[0, 0][plume] == 0:
    #         # get boundary condition for cloud_vapor_mixing_ratio_forced
    #         vapor_boundary_condition = get_cloud_boundary_conditions(
    #             field=vapor_cloud_levels_forced,
    #             scalar_perturbation=0,
    #             p=p_forced,
    #             updraft_origin_level=updraft_origin_level[0, 0][plume],
    #             ocean_fraction=ocean_fraction,
    #             BOUNDARY_CONDITION_METHOD=BOUNDARY_CONDITION_METHOD,
    #             AVERAGE_LAYER_DEPTH=AVERAGE_LAYER_DEPTH,
    #             k_end=k_end,
    #             compute_perturbation=False,
    #             perturbation_field=dummy_field_no_read,
    #         )
    #         cloud_vapor_mixing_ratio_forced = (
    #             vapor_boundary_condition
    #             + vapor_excess
    #             + 0.5 * add_buoyancy / cumulus_parameterization_constants.XLV
    #         )

    # with computation(FORWARD), interval(...):
    #     if error_code[0, 0][plume] == 0:
    #         if K >= start_level + 1 and K <= cloud_top_level[0, 0][plume] + 1:
    #             dz = (
    #                 geopotential_height_cloud_levels_forced
    #                 - geopotential_height_cloud_levels_forced[0, 0, -1]
    #             )

    #             # saturation in cloud, this is what is allowed to be in it
    #             saturation_cloud_liquid_rain_forced = (
    #                 environment_saturation_mixing_ratio_cloud_levels_forced
    #                 + (1 / cumulus_parameterization_constants.XLV)
    #                 * (gamma_cloud_levels_forced / (1.0 + gamma_cloud_levels_forced))
    #                 * buoyancy_forced
    #             )

    #             # 1. steady state plume equation, for what could be in cloud without condensation
    #             denom = (
    #                 normalized_massflux_updraft_forced[0, 0, -1][plume]
    #                 - 0.5 * mass_detrainment_updraft[0, 0, -1]
    #                 + mass_entrainment_updraft[0, 0, -1]
    #             )
    #             if denom > 0.0:
    #                 cloud_vapor_mixing_ratio_forced = (
    #                     cloud_vapor_mixing_ratio_forced[0, 0, -1]
    #                     * normalized_massflux_updraft_forced[0, 0, -1][plume]
    #                     - 0.5 * mass_detrainment_updraft[0, 0, -1] * cloud_vapor_mixing_ratio_forced[0, 0, -1]
    #                     + mass_entrainment_updraft[0, 0, -1] * vapor_forced[0, 0, -1]
    #                 ) / denom
    #                 if K == start_level + 1:
    #                     cloud_vapor_mixing_ratio_forced = (
    #                         cloud_vapor_mixing_ratio_forced
    #                         + vapor_excess * mass_entrainment_updraft[0, 0, -1] / denom
    #                     )
    #             else:
    #                 cloud_vapor_mixing_ratio_forced = cloud_vapor_mixing_ratio_forced[0, 0, -1]

    #             # total condensed water before rainout
    #             cloud_liquid_before_rain_forced = max(
    #                 0.0, cloud_vapor_mixing_ratio_forced - saturation_cloud_liquid_rain_forced
    #             )
    #             # updraft temp
    #             unspecifid_temperature = (1.0 / cumulus_parameterization_constants.CP) * (
    #                 cloud_moist_static_energy_forced
    #                 - constants.MAPL_GRAV * geopotential_height_cloud_levels_forced
    #                 - cumulus_parameterization_constants.XLV * saturation_cloud_liquid_rain_forced
    #             )

    #             # add glaciation effect on the MSE
    #             if MELT_GLAC == True:
    #                 delta_cloud_mse_glac = (
    #                     cloud_liquid_before_rain_forced
    #                     * (
    #                         1.0
    #                         - liquid_fraction(
    #                             unspecifid_temperature, convection_fraction, surface_type, FRAC_MODIS
    #                         )
    #                     )
    #                     * cumulus_parameterization_constants.XLF
    #                 )

    #                 unspecifid_temperature = (
    #                     unspecifid_temperature
    #                     + (1.0 / cumulus_parameterization_constants.CP) * delta_cloud_mse_glac
    #                 )

    #             if C0 < 1.0e-6:
    #                 cx0 = 0.0
    #             else:
    #                 cx0 = (c1d + C0) * dz

    #             cloud_liquid_after_rain_forced[0, 0, 0][plume] = cloud_liquid_before_rain_forced / (1.0 + cx0)
    #             min_liq = ocean_fraction * QRC_CRIT_OCN + (1.0 - ocean_fraction) * QRC_CRIT_LND
    #             precipitable_water_updraft_forced[0, 0, 0][plume] = cx0 * max(
    #                 0.0, cloud_liquid_after_rain_forced[0, 0, 0][plume] - min_liq
    #             )  # units kg[rain]/kg[air]

    #             # convert pw to normalized pw
    #             precipitable_water_updraft_forced[0, 0, 0][plume] = (
    #                 precipitable_water_updraft_forced[0, 0, 0][plume]
    #                 * normalized_massflux_updraft_forced[0, 0, 0][plume]
    #             )

    #             # total water (vapor + condensed) in updraft after the rainout
    #             cloud_vapor_mixing_ratio_forced = cloud_liquid_after_rain_forced[0, 0, 0][plume] + min(
    #                 cloud_vapor_mixing_ratio_forced, saturation_cloud_liquid_rain_forced
    #             )

    # with computation(PARALLEL), interval(...):
    #     if error_code[0, 0][plume] == 0:
    #         # get back water vapor qc
    #         if K <= cloud_top_level[0, 0][plume]:
    #             cloud_vapor_mixing_ratio_forced = (
    #                 cloud_vapor_mixing_ratio_forced - cloud_liquid_after_rain_forced[0, 0, 0][plume]
    #             )


def in_cloud_updraft_air_temperature(
    error_code: IntFieldIJ_Plume,
    plume: Int,
    local_incloud_air_temp_forced: FloatField,
    local_hcdo: FloatField,
    local_geopotential_height_cloud_levels_forced: FloatField,
    local_incloud_water_vapor_mixing_ratio_forced: FloatField,
    local_t_cloud_levels_forced: FloatField,
):
    with computation(PARALLEL), interval(0, -1):
        if error_code[0, 0][plume] == 0:
            local_incloud_air_temp_forced = (1.0 / cumulus_parameterization_constants.CP) * (
                local_hcdo
                - constants.MAPL_GRAV * local_geopotential_height_cloud_levels_forced
                - cumulus_parameterization_constants.XLV * local_incloud_water_vapor_mixing_ratio_forced
            )
    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] != 0:
            local_incloud_air_temp_forced = local_t_cloud_levels_forced


def get_melting_profile(
    error_code: IntFieldIJ_Plume,
    plume: Int,
    local_melting_layer: FloatField,
    local_partition_liquid_ice: FloatField,
    p_cloud_levels_forced: FloatField_Plume,
    precipitable_water_updraft_forced: FloatField_Plume,
    local_melting: FloatField,
):
    from __externals__ import k_end, MELT_GLAC

    with computation(FORWARD), interval(...):
        ktf = k_end - 1
        if MELT_GLAC == True and plume == 2:
            pwo_solid_phase = 0.0
            pwo_eff = 0.0
            local_melting = 0.0

            if error_code[0, 0][plume] > 0:
                local_melting = 0.0

            total_pwo_solid_phase: FloatFieldIJ = 0.0

    with computation(FORWARD), interval(...):
        if MELT_GLAC == True and plume == 2:
            if K <= ktf - 1:
                if error_code[0, 0][plume] == 0:
                    dp = 100.0 * (
                        p_cloud_levels_forced[0, 0, 0][plume] - p_cloud_levels_forced[0, 0, 1][plume]
                    )

                    pwo_eff = 0.5 * (
                        precipitable_water_updraft_forced[0, 0, 0][plume]
                        + precipitable_water_updraft_forced[0, 0, 1][plume]
                    )

                    pwo_solid_phase = (1.0 - local_partition_liquid_ice) * pwo_eff

                    total_pwo_solid_phase = total_pwo_solid_phase + pwo_solid_phase * dp / constants.MAPL_GRAV

    with computation(PARALLEL), interval(...):
        if MELT_GLAC == True and plume == 2:
            if K <= ktf:
                if error_code[0, 0][plume] == 0:
                    local_melting = local_melting_layer * (
                        total_pwo_solid_phase
                        / (
                            100
                            * (
                                p_cloud_levels_forced.at(K=0, ddim=[plume])
                                - p_cloud_levels_forced.at(K=ktf, ddim=[plume])
                            )
                            / constants.MAPL_GRAV
                        )
                    )
        else:
            local_melting = 0.0


class C1DProfile:
    """
    UNFINISHED implementation of the C1D Profile generator. This code is manually diables with a
    "do not enable" note in fortran, so it is not going to be completed until specifically requested.
    """

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
        self._updraft_moisture_light = stencil_factory.from_dims_halo(
            func=updraft_moisture_light,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "BOUNDARY_CONDITION_METHOD": cumulus_parameterization_config.BOUNDARY_CONDITION_METHOD,
                "MELT_GLAC": cumulus_parameterization_config.MELT_GLAC,
                "FRAC_MODIS": cumulus_parameterization_config.FRAC_MODIS,
                "QRC_CRIT_OCN": cumulus_parameterization_config.QRC_CRIT_OCN,
                "QRC_CRIT_LND": cumulus_parameterization_config.QRC_CRIT_LND,
            },
        )

        self._updraft_vertical_velocity = stencil_factory.from_dims_halo(
            func=updraft_vertical_velocity,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ZERO_DIFF": cumulus_parameterization_config.ZERO_DIFF},
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        # need to wrap this into something that is neither a constant (it can be
        # modified inside the gridcomp) nor a local
        USE_C1D = False
        if plume_dependent_constants.PLUME_INDEX == 2 and USE_C1D:
            # NOTE this needs to go in a stencil once NDSL base functions are merged
            locals.c1d.field[:] = self.config.C1

        if self.cumulus_parameterization_config.FIRST_GUESS_W or self.config.AUTOCONV == 4:
            self._updraft_moisture_light(
                start_level=locals.start_level,
                error_code=state.output.error_code,
                geopotential_height_cloud_levels_forced=locals.geopotential_height_cloud_levels_forced,
                cloud_vapor_mixing_ratio_forced=locals.cloud_vapor_mixing_ratio_forced,
                cloud_liquid_after_rain_forced=state.output.cloud_liquid_after_rain_forced,
                precipitable_water_updraft_forced=state.output.precipitable_water_updraft_forced,
                total_normalized_integrated_condensate_forced=state.output.total_normalized_integrated_condensate_forced,
                cloud_moist_static_energy_forced=locals.cloud_moist_static_energy_forced,
                unspecifid_temperature=locals.unspecifid_temperature,
                ocean_fraction=locals.ocean_fraction,
                convection_fraction=state.input.convection_fraction,
                surface_type=state.input.surface_type,
                p_forced=state.input_output.p_forced,
                cloud_top_level=state.output.cloud_top_level,
                buoyancy_forced=locals.buoyancy_forced,
                cloud_liquid_before_rain_forced=locals.cloud_liquid_before_rain_forced,
                t_cloud_levels=locals.t_cloud_levels,
                vapor_forced=locals.vapor_forced,
                gamma_cloud_levels_forced=locals.gamma_cloud_levels_forced,
                normalized_massflux_updraft_forced=state.output.normalized_massflux_updraft_forced,
                environment_saturation_mixing_ratio_cloud_levels_forced=locals.environment_saturation_mixing_ratio_cloud_levels_forced,
                updraft_origin_level=state.output.updraft_origin_level,
                vapor_cloud_levels_forced=locals.vapor_cloud_levels_forced,
                vapor_excess=locals.vapor_excess,
                mass_entrainment_updraft=locals.mass_entrainment_updraft,
                mass_detrainment_updraft=locals.mass_detrainment_updraft,
                psum=locals.psum,
                psumh=locals.psumh,
                c1d=locals.c1d,
                add_buoyancy=locals.add_buoyancy,
                AVERAGE_LAYER_DEPTH=plume_dependent_constants.AVERAGE_LAYER_DEPTH,
                C0=plume_dependent_constants.C0,
                plume=plume_dependent_constants.PLUME_INDEX,
            )

            self._updraft_vertical_velocity(
                vertical_velocity_3d=locals.vertical_velocity_3d,
                vertical_velocity_2d=locals.vertical_velocity_2d,
                convective_scale_velocity=state.input_output.convective_scale_velocity,
                entrainment_rate=state.output.entrainment_rate,
                detrainment_function_updraft=locals.detrainment_function_updraft,
                geopotential_height_cloud_levels_forced=locals.geopotential_height_cloud_levels_forced,
                t_cloud_levels_forced=locals.t_cloud_levels_forced,
                unspecifid_temperature=locals.unspecifid_temperature,
                cloud_vapor_mixing_ratio_forced=locals.cloud_vapor_mixing_ratio_forced,
                cloud_liquid_after_rain_forced=state.output.cloud_liquid_after_rain_forced,
                vapor_forced=locals.vapor_forced,
                updraft_lfc_level=state.output.updraft_lfc_level,
                cloud_top_level=state.output.cloud_top_level,
                error_code=state.output.error_code,
                plume=plume_dependent_constants.PLUME_INDEX,
            )


class MeltingProfile:
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
        self._get_melting_profile = stencil_factory.from_dims_halo(
            func=get_melting_profile,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "MELT_GLAC": cumulus_parameterization_config.MELT_GLAC,
            },
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):

        self._get_melting_profile(
            error_code=state.output.error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
            local_melting_layer=locals.melting_layer,
            local_partition_liquid_ice=locals.partition_liquid_ice,
            p_cloud_levels_forced=state.output.p_cloud_levels_forced,
            precipitable_water_updraft_forced=state.output.precipitable_water_updraft_forced,
            local_melting=locals.melting,
        )


class InCloudTemperature:
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
        self._in_cloud_updraft_air_temperature = stencil_factory.from_dims_halo(
            func=in_cloud_updraft_air_temperature,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        pass
        self._in_cloud_updraft_air_temperature(
            error_code=state.output.error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
            local_incloud_air_temp_forced=locals.incloud_air_temp_forced,
            local_hcdo=locals.hcdo,
            local_geopotential_height_cloud_levels_forced=locals.geopotential_height_cloud_levels_forced,
            local_incloud_water_vapor_mixing_ratio_forced=locals.cloud_vapor_mixing_ratio_forced,
            local_t_cloud_levels_forced=locals.t_cloud_levels_forced,
        )
