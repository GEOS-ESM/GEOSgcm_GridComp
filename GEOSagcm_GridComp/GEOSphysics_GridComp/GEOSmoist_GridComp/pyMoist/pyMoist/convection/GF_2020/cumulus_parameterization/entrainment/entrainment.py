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
from pyMoist.convection.GF_2020.cumulus_parameterization.entrainment.stencils import (
    entrainment_rates,
    downdraft_entrainment_profiles,
    compute_lateral_massflux,
    compute_uc_vc,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.shared_stencils import (
    unknown_find_level,
)
from ndsl.logging import ndsl_log


class EntrainmentRates:
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
        self._entrainment_rates = stencil_factory.from_dims_halo(
            func=entrainment_rates,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ZERO_DIFF": cumulus_parameterization_config.ZERO_DIFF},
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):

        if self.cumulus_parameterization_config.ZERO_DIFF != 0:
            ndsl_log.warning(
                " GF2020 cumulus parameterization called EntrainmentRates with "
                "untested ZERO_DIFF option. Running untested code... proceed with caution"
            )

        self._entrainment_rates(
            vapor=locals.vapor_cloud_levels_forced,
            environment_saturation_mixing_ratio=locals.environment_saturation_mixing_ratio_cloud_levels_forced,
            lcl_level=state.output.lcl_level,
            error_code=state.output.error_code,
            entrainment_rate=state.output.entrainment_rate,
            detrainment_function_updraft=locals.detrainment_function_updraft,
            plume=plume_dependent_constants.PLUME_INDEX,
        )


class DowndraftEntrainmentProfiles:
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
        self._downdraft_entrainment_profiles = stencil_factory.from_dims_halo(
            func=downdraft_entrainment_profiles,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"DOWNDRAFT": cumulus_parameterization_config.DOWNDRAFT},
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):

        self._downdraft_entrainment_profiles(
            lateral_entrainment_rate=state.input.lateral_entrainment_rate,
            entrainment_rate_downdraft=locals.entrainment_rate_downdraft,
            detrainment_function_downdraft=locals.detrainment_function_downdraft,
            scale_dependence_factor_downdraft=locals.scale_dependence_factor_downdraft,
            plume_entrainment_rate=plume_dependent_constants.ENTRAINMENT_RATE,
        )


class StableDetrainment:
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
        self._unknown_find_level = stencil_factory.from_dims_halo(
            func=unknown_find_level,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._unknown_find_level(
            array=locals.environment_saturation_moist_static_energy_cloud_levels_forced,
            start_index=state.output.updraft_lfc_level,
            end_index=locals.kstabm,
            out_index=state.output.kstabi,
            error_code=state.output.error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
        )


class CalculateMassEntrainmentDetrainment:
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
        self._compute_lateral_massflux = stencil_factory.from_dims_halo(
            func=compute_lateral_massflux,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._compute_uc_vc = stencil_factory.from_dims_halo(
            func=compute_uc_vc,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "BOUNDARY_CONDITION_METHOD": cumulus_parameterization_config.BOUNDARY_CONDITION_METHOD
            },
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):

        if self.cumulus_parameterization_config.BOUNDARY_CONDITION_METHOD != 1:
            ndsl_log.warning(
                " GF2020 cumulus parameterization called CalculateMassEntrainmentDetrainment with "
                "untested BOUNDARY_CONDITION_METHOD option. Running untested code... proceed with caution"
            )

        self._compute_lateral_massflux(
            error_code=state.output.error_code,
            cloud_top=state.output.cloud_top,
            geopotential_height=locals.geopotential_height_cloud_levels_forced,
            normalized_massflux_updraft=state.output.normalized_massflux_updraft_forced,
            detrainment_function_updraft=locals.detrainment_function_updraft,
            entrainment_rate=state.output.entrainment_rate,
            p_cloud_levels_forced=state.output.p_cloud_levels_forced,
            mass_entrainment_updraft_forced=state.output.mass_entrainment_updraft_forced,
            mass_detrainment_updraft_forced=state.output.mass_detrainment_updraft_forced,
            mass_entrainment_updraft=locals.mass_entrainment_updraft,
            mass_detrainment_updraft=locals.mass_detrainment_updraft,
            updraft_lfc_level=state.output.updraft_lfc_level,
            updraft_origin_level=state.output.updraft_origin_level,
            pbl_level=state.input_output.pbl_level,
            mass_entrainment_u_updraft=locals.mass_entrainment_u_updraft,
            mass_detrainment_u_updraft=locals.mass_detrainment_u_updraft,
            LAMBDA_DEEP=plume_dependent_constants.LAMBDA_DEEP,
            plume=plume_dependent_constants.PLUME_INDEX,
        )

        self._compute_uc_vc(
            u_c=locals.u_c,
            v_c=locals.v_c,
            cloud_moist_static_energy=locals.cloud_moist_static_energy,
            cloud_moist_static_energy_forced=locals.cloud_moist_static_energy_forced,
            error_code=state.output.error_code,
            start_level=locals.start_level,
            moist_static_energy_origin_level=locals.moist_static_energy_origin_level,
            moist_static_energy_origin_level_forced=locals.moist_static_energy_origin_level_forced,
            u_cloud_levels=locals.u_cloud_levels,
            v_cloud_levels=locals.v_cloud_levels,
            p=state.input_output.p_forced,
            updraft_origin_level=state.output.updraft_origin_level,
            ocean_fraction=state.input.ocean_fraction,
            AVERAGE_LAYER_DEPTH=plume_dependent_constants.AVERAGE_LAYER_DEPTH,
            plume=plume_dependent_constants.PLUME_INDEX,
        )
