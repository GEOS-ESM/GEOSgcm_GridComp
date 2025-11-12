from ndsl import StencilFactory, QuantityFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.state import GF2020CumulusParameterizationState
from pyMoist.convection.GF_2020.cumulus_parameterization.locals import GF2020CumulusParameterizationLocals
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.entrainment.stencils import (
    entrainment_rates,
    downdraft_entrainment_profiles,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.shared_stencils import unknown_find_level


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
            externals={
                "ZERO_DIFF": cumulus_parameterization_config.ZERO_DIFF,
                "DOWNDRAFT": cumulus_parameterization_config.DOWNDRAFT,
            },
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
            externals={
                "ZERO_DIFF": cumulus_parameterization_config.ZERO_DIFF,
                "DOWNDRAFT": cumulus_parameterization_config.DOWNDRAFT,
            },
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
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass
