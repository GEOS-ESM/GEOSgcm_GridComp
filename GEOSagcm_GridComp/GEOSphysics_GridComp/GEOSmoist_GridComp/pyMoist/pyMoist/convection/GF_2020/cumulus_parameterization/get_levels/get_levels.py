from ndsl import StencilFactory, QuantityFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.state import GF2020CumulusParameterizationState
from pyMoist.convection.GF_2020.cumulus_parameterization.locals import GF2020CumulusParameterizationLocals
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.get_levels.stencils import (
    find_maximum_updraft_origin_level,
    find_detrainmet_start_level,
    find_highest_moist_static_energy_level,
)


class MaximumUpdraftOriginLevel:
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
        self._find_maximum_updraft_origin_level = stencil_factory.from_dims_halo(
            func=find_maximum_updraft_origin_level,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._find_maximum_updraft_origin_level(
            geopotential_height=locals.geopotential_height_cloud_levels_forced,
            topography_height_no_negative=state.input_output.topography_height_no_negative,
            error_code=state.output.error_code,
            maximum_updraft_origin_level=locals.maximum_updraft_origin_level,
            MAX_UPDRAFT_ORIGIN_HEIGHT=plume_dependent_constants.MAX_UPDRAFT_ORIGIN_HEIGHT,
            plume=plume_dependent_constants.PLUME_INDEX,
        )


class DowndraftDetrainmentLevel:
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
        self._find_detrainmet_start_level = stencil_factory.from_dims_halo(
            func=find_detrainmet_start_level,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._find_detrainmet_start_level(
            geopotential_height=locals.geopotential_height_cloud_levels_forced,
            topography_height_no_negative=state.input_output.topography_height_no_negative,
            error_code=state.output.error_code,
            detrainment_start_level=locals.detrainment_start_level,
            DETRAINMENT_CRITICAL_DEPTH=plume_dependent_constants.DETRAINMENT_CRITICAL_DEPTH,
            plume=plume_dependent_constants.PLUME_INDEX,
        )


class HighestMoistStaticEnergyLevel:
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
        self._find_highest_moist_static_energy_level = stencil_factory.from_dims_halo(
            func=find_highest_moist_static_energy_level,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._find_highest_moist_static_energy_level(
            moist_static_energy=locals.environment_moist_static_energy_cloud_levels_forced,
            error_code=state.output.error_code,
            maximum_updraft_origin_level=locals.maximum_updraft_origin_level,
            updraft_origin_level=locals.updraft_origin_level,
            plume=plume_dependent_constants.PLUME_INDEX,
        )


class ConvectiveCloudBaseLevel:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass


class GetLCL:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass


class CloudTop:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass


class DowndraftOriginLevel:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass
