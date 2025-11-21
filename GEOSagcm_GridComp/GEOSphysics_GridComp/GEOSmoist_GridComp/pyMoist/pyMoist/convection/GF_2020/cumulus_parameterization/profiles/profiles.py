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
from pyMoist.convection.GF_2020.cumulus_parameterization.profiles.stencils import (
    in_cloud_updraft_air_temperature,
    get_melting_profile,
)


class C1DProfile:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass


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
                "MELT_ICE": cumulus_parameterization_config.MELT_ICE,
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
            # tempcdo=,
            # hcdo=,
            # zo_cup=,
            # qcdo=,
            # tn_cup=,
        )
