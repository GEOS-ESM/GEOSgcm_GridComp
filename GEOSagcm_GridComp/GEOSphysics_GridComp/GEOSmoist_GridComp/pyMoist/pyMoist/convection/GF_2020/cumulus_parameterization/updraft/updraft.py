from ndsl import StencilFactory, QuantityFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.state import GF2020CumulusParameterizationState
from pyMoist.convection.GF_2020.cumulus_parameterization.locals import GF2020CumulusParameterizationLocals
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)


class UpdraftMassFluxProfile:
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
        self._updraft_rates_pdf = stencil_factory.from_dims_halo(
            func=updraft_rates_pdf,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"OVERSHOOT": cumulus_parameterization_config.OVERSHOOT},
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._updraft_rates_pdf(
            entrainment_rate=state.output.entrainment_rate,
            moist_static_energy=locals.environment_moist_static_energy_forced,
            saturation_moist_static_energy=locals.environment_saturation_mixing_ratio_cloud_levels_forced,
            moist_static_energy_origin_level=locals.moist_static_energy_origin_level_forced,
            updraft_lfc_level=state.output.updraft_lfc_level,
            geopotential_height=locals.geopotential_height_cloud_levels_forced,
            cloud_moist_static_energy=locals.cloud_moist_static_energy_forced_transported,
            error_code=state.output.error_code,
            cloud_top=state.output.cloud_top,
            plume=plume_dependent_constants.PLUME_INDEX,
        )


class UpdraftMoistureProfile:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass


class UpdraftMoistStaticEnergyAndMomentumBudget:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass


class UpdraftInCloudUpdraftAirTemperature:
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
        self._in_cloud_updraft_air_temperature(
            error_code=state.output.error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
            # FIRST_GUESS_W=,
            # tempco=,
            # hco=,
            # zo_cup=,
            # qco=,
            # tn_cup=,
        )


class UpdraftInitialWorkfunctions:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass


class UpdraftCIN:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass


class UpdraftUpdateWorkfunctions:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass
