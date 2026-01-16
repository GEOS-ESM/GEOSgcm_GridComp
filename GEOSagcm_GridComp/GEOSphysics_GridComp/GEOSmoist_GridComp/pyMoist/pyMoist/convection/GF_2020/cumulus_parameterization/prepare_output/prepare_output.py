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
from pyMoist.convection.GF_2020.cumulus_parameterization.prepare_output.stencils import (
    prepare_output,
)


class PrepareOutput:
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
        self._prepare_output = stencil_factory.from_dims_halo(
            func=prepare_output,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._prepare_output(
            error_code=state.output.error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
            cloud_base_mass_flux_modified=state.output.cloud_base_mass_flux_modified,
            total_normalized_integrated_condensate_forced=state.output.total_normalized_integrated_condensate_forced,
            total_normalized_integrated_evaporate_forced=state.output.total_normalized_integrated_evaporate_forced,
            normalized_massflux_updraft_forced=state.output.normalized_massflux_updraft_forced,
            normalized_massflux_downdraft_forced=state.output.normalized_massflux_downdraft_forced,
            condensate_to_fall_forced=state.output.condensate_to_fall_forced,
            evaporate_in_downdraft_forced=state.output.evaporate_in_downdraft_forced,
            mass_entrainment_updraft_forced=state.output.mass_entrainment_updraft_forced,
            mass_detrainment_updraft_forced=state.output.mass_detrainment_updraft_forced,
            mass_entrainment_downdraft_forced=state.output.mass_detrainment_downdraft_forced,
            mass_detrainment_downdraft_forced=state.output.mass_entrainment_downdraft_forced,
            environment_massflux=locals.environment_massflux,
            vapor_tendency_from_environmental_subsidence=locals.vapor_tendency_from_environmental_subsidence,
            moist_static_energy_tendency_from_environmental_subsidence=locals.moist_static_energy_tendency_from_environmental_subsidence,
            t_tendency_from_environmental_subsidence=locals.t_tendency_from_environmental_subsidence,
        )
