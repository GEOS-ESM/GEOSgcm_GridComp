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
from pyMoist.convection.GF_2020.cumulus_parameterization.kinetic_energy_to_heating.stencils import (
    ke_to_heating,
)


class KineticEnergyToHeating:
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
        self._ke_to_heating = stencil_factory.from_dims_halo(
            func=ke_to_heating,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._ke_to_heating(
            del_u_cloud_ensemble=locals.del_u_cloud_ensemble,
            del_v_cloud_ensemble=locals.del_v_cloud_ensemble,
            error_code=state.output.error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
            cloud_top_level=state.output.cloud_top_level,
            p_cloud_levels_forced=state.output.p_cloud_levels_forced,
            u=state.input_output.u,
            v=state.input_output.v,
            del_t_cloud_ensemble=locals.del_t_cloud_ensemble,
        )
