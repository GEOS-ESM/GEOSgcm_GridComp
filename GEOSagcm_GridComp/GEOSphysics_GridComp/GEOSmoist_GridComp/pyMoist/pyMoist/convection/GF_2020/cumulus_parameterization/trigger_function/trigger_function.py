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
from pyMoist.convection.GF_2020.cumulus_parameterization.trigger_function.stencils import (
    trigger_function_convection,
)


class TriggerFunctionConvection:
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
        self._trigger_function_convection = stencil_factory.from_dims_halo(
            func=trigger_function_convection,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"DICYCLE": cumulus_parameterization_config.DICYCLE},
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):

        self._trigger_function_convection(
            error_code=state.output.error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
            local_cloud_work_function_0=locals.cloud_work_function_0,
            convective_scale_velosity=state.input_output.convective_scale_velocity,
        )


class TriggerFunctionXie:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass
