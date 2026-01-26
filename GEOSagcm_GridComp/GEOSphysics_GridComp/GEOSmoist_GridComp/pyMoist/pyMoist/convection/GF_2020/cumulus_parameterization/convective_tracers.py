from ndsl.dsl.gt4py import computation, interval, PARALLEL, FORWARD, K
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import (
    IntFieldIJ_Plume,
    FloatField_Tracers,
)
from ndsl.dsl.typing import Int
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl import StencilFactory, QuantityFactory, Quantity
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)


def environment_cloud_levels_chemistry(
    error_code: IntFieldIJ_Plume,
    chemistry_tracers: FloatField_Tracers,
    chemistry_tracers_cloud_levels: FloatField_Tracers,
    plume: Int,
):
    from __externals__ import k_end, NUMBER_OF_TRACERS

    with computation(FORWARD), interval(0, 1):
        # set up internal constants
        # NOTE only cloud_level_option = 2 has been tested
        cloud_level_option = 2

    with computation(FORWARD), interval(1, -1):
        if error_code[0, 0][plume] == 0 and cloud_level_option == 1:
            tracer = 0
            while tracer < NUMBER_OF_TRACERS:
                chemistry_tracers_cloud_levels[0, 0, 0][tracer] = (
                    0.5 * chemistry_tracers[0, 0, -1][tracer] + chemistry_tracers[0, 0, 0][tracer]
                )
                tracer += 1

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0 and cloud_level_option == 1:
            tracer = 0
            while tracer < NUMBER_OF_TRACERS:
                chemistry_tracers_cloud_levels[0, 0, 0][tracer] = chemistry_tracers[0, 0, 0][tracer]
                tracer += 1

    with computation(FORWARD), interval(-1, None):
        if error_code[0, 0][plume] == 0 and cloud_level_option == 1:
            tracer = 0
            while tracer < NUMBER_OF_TRACERS:
                chemistry_tracers_cloud_levels[0, 0, 0][tracer] = chemistry_tracers[0, 0, 0][tracer]
                tracer += 1

    with computation(PARALLEL), interval(0, -1):
        if error_code[0, 0][plume] == 0 and cloud_level_option != 1:
            tracer = 0
            while tracer < NUMBER_OF_TRACERS:
                chemistry_tracers_cloud_levels[0, 0, 0][tracer] = chemistry_tracers[0, 0, 0][tracer]
                tracer += 1


class ColdPoolParameterization:
    def __init__(self, cumulus_parameterization_config: GF2020CumulusParameterizationConfig):
        if cumulus_parameterization_config.CONVECTION_TRACER == 1:
            raise NotImplementedError(
                "The ColdPoolParameterization section has not been implemented. You should have been caught"
                "before getting here by the config checker. Beware, something likely failing in the config"
                "checker as well - you may be unknowingly calling other untested/unimplemented sections."
            )

    def __call__(self, *args, **kwds):
        pass


class AtmosphericComposition:
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
        self._environment_cloud_levels_chemistry = stencil_factory.from_dims_halo(
            func=environment_cloud_levels_chemistry,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"NUMBER_OF_TRACERS": config.NUMBER_OF_TRACERS},
        )

    def __call__(
        self,
        error_code: Quantity,
        chemistry_tracers: Quantity,
        chemistry_tracers_cloud_levels: Quantity,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        pass
        self._environment_cloud_levels_chemistry(
            error_code=error_code,
            chemistry_tracers=chemistry_tracers,
            chemistry_tracers_cloud_levels=chemistry_tracers_cloud_levels,
            plume=plume_dependent_constants.PLUME_INDEX,
        )
