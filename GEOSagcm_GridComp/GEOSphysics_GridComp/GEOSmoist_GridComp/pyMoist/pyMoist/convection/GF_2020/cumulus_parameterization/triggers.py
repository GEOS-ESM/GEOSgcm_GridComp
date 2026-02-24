from ndsl import QuantityFactory, StencilFactory
from ndsl.dsl.gt4py import FORWARD, computation, interval
from ndsl.dsl.typing import FloatFieldIJ, Int
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import IntFieldIJ_Plume
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)


def convection_trigger(
    error_code: IntFieldIJ_Plume,
    convective_scale_velosity: FloatFieldIJ,
    cin_0: FloatFieldIJ,
    plume: Int,
):
    """Check if sufficient velocity exists to overcome present convective inhibition.

    Args:
        error_code (IntFieldIJ_Plume)
        convective_scale_velosity (FloatFieldIJ)
        cin_0 (FloatFieldIJ)
        plume (Int)
    """
    from __externals__ import DICYCLE

    with computation(FORWARD), interval(...):
        if DICYCLE > 1:
            if error_code[0, 0][plume] == 0:
                # think about including the grid scale vertical velocity at KE calculation
                if cin_0 + 0.5 * convective_scale_velosity ** 2 < 0.0:
                    error_code[0, 0][plume] = 19


class XieTriggerFunction:
    def __init__(
        self,
        config: GF2020Config,
        cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
    ):
        self.config = config
        self.cumulus_parameterization_config = cumulus_parameterization_config

    def __call__(self, plume_dependent_constants: GF2020PlumeDependentConstants):

        if self.config.ADV_TRIGGER == 3 and (plume_dependent_constants.PLUME_INDEX in (1, 2)):
            raise NotImplementedError(
                "[NDSL] GF2020-->CumulusParameterization-->XieTriggerFunction this code"
                "has not been impemented. You should have been caught before getting here by the config"
                "checker. Beware, something likely failing in the config checker as well - you may be"
                "unknowingly calling other untested/unimplemented sections."
            )
