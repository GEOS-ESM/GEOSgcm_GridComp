from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import IntFieldIJ_Plume
from ndsl.dsl.typing import FloatFieldIJ, Int
from ndsl.dsl.gt4py import computation, FORWARD, interval


def convection_trigger(
    error_code: IntFieldIJ_Plume,
    convective_scale_velosity: FloatFieldIJ,
    cin_0: FloatFieldIJ,
    plume: Int,
):
    from __externals__ import DICYCLE

    with computation(FORWARD), interval(...):
        if DICYCLE > 1:
            if error_code[0, 0][plume] == 0:
                # think about including the grid scale vertical velocity at KE calculation
                if cin_0 + 0.5 * convective_scale_velosity**2 < 0.0:
                    error_code[0, 0][plume] = 19


class XieTriggerFunction:
    def __init__(self):
        raise NotImplementedError(
            "[NDSL] GF2020-->CumulusParameterization-->XieTriggerFunction this code"
            "has not been impemented. You should have been caught before getting here, but here we are. Please"
            "choose another option for ADV_TRIGGER or implement to continue."
        )

    def __call__(self, *args, **kwds):
        pass
