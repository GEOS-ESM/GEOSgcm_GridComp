from ndsl.dsl.typing import FloatField, IntFieldIJ, Int
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import IntFieldIJ_Plume
from gt4py.cartesian.gtscript import (
    FORWARD,
    computation,
    interval,
)


def unknown_find_level(
    array: FloatField,
    start_index: IntFieldIJ_Plume,
    end_index: IntFieldIJ,
    out_index: IntFieldIJ_Plume,
    error_code: IntFieldIJ_Plume,
    plume: Int,
):
    """
    Details/purpose unknown

    Args:
        array (in): array to be analyzed
        start_index (in): start index for analysis
        end_index (in): end index for analysis
        out_index (out): output of analysis
        error_code (in): field for stopping flow through the scheme and tracking errors
        plume (in): specifies the current plume
    """
    with computation(FORWARD), interval(0, 1):
        out_index = start_index[0, 0][plume]

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            x = array.at(K=start_index[0, 0][plume])
            stop_index = max(start_index[0, 0][plume] + 1, end_index)

            level = start_index[0, 0][plume] + 1
            while level >= start_index[0, 0][plume] + 1 and level <= stop_index:
                if array.at(K=level) < x:
                    x = array.at(K=level)
                    out_index[0, 0][plume] = level
                level += 1
