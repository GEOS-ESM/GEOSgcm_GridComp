from ndsl.dsl.typing import FloatField, FloatFieldIJ, Float, IntFieldIJ, Int
from ndsl.dsl.gt4py import computation, PARALLEL, interval, FORWARD, K, sqrt
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
import pyMoist.constants as constants
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import (
    IntFieldIJ_Plume,
    FloatFieldIJ_Plume,
    FloatField_Plume,
    FloatFieldIJ_Ensemble,
)


def ke_to_heating(
    dellu: FloatField,
    dellv: FloatField,
    error_code: IntFieldIJ_Plume,
    plume: Int,
    ktop: IntFieldIJ,
    po_cup: FloatField,
    us: FloatField,
    vs: FloatField,
    dellat: FloatField,
):
    with computation(FORWARD), interval(...):
        if error_code[0, 0][plume] == 0:
            dts: FloatFieldIJ = 0.0
            fpi: FloatFieldIJ = 0.0
    with computation(FORWARD), interval(...):
        if error_code[0, 0][plume] == 0:
            if K <= ktop - 1:
                dp = (po_cup - po_cup[0, 0, 1]) * 100.0

                dts = dts - (((dellu * us) + (dellv * vs)) * dp / constants.MAPL_GRAV)

                fpi = fpi + (sqrt((dellu * dellu) + (dellv * dellv)) * dp)

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0:
            if fpi > 0.0:
                if K <= ktop - 1:
                    fp = sqrt((dellu * dellu + dellv * dellv)) / fpi

                    dellat = dellat + (fp * dts * constants.MAPL_GRAV / cumulus_parameterization_constants.CP)
