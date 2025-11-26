from ndsl.dsl.typing import FloatField, FloatFieldIJ, Float, IntFieldIJ, Int
from ndsl.dsl.gt4py import computation, PARALLEL, interval, FORWARD, K
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
import pyMoist.constants as constants
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import (
    IntFieldIJ_Plume,
    FloatFieldIJ_Plume,
    FloatField_Plume,
    FloatFieldIJ_Ensemble,
)


def trigger_function_convection(
    error_code: IntFieldIJ_Plume,
    plume: Int,
    local_cloud_work_function_0: FloatFieldIJ,
    convective_scale_velosity: FloatFieldIJ,
):
    from __externals__ import DICYCLE

    with computation(FORWARD), interval(...):
        if DICYCLE > 1:
            if error_code[0, 0][plume] == 0:
                # print*,"cin=",cin0(i),0.5*zws(i)**2, omeg(i,kpbl(i),1)/(-g*rho(i,kpbl(i))) ;call flush(6)
                if (
                    local_cloud_work_function_0 + 0.5 * convective_scale_velosity**2
                    < 0.0
                ):
                    error_code[0, 0][plume] = 19
                    # ierrc(i)="CIN negat"
