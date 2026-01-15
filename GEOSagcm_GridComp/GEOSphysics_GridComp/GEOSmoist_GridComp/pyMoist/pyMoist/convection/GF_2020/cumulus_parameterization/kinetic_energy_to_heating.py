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


def kinetic_energy_to_heating(
    error_code: IntFieldIJ_Plume,
    cloud_top_level: IntFieldIJ_Plume,
    p_cloud_levels_forced: FloatField_Plume,
    u: FloatField,
    v: FloatField,
    del_u_cloud_ensemble: FloatField,
    del_v_cloud_ensemble: FloatField,
    del_t_cloud_ensemble: FloatField,
    plume: Int,
):
    with computation(FORWARD), interval(...):
        if error_code[0, 0][plume] == 0:
            dts: FloatFieldIJ = 0.0
            fpi: FloatFieldIJ = 0.0
    with computation(FORWARD), interval(...):
        if error_code[0, 0][plume] == 0:
            if K <= cloud_top_level[0, 0][plume]:
                dp = (p_cloud_levels_forced[0, 0, 0][plume] - p_cloud_levels_forced[0, 0, 1][plume]) * 100.0

                dts = dts - (
                    ((del_u_cloud_ensemble * u) + (del_v_cloud_ensemble * v)) * dp / constants.MAPL_GRAV
                )

                fpi = fpi + (
                    sqrt(
                        (del_u_cloud_ensemble * del_u_cloud_ensemble)
                        + (del_v_cloud_ensemble * del_v_cloud_ensemble)
                    )
                    * dp
                )

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0:
            if fpi > 0.0:
                if K <= cloud_top_level[0, 0][plume]:
                    fp = (
                        sqrt(
                            (
                                del_u_cloud_ensemble * del_u_cloud_ensemble
                                + del_v_cloud_ensemble * del_v_cloud_ensemble
                            )
                        )
                        / fpi
                    )

                    del_t_cloud_ensemble = del_t_cloud_ensemble + (
                        fp * dts * constants.MAPL_GRAV / cumulus_parameterization_constants.CP
                    )
