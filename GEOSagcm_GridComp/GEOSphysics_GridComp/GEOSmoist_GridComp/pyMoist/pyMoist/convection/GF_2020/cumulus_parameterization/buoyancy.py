from ndsl.dsl.gt4py import PARALLEL, computation, interval, FORWARD
from ndsl.dsl.typing import FloatField, FloatFieldIJ, Int
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import (
    IntFieldIJ_Plume,
)
from ndsl.dsl.gt4py import K


def get_buoyancy(
    lcl_level: IntFieldIJ_Plume,
    updraft_lfc_level: IntFieldIJ_Plume,
    cloud_top_level: IntFieldIJ_Plume,
    cloud_moist_static_energy: FloatField,
    environment_moist_static_energy: FloatField,
    environment_saturation_moist_static_energy: FloatField,
    d_buoyancy: FloatField,
    error_code: IntFieldIJ_Plume,
    plume: Int,
):
    """
    Determine the "d_buoyancy" of a parcel, defined as the difference between the
    moist static energy of the parcel and the environment

    Args:
        lcl_level (in): lcl level of environment
        updraft_lfc_level (in): lfc level of parcel
        cloud_top_level (in): equilibrium level of cloud
        cloud_moist_static_energy (in)
        environment_moist_static_energy (in)
        environment_saturation_moist_static_energy (in)
        d_buoyancy (out)
        error_code (in): field for stopping flow through the scheme and tracking errors
        plume (in): specifies the current plume
    """

    with computation(PARALLEL), interval(...):
        d_buoyancy = 0

        if error_code[0, 0][plume] == 0:
            if K <= lcl_level[0, 0][plume]:
                d_buoyancy = cloud_moist_static_energy - environment_moist_static_energy
            if K > lcl_level[0, 0][plume] and K <= cloud_top_level[0, 0][plume] + 1:
                d_buoyancy = cloud_moist_static_energy - environment_saturation_moist_static_energy
