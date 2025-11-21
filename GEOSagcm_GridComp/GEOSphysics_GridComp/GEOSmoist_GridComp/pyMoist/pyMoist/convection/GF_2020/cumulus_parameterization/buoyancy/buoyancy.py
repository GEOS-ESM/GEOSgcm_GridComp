from ndsl.dsl.gt4py import PARALLEL, computation, interval, FORWARD
from ndsl.dsl.typing import FloatField, Int
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
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import (
    FloatField_Plume,
    IntFieldIJ_Plume,
)
import pyMoist.constants as constants
from ndsl.dsl.gt4py import K


def get_buoyancy(
    lcl_level: IntFieldIJ_Plume,
    updraft_lfc_level: IntFieldIJ_Plume,
    cloud_top: IntFieldIJ_Plume,
    cloud_moist_static_energy: FloatField,
    environment_moist_static_energy: FloatField,
    environment_saturation_moist_static_energy: FloatField,
    buoyancy: FloatField,
    error_code: IntFieldIJ_Plume,
    plume: Int,
):
    """
    Determine the buoyancy of a parcel.

    Args:
        lcl_level (in): lcl level of environment
        updraft_lfc_level (in): lfc level of parcel
        cloud_top (in): equilibrium level of cloud
        cloud_moist_static_energy (in)
        environment_moist_static_energy (in)
        environment_saturation_moist_static_energy (in)
        buoyancy (out)
        error_code (in): field for stopping flow through the scheme and tracking errors
        plume (in): specifies the current plume
    """

    with computation(PARALLEL), interval(...):
        buoyancy = 0

        if error_code[0, 0][plume] == 0:
            if K <= lcl_level[0, 0][plume]:
                buoyancy = cloud_moist_static_energy - environment_moist_static_energy
            if K > lcl_level[0, 0][plume] and K <= cloud_top[0, 0][plume] + 1:
                buoyancy = (
                    cloud_moist_static_energy
                    - environment_saturation_moist_static_energy
                )


class GetBuoyancy:
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
        self._get_buoyancy = stencil_factory.from_dims_halo(
            func=get_buoyancy,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._get_buoyancy(
            lcl_level=state.output.lcl_level,
            updraft_lfc_level=state.output.updraft_lfc_level,
            cloud_top=state.output.cloud_top,
            cloud_moist_static_energy=locals.cloud_moist_static_energy_forced,
            environment_moist_static_energy=locals.environment_moist_static_energy_cloud_levels_forced,
            environment_saturation_moist_static_energy=locals.environment_saturation_moist_static_energy_cloud_levels_forced,
            buoyancy=locals.buoyancy,
            error_code=state.output.error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
        )
