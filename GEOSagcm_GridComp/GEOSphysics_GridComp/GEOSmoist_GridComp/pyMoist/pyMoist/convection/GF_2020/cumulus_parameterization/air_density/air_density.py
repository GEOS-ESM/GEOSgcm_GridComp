from ndsl.dsl.gt4py import PARALLEL, computation, interval, FORWARD
from ndsl.dsl.typing import FloatField, Int
from ndsl import StencilFactory, QuantityFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.state import GF2020CumulusParameterizationState
from pyMoist.convection.GF_2020.cumulus_parameterization.locals import GF2020CumulusParameterizationLocals
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import FloatField_Plume, IntFieldIJ_Plume
import pyMoist.constants as constants


def hydrostatic_air_density(
    p: FloatField_Plume,
    geopotential_height: FloatField,
    error_code: IntFieldIJ_Plume,
    air_density: FloatField,
    plume: Int,
):
    """
    Compute air density, assuming hydrostatic balance

    Args:
        p (in): pressure
        geopotential_height (in): geopotential height
        air_density (out): air density
        plume (in): specifies the current plume
    """
    with computation(PARALLEL), interval(...):
        # prefil with 0
        air_density = 0.0

    with computation(PARALLEL), interval(0, -1):
        if error_code[0, 0][plume] == 0:
            air_density = (
                100.0
                * (p[0, 0, 0][plume] - p[0, 0, 1][plume])
                / (geopotential_height[0, 0, 1] - geopotential_height)
                / constants.MAPL_GRAV
            )


class HydrostaticAirDensity:
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
        self._hydrostatic_air_density = stencil_factory.from_dims_halo(
            func=hydrostatic_air_density,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._hydrostatic_air_density(
            p=state.output.p_cloud_levels_forced,
            geopotential_height=locals.geopotential_height_cloud_levels_forced,
            error_code=state.output.error_code,
            air_density=locals.air_density,
            plume=plume_dependent_constants.PLUME_INDEX,
        )
