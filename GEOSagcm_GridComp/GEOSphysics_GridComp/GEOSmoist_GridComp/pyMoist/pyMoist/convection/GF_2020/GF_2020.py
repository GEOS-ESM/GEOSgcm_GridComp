import copy
from ndsl import StencilFactory, QuantityFactory
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.cumulus_parameterization import GF2020Config
from pyMoist.convection.GF_2020.state import GF2020State
from pyMoist.convection.GF_2020.locals import GF2020Locals
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.convection.GF_2020.cumulus_parameterization.state import GF2020CumulusParameterizationState
from pyMoist.convection_tracers import ConvectionTracers
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
from pyMoist.convection.GF_2020.setup import GF2020Setup
from pyMoist.convection.GF_2020.cumulus_parameterization.cumulus_parameterization import (
    GF2020CumulusParameterization,
)
from pyMoist.convection.GF_2020.finalize import GF2020Finalize


class GF2020:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GF2020Config,
        cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
        saturation_tables: SaturationVaporPressureTable | None,
    ):
        # make saturation tables visible at runtime
        if saturation_tables is None:
            saturation_tables = SaturationVaporPressureTable(stencil_factory.backend)
        else:
            self.saturation_tables = saturation_tables

        # initialize GF2020 locals
        self.locals = GF2020Locals.zeros(
            quantity_factory,
            data_dimensions={
                "plumes": 3,
                "convection_tracers": config.NUMBER_OF_TRACERS,
            },
        )

        # initialize GF2020 CumulusParameterization state
        self.cumulus_parameterization_state = GF2020CumulusParameterizationState.zeros(
            quantity_factory,
            data_dimensions={
                "plumes": cumulus_parameterization_constants.NUMBER_OF_PLUMES,
                "convection_tracers": config.NUMBER_OF_TRACERS,
            },
        )

        # initialize submodules
        self._setup = GF2020Setup(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            saturation_tables=saturation_tables,
        )

        self._cumulus_parameterization_core = GF2020CumulusParameterization(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
            saturation_tables=saturation_tables,
        )

        self._finalize = GF2020Finalize(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
            saturation_tables=saturation_tables,
        )

    def __call__(self, state: GF2020State, convection_tracers: ConvectionTracers):
        # flag to stop convection scheme for single column models
        # this can be set in setup if surface temperature is very near zero Kelvin
        scm_stop = False

        # call the there parts of the scheme
        self._setup(
            state=state,
            locals=self.locals,
            cumulus_parameterization_state=self.cumulus_parameterization_state,
            convection_tracers=convection_tracers,
            scm_stop=scm_stop,
        )

        if scm_stop == False:
            self._cumulus_parameterization_core(
                state=self.cumulus_parameterization_state,
                convection_tracers=convection_tracers,
            )

            self._finalize(
                state=state,
                locals=self.locals,
                cumulus_parameterization_state=self.cumulus_parameterization_state,
                convection_tracers=convection_tracers,
            )
