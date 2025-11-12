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
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.downdraft.stencils import (
    moist_static_energy_and_moisture_budget,
)


class DowndraftNormalizedMassFlux:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass


class DowndraftLateralMassFlux:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass


class DowndraftWetBlub:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass


class DowndraftMoistStaticEnergyAndMoistureBudget:
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
        self._moist_static_energy_and_moisture_budget = stencil_factory.from_dims_halo(
            func=moist_static_energy_and_moisture_budget,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._moist_static_energy_and_moisture_budget(
            error_code=state.output.error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
            # heso_cup=,
            # u_cup=,
            # v_cup=,
            # cumulus=,
            # use_wetbulb=,
            # jmin=,
            # t_wetbulb=,
            # q_wetbulb=,
            # zo_cup=,
            # hc=,
            # zdo=,
            # dd_massdetro=,
            # dd_massentro=,
            # dd_massdetru=,
            # dd_massentru=,
            # us=,
            # vs=,
            # pgcon=,
            # heo=,
            # hcdo=,
        )


class DowndraftMoistureProperties:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass


class DowndraftWindshear:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass
