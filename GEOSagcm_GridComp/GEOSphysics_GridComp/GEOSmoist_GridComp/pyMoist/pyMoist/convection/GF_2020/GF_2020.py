import copy
from pyMoist.convection.GF_2020.config import GF2020Config
from ndsl import StencilFactory, QuantityFactory
from pyMoist.convection.GF_2020.locals import GF2020Locals
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.convection.GF_2020.state import GF2020State
from pyMoist.convection.GF_2020.setup import GF2020Setup


class GF2020:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        GF_2020_config: GF2020Config,
    ):
        # Create extra quantity factories
        self.nmp_quantity_factory = self.make_nmp_quantity_factory(quantity_factory)

        # initialize saturation tables
        self.saturation_tables = SaturationVaporPressureTable(self.stencil_factory.backend)

        # initialize extra quantities
        self.temporaries = GF2020Locals.make(self.quantity_factory)

        # initialize submodules and build stencils
        self.setup = GF2020Setup(stencil_factory, quantity_factory, GF_2020_config)

    @staticmethod
    def make_nmp_quantity_factory(
        quantity_factory: QuantityFactory,
    ):
        nmp_quantity_factory = copy.deepcopy(quantity_factory)
        nmp_quantity_factory.set_extra_dim_lengths(
            **{
                "nmp": 2,
            }
        )
        return nmp_quantity_factory

    def __call__(self, p_interface, GF_2020_state: GF2020State):
        self.setup()
