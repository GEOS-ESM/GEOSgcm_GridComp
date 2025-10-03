from ndsl import Namelist, Quantity, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.convective_parameterization.GF_2020.config import GF2020Config
from pyMoist.convective_parameterization.GF_2020.temporaries import Temporaries
from pyMoist.convective_parameterization.GF_2020.driver.setup import DriverSetup
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.convective_parameterization.GF_2020.state import MixingRatios


class TranslateGF_2020_driver_setup(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, namelist, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory
        self.quantity_factory.set_extra_dim_lengths(
            **{
                "nmp": 2,
                "maxiens": 3,
            }
        )

        # grid.compute_dict is workaround to remove grid halo, which is hardcoded to 3
        self.in_vars["data_vars"] = {}

        self.out_vars = {}

        # Initalize saturation tables
        self.saturation_tables = SaturationVaporPressureTable(self.stencil_factory.backend)

        # Initalize extra quantities
        self.temporaries = Temporaries.make(self.quantity_factory)

    def make_ijk_quantity(self, data, interface: bool = False) -> Quantity:
        if interface is True:
            quantity = self.quantity_factory.empty([X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a")
            quantity.view[:, :, :] = quantity.np.asarray(data[:, :, :])
            return quantity
        else:
            quantity = self.quantity_factory.empty([X_DIM, Y_DIM, Z_DIM], "n/a")
            quantity.view[:, :, :] = quantity.np.asarray(data[:, :, :])
            return quantity

    def make_ij_quantity(self, data) -> Quantity:
        quantity = self.quantity_factory.empty([X_DIM, Y_DIM], "n/a")
        quantity.view[:, :] = quantity.np.asarray(data[:, :])
        return quantity

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF_2020-constants")

    def compute(self, inputs):
        GF_2020_config = GF2020Config(
            STOCHASTIC_CONVECTION=bool(self.constants["STOCHASTIC_CNV"]),
            STOCH_TOP=self.constants["STOCH_TOP"],
            STOCH_BOT=self.constants["STOCH_BOT"],
            GF_MIN_AREA=self.constants["GF_MIN_AREA"],
            GF_ENV_SETTING=int(self.constants["GF_ENV_SETTING"]),
            ENTRVERSION=self.constants["ENTRVERSION"],
            CONVECTION_TRACER=self.constants["CONVECTION_TRACER"],
            C1=self.constants["C1"],
            ADV_TRIGGER=self.constants["ADV_TRIGGER"],
        )

        # Construct stencils
        driver_setup = DriverSetup(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            GF_2020_config=GF_2020_config,
        )

        driver_setup()

        return {}
