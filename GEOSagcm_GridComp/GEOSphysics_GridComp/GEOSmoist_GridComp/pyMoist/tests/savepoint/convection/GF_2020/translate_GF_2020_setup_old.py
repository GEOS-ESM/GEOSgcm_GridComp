from ndsl import Namelist, Quantity, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.convective_parameterization.GF_2020.config import GF2020Config
from pyMoist.convective_parameterization.GF_2020.temporaries import Temporaries
from pyMoist.convective_parameterization.GF_2020.setup import Setup
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.convective_parameterization.GF_2020.state import MixingRatios


class TranslateGF_2020_setup(TranslateFortranData2Py):
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
        self.in_vars["data_vars"] = {
            "PLE": grid.compute_dict(),
            "ZLE": grid.compute_dict(),
            "T": grid.compute_dict(),
            "Q": grid.compute_dict(),
            "AREA": grid.compute_dict(),
            "CNV_FRC": grid.compute_dict(),
        }

        self.out_vars = {
            "TMP2D": grid.compute_dict(),
            "TPWI_star": grid.compute_dict(),
            "TPWI": grid.compute_dict(),
            "ZL0": grid.compute_dict(),
            "TH": grid.compute_dict(),
            "MASS": grid.compute_dict(),
            "ZLE0": grid.compute_dict(),
            "PL": grid.compute_dict(),
            "PK": grid.compute_dict(),
            "SEEDCNV": grid.compute_dict(),
        }

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
            DT_MOIST=self.constants["DT_MOIST"],
            STOCHASTIC_CONVECTION=bool(self.constants["STOCHASTIC_CNV"]),
            STOCH_TOP=self.constants["STOCH_TOP"],
            STOCH_BOT=self.constants["STOCH_BOT"],
            GF_MIN_AREA=self.constants["GF_MIN_AREA"],
            GF_ENV_SETTING=int(self.constants["GF_ENV_SETTING"]),
            ENTRVERSION=self.constants["ENTRVERSION"],
            CONVECTION_TRACER=self.constants["CONVECTION_TRACER"],
            C1=self.constants["C1"],
            ADV_TRIGGER=self.constants["ADV_TRIGGER"],
            AUTOCONV=self.constants["AUTOCONV"],
            USE_TRACER_TRANSP=self.constants["USE_TRACER_TRANSP"],
            SCLM_DEEP=self.constants["SCLM_DEEP"],
            FIX_CNV_CLOUD=bool(self.constants["FIX_CNV_CLOUD"]),
        )

        p_interface = self.make_ijk_quantity(inputs.pop("PLE"), interface=True)
        geopotential_height_interface = self.make_ijk_quantity(inputs.pop("ZLE"), interface=True)
        t = self.make_ijk_quantity(inputs.pop("T"))
        area = self.make_ij_quantity(inputs.pop("AREA"))
        convection_fraction = self.make_ij_quantity(inputs.pop("CNV_FRC"))

        mixing_ratios = MixingRatios(
            self.make_ijk_quantity(inputs.pop("Q")), None, None, None, None, None, None, None
        )

        # Construct stencils
        setup = Setup(
            stencil_factory=self.stencil_factory,
            GF_2020_config=GF_2020_config,
        )

        setup(
            p_interface,
            geopotential_height_interface,
            t,
            mixing_ratios,
            area,
            convection_fraction,
            self.temporaries,
            self.saturation_tables,
        )

        return {
            "TMP2D": self.temporaries.modified_area.field,
            "TPWI": self.temporaries.tpwi.field,
            "TPWI_star": self.temporaries.tpwi_star.field,
            "ZL0": self.temporaries.layer_height_above_surface.field,
            "TH": self.temporaries.th.field,
            "MASS": self.temporaries.mass.field,
            "ZLE0": self.temporaries.edge_height_above_surface.field,
            "PL": self.temporaries.p.field,
            "PK": self.temporaries.p_kappa.field,
            "SEEDCNV": self.temporaries.seed_convection.field,
        }
