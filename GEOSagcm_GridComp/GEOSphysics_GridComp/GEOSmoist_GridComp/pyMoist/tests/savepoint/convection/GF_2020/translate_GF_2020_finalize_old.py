from ndsl import Namelist, Quantity, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.convective_parameterization.GF_2020.config import GF2020Config
from pyMoist.convective_parameterization.GF_2020.temporaries import Temporaries
from pyMoist.convective_parameterization.GF_2020.finalize import Finalize
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.convective_parameterization.GF_2020.state import MixingRatios, CloudFractions


class TranslateGF_2020_finalize(TranslateFortranData2Py):
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
            "U": grid.compute_dict(),
            "V": grid.compute_dict(),
            "Q": grid.compute_dict(),
            "T": grid.compute_dict(),
            "MASS": grid.compute_dict(),
            "PK": grid.compute_dict(),
            "DUDT_DC": grid.compute_dict(),
            "DVDT_DC": grid.compute_dict(),
            "DQVDT_DC": grid.compute_dict(),
            "DTDT_DC": grid.compute_dict(),
            "CNV_DQCDT": grid.compute_dict(),
            "CNV_FRC": grid.compute_dict(),
            "SRF_TYPE": grid.compute_dict(),
            "MFD_DC": grid.compute_dict(),
            "DQLDT_DC": grid.compute_dict(),
            "DQIDT_DC": grid.compute_dict(),
            "DQADT_DC": grid.compute_dict(),
            "REVSU": grid.compute_dict(),
            "RSU_CN": grid.compute_dict(),
            "REV_CN": grid.compute_dict(),
            "PRFIL": grid.compute_dict(),
            "PFI_CN": grid.compute_dict(),
            "PFL_CN": grid.compute_dict(),
            "QLCN": grid.compute_dict(),
            "QICN": grid.compute_dict(),
            "CLCN": grid.compute_dict(),
            "PL": grid.compute_dict(),
        }

        self.out_vars = {
            "U": grid.compute_dict(),
            "V": grid.compute_dict(),
            "Q": grid.compute_dict(),
            "T": grid.compute_dict(),
            "TH": grid.compute_dict(),
            "DQLDT_DC": grid.compute_dict(),
            "DQIDT_DC": grid.compute_dict(),
            "DQADT_DC": grid.compute_dict(),
            "RSU_CN": grid.compute_dict(),
            "REV_CN": grid.compute_dict(),
            "PFI_CN": grid.compute_dict(),
            "PFL_CN": grid.compute_dict(),
            "QLCN": grid.compute_dict(),
            "QICN": grid.compute_dict(),
            "CLCN": grid.compute_dict(),
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

        U = self.make_ijk_quantity(inputs.pop("U"))
        V = self.make_ijk_quantity(inputs.pop("V"))
        T = self.make_ijk_quantity(inputs.pop("T"))
        MASS = self.make_ijk_quantity(inputs.pop("MASS"))
        PK = self.make_ijk_quantity(inputs.pop("PK"))
        DUDT_DC = self.make_ijk_quantity(inputs.pop("DUDT_DC"))
        DVDT_DC = self.make_ijk_quantity(inputs.pop("DVDT_DC"))
        DQVDT_DC = self.make_ijk_quantity(inputs.pop("DQVDT_DC"))
        DTDT_DC = self.make_ijk_quantity(inputs.pop("DTDT_DC"))
        CNV_DQCDT = self.make_ijk_quantity(inputs.pop("CNV_DQCDT"))
        CNV_FRC = self.make_ij_quantity(inputs.pop("CNV_FRC"))
        SRF_TYPE = self.make_ij_quantity(inputs.pop("SRF_TYPE"))
        MFD_DC = self.make_ijk_quantity(inputs.pop("MFD_DC"))
        DQLDT_DC = self.make_ijk_quantity(inputs.pop("DQLDT_DC"))
        DQIDT_DC = self.make_ijk_quantity(inputs.pop("DQIDT_DC"))
        DQADT_DC = self.make_ijk_quantity(inputs.pop("DQADT_DC"))
        REVSU = self.make_ijk_quantity(inputs.pop("REVSU"))
        RSU_CN = self.make_ijk_quantity(inputs.pop("RSU_CN"))
        REV_CN = self.make_ijk_quantity(inputs.pop("REV_CN"))
        PRFIL = self.make_ijk_quantity(inputs.pop("PRFIL"), interface=True)
        PFI_CN = self.make_ijk_quantity(inputs.pop("PFI_CN"), interface=True)
        PFL_CN = self.make_ijk_quantity(inputs.pop("PFL_CN"), interface=True)
        PL = self.make_ijk_quantity(inputs.pop("PL"))

        mixing_ratios = MixingRatios(
            self.make_ijk_quantity(inputs.pop("Q")),
            None,
            None,
            None,
            self.make_ijk_quantity(inputs.pop("QLCN")),
            self.make_ijk_quantity(inputs.pop("QICN")),
            None,
            None,
        )

        cloud_fractions = CloudFractions(self.make_ijk_quantity(inputs.pop("CLCN")), None)

        # Construct stencils
        finalize = Finalize(
            stencil_factory=self.stencil_factory,
            GF_2020_config=GF_2020_config,
        )

        finalize(
            U,
            V,
            T,
            MASS,
            PK,
            DUDT_DC,
            DVDT_DC,
            DQVDT_DC,
            DTDT_DC,
            DQLDT_DC,
            DQIDT_DC,
            DQADT_DC,
            CNV_DQCDT,
            CNV_FRC,
            SRF_TYPE,
            MFD_DC,
            REVSU,
            RSU_CN,
            REV_CN,
            PRFIL,
            PFI_CN,
            PFL_CN,
            PL,
            mixing_ratios,
            cloud_fractions,
            self.temporaries,
            self.saturation_tables,
        )

        return {
            "U": U.field,
            "V": V.field,
            "Q": mixing_ratios.vapor.field,
            "T": T.field,
            "TH": self.temporaries.th.field,
            "DQLDT_DC": DQLDT_DC.field,
            "DQIDT_DC": DQIDT_DC.field,
            "DQADT_DC": DQADT_DC.field,
            "RSU_CN": RSU_CN.field,
            "REV_CN": REV_CN.field,
            "PFI_CN": PFI_CN.field,
            "PFL_CN": PFL_CN.field,
            "QLCN": mixing_ratios.convective_liquid.field,
            "QICN": mixing_ratios.convective_ice.field,
            "CLCN": cloud_fractions.convective.field,
        }
