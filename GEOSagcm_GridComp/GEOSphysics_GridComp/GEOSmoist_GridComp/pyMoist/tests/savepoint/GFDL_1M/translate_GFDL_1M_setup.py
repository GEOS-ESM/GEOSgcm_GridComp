from ndsl import Namelist, StencilFactory, Quantity
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.GFDL_1M.masks import Masks
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.GFDL_1M.temporaries import Temporaries
from pyMoist.GFDL_1M.outputs import Outputs
from pyMoist.GFDL_1M.setup import Setup
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Float
from pyMoist.GFDL_1M.stencils import prepare_tendencies
from ndsl.stencils.testing.savepoint import DataLoader
from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.GFDL_1M.state import (
    MixingRatios,
    CloudFractions,
)


class TranslateGFDL_1M_setup(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, namelist, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        # grid.compute_dict is workaround to remove grid halo, which is hardcoded to 3
        self.in_vars["data_vars"] = {
            "PLE": grid.compute_dict(),
            "ZLE": grid.compute_dict(),
            "T": grid.compute_dict(),
            "U": grid.compute_dict(),
            "V": grid.compute_dict(),
            "Q": grid.compute_dict(),
            "QRAIN": grid.compute_dict(),
            "QSNOW": grid.compute_dict(),
            "QGRAUPEL": grid.compute_dict(),
            "QLLS": grid.compute_dict(),
            "QLCN": grid.compute_dict(),
            "CLCN": grid.compute_dict(),
            "CLLS": grid.compute_dict(),
            "QILS": grid.compute_dict(),
            "QICN": grid.compute_dict(),
            "SHLW_PRC3": grid.compute_dict(),
            "SHLW_SNO3": grid.compute_dict(),
        }

        self.out_vars = {
            "PLEmb": grid.compute_dict(),
            "PLmb": grid.compute_dict(),
            "ZLE0": grid.compute_dict(),
            "ZL0": grid.compute_dict(),
            "DZET": grid.compute_dict(),
            "DP": grid.compute_dict(),
            "MASS": grid.compute_dict(),
            "iMASS": grid.compute_dict(),
            "U0": grid.compute_dict(),
            "V0": grid.compute_dict(),
            "QST3": grid.compute_dict(),
            "DQST3": grid.compute_dict(),
            "KLCL": grid.compute_dict(),
            "LTS": grid.compute_dict(),
            "EIS": grid.compute_dict(),
            "QRAIN": grid.compute_dict(),
            "QSNOW": grid.compute_dict(),
            "DUDT_macro": grid.compute_dict(),
            "DVDT_macro": grid.compute_dict(),
            "DTDT_macro": grid.compute_dict(),
            "DQVDT_macro": grid.compute_dict(),
            "DQLDT_macro": grid.compute_dict(),
            "DQIDT_macro": grid.compute_dict(),
            "DQADT_macro": grid.compute_dict(),
            "DQRDT_macro": grid.compute_dict(),
            "DQSDT_macro": grid.compute_dict(),
            "DQGDT_macro": grid.compute_dict(),
        }

        # Initalize saturation tables
        self.saturation_tables = SaturationVaporPressureTable(self.stencil_factory.backend)

        # Initalize extra quantities
        self.outputs = Outputs.make(self.quantity_factory)
        self.temporaries = Temporaries.make(self.quantity_factory)
        self.masks = Masks.make(self.quantity_factory)

    def make_ijk_quantity(self, data, interface: bool = False) -> Quantity:
        if interface == True:
            quantity = self.quantity_factory.empty([X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a")
            quantity.view[:, :, :] = quantity.np.asarray(data[:, :, :])
            return quantity
        else:
            quantity = self.quantity_factory.empty([X_DIM, Y_DIM, Z_DIM], "n/a")
            quantity.view[:, :, :] = quantity.np.asarray(data[:, :, :])
            return quantity

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GFDL_1M-constants")

    def compute(self, inputs):
        GFDL_1M_config = GFDL1MConfig(
            PHYS_HYDROSTATIC=bool(self.constants["LPHYS_HYDROSTATIC"]),
            HYDROSTATIC=bool(self.constants["LHYDROSTATIC"]),
            DT_MOIST=self.constants["DT_MOIST"],
            MP_TIME=self.constants["MP_TIME"],
            T_MIN=self.constants["T_MIN"],
            T_SUB=self.constants["T_SUB"],
            TAU_R2G=self.constants["TAU_R2G"],
            TAU_SMLT=self.constants["TAU_SMLT"],
            TAU_G2R=self.constants["TAU_G2R"],
            DW_LAND=self.constants["DW_LAND"],
            DW_OCEAN=self.constants["DW_OCEAN"],
            VI_FAC=self.constants["VI_FAC"],
            VR_FAC=self.constants["VR_FAC"],
            VS_FAC=self.constants["VS_FAC"],
            VG_FAC=self.constants["VG_FAC"],
            QL_MLT=self.constants["QL_MLT"],
            DO_QA=bool(self.constants["DO_QA"]),
            FIX_NEGATIVE=bool(self.constants["FIX_NEGATIVE"]),
            VI_MAX=self.constants["VI_MAX"],
            VS_MAX=self.constants["VS_MAX"],
            VG_MAX=self.constants["VG_MAX"],
            VR_MAX=self.constants["VR_MAX"],
            QS_MLT=self.constants["QS_MLT"],
            QS0_CRT=self.constants["QS0_CRT"],
            QI_GEN=self.constants["QI_GEN"],
            QL0_MAX=self.constants["QL0_MAX"],
            QI0_MAX=self.constants["QI0_MAX"],
            QI0_CRT=self.constants["QI0_CRT"],
            QR0_CRT=self.constants["QR0_CRT"],
            FAST_SAT_ADJ=bool(self.constants["FAST_SAT_ADJ"]),
            RH_INC=self.constants["RH_INC"],
            RH_INS=self.constants["RH_INS"],
            RH_INR=self.constants["RH_INR"],
            CONST_VI=bool(self.constants["CONST_VI"]),
            CONST_VS=bool(self.constants["CONST_VS"]),
            CONST_VG=bool(self.constants["CONST_VG"]),
            CONST_VR=bool(self.constants["CONST_VR"]),
            USE_CCN=bool(self.constants["USE_CCN"]),
            RTHRESHU=self.constants["RTHRESHU"],
            RTHRESHS=self.constants["RTHRESHS"],
            CCN_L=self.constants["CCN_L"],
            CCN_O=self.constants["CCN_O"],
            QC_CRT=self.constants["QC_CRT"],
            TAU_G2V=self.constants["TAU_G2V"],
            TAU_V2G=self.constants["TAU_V2G"],
            TAU_S2V=self.constants["TAU_S2V"],
            TAU_V2S=self.constants["TAU_V2S"],
            TAU_REVP=self.constants["TAU_REVP"],
            TAU_FRZ=self.constants["TAU_FRZ"],
            DO_BIGG=bool(self.constants["DO_BIGG"]),
            DO_EVAP=bool(self.constants["DO_EVAP"]),
            DO_SUBL=bool(self.constants["DO_SUBL"]),
            SAT_ADJ0=self.constants["SAT_ADJ0"],
            C_PIACR=self.constants["C_PIACR"],
            TAU_IMLT=self.constants["TAU_IMLT"],
            TAU_V2L=self.constants["TAU_V2L"],
            TAU_L2V=self.constants["TAU_L2V"],
            TAU_I2V=self.constants["TAU_I2V"],
            TAU_I2S=self.constants["TAU_I2S"],
            TAU_L2R=self.constants["TAU_L2R"],
            QI_LIM=self.constants["QI_LIM"],
            QL_GEN=self.constants["QL_GEN"],
            C_PAUT=self.constants["C_PAUT"],
            C_PSACI=self.constants["C_PSACI"],
            C_PGACS=self.constants["C_PGACS"],
            C_PGACI=self.constants["C_PGACI"],
            Z_SLOPE_LIQ=bool(self.constants["Z_SLOPE_LIQ"]),
            Z_SLOPE_ICE=bool(self.constants["Z_SLOPE_ICE"]),
            PROG_CCN=bool(self.constants["PROG_CCN"]),
            C_CRACW=self.constants["C_CRACW"],
            ALIN=self.constants["ALIN"],
            CLIN=self.constants["CLIN"],
            PRECIPRAD=bool(self.constants["PRECIPRAD"]),
            CLD_MIN=self.constants["CLD_MIN"],
            USE_PPM=bool(self.constants["USE_PPM"]),
            MONO_PROF=bool(self.constants["MONO_PROF"]),
            DO_SEDI_HEAT=bool(self.constants["DO_SEDI_HEAT"]),
            SEDI_TRANSPORT=bool(self.constants["SEDI_TRANSPORT"]),
            DO_SEDI_W=bool(self.constants["DO_SEDI_W"]),
            DE_ICE=bool(self.constants["DE_ICE"]),
            ICLOUD_F=self.constants["ICLOUD_F"],
            IRAIN_F=self.constants["IRAIN_F"],
            MP_PRINT=bool(self.constants["MP_PRINT"]),
            MELTFRZ=bool(self.constants["LMELTFRZ"]),
            USE_BERGERON=bool(self.constants["USE_BERGERON"]),
            TURNRHCRIT_PARAM=self.constants["TURNRHCRIT_PARAM"],
            PDF_SHAPE=self.constants["PDFSHAPE"],
            ANV_ICEFALL=self.constants["ANV_ICEFALL"],
            LS_ICEFALL=self.constants["LS_ICEFALL"],
            LIQ_RADII_PARAM=self.constants["LIQ_RADII_PARAM"],
            ICE_RADII_PARAM=self.constants["ICE_RADII_PARAM"],
            FAC_RI=self.constants["FAC_RI"],
            MIN_RI=self.constants["MIN_RI"],
            MAX_RI=self.constants["MAX_RI"],
            FAC_RL=self.constants["FAC_RL"],
            MIN_RL=self.constants["MIN_RL"],
            MAX_RL=self.constants["MAX_RL"],
            CCW_EVAP_EFF=self.constants["CCW_EVAP_EFF"],
            CCI_EVAP_EFF=self.constants["CCI_EVAP_EFF"],
        )
        p_interface = self.make_ijk_quantity(inputs.pop("PLE"), interface=True)
        geopotential_height_interface = self.make_ijk_quantity(inputs.pop("ZLE"), interface=True)
        t = self.make_ijk_quantity(inputs.pop("T"))
        u = self.make_ijk_quantity(inputs.pop("U"))
        v = self.make_ijk_quantity(inputs.pop("V"))
        mixing_ratios = MixingRatios(
            vapor=self.make_ijk_quantity(inputs.pop("Q")),
            rain=self.make_ijk_quantity(inputs.pop("QRAIN")),
            snow=self.make_ijk_quantity(inputs.pop("QSNOW")),
            graupel=self.make_ijk_quantity(inputs.pop("QGRAUPEL")),
            convective_liquid=self.make_ijk_quantity(inputs.pop("QLCN")),
            convective_ice=self.make_ijk_quantity(inputs.pop("QICN")),
            large_scale_liquid=self.make_ijk_quantity(inputs.pop("QLLS")),
            large_scale_ice=self.make_ijk_quantity(inputs.pop("QILS")),
        )
        cloud_fractions = CloudFractions(
            convective=self.make_ijk_quantity(inputs.pop("CLCN")),
            large_scale=self.make_ijk_quantity(inputs.pop("CLLS")),
        )
        shallow_convective_rain = self.make_ijk_quantity(inputs.pop("SHLW_PRC3"))
        shallow_convective_snow = self.make_ijk_quantity(inputs.pop("SHLW_SNO3"))

        # Construct stencils
        self.prepare_tendencies = self.stencil_factory.from_dims_halo(
            func=prepare_tendencies,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        setup = Setup(
            stencil_factory=self.stencil_factory,
            GFDL_1M_config=GFDL_1M_config,
            saturation_tables=self.saturation_tables,
            prepare_tendencies=self.prepare_tendencies,
        )

        setup(
            geopotential_height_interface=geopotential_height_interface,
            p_interface=p_interface,
            t=t,
            u=u,
            v=v,
            shallow_convective_rain=shallow_convective_rain,
            shallow_convective_snow=shallow_convective_snow,
            mixing_ratios=mixing_ratios,
            cloud_fractions=cloud_fractions,
            masks=self.masks,
            outputs=self.outputs,
            temporaries=self.temporaries,
        )

        return {
            "PLEmb": self.temporaries.p_interface_mb.field,
            "PLmb": self.temporaries.p_mb.field,
            "ZLE0": self.temporaries.edge_height_above_surface.field,
            "ZL0": self.temporaries.layer_height_above_surface.field,
            "DZET": self.temporaries.layer_thickness.field,
            "DP": self.temporaries.dp.field,
            "MASS": self.temporaries.mass.field,
            "iMASS": 1 / self.temporaries.mass.field,
            "U0": self.temporaries.u_unmodified.field,
            "V0": self.temporaries.v_unmodified.field,
            "QST3": self.temporaries.qsat.field,
            "DQST3": self.temporaries.dqsat.field,
            "KLCL": self.temporaries.k_lcl.field + 1,  # add 1 b/c python indexing starts at 0
            "LTS": self.outputs.lower_tropospheric_stability.field,
            "EIS": self.outputs.estimated_inversion_strength.field,
            "QRAIN": mixing_ratios.rain.field,
            "QSNOW": mixing_ratios.snow.field,
            "DUDT_macro": self.outputs.du_dt_macro.field,
            "DVDT_macro": self.outputs.dv_dt_macro.field,
            "DTDT_macro": self.outputs.dt_dt_macro.field,
            "DQVDT_macro": self.outputs.dvapor_dt_macro.field,
            "DQLDT_macro": self.outputs.dliquid_dt_macro.field,
            "DQIDT_macro": self.outputs.dice_dt_macro.field,
            "DQADT_macro": self.outputs.dcloud_fraction_dt_macro.field,
            "DQRDT_macro": self.outputs.drain_dt_macro.field,
            "DQSDT_macro": self.outputs.dsnow_dt_macro.field,
            "DQGDT_macro": self.outputs.dgraupel_dt_macro.field,
        }
