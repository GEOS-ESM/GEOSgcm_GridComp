from ndsl import Namelist, StencilFactory, Quantity
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.GFDL_1M.masks import Masks
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.GFDL_1M.temporaries import Temporaries
from pyMoist.GFDL_1M.outputs import Outputs
from pyMoist.GFDL_1M.finalize import Finalize
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Float
from pyMoist.GFDL_1M.stencils import update_tendencies
from ndsl.stencils.testing.savepoint import DataLoader
from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.GFDL_1M.state import (
    MixingRatios,
    CloudFractions,
)
from pyMoist.GFDL_1M.driver.driver import MicrophysicsDriver


class TranslateGFDL_1M_finalize(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        # self.override_input_netcdf_name = "None"
        # self.override_output_netcdf_name = "GFDL_1M_finalize-part-4"
        super().__init__(grid, namelist, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        # grid.compute_dict is workaround to remove grid halo, which is hardcoded to 3
        self.in_vars["data_vars"] = {
            "T": grid.compute_dict(),
            "U": grid.compute_dict(),
            "V": grid.compute_dict(),
            "PLmb": grid.compute_dict(),
            "MASS": grid.compute_dict(),
            "RAD_CF": grid.compute_dict(),
            "RAD_QI": grid.compute_dict(),
            "RAD_QL": grid.compute_dict(),
            "RAD_QV": grid.compute_dict(),
            "RAD_QR": grid.compute_dict(),
            "RAD_QS": grid.compute_dict(),
            "RAD_QG": grid.compute_dict(),
            "CLCN": grid.compute_dict(),
            "CLLS": grid.compute_dict(),
            "QLCN": grid.compute_dict(),
            "QLLS": grid.compute_dict(),
            "NACTL": grid.compute_dict(),
            "CLDREFFL": grid.compute_dict(),
            "QICN": grid.compute_dict(),
            "QILS": grid.compute_dict(),
            "NACTI": grid.compute_dict(),
            "CLDREFFI": grid.compute_dict(),
            "Q": grid.compute_dict(),
            "QRAIN": grid.compute_dict(),
            "QSNOW": grid.compute_dict(),
            "QGRAUPEL": grid.compute_dict(),
            "REV_LS": grid.compute_dict(),
            "RSU_LS": grid.compute_dict(),
            "PRCP_RAIN": grid.compute_dict(),
            "PRCP_SNOW": grid.compute_dict(),
            "PRCP_ICE": grid.compute_dict(),
            "PRCP_GRAUPEL": grid.compute_dict(),
            "PFL_LS": grid.compute_dict(),
            "PFI_LS": grid.compute_dict(),
            "PFL_AN": grid.compute_dict(),
            "PFI_AN": grid.compute_dict(),
            "LS_PRCP": grid.compute_dict(),
            "LS_SNR": grid.compute_dict(),
            "ICE": grid.compute_dict(),
            "FRZR": grid.compute_dict(),
            "DQVDT_micro": grid.compute_dict(),
            "DQLDT_micro": grid.compute_dict(),
            "DQIDT_micro": grid.compute_dict(),
            "DQADT_micro": grid.compute_dict(),
            "DQRDT_micro": grid.compute_dict(),
            "DQSDT_micro": grid.compute_dict(),
            "DQGDT_micro": grid.compute_dict(),
            "DUDT_micro": grid.compute_dict(),
            "DVDT_micro": grid.compute_dict(),
            "DTDT_micro": grid.compute_dict(),
        }

        self.out_vars = self.in_vars["data_vars"].copy()

    def make_ijk_quantity(self, data, interface: bool = False) -> Quantity:
        if interface == True:
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
        self.constants = data_loader.load("GFDL_1M-constants")

    def compute(self, inputs):
        # Initalize GFDL_1M configuration
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

        # Initalize saturation tables
        saturation_tables = SaturationVaporPressureTable(self.stencil_factory.backend)

        # Initalize extra quantities
        outputs = Outputs.make(self.quantity_factory)
        temporaries = Temporaries.make(self.quantity_factory)
        masks = Masks.make(self.quantity_factory)

        # Make input quantities
        t = self.make_ijk_quantity(inputs.pop("T"))
        u = self.make_ijk_quantity(inputs.pop("U"))
        v = self.make_ijk_quantity(inputs.pop("V"))
        temporaries.p_mb = self.make_ijk_quantity(inputs.pop("PLmb"))
        temporaries.mass = self.make_ijk_quantity(inputs.pop("MASS"))
        ice_concentration = self.make_ijk_quantity(inputs.pop("NACTI"))
        liquid_concentration = self.make_ijk_quantity(inputs.pop("NACTL"))
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
        outputs.radiation_cloud_fraction = self.make_ijk_quantity(inputs.pop("RAD_CF"))
        outputs.radiation_ice = self.make_ijk_quantity(inputs.pop("RAD_QI"))
        outputs.radiation_liquid = self.make_ijk_quantity(inputs.pop("RAD_QL"))
        outputs.radiation_vapor = self.make_ijk_quantity(inputs.pop("RAD_QV"))
        outputs.radiation_rain = self.make_ijk_quantity(inputs.pop("RAD_QR"))
        outputs.radiation_snow = self.make_ijk_quantity(inputs.pop("RAD_QS"))
        outputs.radiation_graupel = self.make_ijk_quantity(inputs.pop("RAD_QG"))
        outputs.liquid_radius = self.make_ijk_quantity(inputs.pop("CLDREFFL"))
        outputs.ice_radius = self.make_ijk_quantity(inputs.pop("CLDREFFI"))
        outputs.du_dt_micro = self.make_ijk_quantity(inputs.pop("DUDT_micro"))
        outputs.dv_dt_micro = self.make_ijk_quantity(inputs.pop("DVDT_micro"))
        outputs.dt_dt_micro = self.make_ijk_quantity(inputs.pop("DTDT_micro"))
        outputs.dvapor_dt_micro = self.make_ijk_quantity(inputs.pop("DQVDT_micro"))
        outputs.dliquid_dt_micro = self.make_ijk_quantity(inputs.pop("DQLDT_micro"))
        outputs.dice_dt_micro = self.make_ijk_quantity(inputs.pop("DQIDT_micro"))
        outputs.dcloud_fraction_dt_micro = self.make_ijk_quantity(inputs.pop("DQADT_micro"))
        outputs.drain_dt_micro = self.make_ijk_quantity(inputs.pop("DQRDT_micro"))
        outputs.dsnow_dt_micro = self.make_ijk_quantity(inputs.pop("DQSDT_micro"))
        outputs.dgraupel_dt_micro = self.make_ijk_quantity(inputs.pop("DQGDT_micro"))

        # Spoof driver output with Fortran data
        driver = MicrophysicsDriver(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            GFDL_1M_config=GFDL_1M_config,
        )
        driver.outputs.rain = self.make_ij_quantity(inputs.pop("PRCP_RAIN"))
        driver.outputs.snow = self.make_ij_quantity(inputs.pop("PRCP_SNOW"))
        driver.outputs.ice = self.make_ij_quantity(inputs.pop("PRCP_ICE"))
        driver.outputs.graupel = self.make_ij_quantity(inputs.pop("PRCP_GRAUPEL"))
        driver.outputs.m2_rain = self.make_ijk_quantity(inputs.pop("PFL_LS"), interface=True)
        driver.outputs.m2_sol = self.make_ijk_quantity(inputs.pop("PFI_LS"), interface=True)
        driver.outputs.revap = self.make_ijk_quantity(inputs.pop("REV_LS"))
        driver.outputs.isubl = self.make_ijk_quantity(inputs.pop("RSU_LS"))

        # Construct stencils
        self.update_tendencies = self.stencil_factory.from_dims_halo(
            func=update_tendencies,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "DT_MOIST": GFDL_1M_config.DT_MOIST,
            },
        )

        finalize = Finalize(
            stencil_factory=self.stencil_factory,
            GFDL_1M_config=GFDL_1M_config,
            saturation_tables=saturation_tables,
            update_tendencies=self.update_tendencies,
        )

        # Call tested code
        finalize(
            t=t,
            u=u,
            v=v,
            ice_concentration=ice_concentration,
            liquid_concentration=liquid_concentration,
            mixing_ratios=mixing_ratios,
            cloud_fractions=cloud_fractions,
            masks=masks,
            outputs=outputs,
            temporaries=temporaries,
            driver=driver,
            saturation_tables=saturation_tables,
        )

        return {
            "T": t.field,
            "U": u.field,
            "V": v.field,
            "PLmb": temporaries.p_mb.field,
            "MASS": temporaries.mass.field,
            "RAD_CF": outputs.radiation_cloud_fraction.field,
            "RAD_QI": outputs.radiation_ice.field,
            "RAD_QL": outputs.radiation_liquid.field,
            "RAD_QV": outputs.radiation_vapor.field,
            "RAD_QR": outputs.radiation_rain.field,
            "RAD_QS": outputs.radiation_snow.field,
            "RAD_QG": outputs.radiation_graupel.field,
            "CLCN": cloud_fractions.convective.field,
            "CLLS": cloud_fractions.large_scale.field,
            "QLCN": mixing_ratios.convective_liquid.field,
            "QLLS": mixing_ratios.large_scale_liquid.field,
            "NACTL": liquid_concentration.field,
            "CLDREFFL": outputs.liquid_radius.field,
            "QICN": mixing_ratios.convective_ice.field,
            "QILS": mixing_ratios.large_scale_ice.field,
            "CLDREFFI": outputs.ice_radius.field,
            "NACTI": ice_concentration.field,
            "Q": mixing_ratios.vapor.field,
            "QRAIN": mixing_ratios.rain.field,
            "QSNOW": mixing_ratios.snow.field,
            "QGRAUPEL": mixing_ratios.graupel.field,
            "REV_LS": outputs.large_scale_nonanvil_precipitation_evaporation.field,
            "RSU_LS": outputs.large_scale_nonanvil_precipitation_sublimation.field,
            "PRCP_RAIN": outputs.precipitated_rain.field,
            "PRCP_SNOW": outputs.precipitated_snow.field,
            "PRCP_ICE": outputs.precipitated_ice.field,
            "PRCP_GRAUPEL": outputs.precipitated_graupel.field,
            "PFL_LS": outputs.large_scale_nonanvil_liquid_flux.field,
            "PFI_LS": outputs.large_scale_nonanvil_ice_flux.field,
            "PFL_AN": outputs.anvil_liquid_flux.field,
            "PFI_AN": outputs.anvil_ice_flux.field,
            "LS_PRCP": outputs.large_scale_precip.field,
            "LS_SNR": outputs.large_scale_snow.field,
            "ICE": outputs.icefall.field,
            "FRZR": outputs.freezing_rainfall.field,
            "DQVDT_micro": outputs.dvapor_dt_micro.field,
            "DQLDT_micro": outputs.dliquid_dt_micro.field,
            "DQIDT_micro": outputs.dice_dt_micro.field,
            "DQADT_micro": outputs.dcloud_fraction_dt_micro.field,
            "DQRDT_micro": outputs.drain_dt_micro.field,
            "DQSDT_micro": outputs.dsnow_dt_micro.field,
            "DQGDT_micro": outputs.dgraupel_dt_micro.field,
            "DUDT_micro": outputs.du_dt_micro.field,
            "DVDT_micro": outputs.dv_dt_micro.field,
            "DTDT_micro": outputs.dt_dt_micro.field,
        }
