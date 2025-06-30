from ndsl import Namelist, Quantity, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.radiation_coupling import GFDL1MRadiationCoupling
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable


class TranslateGFDL_1M_radiation_coupling(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory
        # Inputs
        self.in_vars["data_vars"] = {
            "Q": grid.compute_dict(),
            "T": grid.compute_dict(),
            "QLLS": grid.compute_dict(),
            "QILS": grid.compute_dict(),
            "CLLS": grid.compute_dict(),
            "QLCN": grid.compute_dict(),
            "QICN": grid.compute_dict(),
            "CLCN": grid.compute_dict(),
            "PLmb": grid.compute_dict(),
            "QRAIN": grid.compute_dict(),
            "QSNOW": grid.compute_dict(),
            "QGRAUPEL": grid.compute_dict(),
            "NACTL": grid.compute_dict(),
            "NACTI": grid.compute_dict(),
            "RAD_QV": grid.compute_dict(),
            "RAD_QL": grid.compute_dict(),
            "RAD_QI": grid.compute_dict(),
            "RAD_QR": grid.compute_dict(),
            "RAD_QS": grid.compute_dict(),
            "RAD_QG": grid.compute_dict(),
            "RAD_CF": grid.compute_dict(),
            "CLDREFFI": grid.compute_dict(),
            "CLDREFFL": grid.compute_dict(),
            "RHX": grid.compute_dict(),
        }

        # Outputs
        self.out_vars = {
            "RAD_QV": grid.compute_dict(),
            "RAD_QL": grid.compute_dict(),
            "RAD_QI": grid.compute_dict(),
            "RAD_QR": grid.compute_dict(),
            "RAD_QS": grid.compute_dict(),
            "RAD_QG": grid.compute_dict(),
            "RAD_CF": grid.compute_dict(),
            "CLDREFFI": grid.compute_dict(),
            "CLDREFFL": grid.compute_dict(),
            "RHX": grid.compute_dict(),
        }

    def make_ijk_quantity(self, data, interface: bool = False) -> Quantity:
        if interface is True:
            quantity = self.quantity_factory.empty([X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a")
            quantity.view[:, :, :] = quantity.np.asarray(data[:, :, :])
            return quantity
        else:
            quantity = self.quantity_factory.empty([X_DIM, Y_DIM, Z_DIM], "n/a")
            quantity.view[:, :, :] = quantity.np.asarray(data[:, :, :])
            return quantity

    def make_ij_quantity(self, data, interface: bool = False) -> Quantity:
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

        vapor = self.make_ijk_quantity(inputs.pop("Q"))
        temperature = self.make_ijk_quantity(inputs.pop("T"))
        large_scale_liquid = self.make_ijk_quantity(inputs.pop("QLLS"))
        large_scale_ice = self.make_ijk_quantity(inputs.pop("QILS"))
        large_scale_cloud_fraction = self.make_ijk_quantity(inputs.pop("CLLS"))
        convective_liquid = self.make_ijk_quantity(inputs.pop("QLCN"))
        convective_ice = self.make_ijk_quantity(inputs.pop("QICN"))
        convective_cloud_fraction = self.make_ijk_quantity(inputs.pop("CLCN"))
        pressure = self.make_ijk_quantity(inputs.pop("PLmb"))
        rain = self.make_ijk_quantity(inputs.pop("QRAIN"))
        snow = self.make_ijk_quantity(inputs.pop("QSNOW"))
        graupel = self.make_ijk_quantity(inputs.pop("QGRAUPEL"))
        liquid_concentration = self.make_ijk_quantity(inputs.pop("NACTL"))
        ice_concentration = self.make_ijk_quantity(inputs.pop("NACTI"))
        radiation_vapor = self.make_ijk_quantity(inputs.pop("RAD_QV"))
        radiation_liquid = self.make_ijk_quantity(inputs.pop("RAD_QL"))
        radiation_ice = self.make_ijk_quantity(inputs.pop("RAD_QI"))
        radiation_rain = self.make_ijk_quantity(inputs.pop("RAD_QR"))
        radiation_snow = self.make_ijk_quantity(inputs.pop("RAD_QS"))
        radiation_graupel = self.make_ijk_quantity(inputs.pop("RAD_QG"))
        radiation_cloud_fraction = self.make_ijk_quantity(inputs.pop("RAD_CF"))
        liquid_radius = self.make_ijk_quantity(inputs.pop("CLDREFFL"))
        ice_radius = self.make_ijk_quantity(inputs.pop("CLDREFFI"))
        humidity = self.make_ijk_quantity(inputs.pop("RHX"))
        radiation_coupling = GFDL1MRadiationCoupling(  # type: ignore
            stencil_factory=self.stencil_factory,
            DO_QA=GFDL_1M_config.DO_QA,
            FAC_RL=GFDL_1M_config.FAC_RL,
            MIN_RL=GFDL_1M_config.MIN_RL,
            MAX_RL=GFDL_1M_config.MAX_RL,
            FAC_RI=GFDL_1M_config.FAC_RI,
            MIN_RI=GFDL_1M_config.MIN_RI,
            MAX_RI=GFDL_1M_config.MAX_RI,
        )
        radiation_coupling(
            vapor,
            temperature,
            large_scale_liquid,
            large_scale_ice,
            large_scale_cloud_fraction,
            convective_liquid,
            convective_ice,
            convective_cloud_fraction,
            pressure,
            rain,
            snow,
            graupel,
            liquid_concentration,
            ice_concentration,
            radiation_vapor,
            radiation_liquid,
            radiation_ice,
            radiation_rain,
            radiation_snow,
            radiation_graupel,
            radiation_cloud_fraction,
            liquid_radius,
            ice_radius,
            humidity,
            saturation_tables.ese,
            saturation_tables.esx,
        )
        return {
            "Q": vapor.field,
            "T": temperature.field,
            "QLLS": large_scale_liquid.field,
            "QILS": large_scale_ice.field,
            "CLLS": convective_cloud_fraction.field,
            "QLCN": convective_liquid.field,
            "QICN": convective_ice.field,
            "CLCN": convective_cloud_fraction.field,
            "PLmb": pressure.field,
            "QRAIN": rain.field,
            "QSNOW": snow.field,
            "QGRAUPEL": graupel.field,
            "NACTL": liquid_concentration.field,
            "NACTI": ice_concentration.field,
            "RAD_QV": radiation_vapor.field,
            "RAD_QL": radiation_liquid.field,
            "RAD_QI": radiation_ice.field,
            "RAD_QR": radiation_rain.field,
            "RAD_QS": radiation_snow.field,
            "RAD_QG": radiation_graupel.field,
            "RAD_CF": radiation_cloud_fraction.field,
            "CLDREFFL": liquid_radius.field,
            "CLDREFFI": ice_radius.field,
            "RHX": humidity.field,
        }
