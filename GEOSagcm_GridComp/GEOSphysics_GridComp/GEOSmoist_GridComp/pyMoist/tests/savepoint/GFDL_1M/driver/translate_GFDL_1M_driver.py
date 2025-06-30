from ndsl import Namelist, StencilFactory
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.GFDL_1M.driver.driver import MicrophysicsDriver


class TranslateGFDL_1M_driver(TranslateFortranData2Py):
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
            "t": grid.compute_dict() | {"serialname": "T"},
            "w": grid.compute_dict() | {"serialname": "W"},
            "u": grid.compute_dict() | {"serialname": "U"},
            "v": grid.compute_dict() | {"serialname": "V"},
            "dz": grid.compute_dict() | {"serialname": "DZ"},
            "dp": grid.compute_dict() | {"serialname": "DP"},
            "area": grid.compute_dict() | {"serialname": "AREA"},
            "land_fraction": grid.compute_dict() | {"serialname": "FRLAND"},
            "convection_fraction": grid.compute_dict() | {"serialname": "CNV_FRC"},
            "surface_type": grid.compute_dict() | {"serialname": "SRF_TYPE"},
            "estimated_inversion_strength": grid.compute_dict() | {"serialname": "EIS"},
            "rh_crit": grid.compute_dict() | {"serialname": "RHCRIT3D"},
            "vapor": grid.compute_dict() | {"serialname": "RAD_QV"},
            "liquid": grid.compute_dict() | {"serialname": "RAD_QL"},
            "rain": grid.compute_dict() | {"serialname": "RAD_QR"},
            "ice": grid.compute_dict() | {"serialname": "RAD_QI"},
            "snow": grid.compute_dict() | {"serialname": "RAD_QS"},
            "graupel": grid.compute_dict() | {"serialname": "RAD_QG"},
            "cloud_fraction": grid.compute_dict() | {"serialname": "RAD_CF"},
            "ice_concentration": grid.compute_dict() | {"serialname": "NACTI"},
            "liquid_concentration": grid.compute_dict() | {"serialname": "NACTL"},
            "dvapor_dt": grid.compute_dict() | {"serialname": "DQVDTmic"},
            "dliquid_dt": grid.compute_dict() | {"serialname": "DQLDTmic"},
            "drain_dt": grid.compute_dict() | {"serialname": "DQRDTmic"},
            "dice_dt": grid.compute_dict() | {"serialname": "DQIDTmic"},
            "dsnow_dt": grid.compute_dict() | {"serialname": "DQSDTmic"},
            "dgraupel_dt": grid.compute_dict() | {"serialname": "DQGDTmic"},
            "dcloud_fraction_dt": grid.compute_dict() | {"serialname": "DQADTmic"},
            "dt_dt": grid.compute_dict() | {"serialname": "DTDTmic"},
            "du_dt": grid.compute_dict() | {"serialname": "DUDTmic"},
            "dv_dt": grid.compute_dict() | {"serialname": "DVDTmic"},
        }

        self.out_vars = self.in_vars["data_vars"].copy()
        self.out_vars.update(
            {
                "REV_LS": {},
                "RSU_LS": {},
                "PRCP_RAIN": {},
                "PRCP_SNOW": {},
                "PRCP_ICE": {},
                "PRCP_GRAUPEL": {},
                "PFL_LS_driver": {},
                "PFI_LS_driver": {},
            }
        )
        del (
            self.out_vars["t"],
            self.out_vars["w"],
            self.out_vars["u"],
            self.out_vars["v"],
            self.out_vars["dz"],
            self.out_vars["dp"],
            self.out_vars["area"],
            self.out_vars["land_fraction"],
            self.out_vars["convection_fraction"],
            self.out_vars["surface_type"],
            self.out_vars["estimated_inversion_strength"],
            self.out_vars["rh_crit"],
            self.out_vars["vapor"],
            self.out_vars["liquid"],
            self.out_vars["rain"],
            self.out_vars["graupel"],
            self.out_vars["cloud_fraction"],
            self.out_vars["ice_concentration"],
            self.out_vars["liquid_concentration"],
        )

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GFDL_1M-constants")

    def compute_from_storage(self, inputs):
        self.GFDL_1M_config = GFDL1MConfig(
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

        # Initalize object to be tested
        self.driver = MicrophysicsDriver(
            self.stencil_factory,
            self.quantity_factory,
            self.GFDL_1M_config,
        )

        self.driver(
            **inputs,
        )

        extra_outputs = {
            "REV_LS": self.driver.outputs.revap.view[:],
            "RSU_LS": self.driver.outputs.isubl.view[:],
            "PRCP_RAIN": self.driver.outputs.rain.view[:],
            "PRCP_SNOW": self.driver.outputs.snow.view[:],
            "PRCP_ICE": self.driver.outputs.ice.view[:],
            "PRCP_GRAUPEL": self.driver.outputs.graupel.view[:],
            "PFL_LS": self.driver.outputs.m2_rain.view[:],
            "PFI_LS": self.driver.outputs.m2_sol.view[:],
        }
        # breakpoint()
        inputs.update(
            {
                "REV_LS": self.driver.outputs.revap.view[:],
                "RSU_LS": self.driver.outputs.isubl.view[:],
                "PRCP_RAIN": self.driver.outputs.rain.view[:],
                "PRCP_SNOW": self.driver.outputs.snow.view[:],
                "PRCP_ICE": self.driver.outputs.ice.view[:],
                "PRCP_GRAUPEL": self.driver.outputs.graupel.view[:],
                "PFL_LS_driver": self.driver.outputs.m2_rain.view[:],
                "PFI_LS_driver": self.driver.outputs.m2_sol.view[:],
            }
        )

        return inputs
