from ndsl import Namelist, StencilFactory
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.GFDL_1M.driver.config_constants import ConfigConstants
from pyMoist.GFDL_1M.driver.masks import Masks
from pyMoist.GFDL_1M.driver.outputs import Outputs
from pyMoist.GFDL_1M.driver.sat_tables import get_tables
from pyMoist.GFDL_1M.driver.temporaries import Temporaries
from pyMoist.GFDL_1M.driver.warm_rain.main import WarmRain


class TranslateGFDL_1M_warm_rain(TranslateFortranData2Py):
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
            "dz1": grid.compute_dict() | {"serialname": "dz1_warm_rain"},
            "dp1": grid.compute_dict() | {"serialname": "dp1_warm_rain"},
            "t1": grid.compute_dict() | {"serialname": "t1_warm_rain"},
            "qv1": grid.compute_dict() | {"serialname": "qv_warm_rain"},
            "ql1": grid.compute_dict() | {"serialname": "ql_warm_rain"},
            "qr1": grid.compute_dict() | {"serialname": "qr_warm_rain"},
            "qi1": grid.compute_dict() | {"serialname": "qi_warm_rain"},
            "qs1": grid.compute_dict() | {"serialname": "qs_warm_rain"},
            "qg1": grid.compute_dict() | {"serialname": "qg_warm_rain"},
            "qa1": grid.compute_dict() | {"serialname": "qa_warm_rain"},
            "ccn": grid.compute_dict() | {"serialname": "ccn_warm_rain"},
            "den": grid.compute_dict() | {"serialname": "den_warm_rain"},
            "denfac": grid.compute_dict() | {"serialname": "denfac_warm_rain"},
            "c_praut": grid.compute_dict() | {"serialname": "c_praut_warm_rain"},
            "vtr": grid.compute_dict() | {"serialname": "vtr_warm_rain"},
            "evap1": grid.compute_dict() | {"serialname": "evap1_warm_rain"},
            "m1_rain": grid.compute_dict() | {"serialname": "m1_rain_warm_rain"},
            "w1": grid.compute_dict() | {"serialname": "w1_warm_rain"},
            "rh_limited": grid.compute_dict() | {"serialname": "rh_limited_warm_rain"},
            "eis": grid.compute_dict() | {"serialname": "eis_warm_rain"},
            "rain": grid.compute_dict() | {"serialname": "rain_warm_rain"},
            "rain1": grid.compute_dict() | {"serialname": "r1_warm_rain"},
            "m2_sol": grid.compute_dict() | {"serialname": "m2_sol_warm_rain"},
            "m2_rain": grid.compute_dict() | {"serialname": "m2_rain_warm_rain"},
            "revap": grid.compute_dict() | {"serialname": "revap_warm_rain"},
            "m1": grid.compute_dict() | {"serialname": "m1_warm_rain"},
            "m1_sol": grid.compute_dict() | {"serialname": "m1_sol_warm_rain"},
            "onemsig": grid.compute_dict() | {"serialname": "onemsig_warm_rain"},
        }

        self.out_vars = self.in_vars["data_vars"].copy()
        del (
            self.out_vars["dz1"],
            self.out_vars["dp1"],
            self.out_vars["den"],
            self.out_vars["denfac"],
            self.out_vars["c_praut"],
            self.out_vars["ccn"],
            self.out_vars["onemsig"],
            self.out_vars["rh_limited"],
            self.out_vars["eis"],
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

        self.config_dependent_constants = ConfigConstants.make(self.GFDL_1M_config)

        # Initalize saturation tables
        self.sat_tables = get_tables(self.stencil_factory.backend)

        # Initalize extra quantities
        temporaries = Temporaries.make(self.quantity_factory)
        outputs = Outputs.make(self.quantity_factory)
        masks = Masks.make(self.quantity_factory)

        # Initalize object to be tested
        self.warm_rain = WarmRain(
            self.stencil_factory,
            self.quantity_factory,
            self.GFDL_1M_config,
            self.config_dependent_constants,
        )
        self.warm_rain(
            ze=temporaries.ze,
            zt=temporaries.zt,
            precip_fall=masks.precip_fall,
            table1=self.sat_tables.table1,
            table2=self.sat_tables.table2,
            table3=self.sat_tables.table3,
            table4=self.sat_tables.table4,
            des1=self.sat_tables.des1,
            des2=self.sat_tables.des2,
            des3=self.sat_tables.des3,
            des4=self.sat_tables.des4,
            **inputs,
        )

        return inputs
