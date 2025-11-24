from f90nml import Namelist

from ndsl import StencilFactory
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.GFDL_1M.driver.config_constants import ConfigConstants
from pyMoist.GFDL_1M.driver.finish.main import Finish


class TranslateGFDL_1M_driver_finish(TranslateFortranData2Py):
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
            "qv0": grid.compute_dict() | {"serialname": "qv0_driver_finish"},
            "ql0": grid.compute_dict() | {"serialname": "ql0_driver_finish"},
            "qr0": grid.compute_dict() | {"serialname": "qr0_driver_finish"},
            "qi0": grid.compute_dict() | {"serialname": "qi0_driver_finish"},
            "qs0": grid.compute_dict() | {"serialname": "qs0_driver_finish"},
            "qg0": grid.compute_dict() | {"serialname": "qg0_driver_finish"},
            "qa0": grid.compute_dict() | {"serialname": "qa0_driver_finish"},
            "qv1": grid.compute_dict() | {"serialname": "qv1_driver_finish"},
            "ql1": grid.compute_dict() | {"serialname": "ql1_driver_finish"},
            "qr1": grid.compute_dict() | {"serialname": "qr1_driver_finish"},
            "qi1": grid.compute_dict() | {"serialname": "qi1_driver_finish"},
            "qs1": grid.compute_dict() | {"serialname": "qs1_driver_finish"},
            "qg1": grid.compute_dict() | {"serialname": "qg1_driver_finish"},
            "qv_dt": grid.compute_dict() | {"serialname": "qv_dt_driver_finish"},
            "ql_dt": grid.compute_dict() | {"serialname": "ql_dt_driver_finish"},
            "qr_dt": grid.compute_dict() | {"serialname": "qr_dt_driver_finish"},
            "qi_dt": grid.compute_dict() | {"serialname": "qi_dt_driver_finish"},
            "qs_dt": grid.compute_dict() | {"serialname": "qs_dt_driver_finish"},
            "qg_dt": grid.compute_dict() | {"serialname": "qg_dt_driver_finish"},
            "qa_dt": grid.compute_dict() | {"serialname": "qa_dt_driver_finish"},
            "t": grid.compute_dict() | {"serialname": "t_driver_finish"},
            "t1": grid.compute_dict() | {"serialname": "t1_driver_finish"},
            "t_dt": grid.compute_dict() | {"serialname": "t_dt_driver_finish"},
            "w": grid.compute_dict() | {"serialname": "w_driver_finish"},
            "w1": grid.compute_dict() | {"serialname": "w1_driver_finish"},
            "u": grid.compute_dict() | {"serialname": "u_driver_finish"},
            "u1": grid.compute_dict() | {"serialname": "u1_driver_finish"},
            "udt": grid.compute_dict() | {"serialname": "udt_driver_finish"},
            "v": grid.compute_dict() | {"serialname": "v_driver_finish"},
            "v1": grid.compute_dict() | {"serialname": "v1_driver_finish"},
            "vdt": grid.compute_dict() | {"serialname": "vdt_driver_finish"},
            "dp": grid.compute_dict() | {"serialname": "dp_driver_finish"},
            "dp1": grid.compute_dict() | {"serialname": "dp1_driver_finish"},
            "m1": grid.compute_dict() | {"serialname": "m1_driver_finish"},
            "rain": grid.compute_dict() | {"serialname": "rain_driver_finish"},
            "snow": grid.compute_dict() | {"serialname": "snow_driver_finish"},
            "ice": grid.compute_dict() | {"serialname": "ice_driver_finish"},
            "graupel": grid.compute_dict() | {"serialname": "graupel_driver_finish"},
        }

        self.out_vars = self.in_vars["data_vars"].copy()
        del (
            self.out_vars["rain"],
            self.out_vars["snow"],
            self.out_vars["ice"],
            self.out_vars["graupel"],
        )

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GFDL_1M-constants")

    def compute_from_storage(self, inputs):
        self.GFDL_1M_config = GFDL1MConfig(
            LPHYS_HYDROSTATIC=bool(self.constants["LPHYS_HYDROSTATIC"]),
            LHYDROSTATIC=bool(self.constants["LHYDROSTATIC"]),
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
            LMELTFRZ=bool(self.constants["LMELTFRZ"]),
            USE_BERGERON=bool(self.constants["USE_BERGERON"]),
            TURNRHCRIT_PARAM=self.constants["TURNRHCRIT_PARAM"],
            PDFSHAPE=self.constants["PDFSHAPE"],
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

        config_dependent_constants = ConfigConstants.make(self.GFDL_1M_config)

        # Initalize object to be tested
        finish = Finish(self.stencil_factory, self.GFDL_1M_config, config_dependent_constants)

        finish(**inputs)

        return inputs
