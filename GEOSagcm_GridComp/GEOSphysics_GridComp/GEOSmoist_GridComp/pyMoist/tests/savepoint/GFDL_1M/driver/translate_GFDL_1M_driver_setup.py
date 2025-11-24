from f90nml import Namelist

from ndsl import StencilFactory
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.GFDL_1M.driver.config_constants import ConfigConstants
from pyMoist.GFDL_1M.driver.setup.main import Setup


class TranslateGFDL_1M_driver_setup(TranslateFortranData2Py):
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
            "t": grid.compute_dict() | {"serialname": "t_driver_setup"},
            "dp": grid.compute_dict() | {"serialname": "dp_driver_setup"},
            "rhcrit3d": grid.compute_dict() | {"serialname": "rhcrit3d_driver_setup"},
            "qv": grid.compute_dict() | {"serialname": "qv_driver_setup"},
            "ql": grid.compute_dict() | {"serialname": "ql_driver_setup"},
            "qi": grid.compute_dict() | {"serialname": "qi_driver_setup"},
            "qr": grid.compute_dict() | {"serialname": "qr_driver_setup"},
            "qs": grid.compute_dict() | {"serialname": "qs_driver_setup"},
            "qg": grid.compute_dict() | {"serialname": "qg_driver_setup"},
            "qa": grid.compute_dict() | {"serialname": "qa_driver_setup"},
            "ice_concentration": grid.compute_dict() | {"serialname": "ice_concentration_driver_setup"},
            "liquid_concentration": grid.compute_dict() | {"serialname": "liquid_concentration_driver_setup"},
            "qv0": grid.compute_dict() | {"serialname": "qv0_driver_setup"},
            "ql0": grid.compute_dict() | {"serialname": "ql0_driver_setup"},
            "qr0": grid.compute_dict() | {"serialname": "qr0_driver_setup"},
            "qi0": grid.compute_dict() | {"serialname": "qi0_driver_setup"},
            "qs0": grid.compute_dict() | {"serialname": "qs0_driver_setup"},
            "qg0": grid.compute_dict() | {"serialname": "qg0_driver_setup"},
            "qa0": grid.compute_dict() | {"serialname": "qa0_driver_setup"},
            "qv1": grid.compute_dict() | {"serialname": "qv1_driver_setup"},
            "ql1": grid.compute_dict() | {"serialname": "ql1_driver_setup"},
            "qr1": grid.compute_dict() | {"serialname": "qr1_driver_setup"},
            "qi1": grid.compute_dict() | {"serialname": "qi1_driver_setup"},
            "qs1": grid.compute_dict() | {"serialname": "qs1_driver_setup"},
            "qg1": grid.compute_dict() | {"serialname": "qg1_driver_setup"},
            "qa1": grid.compute_dict() | {"serialname": "qa1_driver_setup"},
            "dz": grid.compute_dict() | {"serialname": "dz_driver_setup"},
            "u": grid.compute_dict() | {"serialname": "u_driver_setup"},
            "v": grid.compute_dict() | {"serialname": "v_driver_setup"},
            "w": grid.compute_dict() | {"serialname": "w_driver_setup"},
            "area": grid.compute_dict() | {"serialname": "area_driver_setup"},
            "t1": grid.compute_dict() | {"serialname": "t1_driver_setup"},
            "dp1": grid.compute_dict() | {"serialname": "dp1_driver_setup"},
            "omq": grid.compute_dict() | {"serialname": "omq_driver_setup"},
            "den": grid.compute_dict() | {"serialname": "den_driver_setup"},
            "p_dry": grid.compute_dict() | {"serialname": "p_dry_driver_setup"},
            "m1": grid.compute_dict() | {"serialname": "m1_driver_setup"},
            "u1": grid.compute_dict() | {"serialname": "u1_driver_setup"},
            "v1": grid.compute_dict() | {"serialname": "v1_driver_setup"},
            "w1": grid.compute_dict() | {"serialname": "w1_driver_setup"},
            "onemsig": grid.compute_dict() | {"serialname": "onemsig_driver_setup"},
            "ccn": grid.compute_dict() | {"serialname": "ccn_driver_setup"},
            "c_praut": grid.compute_dict() | {"serialname": "c_praut_driver_setup"},
            "rh_limited": grid.compute_dict() | {"serialname": "rh_limited_driver_setup"},
            "rain": grid.compute_dict() | {"serialname": "rain_driver_setup"},
            "snow": grid.compute_dict() | {"serialname": "snow_driver_setup"},
            "graupel": grid.compute_dict() | {"serialname": "graupel_driver_setup"},
            "ice": grid.compute_dict() | {"serialname": "ice_driver_setup"},
            "m2_rain": grid.compute_dict() | {"serialname": "m2_rain_driver_setup"},
            "m2_sol": grid.compute_dict() | {"serialname": "m2_sol_driver_setup"},
            "revap": grid.compute_dict() | {"serialname": "revap_driver_setup"},
            "isubl": grid.compute_dict() | {"serialname": "isubl_driver_setup"},
        }

        self.out_vars = self.in_vars["data_vars"].copy()

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
        setup = Setup(self.stencil_factory, self.GFDL_1M_config, config_dependent_constants)

        setup(**inputs)

        return inputs
