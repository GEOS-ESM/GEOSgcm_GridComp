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
            "qv": grid.compute_dict() | {"serialname": "RAD_QV"},
            "ql": grid.compute_dict() | {"serialname": "RAD_QL"},
            "qr": grid.compute_dict() | {"serialname": "RAD_QR"},
            "qi": grid.compute_dict() | {"serialname": "RAD_QI"},
            "qs": grid.compute_dict() | {"serialname": "RAD_QS"},
            "qg": grid.compute_dict() | {"serialname": "RAD_QG"},
            "qa": grid.compute_dict() | {"serialname": "RAD_CF"},
            "qn": grid.compute_dict() | {"serialname": "NACTLI"},
            "qv_dt": grid.compute_dict() | {"serialname": "DQVDTmic"},
            "ql_dt": grid.compute_dict() | {"serialname": "DQLDTmic"},
            "qr_dt": grid.compute_dict() | {"serialname": "DQRDTmic"},
            "qi_dt": grid.compute_dict() | {"serialname": "DQIDTmic"},
            "qs_dt": grid.compute_dict() | {"serialname": "DQSDTmic"},
            "qg_dt": grid.compute_dict() | {"serialname": "DQGDTmic"},
            "qa_dt": grid.compute_dict() | {"serialname": "DQADTmic"},
            "t_dt": grid.compute_dict() | {"serialname": "DTDTmic"},
            "t": grid.compute_dict() | {"serialname": "T"},
            "w": grid.compute_dict() | {"serialname": "W"},
            "u": grid.compute_dict() | {"serialname": "U"},
            "v": grid.compute_dict() | {"serialname": "V"},
            "u_dt": grid.compute_dict() | {"serialname": "DUDTmic"},
            "v_dt": grid.compute_dict() | {"serialname": "DVDTmic"},
            "dz": grid.compute_dict() | {"serialname": "DZ"},
            "dp": grid.compute_dict() | {"serialname": "DP"},
            "area": grid.compute_dict() | {"serialname": "AREA"},
            "fr_land": grid.compute_dict() | {"serialname": "FRLAND"},
            "cnv_frc": grid.compute_dict() | {"serialname": "CNV_FRC"},
            "srf_type": grid.compute_dict() | {"serialname": "SRF_TYPE"},
            "eis": grid.compute_dict() | {"serialname": "EIS"},
            "rhcrit3d": grid.compute_dict() | {"serialname": "RHCRIT3D"},
        }
        self.in_vars["parameters"] = [
            "ANV_ICEFALL",
            "LS_ICEFALL",
        ]

        self.out_vars = self.in_vars["data_vars"].copy()
        self.out_vars.update(
            {
                "REV_LS": {},
                "RSU_LS": {},
                "PRCP_RAIN": {},
                "PRCP_SNOW": {},
                "PRCP_ICE": {},
                "PRCP_GRAUPEL": {},
                "PFL_LS": {},
                "PFI_LS": {},
            }
        )
        del (
            self.out_vars["qv"],
            self.out_vars["ql"],
            self.out_vars["qr"],
            self.out_vars["qg"],
            self.out_vars["qa"],
            self.out_vars["qn"],
            self.out_vars["t"],
            self.out_vars["w"],
            self.out_vars["u"],
            self.out_vars["v"],
            self.out_vars["dz"],
            self.out_vars["dp"],
            self.out_vars["area"],
            self.out_vars["fr_land"],
            self.out_vars["cnv_frc"],
            self.out_vars["srf_type"],
            self.out_vars["eis"],
            self.out_vars["rhcrit3d"],
        )

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GFDL_1M_driver-constants")

    def compute_from_storage(self, inputs):
        self.GFDL_1M_config = GFDL1MConfig(
            bool(self.constants["PHYS_HYDROSTATIC"]),
            bool(self.constants["HYDROSTATIC"]),
            self.constants["DT_MOIST"],
            self.constants["MP_TIME"],
            self.constants["T_MIN"],
            self.constants["T_SUB"],
            self.constants["TAU_R2G"],
            self.constants["TAU_SMLT"],
            self.constants["TAU_G2R"],
            self.constants["DW_LAND"],
            self.constants["DW_OCEAN"],
            self.constants["VI_FAC"],
            self.constants["VR_FAC"],
            self.constants["VS_FAC"],
            self.constants["VG_FAC"],
            self.constants["QL_MLT"],
            bool(self.constants["DO_QA"]),
            bool(self.constants["FIX_NEGATIVE"]),
            self.constants["VI_MAX"],
            self.constants["VS_MAX"],
            self.constants["VG_MAX"],
            self.constants["VR_MAX"],
            self.constants["QS_MLT"],
            self.constants["QS0_CRT"],
            self.constants["QI_GEN"],
            self.constants["QL0_MAX"],
            self.constants["QI0_MAX"],
            self.constants["QI0_CRT"],
            self.constants["QR0_CRT"],
            bool(self.constants["FAST_SAT_ADJ"]),
            self.constants["RH_INC"],
            self.constants["RH_INS"],
            self.constants["RH_INR"],
            bool(self.constants["CONST_VI"]),
            bool(self.constants["CONST_VS"]),
            bool(self.constants["CONST_VG"]),
            bool(self.constants["CONST_VR"]),
            bool(self.constants["USE_CCN"]),
            self.constants["RTHRESHU"],
            self.constants["RTHRESHS"],
            self.constants["CCN_L"],
            self.constants["CCN_O"],
            self.constants["QC_CRT"],
            self.constants["TAU_G2V"],
            self.constants["TAU_V2G"],
            self.constants["TAU_S2V"],
            self.constants["TAU_V2S"],
            self.constants["TAU_REVP"],
            self.constants["TAU_FRZ"],
            bool(self.constants["DO_BIGG"]),
            bool(self.constants["DO_EVAP"]),
            bool(self.constants["DO_SUBL"]),
            self.constants["SAT_ADJ0"],
            self.constants["C_PIACR"],
            self.constants["TAU_IMLT"],
            self.constants["TAU_V2L"],
            self.constants["TAU_L2V"],
            self.constants["TAU_I2V"],
            self.constants["TAU_I2S"],
            self.constants["TAU_L2R"],
            self.constants["QI_LIM"],
            self.constants["QL_GEN"],
            self.constants["C_PAUT"],
            self.constants["C_PSACI"],
            self.constants["C_PGACS"],
            self.constants["C_PGACI"],
            bool(self.constants["Z_SLOPE_LIQ"]),
            bool(self.constants["Z_SLOPE_ICE"]),
            bool(self.constants["PROG_CCN"]),
            self.constants["C_CRACW"],
            self.constants["ALIN"],
            self.constants["CLIN"],
            bool(self.constants["PRECIPRAD"]),
            self.constants["CLD_MIN"],
            bool(self.constants["USE_PPM"]),
            bool(self.constants["MONO_PROF"]),
            bool(self.constants["DO_SEDI_HEAT"]),
            bool(self.constants["SEDI_TRANSPORT"]),
            bool(self.constants["DO_SEDI_W"]),
            bool(self.constants["DE_ICE"]),
            self.constants["ICLOUD_F"],
            self.constants["IRAIN_F"],
            bool(self.constants["MP_PRINT"]),
        )

        # Initalize object to be tested
        self.driver = MicrophysicsDriver(
            self.stencil_factory,
            self.quantity_factory,
            self.GFDL_1M_config,
        )

        self.driver(
            anv_icefall=inputs.pop("ANV_ICEFALL"),
            ls_icefall=inputs.pop("LS_ICEFALL"),
            GFDL_1M_config=self.GFDL_1M_config,
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
                "PFL_LS": self.driver.outputs.m2_rain.view[:],
                "PFI_LS": self.driver.outputs.m2_sol.view[:],
            }
        )

        return inputs
