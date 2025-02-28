from ndsl import Namelist, StencilFactory
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.GFDL_1M.driver.config import config
from ndsl.stencils.testing.savepoint import DataLoader
from pyMoist.GFDL_1M.driver.config import config
from pyMoist.GFDL_1M.driver.config_constants import ConfigConstants
from pyMoist.GFDL_1M.driver.temporaries import Temporaries
from pyMoist.GFDL_1M.driver.outputs import Outputs
from pyMoist.GFDL_1M.driver.masks import Masks
from pyMoist.GFDL_1M.driver.sat_tables import get_tables
from pyMoist.GFDL_1M.driver.fall_speed.main import FallSpeed
from ndsl.stencils.testing.grid import Grid


class Translatefall_speed(TranslateFortranData2Py):
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
            "p_dry": grid.compute_dict() | {"serialname": "p1_fall_speed"},
            "qs1": grid.compute_dict() | {"serialname": "qs_fall_speed"},
            "qi1": grid.compute_dict() | {"serialname": "qi_fall_speed"},
            "qg1": grid.compute_dict() | {"serialname": "qg_fall_speed"},
            "ql1": grid.compute_dict() | {"serialname": "ql_fall_speed"},
            "t": grid.compute_dict() | {"serialname": "t_fall_speed"},
            "t1": grid.compute_dict() | {"serialname": "t1_fall_speed"},
            "vts": grid.compute_dict() | {"serialname": "vts_fall_speed"},
            "vti": grid.compute_dict() | {"serialname": "vti_fall_speed"},
            "vtg": grid.compute_dict() | {"serialname": "vtg_fall_speed"},
            "den": grid.compute_dict() | {"serialname": "den_fall_speed"},
            "den1": grid.compute_dict() | {"serialname": "den1_fall_speed"},
            "denfac": grid.compute_dict() | {"serialname": "denfac_fall_speed"},
            "dz": grid.compute_dict() | {"serialname": "dz_fall_speed"},
            "dz1": grid.compute_dict() | {"serialname": "dz1_fall_speed"},
            "cnv_frc": grid.compute_dict() | {"serialname": "cnv_frc_fall_speed"},
        }
        self.in_vars["parameters"] = [
            "anv_icefall_fall_speed",
            "lsc_icefall_fall_speed",
        ]

        self.out_vars = self.in_vars["data_vars"].copy()
        del (
            self.out_vars["p_dry"],
            self.out_vars["den"],
            self.out_vars["dz"],
            self.out_vars["qs1"],
            self.out_vars["qi1"],
            self.out_vars["qg1"],
            self.out_vars["ql1"],
            self.out_vars["t"],
            self.out_vars["t1"],
            self.out_vars["cnv_frc"],
        )

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GFDL_1M_driver-constants")

    def compute_from_storage(self, inputs):
        self.GFDL_1M_config = config(
            self.constants["PHYS_HYDROSTATIC"],
            self.constants["HYDROSTATIC"],
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
            self.constants["DO_QA"],
            self.constants["FIX_NEGATIVE"],
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
            self.constants["FAST_SAT_ADJ"],
            self.constants["RH_INC"],
            self.constants["RH_INS"],
            self.constants["RH_INR"],
            self.constants["CONST_VI"],
            self.constants["CONST_VS"],
            self.constants["CONST_VG"],
            self.constants["CONST_VR"],
            self.constants["USE_CCN"],
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
            self.constants["DO_BIGG"],
            self.constants["DO_EVAP"],
            self.constants["DO_SUBL"],
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
            self.constants["Z_SLOPE_LIQ"],
            self.constants["Z_SLOPE_ICE"],
            self.constants["PROG_CCN"],
            self.constants["C_CRACW"],
            self.constants["ALIN"],
            self.constants["CLIN"],
            self.constants["PRECIPRAD"],
            self.constants["CLD_MIN"],
            self.constants["USE_PPM"],
            self.constants["MONO_PROF"],
            self.constants["DO_SEDI_HEAT"],
            self.constants["SEDI_TRANSPORT"],
            self.constants["DO_SEDI_W"],
            self.constants["DE_ICE"],
            self.constants["ICLOUD_F"],
            self.constants["IRAIN_F"],
            self.constants["MP_PRINT"],
        )

        self.config_dependent_constants = ConfigConstants(self.GFDL_1M_config)

        # Initalize saturation tables
        self.sat_tables = get_tables(self.stencil_factory.backend)

        # Initalize extra quantities
        temporaries = Temporaries(self.quantity_factory)
        outputs = Outputs(self.quantity_factory)
        masks = Masks(self.quantity_factory)

        # Initalize object to be tested
        self.fall_speed = FallSpeed(
            self.stencil_factory,
            self.quantity_factory,
            self.GFDL_1M_config,
            self.config_dependent_constants,
        )

        self.fall_speed(
            anv_icefall=inputs.pop("anv_icefall_fall_speed"),
            ls_icefall=inputs.pop("lsc_icefall_fall_speed"),
            **inputs
        )
        return inputs
