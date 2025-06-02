from ndsl import Namelist, StencilFactory
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.GFDL_1M.driver.config_constants import ConfigConstants
from pyMoist.GFDL_1M.driver.ice_cloud.main import IceCloud
from pyMoist.GFDL_1M.driver.masks import Masks
from pyMoist.GFDL_1M.driver.outputs import Outputs
from pyMoist.GFDL_1M.driver.sat_tables import get_tables
from pyMoist.GFDL_1M.driver.temporaries import Temporaries


class Translateice_cloud(TranslateFortranData2Py):
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
            "t1": grid.compute_dict() | {"serialname": "t1_icloud"},
            "p_dry": grid.compute_dict() | {"serialname": "p1_icloud"},
            "dp1": grid.compute_dict() | {"serialname": "dp1_icloud"},
            "qv1": grid.compute_dict() | {"serialname": "qv_icloud"},
            "ql1": grid.compute_dict() | {"serialname": "ql_icloud"},
            "qr1": grid.compute_dict() | {"serialname": "qr_icloud"},
            "qi1": grid.compute_dict() | {"serialname": "qi_icloud"},
            "qs1": grid.compute_dict() | {"serialname": "qs_icloud"},
            "qg1": grid.compute_dict() | {"serialname": "qg_icloud"},
            "qa1": grid.compute_dict() | {"serialname": "qa_icloud"},
            "vts": grid.compute_dict() | {"serialname": "vts_icloud"},
            "vtg": grid.compute_dict() | {"serialname": "vtg_icloud"},
            "vtr": grid.compute_dict() | {"serialname": "vtr_icloud"},
            "den1": grid.compute_dict() | {"serialname": "den1_icloud"},
            "denfac": grid.compute_dict() | {"serialname": "denfac_icloud"},
            "subl1": grid.compute_dict() | {"serialname": "subl1_icloud"},
            "rh_limited": grid.compute_dict() | {"serialname": "rh_limited_icloud"},
            "ccn": grid.compute_dict() | {"serialname": "ccn_icloud"},
            "cnv_frc": grid.compute_dict() | {"serialname": "cnv_frc_icloud"},
            "srf_type": grid.compute_dict() | {"serialname": "srf_type_icloud"},
            "isubl": grid.compute_dict() | {"serialname": "isubl_icloud"},
        }

        self.out_vars = self.in_vars["data_vars"].copy()
        del (
            self.out_vars["p_dry"],
            self.out_vars["dp1"],
            self.out_vars["vts"],
            self.out_vars["vtg"],
            self.out_vars["vtr"],
            self.out_vars["den1"],
            self.out_vars["denfac"],
            self.out_vars["subl1"],
            self.out_vars["rh_limited"],
            self.out_vars["ccn"],
            self.out_vars["cnv_frc"],
            self.out_vars["srf_type"],
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

        self.config_dependent_constants = ConfigConstants.make(self.GFDL_1M_config)

        # Initalize saturation tables
        self.sat_tables = get_tables(self.stencil_factory.backend)

        # Initalize extra quantities
        temporaries = Temporaries.make(self.quantity_factory)
        outputs = Outputs.make(self.quantity_factory)
        masks = Masks.make(self.quantity_factory)

        # Initalize object to be tested
        self.ice_cloud = IceCloud(
            self.stencil_factory,
            self.GFDL_1M_config,
            self.config_dependent_constants,
        )

        self.ice_cloud(
            table2=self.sat_tables.table2,
            table3=self.sat_tables.table3,
            des2=self.sat_tables.des2,
            des3=self.sat_tables.des3,
            **inputs,
        )
        return inputs
