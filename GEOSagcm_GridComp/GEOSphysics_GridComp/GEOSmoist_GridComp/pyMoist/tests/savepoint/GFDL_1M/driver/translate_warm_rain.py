from ndsl import Namelist, StencilFactory
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.GFDL_1M.driver.config import MicrophysicsConfiguration
from pyMoist.GFDL_1M.driver.config_constants import ConfigConstants
from pyMoist.GFDL_1M.driver.masks import Masks
from pyMoist.GFDL_1M.driver.outputs import Outputs
from pyMoist.GFDL_1M.driver.sat_tables import get_tables
from pyMoist.GFDL_1M.driver.temporaries import Temporaries
from pyMoist.GFDL_1M.driver.warm_rain.main import WarmRain


class Translatewarm_rain(TranslateFortranData2Py):
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
        self.constants = data_loader.load("GFDL_1M_driver-constants")

    def compute_from_storage(self, inputs):
        self.GFDL_1M_config = MicrophysicsConfiguration(
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
            **inputs
        )
        return inputs
