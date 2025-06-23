from ndsl import Namelist, StencilFactory, Quantity
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.stencils.testing.savepoint import DataLoader
from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.GFDL_1M.state import (
    LiquidWaterStaticEnergy,
    TotalWater,
    VericalMotion,
    MixingRatios,
    CloudFractions,
)
from pyMoist.GFDL_1M.GFDL_1M import GFDL1M


class TranslateGFDL_1M(TranslateFortranData2Py):
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
            "Q": grid.compute_dict(),
            "QRAIN": grid.compute_dict(),
            "QSNOW": grid.compute_dict(),
            "QGRAUPEL": grid.compute_dict(),
            "QLCN": grid.compute_dict(),
            "QICN": grid.compute_dict(),
            "QLLS": grid.compute_dict(),
            "QILS": grid.compute_dict(),
            "CLCN": grid.compute_dict(),
            "CLLS": grid.compute_dict(),
            "NACTL": grid.compute_dict(),
            "NACTI": grid.compute_dict(),
            "AREA": grid.compute_dict(),
            "ZLE": grid.compute_dict(),
            "PLE": grid.compute_dict(),
            "T": grid.compute_dict(),
            "U": grid.compute_dict(),
            "V": grid.compute_dict(),
            "FRLAND": grid.compute_dict(),
            "W": grid.compute_dict(),
            "W2": grid.compute_dict(),
            "W3": grid.compute_dict(),
            "WSL": grid.compute_dict(),
            "SL2": grid.compute_dict(),
            "SL3": grid.compute_dict(),
            "WQT": grid.compute_dict(),
            "QT2": grid.compute_dict(),
            "QT3": grid.compute_dict(),
            "CNV_FRC": grid.compute_dict(),
            "SRF_TYPE": grid.compute_dict(),
            "SHLW_PRC3": grid.compute_dict(),
            "SHLW_SNO3": grid.compute_dict(),
            "RHCRIT": grid.compute_dict(),
            "RL": grid.compute_dict(),
            "RI": grid.compute_dict(),
            "EVAPC": grid.compute_dict(),
            "SUBLC": grid.compute_dict(),
            "PRCP_RAIN": grid.compute_dict(),
            "PRCP_SNOW": grid.compute_dict(),
            "PRCP_ICE": grid.compute_dict(),
            "PRCP_GRAUPEL": grid.compute_dict(),
            "FCLD": grid.compute_dict(),
            "QV": grid.compute_dict(),
            "QL": grid.compute_dict(),
            "QI": grid.compute_dict(),
            "QR": grid.compute_dict(),
            "QS": grid.compute_dict(),
            "QG": grid.compute_dict(),
            "LTS": grid.compute_dict(),
            "EIS": grid.compute_dict(),
            "ZLCL": grid.compute_dict(),
            "DUDT_macro": grid.compute_dict(),
            "DVDT_macro": grid.compute_dict(),
            "DTDT_macro": grid.compute_dict(),
            "DQVDT_macro": grid.compute_dict(),
            "DQLDT_macro": grid.compute_dict(),
            "DQIDT_macro": grid.compute_dict(),
            "DQADT_macro": grid.compute_dict(),
            "DQRDT_macro": grid.compute_dict(),
            "DQSDT_macro": grid.compute_dict(),
            "DQGDT_macro": grid.compute_dict(),
            "DUDT_micro": grid.compute_dict(),
            "DVDT_micro": grid.compute_dict(),
            "DTDT_micro": grid.compute_dict(),
            "DQVDT_micro": grid.compute_dict(),
            "DQLDT_micro": grid.compute_dict(),
            "DQIDT_micro": grid.compute_dict(),
            "DQADT_micro": grid.compute_dict(),
            "DQRDT_micro": grid.compute_dict(),
            "DQSDT_micro": grid.compute_dict(),
            "DQGDT_micro": grid.compute_dict(),
            "LS_PRCP": grid.compute_dict(),
            "LS_SNR": grid.compute_dict(),
            "ICE": grid.compute_dict(),
            "FRZR": grid.compute_dict(),
            "RHX": grid.compute_dict(),
            "REV_LS": grid.compute_dict(),
            "RSU_LS": grid.compute_dict(),
            "PFL_LS": grid.compute_dict(),
            "PFI_LS": grid.compute_dict(),
            "PFL_AN": grid.compute_dict(),
            "PFI_AN": grid.compute_dict(),
            "DQRL": grid.compute_dict(),
            "DBZ": grid.compute_dict(),
            "DBZ_MAX": grid.compute_dict(),
            "DBZ_1KM": grid.compute_dict(),
            "DBZ_TOP": grid.compute_dict(),
            "DBZ_M10C": grid.compute_dict(),
            "CN_PRCP": grid.compute_dict(),
            "AN_PRCP": grid.compute_dict(),
            "SC_PRCP": grid.compute_dict(),
            "CN_SNR": grid.compute_dict(),
            "AN_SNR": grid.compute_dict(),
            "SC_SNR": grid.compute_dict(),
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

        # Initalize the module
        gfdl_1m = GFDL1M(self.stencil_factory, self.quantity_factory, GFDL_1M_config)

        ##### Being in fields from netcdf. This replicated the call to fortran memory. #####
        # Get model state
        gfdl_1m.mixing_ratios = MixingRatios(
            vapor=self.make_ijk_quantity(inputs.pop("Q")),
            rain=self.make_ijk_quantity(inputs.pop("QRAIN")),
            snow=self.make_ijk_quantity(inputs.pop("QSNOW")),
            graupel=self.make_ijk_quantity(inputs.pop("QGRAUPEL")),
            convective_liquid=self.make_ijk_quantity(inputs.pop("QLCN")),
            convective_ice=self.make_ijk_quantity(inputs.pop("QICN")),
            large_scale_liquid=self.make_ijk_quantity(inputs.pop("QLLS")),
            large_scale_ice=self.make_ijk_quantity(inputs.pop("QILS")),
        )
        gfdl_1m.cloud_fractions = CloudFractions(
            convective=self.make_ijk_quantity(inputs.pop("CLCN")),
            large_scale=self.make_ijk_quantity(inputs.pop("CLLS")),
        )
        gfdl_1m.liquid_concentration = self.make_ijk_quantity(inputs.pop("NACTL"))
        gfdl_1m.ice_concentration = self.make_ijk_quantity(inputs.pop("NACTI"))
        gfdl_1m.area = self.make_ij_quantity(inputs.pop("AREA"))
        gfdl_1m.geopotential_height_interface = self.make_ijk_quantity(inputs.pop("ZLE"), interface=True)
        gfdl_1m.p_interface = self.make_ijk_quantity(inputs.pop("PLE"), interface=True)
        gfdl_1m.t = self.make_ijk_quantity(inputs.pop("T"))
        gfdl_1m.u = self.make_ijk_quantity(inputs.pop("U"))
        gfdl_1m.v = self.make_ijk_quantity(inputs.pop("V"))
        gfdl_1m.land_fraction = self.make_ij_quantity(inputs.pop("FRLAND"))
        gfdl_1m.vertical_motion = VericalMotion(
            velocity=self.make_ijk_quantity(inputs.pop("W")),
            variance=self.make_ijk_quantity(inputs.pop("W2")),
            third_moment=self.make_ijk_quantity(inputs.pop("W3")),
        )
        gfdl_1m.liquid_water_static_energy = LiquidWaterStaticEnergy(
            flux=self.make_ijk_quantity(inputs.pop("WSL")),
            variance=self.make_ijk_quantity(inputs.pop("SL2")),
            third_moment=self.make_ijk_quantity(inputs.pop("SL3")),
        )
        gfdl_1m.total_water = TotalWater(
            flux=self.make_ijk_quantity(inputs.pop("WQT")),
            variance=self.make_ijk_quantity(inputs.pop("QT2")),
            third_moment=self.make_ijk_quantity(inputs.pop("QT3")),
        )
        gfdl_1m.convection_fraction = self.make_ij_quantity(inputs.pop("CNV_FRC"))
        gfdl_1m.surface_type = self.make_ij_quantity(inputs.pop("SRF_TYPE"))
        gfdl_1m.shallow_convective_rain = self.make_ijk_quantity(inputs.pop("SHLW_PRC3"))
        gfdl_1m.shallow_convective_snow = self.make_ijk_quantity(inputs.pop("SHLW_SNO3"))
        gfdl_1m.rh_crit = self.make_ijk_quantity(inputs.pop("RHCRIT"))

        # Outputs: model fields originating from within GFDL
        gfdl_1m.outputs.liquid_radius = self.make_ijk_quantity(inputs.pop("RL"))
        gfdl_1m.outputs.ice_radius = self.make_ijk_quantity(inputs.pop("RI"))
        gfdl_1m.outputs.large_scale_nonanvil_precipitation_evaporation = self.make_ijk_quantity(
            inputs.pop("EVAPC")
        )
        gfdl_1m.outputs.large_scale_nonanvil_precipitation_sublimation = self.make_ijk_quantity(
            inputs.pop("SUBLC")
        )
        gfdl_1m.outputs.precipitated_rain = self.make_ij_quantity(inputs.pop("PRCP_RAIN"))
        gfdl_1m.outputs.precipitated_snow = self.make_ij_quantity(inputs.pop("PRCP_SNOW"))
        gfdl_1m.outputs.precipitated_ice = self.make_ij_quantity(inputs.pop("PRCP_ICE"))
        gfdl_1m.outputs.precipitated_graupel = self.make_ij_quantity(inputs.pop("PRCP_GRAUPEL"))

        # Outputs: model fields originating from within GFDL; radiation fields
        gfdl_1m.outputs.radiation_cloud_fraction = self.make_ijk_quantity(inputs.pop("FCLD"))
        gfdl_1m.outputs.radiation_vapor = self.make_ijk_quantity(inputs.pop("QV"))
        gfdl_1m.outputs.radiation_liquid = self.make_ijk_quantity(inputs.pop("QL"))
        gfdl_1m.outputs.radiation_ice = self.make_ijk_quantity(inputs.pop("QI"))
        gfdl_1m.outputs.radiation_rain = self.make_ijk_quantity(inputs.pop("QR"))
        gfdl_1m.outputs.radiation_snow = self.make_ijk_quantity(inputs.pop("QS"))
        gfdl_1m.outputs.radiation_graupel = self.make_ijk_quantity(inputs.pop("QG"))
        gfdl_1m.outputs.lower_tropospheric_stability = self.make_ijk_quantity(inputs.pop("LTS"))
        gfdl_1m.outputs.estimated_inversion_strength = self.make_ijk_quantity(inputs.pop("EIS"))
        gfdl_1m.outputs.z_lcl = self.make_ij_quantity(inputs.pop("ZLCL"))

        # Outputs: model fields originating from within GFDL; macrophysics/microphysics tendencies
        gfdl_1m.outputs.du_dt_macro = self.make_ijk_quantity(inputs.pop("DUDT_macro"))
        gfdl_1m.outputs.dv_dt_macro = self.make_ijk_quantity(inputs.pop("DVDT_macro"))
        gfdl_1m.outputs.dt_dt_macro = self.make_ijk_quantity(inputs.pop("DTDT_macro"))
        gfdl_1m.outputs.dvapor_dt_macro = self.make_ijk_quantity(inputs.pop("DQVDT_macro"))
        gfdl_1m.outputs.dliquid_dt_macro = self.make_ijk_quantity(inputs.pop("DQLDT_macro"))
        gfdl_1m.outputs.dice_dt_macro = self.make_ijk_quantity(inputs.pop("DQIDT_macro"))
        gfdl_1m.outputs.dcloud_fraction_dt_macro = self.make_ijk_quantity(inputs.pop("DQADT_macro"))
        gfdl_1m.outputs.drain_dt_macro = self.make_ijk_quantity(inputs.pop("DQRDT_macro"))
        gfdl_1m.outputs.dsnow_dt_macro = self.make_ijk_quantity(inputs.pop("DQSDT_macro"))
        gfdl_1m.outputs.dgraupel_dt_macro = self.make_ijk_quantity(inputs.pop("DQGDT_macro"))
        gfdl_1m.outputs.du_dt_micro = self.make_ijk_quantity(inputs.pop("DUDT_micro"))
        gfdl_1m.outputs.dv_dt_micro = self.make_ijk_quantity(inputs.pop("DVDT_micro"))
        gfdl_1m.outputs.dt_dt_micro = self.make_ijk_quantity(inputs.pop("DTDT_micro"))
        gfdl_1m.outputs.dvapor_dt_micro = self.make_ijk_quantity(inputs.pop("DQVDT_micro"))
        gfdl_1m.outputs.dliquid_dt_micro = self.make_ijk_quantity(inputs.pop("DQLDT_micro"))
        gfdl_1m.outputs.dice_dt_micro = self.make_ijk_quantity(inputs.pop("DQIDT_micro"))
        gfdl_1m.outputs.dcloud_fraction_dt_micro = self.make_ijk_quantity(inputs.pop("DQADT_micro"))
        gfdl_1m.outputs.drain_dt_micro = self.make_ijk_quantity(inputs.pop("DQRDT_micro"))
        gfdl_1m.outputs.dsnow_dt_micro = self.make_ijk_quantity(inputs.pop("DQSDT_micro"))
        gfdl_1m.outputs.dgraupel_dt_micro = self.make_ijk_quantity(inputs.pop("DQGDT_micro"))
        # Outputs: Exports to be filled
        gfdl_1m.outputs.large_scale_precip = self.make_ij_quantity(inputs.pop("LS_PRCP"))
        gfdl_1m.outputs.large_scale_snow = self.make_ij_quantity(inputs.pop("LS_SNR"))
        gfdl_1m.outputs.icefall = self.make_ij_quantity(inputs.pop("ICE"))
        gfdl_1m.outputs.freezing_rainfall = self.make_ij_quantity(inputs.pop("FRZR"))
        gfdl_1m.outputs.relative_humidity_after_pdf = self.make_ijk_quantity(inputs.pop("RHX"))
        gfdl_1m.outputs.large_scale_nonanvil_precipitation_evaporation = self.make_ijk_quantity(
            inputs.pop("REV_LS")
        )
        gfdl_1m.outputs.large_scale_nonanvil_precipitation_sublimation = self.make_ijk_quantity(
            inputs.pop("RSU_LS")
        )
        gfdl_1m.outputs.large_scale_nonanvil_liquid_flux = self.make_ijk_quantity(
            inputs.pop("PFL_LS"), interface=True
        )
        gfdl_1m.outputs.large_scale_nonanvil_ice_flux = self.make_ijk_quantity(
            inputs.pop("PFI_LS"), interface=True
        )
        gfdl_1m.outputs.anvil_liquid_flux = self.make_ijk_quantity(inputs.pop("PFL_AN"), interface=True)
        gfdl_1m.outputs.anvil_ice_flux = self.make_ijk_quantity(inputs.pop("PFI_AN"), interface=True)
        gfdl_1m.outputs.large_scale_rainwater_source = self.make_ijk_quantity(inputs.pop("DQRL"))
        gfdl_1m.outputs.simulated_reflectivity = self.make_ijk_quantity(inputs.pop("DBZ"))
        gfdl_1m.outputs.maximum_reflectivity = self.make_ij_quantity(inputs.pop("DBZ_MAX"))
        gfdl_1m.outputs.one_km_agl_reflectivity = self.make_ij_quantity(inputs.pop("DBZ_1KM"))
        gfdl_1m.outputs.echo_top_reflectivity = self.make_ij_quantity(inputs.pop("DBZ_TOP"))
        gfdl_1m.outputs.minus_10c_reflectivity = self.make_ij_quantity(inputs.pop("DBZ_M10C"))
        # Unused fields, force to zero
        gfdl_1m.temporaries.all_zeros_3d = self.make_ijk_quantity(inputs.pop("CN_PRCP"))
        gfdl_1m.temporaries.all_zeros_3d = self.make_ijk_quantity(inputs.pop("AN_PRCP"))
        gfdl_1m.temporaries.all_zeros_3d = self.make_ijk_quantity(inputs.pop("SC_PRCP"))
        gfdl_1m.temporaries.all_zeros_3d = self.make_ijk_quantity(inputs.pop("CN_SNR"))
        gfdl_1m.temporaries.all_zeros_3d = self.make_ijk_quantity(inputs.pop("AN_SNR"))
        gfdl_1m.temporaries.all_zeros_3d = self.make_ijk_quantity(inputs.pop("SC_SNR"))

        # Initalize the module

        return {
            "Q": gfdl_1m.mixing_ratios.vapor.field,
            "QRAIN": gfdl_1m.mixing_ratios.rain.field,
            "QSNOW": gfdl_1m.mixing_ratios.snow.field,
            "QGRAUPEL": gfdl_1m.mixing_ratios.graupel.field,
            "QLCN": gfdl_1m.mixing_ratios.convective_liquid.field,
            "QICN": gfdl_1m.mixing_ratios.convective_ice.field,
            "QLLS": gfdl_1m.mixing_ratios.large_scale_liquid.field,
            "QILS": gfdl_1m.mixing_ratios.large_scale_ice.field,
            "CLCN": gfdl_1m.cloud_fractions.convective.field,
            "CLLS": gfdl_1m.cloud_fractions.large_scale.field,
            "NACTL": gfdl_1m.liquid_concentration.field,
            "NACTI": gfdl_1m.ice_concentration.field,
            "AREA": gfdl_1m.area.field,
            "ZLE": gfdl_1m.geopotential_height_interface.field,
            "PLE": gfdl_1m.p_interface.field,
            "T": gfdl_1m.t.field,
            "U": gfdl_1m.u.field,
            "V": gfdl_1m.v.field,
            "FRLAND": gfdl_1m.land_fraction.field,
            "W": gfdl_1m.vertical_motion.velocity.field,
            "W2":gfdl_1m.vertical_motion.variance.field, 
            "W3":gfdl_1m.vertical_motion.third_moment.field, 
            "WSL": gfdl_1m.liquid_water_static_energy.flux.field,
            "SL2": gfdl_1m.liquid_water_static_energy.variance.field,
            "SL3": gfdl_1m.liquid_water_static_energy.third_moment.field,
            "WQT": gfdl_1m.total_water.flux.field,
            "QT2": gfdl_1m.total_water.variance.field,
            "QT3": gfdl_1m.total_water.third_moment.field,
            "CNV_FRC": gfdl_1m.convection_fraction.field,
            "SRF_TYPE": gfdl_1m.surface_type.field,
            "SHLW_PRC3": gfdl_1m.shallow_convective_rain.field,
            "SHLW_SNO3": gfdl_1m.shallow_convective_snow.field,
            "RHCRIT": gfdl_1m.rh_crit.field,
            "RL": gfdl_1m.outputs.liquid_radius.field,
            "RI": gfdl_1m.outputs.ice_radius.field,
            "EVAPC": gfdl_1m.outputs.large_scale_nonanvil_precipitation_evaporation.field,
            "SUBLC": gfdl_1m.outputs.large_scale_nonanvil_precipitation_sublimation.field,
            "PRCP_RAIN": gfdl_1m.outputs.precipitated_rain.field,
            "PRCP_SNOW": gfdl_1m.outputs.precipitated_snow.field,
            "PRCP_ICE": gfdl_1m.outputs.precipitated_ice.field,
            "PRCP_GRAUPEL": gfdl_1m.outputs.precipitated_graupel.field,
            "FCLD": gfdl_1m.outputs.radiation_cloud_fraction.field,
            "QV": gfdl_1m.outputs.radiation_vapor.field,
            "QL": gfdl_1m.outputs.radiation_liquid.field,
            "QI": gfdl_1m.outputs.radiation_ice.field,
            "QR": gfdl_1m.outputs.radiation_rain.field,
            "QS": gfdl_1m.outputs.radiation_snow.field,
            "QG": gfdl_1m.outputs.radiation_graupel.field,
            "LTS": gfdl_1m.outputs.lower_tropospheric_stability.field,
            "EIS": gfdl_1m.outputs.estimated_inversion_strength.field,
            "ZLCL": gfdl_1m.outputs.z_lcl.field,
            "DUDT_macro": gfdl_1m.outputs.du_dt_macro.field,
            "DVDT_macro": gfdl_1m.outputs.dv_dt_macro.field,
            "DTDT_macro": gfdl_1m.outputs.dt_dt_macro.field,
            "DQVDT_macro": gfdl_1m.outputs.dvapor_dt_macro.field,
            "DQLDT_macro": gfdl_1m.outputs.dliquid_dt_macro.field,
            "DQIDT_macro": gfdl_1m.outputs.dice_dt_macro.field,
            "DQADT_macro": gfdl_1m.outputs.dcloud_fraction_dt_macro.field,
            "DQRDT_macro": gfdl_1m.outputs.drain_dt_macro.field,
            "DQSDT_macro": gfdl_1m.outputs.dsnow_dt_macro.field,
            "DQGDT_macro": gfdl_1m.outputs.dgraupel_dt_macro.field,
            "DUDT_micro": gfdl_1m.outputs.du_dt_micro.field,
            "DVDT_micro": gfdl_1m.outputs.dv_dt_micro.field,
            "DTDT_micro": gfdl_1m.outputs.dt_dt_micro.field,
            "DQVDT_micro": gfdl_1m.outputs.dvapor_dt_micro.field,
            "DQLDT_micro": gfdl_1m.outputs.dliquid_dt_micro.field,
            "DQIDT_micro": gfdl_1m.outputs.dice_dt_micro.field,
            "DQADT_micro": gfdl_1m.outputs.dcloud_fraction_dt_micro.field,
            "DQRDT_micro": gfdl_1m.outputs.drain_dt_micro.field,
            "DQSDT_micro": gfdl_1m.outputs.dsnow_dt_micro.field,
            "DQGDT_micro": gfdl_1m.outputs.dgraupel_dt_micro.field,
            "LS_PRCP": gfdl_1m.outputs.large_scale_precip.field,
            "LS_SNR": gfdl_1m.outputs.large_scale_snow.field,
            "ICE": gfdl_1m.outputs.icefall.field,
            "FRZR": gfdl_1m.outputs.freezing_rainfall.field,
            "RHX": gfdl_1m.outputs.relative_humidity_after_pdf.field,
            "REV_LS": gfdl_1m.outputs.large_scale_nonanvil_precipitation_evaporation.field,
            "RSU_LS": gfdl_1m.outputs.large_scale_nonanvil_precipitation_sublimation.field,
            "PFL_LS": gfdl_1m.outputs.large_scale_nonanvil_liquid_flux.field,
            "PFI_LS": gfdl_1m.outputs.large_scale_nonanvil_ice_flux.field,
            "PFL_AN": gfdl_1m.outputs.anvil_liquid_flux.field,
            "PFI_AN": gfdl_1m.outputs.anvil_ice_flux.field,
            "DQRL": gfdl_1m.outputs.large_scale_rainwater_source.field,
            "DBZ": gfdl_1m.outputs.simulated_reflectivity.field,
            "DBZ_MAX": gfdl_1m.outputs.maximum_reflectivity.field,
            "DBZ_1KM": gfdl_1m.outputs.one_km_agl_reflectivity.field,
            "DBZ_TOP": gfdl_1m.outputs.echo_top_reflectivity.field,
            "DBZ_M10C": gfdl_1m.outputs.minus_10c_reflectivity.field,
            "CN_PRCP": gfdl_1m.temporaries.all_zeros_3d.field,
            "AN_PRCP": gfdl_1m.temporaries.all_zeros_3d.field,
            "SC_PRCP": gfdl_1m.temporaries.all_zeros_3d.field,
            "CN_SNR": gfdl_1m.temporaries.all_zeros_3d.field,
            "AN_SNR": gfdl_1m.temporaries.all_zeros_3d.field,
            "SC_SNR": gfdl_1m.temporaries.all_zeros_3d.field,
        }
