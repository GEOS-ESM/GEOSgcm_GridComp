from pyMoist.GFDL_1M.driver.driver import MicrophysicsDriver
from pyMoist.GFDL_1M.PhaseChange.phase_change import PhaseChange
from ndsl import QuantityFactory, StencilFactory
from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.GFDL_1M.getters_temporary import (
    mapl_get_resource_placeholder,
    esmf_placeholder,
    mapl_get_pointer_placeholder,
    mapl_get_pointer_placeholder,
    associated_checker,
)
from ndsl.dsl.typing import Float
from pyMoist.GFDL_1M.state import (
    LiquidWaterStaticEnergy,
    TotalWater,
    VericalMotion,
    MixingRatios,
    CloudFractions,
)
from pyMoist.GFDL_1M.stencils import (
    calculate_derived_states,
    find_klcl,
    vertical_interpolation,
    find_eis,
    prepare_tendencies,
    prepare_radiation_quantities,
    update_tendencies,
    apply_driver_tendencies,
)
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.GFDL_1M.outputs import Outputs
from pyMoist.GFDL_1M.temporaries import Temporaries
from pyMoist.GFDL_1M.masks import Masks
from pyMoist.GFDL_1M.PhaseChange.config import PhaseChangeConfiguration
from pyMoist.redistribute_clouds import RedistributeClouds


class GFDL1M:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        MP_TIME: Float,
        T_MIN: Float,
        T_SUB: Float,
        TAU_R2G: Float,
        TAU_SMLT: Float,
        TAU_G2R: Float,
        DW_LAND: Float,
        DW_OCEAN: Float,
        VI_FAC: Float,
        VR_FAC: Float,
        VS_FAC: Float,
        VG_FAC: Float,
        QL_MLT: Float,
        DO_QA: Float,
        FIX_NEGATIVE: bool,
        VI_MAX: Float,
        VS_MAX: Float,
        VG_MAX: Float,
        VR_MAX: Float,
        QS_MLT: Float,
        QS0_CRT: Float,
        QI_GEN: Float,
        QL0_MAX: Float,
        QI0_MAX: Float,
        QI0_CRT: Float,
        QR0_CRT: Float,
        FAST_SAT_ADJ: bool,
        RH_INC: Float,
        RH_INS: Float,
        RH_INR: Float,
        CONST_VI: bool,
        CONST_VS: bool,
        CONST_VG: bool,
        CONST_VR: bool,
        USE_CCN: bool,
        RTHRESHU: Float,
        RTHRESHS: Float,
        CCN_L: Float,
        CCN_O: Float,
        QC_CRT: Float,
        TAU_G2V: Float,
        TAU_V2G: Float,
        TAU_S2V: Float,
        TAU_V2S: Float,
        TAU_REVP: Float,
        TAU_FRZ: Float,
        DO_BIGG: bool,
        DO_EVAP: bool,
        DO_SUBL: bool,
        SAT_ADJ0: Float,
        C_PIACR: Float,
        TAU_IMLT: Float,
        TAU_V2L: Float,
        TAU_L2V: Float,
        TAU_I2V: Float,
        TAU_I2S: Float,
        TAU_L2R: Float,
        QI_LIM: Float,
        QL_GEN: Float,
        C_PAUT: Float,
        C_PSACI: Float,
        C_PGACS: Float,
        C_PGACI: Float,
        Z_SLOPE_LIQ: bool,
        Z_SLOPE_ICE: bool,
        PROG_CCN: bool,
        C_CRACW: Float,
        ALIN: Float,
        CLIN: Float,
        PRECIPRAD: bool,
        CLD_MIN: Float,
        USE_PPM: bool,
        MONO_PROF: bool,
        DO_SEDI_HEAT: bool,
        SEDI_TRANSPORT: bool,
        DO_SEDI_W: bool,
        DE_ICE: bool,
        ICLOUD_F: Float,
        IRAIN_F: Float,
        MP_PRINT: bool,
        LMELTFRZ: bool,
        USE_BERGERON: bool,
    ):
        self.stencil_factory = stencil_factory
        if self.stencil_factory.grid_indexing.n_halo != 0:
            raise ValueError("halo needs to be zero for GFDL Single Moment microphysics")
        self.quantity_factory = quantity_factory

        self.PHYS_HYDROSTATIC = mapl_get_resource_placeholder("LPHYS_HYDROSTATIC")
        self.HYDROSTATIC = mapl_get_resource_placeholder("LHYDROSTATIC")
        self.DT_MOIST = esmf_placeholder("DT_R8")
        self.MELTFRZ = mapl_get_resource_placeholder("MELTFRZ")
        self.TURNRHCRIT_PARAM = mapl_get_resource_placeholder("TURNRHCRIT")
        self.PDF_SHAPE = mapl_get_resource_placeholder("PDFSHAPE")
        self.ANV_ICEFALL = mapl_get_resource_placeholder("ANV_ICEFALL")
        self.LS_ICEFALL = mapl_get_resource_placeholder("LS_ICEFALL")
        self.LIQ_RADII_PARAM = mapl_get_resource_placeholder("LIQ_RADII_PARAM")
        self.ICE_RADII_PARAM = mapl_get_resource_placeholder("ICE_RADII_PARAM")
        self.FAC_RI = mapl_get_resource_placeholder("FAC_RI")
        self.MIN_RI = mapl_get_resource_placeholder("MIN_RI")
        self.MAX_RI = mapl_get_resource_placeholder("MAX_RI")
        self.FAC_RL = mapl_get_resource_placeholder("FAC_RL")
        self.MIN_RL = mapl_get_resource_placeholder("MIN_RL")
        self.MAX_RL = mapl_get_resource_placeholder("MAX_RL")
        self.CCW_EVAP_EFF = mapl_get_resource_placeholder("CCW_EVAP_EFF")
        self.CCI_EVAP_EFF = mapl_get_resource_placeholder("CCI_EVAP_EFF")

        self.GFDL_1M_config = GFDL1MConfig(
            PHYS_HYDROSTATIC=self.PHYS_HYDROSTATIC,
            HYDROSTATIC=self.HYDROSTATIC,
            DT_MOIST=self.DT_MOIST,
            MP_TIME=MP_TIME,
            T_MIN=T_MIN,
            T_SUB=T_SUB,
            TAU_R2G=TAU_R2G,
            TAU_SMLT=TAU_SMLT,
            TAU_G2R=TAU_G2R,
            DW_LAND=DW_LAND,
            DW_OCEAN=DW_OCEAN,
            VI_FAC=VI_FAC,
            VR_FAC=VR_FAC,
            VS_FAC=VS_FAC,
            VG_FAC=VG_FAC,
            QL_MLT=QL_MLT,
            DO_QA=DO_QA,
            FIX_NEGATIVE=FIX_NEGATIVE,
            VI_MAX=VI_MAX,
            VS_MAX=VS_MAX,
            VG_MAX=VG_MAX,
            VR_MAX=VR_MAX,
            QS_MLT=QS_MLT,
            QS0_CRT=QS0_CRT,
            QI_GEN=QI_GEN,
            QL0_MAX=QL0_MAX,
            QI0_MAX=QI0_MAX,
            QI0_CRT=QI0_CRT,
            QR0_CRT=QR0_CRT,
            FAST_SAT_ADJ=FAST_SAT_ADJ,
            RH_INC=RH_INC,
            RH_INS=RH_INS,
            RH_INR=RH_INR,
            CONST_VI=CONST_VI,
            CONST_VS=CONST_VS,
            CONST_VG=CONST_VG,
            CONST_VR=CONST_VR,
            USE_CCN=USE_CCN,
            RTHRESHU=RTHRESHU,
            RTHRESHS=RTHRESHS,
            CCN_L=CCN_L,
            CCN_O=CCN_O,
            QC_CRT=QC_CRT,
            TAU_G2V=TAU_G2V,
            TAU_V2G=TAU_V2G,
            TAU_S2V=TAU_S2V,
            TAU_V2S=TAU_V2S,
            TAU_REVP=TAU_REVP,
            TAU_FRZ=TAU_FRZ,
            DO_BIGG=DO_BIGG,
            DO_EVAP=DO_EVAP,
            DO_SUBL=DO_SUBL,
            SAT_ADJ0=SAT_ADJ0,
            C_PIACR=C_PIACR,
            TAU_IMLT=TAU_IMLT,
            TAU_V2L=TAU_V2L,
            TAU_L2V=TAU_L2V,
            TAU_I2V=TAU_I2V,
            TAU_I2S=TAU_I2S,
            TAU_L2R=TAU_L2R,
            QI_LIM=QI_LIM,
            QL_GEN=QL_GEN,
            C_PAUT=C_PAUT,
            C_PSACI=C_PSACI,
            C_PGACS=C_PGACS,
            C_PGACI=C_PGACI,
            Z_SLOPE_LIQ=Z_SLOPE_LIQ,
            Z_SLOPE_ICE=Z_SLOPE_ICE,
            PROG_CCN=PROG_CCN,
            C_CRACW=C_CRACW,
            ALIN=ALIN,
            CLIN=CLIN,
            PRECIPRAD=PRECIPRAD,
            CLD_MIN=CLD_MIN,
            USE_PPM=USE_PPM,
            MONO_PROF=MONO_PROF,
            DO_SEDI_HEAT=DO_SEDI_HEAT,
            SEDI_TRANSPORT=SEDI_TRANSPORT,
            DO_SEDI_W=DO_SEDI_W,
            DE_ICE=DE_ICE,
            ICLOUD_F=ICLOUD_F,
            IRAIN_F=IRAIN_F,
            MP_PRINT=MP_PRINT,
        )

        # Initalize saturation tables
        self.saturation_tables = SaturationVaporPressureTable(self.stencil_factory.backend)

        # Initalize internal fields
        self.outputs = Outputs.make(quantity_factory=quantity_factory)
        self.temporaries = Temporaries.make(quantity_factory=quantity_factory)
        self.masks = Masks.make(quantity_factory=quantity_factory)

        # Construct stencils
        self.calculate_derived_states = stencil_factory.from_dims_halo(
            func=calculate_derived_states,
            compute_dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM],
        )

        self.find_klcl = stencil_factory.from_dims_halo(
            func=find_klcl,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self.vertical_interpolation = stencil_factory.from_dims_halo(
            func=vertical_interpolation,
            compute_dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM],
        )

        self.find_eis = stencil_factory.from_dims_halo(
            func=find_eis,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self.phase_change_config = PhaseChangeConfiguration(
            DT_MOIST=self.DT_MOIST,
            PDF_SHAPE=self.PDF_SHAPE,
            CCW_EVAP_EFF=self.CCW_EVAP_EFF,
            CCI_EVAP_EFF=self.CCI_EVAP_EFF,
            TURNRHCRIT_PARAM=self.TURNRHCRIT_PARAM,
            DW_LAND=self.GFDL_1M_config.DW_LAND,
            DW_OCEAN=self.GFDL_1M_config.DW_OCEAN,
            DO_QA=self.GFDL_1M_config.DO_QA,
            DO_MELT_FREEZE=self.MELTFRZ,
            USE_BERGERON=USE_BERGERON,
        )

        self.phase_change = PhaseChange(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            phase_change_config=self.phase_change_config,
        )

        self.update_macrophysics_tendencies = stencil_factory.from_dims_halo(
            func=update_tendencies,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "DT_MOIST": self.DT_MOIST,
            },
        )

        self.prepare_tendencies = stencil_factory.from_dims_halo(
            func=prepare_tendencies,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self.prepare_radiation_quantities = stencil_factory.from_dims_halo(
            func=prepare_radiation_quantities,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self.driver = MicrophysicsDriver(
            self.stencil_factory,
            self.quantity_factory,
            self.GFDL_1M_config,
        )

        self.apply_driver_tendencies = stencil_factory.from_dims_halo(
            func=apply_driver_tendencies,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "DT_MOIST": self.DT_MOIST,
            },
        )

        self.redistribute_clouds = RedistributeClouds(
            self.stencil_factory,
            self.quantity_factory,
        )

    def __call__(self, *args, **kwds):

        ##### Link FORTRAN memory to Python memory #####
        ##### All links are live, so changes in Python will be reflected in FORTRAN #####
        ##### Not all "input" (model state) fields are modified #####

        # Get model state
        mixing_ratios = MixingRatios(
            vapor=mapl_get_pointer_placeholder("Q"),
            rain=mapl_get_pointer_placeholder("QRAIN"),
            snow=mapl_get_pointer_placeholder("QSNOW"),
            graupel=mapl_get_pointer_placeholder("QGRAUPEL"),
            convective_liquid=mapl_get_pointer_placeholder("QLCN"),
            convective_ice=mapl_get_pointer_placeholder("QICN"),
            large_scale_liquid=mapl_get_pointer_placeholder("QLLS"),
            large_scale_ice=mapl_get_pointer_placeholder("QILS"),
        )
        cloud_fractions = CloudFractions(
            convective=mapl_get_pointer_placeholder("CLCN"),
            large_scale=mapl_get_pointer_placeholder("CLLS"),
        )
        liquid_concentration = mapl_get_pointer_placeholder("NACTL")
        ice_concentration = mapl_get_pointer_placeholder("NACTI")
        area = mapl_get_pointer_placeholder("AREA")
        geopotential_height_interface = mapl_get_pointer_placeholder("ZLE")
        p_interface = mapl_get_pointer_placeholder("PLE")
        t = mapl_get_pointer_placeholder("T")
        u = mapl_get_pointer_placeholder("U")
        v = mapl_get_pointer_placeholder("V")
        land_fraction = mapl_get_pointer_placeholder("FRLAND")
        vertical_motion = VericalMotion(
            velocity=mapl_get_pointer_placeholder("W"),
            variance=mapl_get_pointer_placeholder("W2"),
            third_moment=mapl_get_pointer_placeholder("W3"),
        )
        liquid_water_static_energy = LiquidWaterStaticEnergy(
            flux=mapl_get_pointer_placeholder("WSL"),
            variance=mapl_get_pointer_placeholder("SL2"),
            third_moment=mapl_get_pointer_placeholder("SL3"),
        )
        total_water = TotalWater(
            flux=mapl_get_pointer_placeholder("WQT"),
            variance=mapl_get_pointer_placeholder("QT2"),
            third_moment=mapl_get_pointer_placeholder("QT3"),
        )
        convection_fraction = mapl_get_pointer_placeholder("CNV_FRC")
        surface_type = mapl_get_pointer_placeholder("SRF_TYPE")
        shallow_convective_rain = mapl_get_pointer_placeholder("SHLW_PRC3")
        shallow_convective_snow = mapl_get_pointer_placeholder("SHLW_SNO3")
        rh_crit = mapl_get_pointer_placeholder("RHCRIT")

        # Unused fields, force to zero
        self.temporaries.all_zeros_3d = mapl_get_pointer_placeholder("CN_PRCP")
        self.temporaries.all_zeros_3d = mapl_get_pointer_placeholder("AN_PRCP")
        self.temporaries.all_zeros_3d = mapl_get_pointer_placeholder("SC_PRCP")
        self.temporaries.all_zeros_3d = mapl_get_pointer_placeholder("CN_SNR")
        self.temporaries.all_zeros_3d = mapl_get_pointer_placeholder("AN_SNR")
        self.temporaries.all_zeros_3d = mapl_get_pointer_placeholder("SC_SNR")

        # Outputs: model fields originating from within GFDL
        self.outputs.liquid_radius = mapl_get_pointer_placeholder("RL")
        self.outputs.ice_radius = mapl_get_pointer_placeholder("RI")
        self.outputs.evaporation = mapl_get_resource_placeholder("EVAPC")
        self.outputs.sublimation = mapl_get_resource_placeholder("SUBLC")
        self.outputs.precipitated_rain = mapl_get_resource_placeholder("PRCP_RAIN")
        self.outputs.precipitated_snow = mapl_get_resource_placeholder("PRCP_SNOW")
        self.outputs.precipitated_ice = mapl_get_resource_placeholder("PRCP_ICE")
        self.outputs.precipitated_graupel = mapl_get_resource_placeholder("PRCP_GRAUPEL")

        # Outputs: model fields originating from within GFDL; radiation fields
        self.outputs.radiation_cloud_fraction = mapl_get_pointer_placeholder("FCLD")
        self.outputs.radiation_vapor = mapl_get_pointer_placeholder("QV")
        self.outputs.radiation_liquid = mapl_get_pointer_placeholder("QL")
        self.outputs.radiation_ice = mapl_get_pointer_placeholder("QI")
        self.outputs.radiation_rain = mapl_get_pointer_placeholder("QR")
        self.outputs.radiation_snow = mapl_get_pointer_placeholder("QS")
        self.outputs.radiation_graupel = mapl_get_pointer_placeholder("QG")
        self.outputs.lower_tropospheric_stability = mapl_get_pointer_placeholder("LTS")
        self.outputs.estimated_inversion_strength = mapl_get_pointer_placeholder("EIS")
        self.outputs.z_lcl = mapl_get_resource_placeholder("ZLCL")

        # Outputs: model fields originating from within GFDL; macrophysics/microphysics tendencies
        self.outputs.du_dt_macro = mapl_get_pointer_placeholder("DUDT_macro")
        self.outputs.dv_dt_macro = mapl_get_pointer_placeholder("DVDT_macro")
        self.outputs.dt_dt_macro = mapl_get_pointer_placeholder("DTDT_macro")
        self.outputs.dvapor_dt_macro = mapl_get_pointer_placeholder("DQVDT_macro")
        self.outputs.dliquid_dt_macro = mapl_get_pointer_placeholder("DQLDT_macro")
        self.outputs.dice_dt_macro = mapl_get_pointer_placeholder("DQIDT_macro")
        self.outputs.dcloud_fraction_dt_macro = mapl_get_pointer_placeholder("DQADT_macro")
        self.outputs.drain_dt_macro = mapl_get_pointer_placeholder("DQRDT_macro")
        self.outputs.dsnow_dt_macro = mapl_get_pointer_placeholder("DQSDT_macro")
        self.outputs.dgraupel_dt_macro = mapl_get_pointer_placeholder("DQGDT_macro")
        self.outputs.du_dt_micro = mapl_get_pointer_placeholder("DUDT_micro")
        self.outputs.dv_dt_micro = mapl_get_pointer_placeholder("DVDT_micro")
        self.outputs.dt_dt_micro = mapl_get_pointer_placeholder("DTDT_micro")
        self.outputs.dvapor_dt_micro = mapl_get_pointer_placeholder("DQVDT_micro")
        self.outputs.dliquid_dt_micro = mapl_get_pointer_placeholder("DQLDT_micro")
        self.outputs.dice_dt_micro = mapl_get_pointer_placeholder("DQIDT_micro")
        self.outputs.dcloud_fraction_dt_micro = mapl_get_pointer_placeholder("DQADT_micro")
        self.outputs.drain_dt_micro = mapl_get_pointer_placeholder("DQRDT_micro")
        self.outputs.dsnow_dt_micro = mapl_get_pointer_placeholder("DQSDT_micro")
        self.outputs.dgraupel_dt_micro = mapl_get_pointer_placeholder("DQGDT_micro")

        # prepare macrophysics tendencies
        self.prepare_tendencies(
            u=u,
            v=v,
            t=t,
            vapor=mixing_ratios.vapor,
            rain=mixing_ratios.rain,
            snow=mixing_ratios.snow,
            graupel=mixing_ratios.graupel,
            convective_liquid=mixing_ratios.convective_liquid,
            convective_ice=mixing_ratios.convective_ice,
            large_scale_liquid=mixing_ratios.large_scale_liquid,
            large_scale_ice=mixing_ratios.large_scale_ice,
            convective_cloud_fraction=cloud_fractions.convective,
            large_scale_cloud_fraction=cloud_fractions.large_scale,
            du_dt=self.outputs.du_dt_macro,
            dv_dt=self.outputs.dv_dt_macro,
            dt_dt=self.outputs.dt_dt_macro,
            dvapor_dt=self.outputs.dvapor_dt_macro,
            dliquid_dt=self.outputs.dliquid_dt_macro,
            dice_dt=self.outputs.dice_dt_macro,
            dcloud_fraction_dt=self.outputs.dcloud_fraction_dt_macro,
            drain_dt=self.outputs.drain_dt_macro,
            dsnow_dt=self.outputs.dsnow_dt_macro,
            dgraupel_dt=self.outputs.dgraupel_dt_macro,
        )

        # prepare microphysics tendencies
        self.prepare_tendencies(
            u=u,
            v=v,
            t=t,
            vapor=mixing_ratios.vapor,
            rain=mixing_ratios.rain,
            snow=mixing_ratios.snow,
            graupel=mixing_ratios.graupel,
            convective_liquid=mixing_ratios.convective_liquid,
            convective_ice=mixing_ratios.convective_ice,
            large_scale_liquid=mixing_ratios.large_scale_liquid,
            large_scale_ice=mixing_ratios.large_scale_ice,
            convective_cloud_fraction=cloud_fractions.convective,
            large_scale_cloud_fraction=cloud_fractions.large_scale,
            du_dt=self.outputs.du_dt_micro,
            dv_dt=self.outputs.dv_dt_micro,
            dt_dt=self.outputs.dt_dt_micro,
            dvapor_dt=self.outputs.dvapor_dt_micro,
            dliquid_dt=self.outputs.dliquid_dt_micro,
            dice_dt=self.outputs.dice_dt_micro,
            dcloud_fraction_dt=self.outputs.dcloud_fraction_dt_micro,
            drain_dt=self.outputs.drain_dt_micro,
            dsnow_dt=self.outputs.dsnow_dt_micro,
            dgraupel_dt=self.outputs.dgraupel_dt_micro,
        )

        self.calculate_derived_states(
            p_interface=p_interface,
            p_interface_mb=self.temporaries.p_interface_mb,
            p_mb=self.temporaries.p_mb,
            geopotential_height_interface=geopotential_height_interface,
            edge_height_above_surface=self.temporaries.edge_height_above_surface,
            layer_height_above_surface=self.temporaries.layer_height_above_surface,
            layer_thickness=self.temporaries.layer_thickness,
            dp=self.temporaries.dp,
            mass=self.temporaries.mass,
            t=t,
            ese=self.saturation_tables.ese,
            esx=self.saturation_tables.esx,
            qsat=self.temporaries.qsat,
            dqsat=self.temporaries.dqsat,
            u=u,
            u_unmodified=self.temporaries.u_unmodified,
            v=v,
            v_unmodified=self.temporaries.v_unmodified,
        )

        self.find_klcl(
            t=t,
            p_mb=self.temporaries.p_mb,
            vapor=mixing_ratios.vapor,
            ese=self.saturation_tables.ese,
            esx=self.saturation_tables.esx,
            found_level=self.masks.boolean_2d_mask,
            k_lcl=self.temporaries.k_lcl,
        )

        if associated_checker("ZLCL") == True:
            self.outputs.z_lcl = self.temporaries.layer_height_above_surface.field[
                :, :, self.temporaries.k_lcl.field[:]
            ]  # bad but gets the point across: export height at lcl (should not execute in our test case)

        self.vertical_interpolation(
            field=self.temporaries.th,
            interpolated_field=self.temporaries.th700,
            p_interface_mb=self.temporaries.p_interface_mb,
            target_pressure=Float(70000.0),
            pb=self.temporaries.temporary_2d_1,
            pt=self.temporaries.temporary_2d_2,
            boolean_2d_mask=self.masks.boolean_2d_mask,
        )

        self.vertical_interpolation(
            field=t,
            interpolated_field=self.temporaries.t700,
            p_interface_mb=self.temporaries.p_interface_mb,
            target_pressure=Float(70000.0),
            pb=self.temporaries.temporary_2d_1,
            pt=self.temporaries.temporary_2d_2,
            boolean_2d_mask=self.masks.boolean_2d_mask,
        )
        self.vertical_interpolation(
            field=self.temporaries.layer_height_above_surface,
            interpolated_field=self.temporaries.z700,
            p_interface_mb=self.temporaries.p_interface_mb,
            target_pressure=Float(70000.0),
            pb=self.temporaries.temporary_2d_1,
            pt=self.temporaries.temporary_2d_2,
            boolean_2d_mask=self.masks.boolean_2d_mask,
        )

        self.find_eis(
            t=t,
            th=self.temporaries.th,
            layer_height_above_surface=self.temporaries.layer_height_above_surface,
            t700=self.temporaries.t700,
            th700=self.temporaries.th700,
            z700=self.temporaries.z700,
            k_lcl=self.temporaries.k_lcl,
            ese=self.saturation_tables.ese,
            esx=self.saturation_tables.esx,
            lower_tropospheric_stability=self.outputs.lower_tropospheric_stability,
            estimated_inversion_strength=self.outputs.estimated_inversion_strength,
        )

        if associated_checker("SHLW_PRC3") == True:
            mixing_ratios.rain = mixing_ratios.rain + shallow_convective_rain * self.DT_MOIST

        if associated_checker("SHLW_SNO3") == True:
            mixing_ratios.snow = mixing_ratios.snow + shallow_convective_snow * self.DT_MOIST

        self.phase_change(
            phase_change_config=self.phase_change_config,
            estimated_inversion_strength=self.outputs.estimated_inversion_strength,
            p_mb=self.temporaries.p_mb,
            k_lcl=self.temporaries.k_lcl,
            p_interface_mb=self.temporaries.p_interface_mb,
            area=area,
            convection_fraction=convection_fraction,
            surface_type=surface_type,
            t=t,
            convective_liquid=mixing_ratios.convective_liquid,
            convective_ice=mixing_ratios.convective_ice,
            large_scale_liquid=mixing_ratios.large_scale_liquid,
            large_scale_ice=mixing_ratios.large_scale_ice,
            vapor=mixing_ratios.vapor,
            large_scale_cloud_fraction=cloud_fractions.large_scale,
            convective_cloud_fraction=cloud_fractions.convective,
            nactl=liquid_concentration,
            nacti=ice_concentration,
            qsat=self.temporaries.qsat,
        )

        # pull rh_crit out of phase change component and send it back to the rest of the model
        rh_crit = self.phase_change.outputs.rh_crit

        self.update_macrophysics_tendencies(
            u=u,
            v=v,
            t=t,
            vapor=mixing_ratios.vapor,
            rain=mixing_ratios.rain,
            snow=mixing_ratios.snow,
            graupel=mixing_ratios.graupel,
            convective_liquid=mixing_ratios.convective_liquid,
            convective_ice=mixing_ratios.convective_ice,
            large_scale_liquid=mixing_ratios.large_scale_liquid,
            large_scale_ice=mixing_ratios.large_scale_ice,
            convective_cloud_fraction=cloud_fractions.convective,
            large_scale_cloud_fraction=cloud_fractions.large_scale,
            du_dt=self.outputs.du_dt_macro,
            dv_dt=self.outputs.dv_dt_macro,
            dt_dt=self.outputs.dt_dt_macro,
            dvapor_dt=self.outputs.dvapor_dt_macro,
            dliquid_dt=self.outputs.dliquid_dt_macro,
            dice_dt=self.outputs.dice_dt_macro,
            dcloud_fraction_dt=self.outputs.dcloud_fraction_dt_macro,
            drain_dt=self.outputs.drain_dt_macro,
            dsnow_dt=self.outputs.dsnow_dt_macro,
            dgraupel_dt=self.outputs.dgraupel_dt_macro,
        )

        self.prepare_radiation_quantities(
            convective_cloud_fraction=cloud_fractions.convective,
            large_scale_cloud_fraction=cloud_fractions.large_scale,
            radiation_cloud_fraction=self.outputs.radiation_cloud_fraction,
            convective_liquid=mixing_ratios.convective_liquid,
            large_scale_liquid=mixing_ratios.large_scale_liquid,
            radiation_liquid=self.outputs.radiation_liquid,
            convective_ice=mixing_ratios.convective_ice,
            large_scale_ice=mixing_ratios.large_scale_ice,
            radiation_ice=self.outputs.radiation_ice,
            vapor=mixing_ratios.vapor,
            radiation_vapor=self.outputs.radiation_vapor,
            rain=mixing_ratios.rain,
            radiation_rain=self.outputs.radiation_rain,
            snow=mixing_ratios.snow,
            radiation_snow=self.outputs.radiation_snow,
            graupel=mixing_ratios.graupel,
            radiation_graupel=self.outputs.radiation_graupel,
        )

        self.driver(
            GFDL_1M_config=self.GFDL_1M_config,
            t=t,
            w=vertical_motion.velocity,
            u=u,
            v=v,
            dz=self.temporaries.layer_thickness_negative,
            dp=self.temporaries.dp,
            area=area,
            land_fraction=land_fraction,
            convection_fraction=convection_fraction,
            surface_type=surface_type,
            estimated_inversion_strength=self.outputs.estimated_inversion_strength,
            rh_crit=rh_crit,
            vapor=self.outputs.radiation_vapor,
            liquid=self.outputs.radiation_liquid,
            rain=self.outputs.radiation_rain,
            ice=self.outputs.radiation_ice,
            snow=self.outputs.radiation_snow,
            graupel=self.outputs.radiation_graupel,
            cloud_fraction=self.outputs.radiation_cloud_fraction,
            ice_concentration=ice_concentration,
            liquid_concentration=liquid_concentration,
            dvapor_dt=self.temporaries.dvapor_dt,
            dliquid_dt=self.temporaries.dliquid_dt,
            drain_dt=self.temporaries.drain_dt,
            dice_dt=self.temporaries.dice_dt,
            dsnow_dt=self.temporaries.dsnow_dt,
            dgraupel_dt=self.temporaries.dgraupel_dt,
            dcloud_fraction_dt=self.temporaries.dcloud_fraction_dt,
            dt_dt=self.temporaries.dt_dt,
            du_dt=self.temporaries.du_dt,
            dv_dt=self.temporaries.dv_dt,
            anv_icefall=self.ANV_ICEFALL,
            ls_icefall=self.LS_ICEFALL,
        )

        self.apply_driver_tendencies(
            t=t,
            u=u,
            v=v,
            radiation_cloud_fraction=self.outputs.radiation_cloud_fraction,
            radiation_ice=self.outputs.radiation_ice,
            radiation_liquid=self.outputs.radiation_liquid,
            radiation_vapor=self.outputs.radiation_vapor,
            radiation_rain=self.outputs.radiation_rain,
            radiation_snow=self.outputs.radiation_snow,
            radiation_graupel=self.outputs.radiation_graupel,
            dcloud_fraction_dt=self.temporaries.dcloud_fraction_dt,
            dt_dt=self.temporaries.dt_dt,
            du_dt=self.temporaries.du_dt,
            dv_dt=self.temporaries.dv_dt,
            dice_dt=self.temporaries.dice_dt,
            dliquid_dt=self.temporaries.dliquid_dt,
            dvapor_dt=self.temporaries.dvapor_dt,
            drain_dt=self.temporaries.drain_dt,
            dsnow_dt=self.temporaries.dsnow_dt,
            dgraupel_dt=self.temporaries.dgraupel_dt,
        )

        self.redistribute_clouds(
            cloud_fraction=self.outputs.radiation_cloud_fraction,
            convective_cloud_fraction=cloud_fractions.convective,
            large_scale_cloud_fraction=cloud_fractions.large_scale,
            liquid=self.outputs.radiation_liquid,
            convective_liquid=mixing_ratios.convective_liquid,
            large_scale_liquid=mixing_ratios.large_scale_liquid,
            ice=self.outputs.radiation_ice,
            convective_ice=mixing_ratios.convective_ice,
            large_scale_ice=mixing_ratios.large_scale_ice,
            vapor=mixing_ratios.vapor,
            temperature=t,
        )
