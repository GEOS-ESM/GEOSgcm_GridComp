from pyMoist.GFDL_1M.driver.driver import MicrophysicsDriver
from pyMoist.GFDL_1M.PhaseChange.phase_change import PhaseChange
from ndsl import QuantityFactory, StencilFactory
from pyMoist.GFDL_1M.driver.config import MicrophysicsConfiguration
from pyMoist.GFDL_1M.getters_temporary import mapl_placeholder, esmf_placeholder
from ndsl.dsl.typing import Float
from pyMoist.GFDL_1M.state import (
    LiquidWaterStaticEnergy,
    TotalWater,
    VericalMotion,
    MixingRatios,
    CloudFractions,
)


class GFDL1M:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        GFDL_1M_config: MicrophysicsConfiguration,
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
    ):
        PHYS_HYDROSTATIC = mapl_placeholder("LPHYS_HYDROSTATIC")
        HYDROSTATIC = mapl_placeholder("LHYDROSTATIC")
        DT_MOIST = esmf_placeholder("DT_R8")
        MELTFRZ = mapl_placeholder("MELTFRZ")

        self.GFDL_1M_config = MicrophysicsConfiguration(
            PHYS_HYDROSTATIC=PHYS_HYDROSTATIC,
            HYDROSTATIC=HYDROSTATIC,
            DT_MOIST=DT_MOIST,
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

    def __call__(self, *args, **kwds):
        
        # Get model state
        mixing_ratios = MixingRatios(
            vapor=mapl_placeholder("Q"),
            rain=mapl_placeholder("QRAIN"),
            snow=mapl_placeholder("QSNOW"),
            graupel=mapl_placeholder("QGRAUPEL"),
            convective_liquid=mapl_placeholder("QLLS"),
            convective_ice=mapl_placeholder("QILS"),
            large_scale_liquid=mapl_placeholder("QLCN"),
            large_scale_ice=mapl_placeholder("QICN"),
        )
        cloud_fractions = CloudFractions(
            convective=mapl_placeholder("CLCN"), large_scale=mapl_placeholder("CLLS")
        )
        nactl = mapl_placeholder("NACTL")
        nacti = mapl_placeholder("NACTI")


        area = mapl_placeholder("AREA")
        zle = mapl_placeholder("ZLE")
        ple = mapl_placeholder("PLE")
        t = mapl_placeholder("T")
        u = mapl_placeholder("U")
        v = mapl_placeholder("V")
        land_fraction = mapl_placeholder("FRLAND")
        vertical_motion = VericalMotion(
            velocity=mapl_placeholder("W"),
            variance=mapl_placeholder("W2"),
            third_moment=mapl_placeholder("W3"),
        )
        liquid_water_energy_flux = LiquidWaterStaticEnergy(
            mapl_placeholder("WSL"),
            mapl_placeholder("SL2"),
            mapl_placeholder("SL3"),
        )
        total_water = TotalWater(
            mapl_placeholder("WQT"),
            mapl_placeholder("QT2"),
            mapl_placeholder("QT3"),
        )
        EVAP = mapl_placeholder("EVAP")
        OMEGA = mapl_placeholder("OMEGA")
