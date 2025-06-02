from pyMoist.GFDL_1M.driver.driver import MicrophysicsDriver
from pyMoist.GFDL_1M.PhaseChange.phase_change import PhaseChange
from ndsl import QuantityFactory, StencilFactory
from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.GFDL_1M.getters_temporary import (
    mapl_get_resource_placeholder,
    esmf_placeholder,
)
from ndsl.dsl.typing import Float
from pyMoist.GFDL_1M.state import (
    LiquidWaterStaticEnergy,
    TotalWater,
    VericalMotion,
    MixingRatios,
    CloudFractions,
)
from pyMoist.GFDL_1M.stencils import calculate_derived_states, find_klcl
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.GFDL_1M.temporaries import Temporaries
from pyMoist.GFDL_1M.masks import Masks


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
    ):
        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory

        PHYS_HYDROSTATIC = mapl_get_resource_placeholder("LPHYS_HYDROSTATIC")
        HYDROSTATIC = mapl_get_resource_placeholder("LHYDROSTATIC")
        DT_MOIST = esmf_placeholder("DT_R8")
        MELTFRZ = mapl_get_resource_placeholder("MELTFRZ")

        self.GFDL_1M_config = GFDL1MConfig(
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

        # Initalize internal fields
        self.temporaries = Temporaries.make(quantity_factory=quantity_factory)
        self.masks = Masks.make(quantity_factory=quantity_factory)

        # Construct stencils
        self.calculate_derived_states = stencil_factory.from_dims_halo(
            func=calculate_derived_states,
            compute_dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM],
        )

        self.find_klcl = stencil_factory.from_dims_halo(
            func=find_klcl,
            compute_dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM],
        )

    def __call__(self, *args, **kwds):

        # Get model state
        mixing_ratios = MixingRatios(
            vapor=mapl_get_resource_placeholder("Q"),
            rain=mapl_get_resource_placeholder("QRAIN"),
            snow=mapl_get_resource_placeholder("QSNOW"),
            graupel=mapl_get_resource_placeholder("QGRAUPEL"),
            convective_liquid=mapl_get_resource_placeholder("QLLS"),
            convective_ice=mapl_get_resource_placeholder("QILS"),
            large_scale_liquid=mapl_get_resource_placeholder("QLCN"),
            large_scale_ice=mapl_get_resource_placeholder("QICN"),
        )
        cloud_fractions = CloudFractions(
            convective=mapl_get_resource_placeholder("CLCN"),
            large_scale=mapl_get_resource_placeholder("CLLS"),
        )
        nactl = mapl_get_resource_placeholder("NACTL")
        nacti = mapl_get_resource_placeholder("NACTI")
        area = mapl_get_resource_placeholder("AREA")
        geopotential_height_interface = mapl_get_resource_placeholder("ZLE")
        p_interface = mapl_get_resource_placeholder("PLE")
        t = mapl_get_resource_placeholder("T")
        u = mapl_get_resource_placeholder("U")
        v = mapl_get_resource_placeholder("V")
        land_fraction = mapl_get_resource_placeholder("FRLAND")
        vertical_motion = VericalMotion(
            velocity=mapl_get_resource_placeholder("W"),
            variance=mapl_get_resource_placeholder("W2"),
            third_moment=mapl_get_resource_placeholder("W3"),
        )
        liquid_water_static_energy = LiquidWaterStaticEnergy(
            flux=mapl_get_resource_placeholder("WSL"),
            variance=mapl_get_resource_placeholder("SL2"),
            third_moment=mapl_get_resource_placeholder("SL3"),
        )
        total_water = TotalWater(
            flux=mapl_get_resource_placeholder("WQT"),
            variance=mapl_get_resource_placeholder("QT2"),
            third_moment=mapl_get_resource_placeholder("QT3"),
        )

        # Initalize saturation tables
        saturation_tables = SaturationVaporPressureTable(self.stencil_factory.backend)

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
            ese=saturation_tables.ese,
            esx=saturation_tables.esx,
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
            ese=saturation_tables.ese,
            esx=saturation_tables.esx,
            found_level=self.masks.boolean_2d_mask,
            k_lcl=self.temporaries.k_lcl,
        )
