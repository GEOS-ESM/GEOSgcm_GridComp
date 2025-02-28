from ndsl.dsl.typing import Float


class MicrophysicsConfiguration:
    def __init__(
        self,
        PHYS_HYDROSTATIC: bool,
        HYDROSTATIC: bool,
        DT_MOIST: Float,
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
        DO_QA: bool,
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
        self.PHYS_HYDROSTATIC = PHYS_HYDROSTATIC
        self.HYDROSTATIC = HYDROSTATIC
        self.DT_MOIST = DT_MOIST
        self.MP_TIME = MP_TIME
        self.T_MIN = T_MIN
        self.T_SUB = T_SUB
        self.TAU_R2G = TAU_R2G
        self.TAU_SMLT = TAU_SMLT
        self.TAU_G2R = TAU_G2R
        self.DW_LAND = DW_LAND
        self.DW_OCEAN = DW_OCEAN
        self.VI_FAC = VI_FAC
        self.VR_FAC = VR_FAC
        self.VS_FAC = VS_FAC
        self.VG_FAC = VG_FAC
        self.QL_MLT = QL_MLT
        self.DO_QA = DO_QA
        self.FIX_NEGATIVE = FIX_NEGATIVE
        self.VI_MAX = VI_MAX
        self.VS_MAX = VS_MAX
        self.VG_MAX = VG_MAX
        self.VR_MAX = VR_MAX
        self.QS_MLT = QS_MLT
        self.QS0_CRT = QS0_CRT
        self.QI_GEN = QI_GEN
        self.QL0_MAX = QL0_MAX
        self.QI0_MAX = QI0_MAX
        self.QI0_CRT = QI0_CRT
        self.QR0_CRT = QR0_CRT
        self.FAST_SAT_ADJ = FAST_SAT_ADJ
        self.RH_INC = RH_INC
        self.RH_INS = RH_INS
        self.RH_INR = RH_INR
        self.CONST_VI = CONST_VI
        self.CONST_VS = CONST_VS
        self.CONST_VG = CONST_VG
        self.CONST_VR = CONST_VR
        self.USE_CCN = USE_CCN
        self.RTHRESHU = RTHRESHU
        self.RTHRESHS = RTHRESHS
        self.CCN_L = CCN_L
        self.CCN_O = CCN_O
        self.QC_CRT = QC_CRT
        self.TAU_G2V = TAU_G2V
        self.TAU_V2G = TAU_V2G
        self.TAU_S2V = TAU_S2V
        self.TAU_V2S = TAU_V2S
        self.TAU_REVP = TAU_REVP
        self.TAU_FRZ = TAU_FRZ
        self.DO_BIGG = DO_BIGG
        self.DO_EVAP = DO_EVAP
        self.DO_SUBL = DO_SUBL
        self.SAT_ADJ0 = SAT_ADJ0
        self.C_PIACR = C_PIACR
        self.TAU_IMLT = TAU_IMLT
        self.TAU_V2L = TAU_V2L
        self.TAU_L2V = TAU_L2V
        self.TAU_I2V = TAU_I2V
        self.TAU_I2S = TAU_I2S
        self.TAU_L2R = TAU_L2R
        self.QI_LIM = QI_LIM
        self.QL_GEN = QL_GEN
        self.C_PAUT = C_PAUT
        self.C_PSACI = C_PSACI
        self.C_PGACS = C_PGACS
        self.C_PGACI = C_PGACI
        self.Z_SLOPE_LIQ = Z_SLOPE_LIQ
        self.Z_SLOPE_ICE = Z_SLOPE_ICE
        self.PROG_CCN = PROG_CCN
        self.C_CRACW = C_CRACW
        self.ALIN = ALIN
        self.CLIN = CLIN
        self.PRECIPRAD = PRECIPRAD
        self.CLD_MIN = CLD_MIN
        self.USE_PPM = USE_PPM
        self.MONO_PROF = MONO_PROF
        self.DO_SEDI_HEAT = DO_SEDI_HEAT
        self.SEDI_TRANSPORT = SEDI_TRANSPORT
        self.DO_SEDI_W = DO_SEDI_W
        self.DE_ICE = DE_ICE
        self.ICLOUD_F = ICLOUD_F
        self.IRAIN_F = IRAIN_F
        self.MP_PRINT = MP_PRINT
