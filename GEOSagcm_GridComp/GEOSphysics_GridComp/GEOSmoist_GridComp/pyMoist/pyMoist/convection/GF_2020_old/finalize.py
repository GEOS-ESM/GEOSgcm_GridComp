from ndsl import StencilFactory
from ndsl.dsl.gt4py import PARALLEL, interval, computation, FORWARD, sqrt, max, min, abs, floor
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import FloatField, FloatFieldIJ
from pyMoist.convective_parameterization.GF_2020.config import GF2020Config
import pyMoist.constants as constants
from pyMoist.convective_parameterization.GF_2020.temporaries import Temporaries
from pyMoist.convective_parameterization.GF_2020.state import MixingRatios, CloudFractions
from pyMoist.saturation_tables.qsat_functions import saturation_specific_humidity
from pyMoist.saturation_tables.types import GlobalTable_saturation_tables
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.shared_incloud_processes import ice_fraction


def finalize(
    U: FloatField,
    V: FloatField,
    Q: FloatField,
    T: FloatField,
    MASS: FloatField,
    TH: FloatField,
    PK: FloatField,
    DUDT_DC: FloatField,
    DVDT_DC: FloatField,
    DQVDT_DC: FloatField,
    DTDT_DC: FloatField,
    DQLDT_DC: FloatField,
    DQIDT_DC: FloatField,
    DQADT_DC: FloatField,
    CNV_DQCDT: FloatField,
    CNV_FRC: FloatFieldIJ,
    SRF_TYPE: FloatFieldIJ,
    MFD_DC: FloatField,
    REVSU: FloatField,
    RSU_CN: FloatField,
    REV_CN: FloatField,
    PRFIL: FloatField,
    PFI_CN: FloatField,
    PFL_CN: FloatField,
    QLCN: FloatField,
    QICN: FloatField,
    CLCN: FloatField,
    PL: FloatField,
    ese: GlobalTable_saturation_tables,
    esx: GlobalTable_saturation_tables,
):
    """
    This stencil MUST be built using Z_INTERFACE_DIM to function properly.
    """
    from __externals__ import DT_MOIST, SCLM_DEEP, FIX_CNV_CLOUD

    with computation(PARALLEL), interval(0, -1):
        # add tendencies to the moist import state
        U = U + DUDT_DC * DT_MOIST
        V = V + DVDT_DC * DT_MOIST
        Q = Q + DQVDT_DC * DT_MOIST
        T = T + DTDT_DC * DT_MOIST
        TH = T / PK
        # update DeepCu QL/QI/CF tendencies
        fQi = ice_fraction(T, CNV_FRC, SRF_TYPE)
        TMP3D = CNV_DQCDT / MASS
        DQLDT_DC = (1.0 - fQi) * TMP3D
        DQIDT_DC = fQi * TMP3D
        DQADT_DC = MFD_DC * SCLM_DEEP / MASS
        # sublimation/evaporation tendencies (kg/kg/s)
        RSU_CN = REVSU * fQi
        REV_CN = REVSU * (1.0 - fQi)
        # # add QI/QL/CL tendencies
        QLCN = QLCN + DQLDT_DC * DT_MOIST
        QICN = QICN + DQIDT_DC * DT_MOIST
        CLCN = max(min(CLCN + DQADT_DC * DT_MOIST, 1.0), 0.0)

    with computation(PARALLEL), interval(1, None):
        # preciptation fluxes (kg/kg/s)
        PFI_CN = PRFIL * fQi[0, 0, -1]
        PFL_CN = PRFIL * (1.0 - fQi[0, 0, -1])

    with computation(PARALLEL), interval(0, -1):
        # fix 'convective' cloud fraction
        if FIX_CNV_CLOUD == True:
            # fix convective cloud
            QST3, _ = saturation_specific_humidity(T, PL, ese, esx)
            TMP3D = QST3
            if CLCN < 1.0:
                TMP3D = (Q - QST3 * CLCN) / (1.0 - CLCN)
            minrhx = 0.001
            if ((TMP3D - minrhx * QST3) < 0.0) and (CLCN > 0.0):
                CLCN = (Q - minrhx * QST3) / (QST3 * (1.0 - minrhx))
            # If still cant make suitable env RH then destroy anvil
            if CLCN < 0.0:
                CLCN = 0.0
                DQLDT_DC = DQLDT_DC - (QLCN) / DT_MOIST
                DQIDT_DC = DQIDT_DC - (QICN) / DT_MOIST
                DQVDT_DC = DQVDT_DC + (QLCN + QICN) / DT_MOIST
                Q = Q + (QLCN + QICN)
                TMP3D = (constants.MAPL_ALHL * QLCN + constants.MAPL_ALHS * QICN) / constants.MAPL_CP
                DTDT_DC = DTDT_DC - TMP3D / DT_MOIST
                T = T - TMP3D
                TH = T / PK
                QLCN = 0.0
                QICN = 0.0

    # NOTE need to figure out how to do this with the proper associated checks
    #   ! Export
    #   call MAPL_GetPointer(EXPORT, PTR3D, 'CNV_FICE', RC=STATUS); VERIFY_(STATUS)
    #   if (associated(PTR3D)) PTR3D = fQi
    #   !--------------------------------------------------------------
    #   !  For Now add DeepCu contribution to total/detraining mass flux exports
    #   !--------------------------------------------------------------
    #   call MAPL_GetPointer(EXPORT, PTR3D, 'CNV_MFC', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    #   PTR3D = PTR3D + UMF_DC
    #   call MAPL_GetPointer(EXPORT, PTR3D, 'CNV_MFD', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    #   PTR3D = PTR3D + MFD_DC

    #   call MAPL_GetPointer(EXPORT, PTR3D, 'DQRC', RC=STATUS); VERIFY_(STATUS)
    #   if(associated(PTR3D)) PTR3D = CNV_PRC3 / GF_DT


class Finalize:
    def __init__(self, stencil_factory: StencilFactory, GF_2020_config: GF2020Config):
        self.finalize = stencil_factory.from_dims_halo(
            func=finalize,
            compute_dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM],
            externals={
                "DT_MOIST": GF_2020_config.DT_MOIST,
                "SCLM_DEEP": GF_2020_config.SCLM_DEEP,
                "FIX_CNV_CLOUD": GF_2020_config.FIX_CNV_CLOUD,
            },
        )

    def __call__(
        self,
        U,
        V,
        T,
        MASS,
        PK,
        DUDT_DC,
        DVDT_DC,
        DQVDT_DC,
        DTDT_DC,
        DQLDT_DC,
        DQIDT_DC,
        DQADT_DC,
        CNV_DQCDT,
        CNV_FRC,
        SRF_TYPE,
        MFD_DC,
        REVSU,
        RSU_CN,
        REV_CN,
        PRFIL,
        PFI_CN,
        PFL_CN,
        PL,
        mixing_ratios: MixingRatios,
        cloud_fractions: CloudFractions,
        temporaries: Temporaries,
        saturation_tables: SaturationVaporPressureTable,
    ):
        self.finalize(
            U,
            V,
            mixing_ratios.vapor,
            T,
            MASS,
            temporaries.th,
            PK,
            DUDT_DC,
            DVDT_DC,
            DQVDT_DC,
            DTDT_DC,
            DQLDT_DC,
            DQIDT_DC,
            DQADT_DC,
            CNV_DQCDT,
            CNV_FRC,
            SRF_TYPE,
            MFD_DC,
            REVSU,
            RSU_CN,
            REV_CN,
            PRFIL,
            PFI_CN,
            PFL_CN,
            mixing_ratios.convective_liquid,
            mixing_ratios.convective_ice,
            cloud_fractions.convective,
            PL,
            saturation_tables.ese,
            saturation_tables.esx,
        )
