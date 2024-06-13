from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField, Int, Float 
from ndsl import Quantity, QuantityFactory, StencilFactory, orchestrate
import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import computation, interval, PARALLEL, log10, exp  # type: ignore
import pyMoist.constants as radconstants

#Calculate air density [kg/m^3]
@gtscript.function
def RHO(PL, TE):
    RHO = (100.0 * PL) / (radconstants.MAPL_RGAS * TE)
    return RHO

#Calculate ice cloud effective radius [m]
@gtscript.function
def LDRADIUS4_ICE(
    PL: Float,
    TE: Float, 
    QC: Float, 
    NNL: Float,
    NNI: Float,
)-> Float:
    #Calculate ice water content
    WC = 1.0e3 * RHO(PL, TE) * QC #air density [g/m3] * ice cloud mixing ratio [kg/kg]
    #Calculate radius in meters [m]
    if radconstants.ICE_RADII_PARAM == 1:
        #Ice cloud effective radius -- [klaus wyser, 1998]
        if TE > radconstants.MAPL_TICE or QC <= 0.0:
            BB = -2.0
        else:
            BB = -2.0 + log10(WC / 50.0) * (1.0e-3 * (radconstants.MAPL_TICE - TE)**1.5)
        BB = min(max(BB, -6.0), -2.0)
        RADIUS = 377.4 + 203.3 * BB + 37.91 * BB**2 + 2.3696 * BB**3
        RADIUS = min(150.0e-6, max(5.0e-6, 1.0e-6 * RADIUS))
    else: 
        #Ice cloud effective radius ----- [Sun, 2001]
        TC = TE - radconstants.MAPL_TICE
        ZFSR = 1.2351 + 0.0105 * TC
        AA = 45.8966 * (WC**0.2214)
        BB = 0.79570 * (WC**0.2535)
        RADIUS = ZFSR * (AA + BB * (TE - 83.15))
        RADIUS = min(150.0e-6, max(5.0e-6, 1.0e-6 * RADIUS * 0.64952))
    return RADIUS

#Calculate liquid cloud effective radius [m]
@gtscript.function
def LDRADIUS4_LIQUID(
    PL: Float,
    TE: Float, 
    QC: Float, 
    NNL: Float,
    NNI: Float,
)-> Float:
    #Calculate liquid water content
    WC = 1.0e3 * RHO(PL,TE) * QC #air density [g/m3] * liquid cloud mixing ratio [kg/kg]
    #Calculate cloud drop number concentration from the aerosol model + ....
    NNX = max(NNL * 1.0e-6, 10.0)
    #Calculate Radius in meters [m]
    if radconstants.LIQ_RADII_PARAM == 1:
        #Jason Version
        RADIUS = min(60.0e-6, max(2.5e-6, 1.0e-6 * radconstants.bx * (WC / NNX)**radconstants.r13bbeta * radconstants.abeta * 6.92))
    else:
        #[liu&daum, 2000 and 2005. liu et al 2008]
        RADIUS = min(60.0e-6, max(2.5e-6, 1.0e-6 * radconstants.Lbx * (WC / NNX)**radconstants.Lbe))
    return RADIUS

def _fix_up_clouds_stencil(
    QV: FloatField,
    TE: FloatField,
    QLC: FloatField, 
    QIC: FloatField,
    CF: FloatField, 
    QLA: FloatField, 
    QIA: FloatField,
    AF: FloatField, 
):
    with computation(PARALLEL), interval(...):
        #Fix if Anvil cloud fraction too small
        if AF < 1.E-5:
            QV = QV + QLA + QIA
            TE = TE - (radconstants.alhlbcp) * QLA - (radconstants.alhsbcp) * QIA
            AF = 0.0
            QLA = 0.0
            QIA = 0.0
        #Fix if LS cloud fraction too small
        if CF < 1.E-5:
            QV = QV + QLC + QIC
            TE = TE - (radconstants.alhlbcp) * QLC - (radconstants.alhsbcp) * QIC
            CF = 0.0
            QLC = 0.0
            QIC = 0.0
        #LS LIQUID too small
        if QLC < 1.E-8:
            QV = QV + QLC
            TE = TE - (radconstants.alhlbcp) * QLC
            QLC = 0.0
        #LS ICE too small
        if QIC < 1.E-8:
            QV = QV + QIC
            TE = TE - (radconstants.alhsbcp) * QIC
            QIC = 0.0
        #Anvil LIQUID too small
        if QLA < 1.E-8:
            QV = QV + QLA
            TE = TE - (radconstants.alhlbcp) * QLA
            QLA = 0.0
        #Anvil ICE too small
        if QIA < 1.E-8:
            QV = QV + QIA
            TE = TE - (radconstants.alhsbcp) * QIA
            QIA = 0.0
        #Fix ALL cloud quants if Anvil cloud LIQUID+ICE too small
        if (QLA + QIA) < 1.E-8:
            QV = QV + QLA + QIA
            TE = TE - (radconstants.alhlbcp) * QLA - (radconstants.alhsbcp) * QIA
            AF = 0.0
            QLA = 0.0
            QIA = 0.0
        #Fix ALL cloud quants if LS cloud LIQUID+ICE too small
        if (QLC + QIC) < 1.E-8:
            QV = QV + QLC + QIC
            TE = TE - (radconstants.alhlbcp) * QLC - (radconstants.alhsbcp) * QIC
            CF = 0.0
            QLC = 0.0
            QIC = 0.0

def _radcouple_stencil(
        TE: FloatField,
        PL: FloatField,
        CF: FloatField, 
        AF: FloatField, 
        QV: FloatField, 
        QClLS: FloatField, 
        QCiLS: FloatField, 
        QClAN: FloatField, 
        QCiAN: FloatField, 
        QRN_ALL: FloatField, 
        QSN_ALL: FloatField, 
        QGR_ALL: FloatField, 
        NL: FloatField, 
        NI: FloatField, 
        RAD_QV: FloatField, 
        RAD_QL: FloatField, 
        RAD_QI: FloatField, 
        RAD_QR: FloatField, 
        RAD_QS: FloatField, 
        RAD_QG: FloatField, 
        RAD_CF: FloatField, 
        RAD_RL: FloatField, 
        RAD_RI: FloatField, 
        FAC_RL: Float,
        MIN_RL: Float,
        MAX_RL: Float,
        FAC_RI: Float,
        MIN_RI: Float,
        MAX_RI: Float,
):
    with computation(PARALLEL), interval(...):
        #water vapor
        RAD_QV = QV

        #total cloud fraction
        RAD_CF = max(min(CF + AF, 1.0), 0.0)
        if RAD_CF >= 1.e-5:
            RAD_QL = (QClLS + QClAN) / RAD_CF if (QClLS + QClAN) >= 1.e-8 else 0.0
            RAD_QI = (QCiLS + QCiAN) / RAD_CF if (QCiLS + QCiAN) >= 1.e-8 else 0.0
            RAD_QR = QRN_ALL / RAD_CF if QRN_ALL >= 1.e-8 else 0.0
            RAD_QS = QSN_ALL / RAD_CF if QSN_ALL >= 1.e-8 else 0.0
            RAD_QG = QGR_ALL / RAD_CF if QGR_ALL >= 1.e-8 else 0.0
        else:
            RAD_CF = 0.0
            RAD_QL = 0.0
            RAD_QI = 0.0
            RAD_QR = 0.0
            RAD_QS = 0.0
            RAD_QG = 0.0

        #Cap the high end of condensates
        RAD_QL = min(RAD_QL, 0.01)
        RAD_QI = min(RAD_QI, 0.01)
        RAD_QR = min(RAD_QR, 0.01)
        RAD_QS = min(RAD_QS, 0.01)
        RAD_QG = min(RAD_QG, 0.01)

        #Liquid radii - Brams formulation with limits
        RAD_RL = max(MIN_RL, min(LDRADIUS4_LIQUID(PL, TE, RAD_QL, NL, NI) * FAC_RL, MAX_RL))
        #Ice radii - Brams formulation with limits
        RAD_RI = max(MIN_RI, min(LDRADIUS4_ICE(PL, TE, RAD_QI, NL, NI) * FAC_RI, MAX_RI))

class RadiationCoupling:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        do_qa: bool,
    ):
        self._fix_up_clouds = stencil_factory.from_dims_halo(
            func=_fix_up_clouds_stencil, 
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )
        self._radcouple = stencil_factory.from_dims_halo(
            func=_radcouple_stencil, 
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )
        self.do_qa = do_qa

    def __call__(
        self,
        Q: FloatField,
        T: FloatField,
        QLLS: FloatField,
        QILS: FloatField,
        CLLS: FloatField,
        QLCN: FloatField,
        QICN: FloatField,
        CLCN: FloatField,
        PLmb: FloatField,
        QRAIN: FloatField,
        QSNOW: FloatField,
        QGRAUPEL: FloatField,
        NACTL: FloatField,
        NACTI: FloatField,
        RAD_QV: FloatField,
        RAD_QL: FloatField,
        RAD_QI: FloatField,
        RAD_QR: FloatField,
        RAD_QS: FloatField,
        RAD_QG: FloatField,
        RAD_CF: FloatField,
        CLDREFFL: FloatField,
        CLDREFFI: FloatField,
        FAC_RL: Float,
        MIN_RL: Float,
        MAX_RL: Float,
        FAC_RI: Float,
        MIN_RI: Float,
        MAX_RI: Float,
    ):
        self._fix_up_clouds(QV = Q, 
                            TE = T, 
                            QLC = QLLS, 
                            QIC = QILS,
                            CF = CLLS, 
                            QLA = QLCN,  
                            QIA = QICN, 
                            AF = CLCN,
                            )
        self._radcouple(TE = T, 
                        PL = PLmb,
                        CF = CLLS, 
                        AF = CLCN, 
                        QV = Q, 
                        QClLS = QLLS, 
                        QCiLS = QILS, 
                        QClAN = QLCN, 
                        QCiAN = QICN, 
                        QRN_ALL = QRAIN, 
                        QSN_ALL = QSNOW, 
                        QGR_ALL = QGRAUPEL, 
                        NL = NACTL, 
                        NI = NACTI, 
                        RAD_QV = RAD_QV, 
                        RAD_QL = RAD_QL, 
                        RAD_QI = RAD_QI, 
                        RAD_QR = RAD_QR, 
                        RAD_QS = RAD_QS, 
                        RAD_QG = RAD_QG, 
                        RAD_CF = RAD_CF, 
                        RAD_RL = CLDREFFL, 
                        RAD_RI = CLDREFFI, 
                        FAC_RL = FAC_RL, 
                        MIN_RL = MIN_RL, 
                        MAX_RL = MAX_RL, 
                        FAC_RI = FAC_RI, 
                        MIN_RI = MIN_RI, 
                        MAX_RI = MAX_RI,
                        )
        '''
        From GEOS_GFDL_1M_InterfaceMod.F90 Line 866 not implemented
        Implement QSAT0 and QSAT3 logic if diagnostics are needed. 
        QSAT0 and QSAT3 are found in GEOS/src/Shared/@GMAO_Shared/GEOS_Shared/GEOS_Utilities.F90 
        
        Variables: RHX - FloatField
        '''
        if self.do_qa:
            #RHX = Q/GEOS_QSAT( T, PLmb)
            raise NotImplementedError("Module procedures QSAT0 and QSAT3 that are needed for GEOS_QSAT are not implemented")