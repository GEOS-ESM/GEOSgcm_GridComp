from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField, Int 
from ndsl import Quantity, QuantityFactory, StencilFactory, orchestrate
import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import computation, interval, PARALLEL, log10, exp  # type: ignore
import pyMoist.constants as radconstants

# GEOS_QSAT Function --This might be wrong
def GEOS_QSAT(
    T: FloatField, 
    P: FloatField,
):
    Rv = 461.5  # Gas constant for water vapor [J/(kg·K)]
    Lv = 2.5e6  # Latent heat of vaporization [J/kg]
    es0 = 610.78  # Saturation vapor pressure at 0°C [Pa]
    
    # Saturation vapor pressure formula (Clausius-Clapeyron equation)
    es = es0 * exp(Lv / Rv * (1/273.15 - 1/T))
    
    # Saturation mixing ratio
    qsat = 0.622 * es / (P - es)
    
    return qsat

#@gtscript.function
def LDRADIUS4(
    PL: FloatField,
    TE: FloatField, 
    QC: FloatField, 
    NNL: FloatField,
    NNI: FloatField,
    ITYPE: Int,
)-> FloatField:
    RHO = (100.0 * PL) / (radconstants.MAPL_RGAS * TE)
    if ITYPE == 1:  # Liquid
        WC = 1.0e3 * RHO * QC
        NNX = max(NNL * 1.0e-6, 10.0)
        if radconstants.LIQ_RADII_PARAM == 1:
            RADIUS = min(60.0e-6, max(2.5e-6, 1.0e-6 * radconstants.bx * (WC / NNX)**radconstants.r13bbeta * radconstants.abeta * 6.92))
        else:
            RADIUS = min(60.0e-6, max(2.5e-6, 1.0e-6 * radconstants.Lbx * (WC / NNX)**radconstants.Lbe))
    elif ITYPE == 2:  # Ice
        WC = 1.0e3 * RHO * QC
        if radconstants.ICE_RADII_PARAM == 1:
            if TE > radconstants.MAPL_TICE or QC <= 0.0:
                BB = -2.0
            else:
                BB = -2.0 + log10(WC / 50.0) * (1.0e-3 * (radconstants.MAPL_TICE - TE)**1.5)
            BB = min(max(BB, -6.0), -2.0)
            RADIUS = 377.4 + 203.3 * BB + 37.91 * BB**2 + 2.3696 * BB**3
            RADIUS = min(150.0e-6, max(5.0e-6, 1.0e-6 * RADIUS))
        else:
            TC = TE - radconstants.MAPL_TICE
            ZFSR = 1.2351 + 0.0105 * TC
            AA = 45.8966 * (WC**0.2214)
            BB = 0.79570 * (WC**0.2535)
            RADIUS = ZFSR * (AA + BB * (TE - 83.15))
            RADIUS = min(150.0e-6, max(5.0e-6, 1.0e-6 * RADIUS * 0.64952))
    else:
        raise ValueError("WRONG HYDROMETEOR type: CLOUD = 1 OR ICE = 2")
    return RADIUS

def _fix_up_clouds_stencil(
    QV: FloatField,
    TE: FloatField,
    QLC: FloatField, 
    CF: FloatField, 
    QLA: FloatField, 
    AF: FloatField, 
    QIC: FloatField, 
    QIA: FloatField,
):
    '''
    from __externals__ import (
    )
    '''
    with computation(PARALLEL), interval(...):
        if AF < 1.E-5:
            QV = QV + QLA + QIA
            TE = TE - (radconstants.alhlbcp) * QLA - (radconstants.alhsbcp) * QIA
            AF = 0.0
            QLA = 0.0
            QIA = 0.0

        if CF < 1.E-5:
            QV = QV + QLC + QIC
            TE = TE - (radconstants.alhlbcp) * QLC - (radconstants.alhsbcp) * QIC
            CF = 0.0
            QLC = 0.0
            QIC = 0.0

        if QLC < 1.E-8:
            QV = QV + QLC
            TE = TE - (radconstants.alhlbcp) * QLC
            QLC = 0.0

        if QIC < 1.E-8:
            QV = QV + QIC
            TE = TE - (radconstants.alhsbcp) * QIC
            QIC = 0.0

        if QLA < 1.E-8:
            QV = QV + QLA
            TE = TE - (radconstants.alhlbcp) * QLA
            QLA = 0.0

        if QIA < 1.E-8:
            QV = QV + QIA
            TE = TE - (radconstants.alhsbcp) * QIA
            QIA = 0.0

        if (QLA + QIA) < 1.E-8:
            QV = QV + QLA + QIA
            TE = TE - (radconstants.alhlbcp) * QLA - (radconstants.alhsbcp) * QIA
            AF = 0.0
            QLA = 0.0
            QIA = 0.0

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
        FAC_RL: FloatField, 
        MIN_RL: FloatField, 
        MAX_RL: FloatField, 
        FAC_RI: FloatField, 
        MIN_RI: FloatField, 
        MAX_RI: FloatField,
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
        #RAD_RL = max(MIN_RL, min(LDRADIUS4(PL, TE, RAD_QL, NL, NI, 1) * FAC_RL, MAX_RL))
        #Ice radii - Brams formulation with limits
        #RAD_RI = max(MIN_RI, min(LDRADIUS4(PL, TE, RAD_QI, NL, NI, 2) * FAC_RI, MAX_RI))

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
            externals = {

            }
        )
        self._radcouple = stencil_factory.from_dims_halo(
            func=_radcouple_stencil, 
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals = {

            }
        )
        self.do_qa = do_qa

        '''
        self._dx = grid_data.dx
        self._dy = grid_data.dy
        self._a11 = grid_data.a11
        self._a12 = grid_data.a12
        self._a21 = grid_data.a21
        self._a22 = grid_data.a22
        self.comm = comm
        '''

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
        FAC_RL: FloatField,
        MIN_RL: FloatField,
        MAX_RL: FloatField,
        FAC_RI: FloatField,
        MIN_RI: FloatField,
        MAX_RI: FloatField,
    ):
        self._fix_up_clouds(Q, 
                            T, 
                            QLLS, 
                            CLLS, 
                            QLCN, 
                            CLCN, 
                            QICN, 
                            QGRAUPEL,
                            )
        #self._fix_up_clouds(Q)
        self._radcouple(T, PLmb, CLLS, CLCN, Q, QLLS, QILS, QLCN, QICN, QRAIN, QSNOW, QGRAUPEL, NACTL, NACTI, RAD_QV, RAD_QL, RAD_QI, RAD_QR, RAD_QS, RAD_QG, RAD_CF, CLDREFFL, CLDREFFI, FAC_RL, MIN_RL, MAX_RL, FAC_RI, MIN_RI, MAX_RI)
        #self._radcouple(Q)
        #if self.do_qa:
            # Implement QA logic if needed
          #  RHX = Q / GEOS_QSAT(T, PLmb)
