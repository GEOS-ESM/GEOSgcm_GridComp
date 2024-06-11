from ndsl import QuantityFactory, StencilFactory  # type: ignore
from ndsl.constants import X_DIM, Y_DIM, Z_DIM  # type: ignore
from ndsl.comm.communicator import Communicator  # type: ignore
from ndsl.dsl.typing import FloatField  # type: ignore
from gt4py.cartesian.gtscript import computation, interval, PARALLEL, function, log10,   # type: ignore
import numpy as np

#Math constants
MAPL_PI_R8 = 3.14159265358979323846e0
MAPL_PI = MAPL_PI_R8
MAPL_DEGREES_TO_RADIANS_R8 = MAPL_PI_R8 / 180.e0
MAPL_DEGREES_TO_RADIANS = MAPL_PI / 180.0
MAPL_RADIANS_TO_DEGREES = 180.e0 / MAPL_PI_R8

#Define whether or not to use CODATA 2018 constants
CODATA_2018_CONSTANTS = True

# Universal Constants
if CODATA_2018_CONSTANTS:
    MAPL_STFBOL = 5.670374419E-8  # W/(m^2 K^4)
    MAPL_AVOGAD = 6.02214076E26   # 1/kmol
    MAPL_RUNIV = 8314.462618      # J/(Kmole K)
else: 
    MAPL_STFBOL = 5.6734E-8  # W/(m^2 K^4)
    MAPL_AVOGAD = 6.023E26   # 1/kmol
    MAPL_RUNIV = 8314.47    # J/(Kmole K)

# Earth Constants
MAPL_PSDRY = 98305.0  # Pa
MAPL_SECONDS_PER_SIDEREAL_DAY = 86164.0  # s
MAPL_GRAV = 9.80665  # m^2/s
MAPL_RADIUS = 6371.0E3  # m
MAPL_OMEGA_R8 = 2.0 * MAPL_PI_R8 / MAPL_SECONDS_PER_SIDEREAL_DAY  # 1/s
MAPL_OMEGA = 2.0 * MAPL_PI / MAPL_SECONDS_PER_SIDEREAL_DAY  # 1/s
MAPL_EARTH_ECCENTRICITY = 8.181919084262200e-2  # --
MAPL_EARTH_SEMIMAJOR_AXIS = 6378137.0  # m
MAPL_KM_PER_DEG = (1.0 / (MAPL_RADIUS / 1000.0)) * MAPL_RADIANS_TO_DEGREES
MAPL_DEG_PER_KM = (MAPL_RADIUS / 1000.0) * MAPL_DEGREES_TO_RADIANS_R8

# Physical properties
MAPL_H2OMW = 18.015  # kg/Kmole
MAPL_O3MW = 47.9982  # kg/Kmole
MAPL_LATENT_HEAT_VAPORIZATION = 2.4665E6  # J/kg @15C @1atm
MAPL_ALHL = MAPL_LATENT_HEAT_VAPORIZATION  # J/kg
MAPL_LATENT_HEAT_FUSION = 3.3370E5  # J/kg @1atm
MAPL_ALHF = MAPL_LATENT_HEAT_FUSION  # J/kg
MAPL_LATENT_HEAT_SUBLIMATION = MAPL_ALHL + MAPL_ALHF  # J/kg
MAPL_ALHS = MAPL_LATENT_HEAT_SUBLIMATION  # J/kg

# Earth Specific Chemistry and Thermodynamic Constants
MAPL_AIRMW = 28.965  # kg/Kmole
MAPL_RDRY = MAPL_RUNIV / MAPL_AIRMW  # J/(kg K)
MAPL_CPDRY = 3.5 * MAPL_RDRY  # J/(kg K)
MAPL_CVDRY = MAPL_CPDRY - MAPL_RDRY  # J/(kg K)
MAPL_RVAP = MAPL_RUNIV / MAPL_H2OMW  # J/(kg K)
MAPL_CPVAP = 4.0 * MAPL_RVAP  # J/(kg K)
MAPL_CVVAP = MAPL_CPVAP - MAPL_RVAP  # J/(kg K)
MAPL_KAPPA = MAPL_RDRY / MAPL_CPDRY  # (2.0/7.0)
MAPL_EPSILON = MAPL_H2OMW / MAPL_AIRMW  # --
MAPL_DELTAP = MAPL_CPVAP / MAPL_CPDRY  # --
MAPL_DELTAV = MAPL_CVVAP / MAPL_CVDRY  # --
MAPL_GAMMAD = MAPL_CPDRY / MAPL_CVDRY  # --
MAPL_RGAS = MAPL_RDRY  # J/(kg K) (DEPRECATED)
MAPL_CP = MAPL_RGAS / MAPL_KAPPA  # J/(kg K) (DEPRECATED)
MAPL_VIREPS = 1.0 / MAPL_EPSILON - 1.0  # (DEPRECATED)
MAPL_P00 = 100000.0  # Pa
MAPL_CAPICE = 2000.0  # J/(K kg)
MAPL_CAPWTR = 4218.0  # J/(K kg)
MAPL_RHOWTR = 1000.0  # kg/m^3
MAPL_NUAIR = 1.533E-5  # m^2/S (@ 18C)
MAPL_TICE = 273.16  # K
MAPL_SRFPRS = 98470  # Pa
MAPL_KARMAN = 0.40  # --
MAPL_USMIN = 1.00  # m/s
MAPL_RHO_SEAWATER = 1026.0  # sea water density [kg/m^3]
MAPL_RHO_SEAICE = 917.0  # sea ice density [kg/m^3]
MAPL_RHO_SNOW = 330.0  # snow density [kg/m^3]
MAPL_CELSIUS_TO_KELVIN = 273.15  # K

# Define constants
alhlbcp = MAPL_ALHL / MAPL_CP
alhsbcp = MAPL_ALHS / MAPL_CP

#for LDRADIUS4
MAPL_RGAS = 287.0
MAPL_TICE = 273.15
LIQ_RADII_PARAM = 1
ICE_RADII_PARAM = 1

bx = 100.0 * (3.0 / (4.0 * MAPL_PI)) ** (1.0 / 3.0)
r13bbeta = 1.0 / 3.0 - 0.14
abeta = 0.07

RHO_W = 1000.0
Ldiss = 0.07
Lk = 0.75
Lbx = Ldiss * 1.0e3 * (3.0 / (4.0 * MAPL_PI * Lk * RHO_W * 1.0e-3)) ** (1.0 / 3.0)
Lbe = 1.0 / 3.0 - 0.14

#Taken from GEOS_GFDL_1M_InterfaceMod.F90:269-274
#Uses MAPL_GetResource to get variables. These do not change in the radiation coupling unless stated otherwise in the code.
#These might be wrong... Checking on it atm.
MIN_RL = 2.5e-6
MAX_RL = 60.0e-6
FAC_RL = 1.0
MIN_RI = 5.0e-6
MAX_RI = 100.0e-6
FAC_RI = 1.0


# GEOS_QSAT Function --This might be wrong
def GEOS_QSAT(T: float, P: float) -> float:
    Rv = 461.5  # Gas constant for water vapor [J/(kg·K)]
    Lv = 2.5e6  # Latent heat of vaporization [J/kg]
    es0 = 610.78  # Saturation vapor pressure at 0°C [Pa]
    
    # Saturation vapor pressure formula (Clausius-Clapeyron equation)
    es = es0 * np.exp(Lv / Rv * (1/273.15 - 1/T))
    
    # Saturation mixing ratio
    qsat = 0.622 * es / (P - es)
    
    return qsat

#@function
def LDRADIUS4(PL: float, TE: float, QC: float, NNL: float, NNI: float, ITYPE: int) -> float:
    RHO = (100.0 * PL) / (MAPL_RGAS * TE)
    if ITYPE == 1:  # Liquid
        WC = 1.0e3 * RHO * QC
        NNX = max(NNL * 1.0e-6, 10.0)
        if LIQ_RADII_PARAM == 1:
            RADIUS = min(60.0e-6, max(2.5e-6, 1.0e-6 * bx * (WC / NNX)**r13bbeta * abeta * 6.92))
        else:
            RADIUS = min(60.0e-6, max(2.5e-6, 1.0e-6 * Lbx * (WC / NNX)**Lbe))
    elif ITYPE == 2:  # Ice
        WC = 1.0e3 * RHO * QC
        if ICE_RADII_PARAM == 1:
            if TE > MAPL_TICE or QC <= 0.0:
                BB = -2.0
            else:
                BB = -2.0 + np.log10(WC / 50.0) * (1.0e-3 * (MAPL_TICE - TE)**1.5)
            BB = min(max(BB, -6.0), -2.0)
            RADIUS = 377.4 + 203.3 * BB + 37.91 * BB**2 + 2.3696 * BB**3
            RADIUS = min(150.0e-6, max(5.0e-6, 1.0e-6 * RADIUS))
        else:
            TC = TE - MAPL_TICE
            ZFSR = 1.2351 + 0.0105 * TC
            AA = 45.8966 * (WC**0.2214)
            BB = 0.79570 * (WC**0.2535)
            RADIUS = ZFSR * (AA + BB * (TE - 83.15))
            RADIUS = min(150.0e-6, max(5.0e-6, 1.0e-6 * RADIUS * 0.64952))
    else:
        raise ValueError("WRONG HYDROMETEOR type: CLOUD = 1 OR ICE = 2")
    return RADIUS

def _fix_up_clouds_stencil(QV: FloatField, TE: FloatField, QLC: FloatField, CF: FloatField, QLA: FloatField, AF: FloatField, QIC: FloatField, QIA: FloatField) -> None:
    with computation(PARALLEL), interval(...):
        if AF < 1.E-5:
            QV = QV + QLA + QIA
            TE = TE - (alhlbcp) * QLA - (alhsbcp) * QIA
            AF = 0.
            QLA = 0.
            QIA = 0.

        if CF < 1.E-5:
            QV = QV + QLC + QIC
            TE = TE - (alhlbcp) * QLC - (alhsbcp) * QIC
            CF = 0.
            QLC = 0.
            QIC = 0.

        if QLC < 1.E-8:
            QV = QV + QLC
            TE = TE - (alhlbcp) * QLC
            QLC = 0.

        if QIC < 1.E-8:
            QV = QV + QIC
            TE = TE - (alhsbcp) * QIC
            QIC = 0.

        if QLA < 1.E-8:
            QV = QV + QLA
            TE = TE - (alhlbcp) * QLA
            QLA = 0.

        if QIA < 1.E-8:
            QV = QV + QIA
            TE = TE - (alhsbcp) * QIA
            QIA = 0.

        if (QLA + QIA) < 1.E-8:
            QV = QV + QLA + QIA
            TE = TE - (alhlbcp) * QLA - (alhsbcp) * QIA
            AF = 0.
            QLA = 0.
            QIA = 0.

        if (QLC + QIC) < 1.E-8:
            QV = QV + QLC + QIC
            TE = TE - (alhlbcp) * QLC - (alhsbcp) * QIC
            CF = 0.
            QLC = 0.
            QIC = 0.

def _radcouple_stencil(TE: FloatField, PL: FloatField, CF: FloatField, AF: FloatField, QV: FloatField, QClLS: FloatField, QCiLS: FloatField, QClAN: FloatField, QCiAN: FloatField, QRN_ALL: FloatField, QSN_ALL: FloatField, QGR_ALL: FloatField, NL: FloatField, NI: FloatField, RAD_QV: FloatField, RAD_QL: FloatField, RAD_QI: FloatField, RAD_QR: FloatField, RAD_QS: FloatField, RAD_QG: FloatField, RAD_CF: FloatField, RAD_RL: FloatField, RAD_RI: FloatField, FAC_RL: FloatField, MIN_RL: FloatField, MAX_RL: FloatField, FAC_RI: FloatField, MIN_RI: FloatField, MAX_RI: FloatField) -> None:
    with computation(PARALLEL), interval(...):
        RAD_QV = QV
        RAD_CF = max(min(CF + AF, 1.0), 0.0)

        if RAD_CF >= 1.e-5:
            RAD_QL = (QClLS + QClAN) / RAD_CF if (QClLS + QClAN) >= 1.e-8 else 0.0
            RAD_QI = (QCiLS + QCiAN) / RAD_CF if (QCiLS + QCiAN) >= 1.e-8 else 0.0
            RAD_QR = QRN_ALL / RAD_CF if QRN_ALL >= 1.e-8 else 0.0
            RAD_QS = QSN_ALL / RAD_CF if QSN_ALL >= 1.e-8 else 0.0
            RAD_QG = QGR_ALL / RAD_CF if QGR_ALL >= 1.e-8 else 0.0
        else:
            RAD_CF = RAD_QL = RAD_QI = RAD_QR = RAD_QS = RAD_QG = 0.0

        RAD_QL = min(RAD_QL, 0.01)
        RAD_QI = min(RAD_QI, 0.01)
        RAD_QR = min(RAD_QR, 0.01)
        RAD_QS = min(RAD_QS, 0.01)
        RAD_QG = min(RAD_QG, 0.01)

        RAD_RL = max(MIN_RL, min(LDRADIUS4(PL, TE, RAD_QL, NL, NI, 1) * FAC_RL, MAX_RL))
        RAD_RI = max(MIN_RI, min(LDRADIUS4(PL, TE, RAD_QI, NL, NI, 2) * FAC_RI, MAX_RI))

class RadiationCoupling:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        do_qa: bool,
        #grid_data: GridData,
        #comm: Communicator,
    ):
        self._fix_up_clouds = stencil_factory.from_dims_halo(
            func=_fix_up_clouds_stencil, compute_dims=[X_DIM, Y_DIM, Z_DIM]
        )
        self._radcouple = stencil_factory.from_dims_halo(
            func=_radcouple_stencil, compute_dims=[X_DIM, Y_DIM, Z_DIM]
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
        self._fix_up_clouds(Q, T, QLLS, CLLS, QLCN, CLCN, QICN, QGRAUPEL)
        self._radcouple(T, PLmb, CLLS, CLCN, Q, QLLS, QILS, QLCN, QICN, QRAIN, QSNOW, QGRAUPEL, NACTL, NACTI, RAD_QV, RAD_QL, RAD_QI, RAD_QR, RAD_QS, RAD_QG, RAD_CF, CLDREFFL, CLDREFFI, FAC_RL, MIN_RL, MAX_RL, FAC_RI, MIN_RI, MAX_RI)
        if self.do_qa:
            # Implement QA logic if needed
            RHX = Q / GEOS_QSAT(T, PLmb)
