import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import computation, interval, PARALLEL, FORWARD, atan, sin, tan, sqrt, tanh, exp, log10
from ndsl.boilerplate import get_factories_single_tile_numpy
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField, FloatFieldIJ, Float, IntField, IntFieldIJ, Int
from ndsl import StencilFactory, QuantityFactory, orchestrate
import numpy as np
import pyMoist.pyMoist_constants as constants
from pyMoist.saturation.qsat import QSat_Float


@gtscript.function
def air_density(PL: Float, TE: Float) -> Float:
    """
    Calculate air density [kg/m^3]

    Parameters:
    PL (Float): Pressure level.
    TE (Float): Temperature.

    Returns:
    Float: Calculated air density.
    """
    air_density = (100.0 * PL) / (constants.rdry * TE)
    return air_density

@gtscript.function
def ice_fraction_modis(
    temp: Float,
    cnv_frc: Float,
    srf_type: Float,
):
    # Use MODIS polynomial from Hu et al, DOI: (10.1029/2009JD012384)
    tc = max(-46.0, min(temp-constants.t_ice, 46.0)) # convert to celcius and limit range from -46:46 C
    ptc = 7.6725 + 1.0118*tc + 0.1422*tc**2 + 0.0106*tc**3 + 0.000339*tc**4 + 0.00000395*tc**5
    ice_frct = 1.0 - (1.0/(1.0 + exp(-1*ptc)))
    return ice_frct

@gtscript.function
def ice_fraction(
    temp: Float,
    cnv_frc: Float,
    srf_type: Float,
):
    # Anvil clouds
    # Anvil-Convective sigmoidal function like figure 6(right)
    # Sigmoidal functions Hu et al 2010, doi:10.1029/2009JD012384
    if temp <= constants.JaT_ICE_ALL:
        icefrct_c = 1.000
    elif temp > constants.JaT_ICE_ALL and temp <= constants.JaT_ICE_MAX:
        icefrct_c = sin(0.5 * constants.PI * (1.00 - (temp - constants.JaT_ICE_ALL) / (constants.JaT_ICE_MAX - constants.JaT_ICE_ALL)))
    else:
        icefrct_c = 0.00
    icefrct_c = max(min(icefrct_c,1.00),0.00) ** constants.aICEFRPWR
    # Sigmoidal functions like figure 6b/6c of Hu et al 2010, doi:10.1029/2009JD012384
    if srf_type == 2.0:
        if temp <= constants.JiT_ICE_ALL:
            icefrct_m = 1.000
        elif temp > constants.JiT_ICE_ALL and temp <= constants.JiT_ICE_MAX:
            icefrct_m = 1.00 - (temp - constants.JiT_ICE_ALL) / (constants.JiT_ICE_MAX - constants.JiT_ICE_ALL)
        else:
            icefrct_m = 0.00
        icefrct_m = max(min(icefrct_m, 1.00), 0.00) ** constants.iICEFRPWR
    elif srf_type > 1.0:
        if temp <= constants.lT_ICE_ALL:
           icefrct_m = 1.000
        elif temp > constants.lT_ICE_ALL and temp <= constants.lT_ICE_MAX:
           icefrct_m = sin(0.5 * constants.PI * (1.00 - (temp - constants.lT_ICE_ALL) / (constants.lT_ICE_MAX - constants.lT_ICE_ALL)))
        else:
            icefrct_m = 0.00
        icefrct_m = max(min(icefrct_m, 1.00), 0.00) ** constants.lICEFRPWR
    else:
        if temp <= constants.oT_ICE_ALL:
           icefrct_m = 1.000
        elif temp > constants.oT_ICE_ALL and temp <= constants.oT_ICE_MAX:
           icefrct_m = sin(0.5 * constants.PI * (1.00 - (temp - constants.oT_ICE_ALL) / (constants.oT_ICE_MAX - constants.oT_ICE_ALL)))
        else:
            icefrct_m = 0.00
        icefrct_m = max(min(icefrct_m, 1.00), 0.00) ** constants.oICEFRPWR
    ice_frac = icefrct_m * (1.0-cnv_frc) + icefrct_c * cnv_frc
    return ice_frac

@gtscript.function
def cloud_effective_radius_liquid(
    PL: Float,
    TE: Float,
    QC: Float,
    NNL: Float,
    NNI: Float,
) -> Float:
    """
    Calculate the effective radius of liquid clouds [m]
    Implementation of LDRADIUS4 for liquid clouds

    Parameters:
    PL (Float): Pressure level.
    TE (Float): Temperature.
    QC (Float): Liquid cloud mixing ratio.
    NNL (Float): Number concentration of liquid cloud droplets.
    NNI (Float): Number concentration of ice cloud crystals. Not used in function body.

    Returns:
    Float: Effective radius of liquid clouds.
    """
    # Calculate liquid water content
    WC = (
        1.0e3 * air_density(PL, TE) * QC
    )  # air density [g/m3] * liquid cloud mixing ratio [kg/kg]
    # Calculate cloud drop number concentration from the aerosol model + ....
    NNX = max(NNL * 1.0e-6, 10.0)
    # Calculate Radius in meters [m]
    if constants.LIQ_RADII_PARAM == 1:
        # Jason Version
        RADIUS = min(
            60.0e-6,
            max(
                2.5e-6,
                1.0e-6
                * constants.bx
                * (WC / NNX) ** constants.r13bbeta
                * constants.abeta
                * 6.92,
            ),
        )
    else:
        # [liu&daum, 2000 and 2005. liu et al 2008]
        RADIUS = min(
            60.0e-6,
            max(2.5e-6, 1.0e-6 * constants.Lbx * (WC / NNX) ** constants.Lbe),
        )
    return RADIUS

@gtscript.function
def cloud_effective_radius_ice(
    PL: Float,
    TE: Float,
    QC: Float,
    NNL: Float,
    NNI: Float,
) -> Float:
    """
    Calculate the effective radius of ice clouds [m]
    Implementation of LDRADIUS4 for Ice clouds

    Parameters:
    PL (Float): Pressure level.
    TE (Float): Temperature.
    QC (Float): Ice cloud mixing ratio.
    NNL (Float): Number concentration of liquid cloud droplets. Not used in function body, but included in the Fortran code.
    NNI (Float): Number concentration of ice cloud crystals. Not used in function body, but included in the Fortran code.

    Returns:
    Float: Effective radius of ice clouds.
    """
    # Calculate ice water content
    WC = (
        1.0e3 * air_density(PL, TE) * QC
    )  # air density [g/m3] * ice cloud mixing ratio [kg/kg]
    # Calculate radius in meters [m]
    if constants.ICE_RADII_PARAM == 1:
        # Ice cloud effective radius -- [klaus wyser, 1998]
        if TE > constants.t_ice or QC <= 0.0:
            BB = -2.0
        else:
            BB = -2.0 + log10(WC / 50.0) * (
                1.0e-3 * (constants.t_ice - TE) ** 1.5
            )
        BB = min(max(BB, -6.0), -2.0)
        RADIUS = 377.4 + 203.3 * BB + 37.91 * BB**2 + 2.3696 * BB**3
        RADIUS = min(150.0e-6, max(5.0e-6, 1.0e-6 * RADIUS))
    else:
        # Ice cloud effective radius ----- [Sun, 2001]
        TC = TE - constants.t_ice
        ZFSR = 1.2351 + 0.0105 * TC
        AA = 45.8966 * (WC**0.2214)
        BB = 0.79570 * (WC**0.2535)
        RADIUS = ZFSR * (AA + BB * (TE - 83.15))
        RADIUS = min(150.0e-6, max(5.0e-6, 1.0e-6 * RADIUS * 0.64952))
    return RADIUS

@gtscript.function
def pdffrac (
    flag: Int,
    qtmean: Float,
    sigmaqt1: Float,
    sigmaqt2: Float,
    qstar: Float,
    clfrac: Float,
):
    if flag == 1:
        if (qtmean + sigmaqt1) < qstar:
            clfrac = 0.
        else:
            if sigmaqt1 > 0.:
                clfrac = min((qtmean + sigmaqt1 - qstar),2.*sigmaqt1)/(2.*sigmaqt1)
            else:
                clfrac = 1.
    elif flag == 2:
        qtmode =  qtmean + (sigmaqt1-sigmaqt2)/3.
        qtmin = max(qtmode-sigmaqt1,0.)
        qtmax = qtmode + sigmaqt2
        if qtmax <= qstar:
            clfrac = 0.
        elif (qtmode <= qstar) and (qstar <= qtmax):
            clfrac = (qtmax - qstar) * (qtmax - qstar) / ((qtmax - qtmin) * (qtmax - qtmode))
        elif (qtmin <= qstar)and (qstar < qtmode):
            clfrac = 1. - ( (qstar-qtmin)*(qstar-qtmin) / ( (qtmax-qtmin)*(qtmode-qtmin) ) )
        elif qstar <= qtmin:
            clfrac = 1.

    return clfrac

@gtscript.function
def pdfcondensate(
    flag: Int,
    qtmean: Float,
    sigmaqt1: Float,
    sigmaqt2: Float,
    qstar: Float,
    condensate: Float,
):
    if flag ==1:
        if (qtmean + sigmaqt1) < qstar:
            condensate = 0.
        elif qstar > (qtmean-sigmaqt1):
            if sigmaqt1 > 0.:
                condensate = (min(qtmean + sigmaqt1 - qstar,2.d0*sigmaqt1)**2)/ (4.d0*sigmaqt1)
            else:
                condensate = qtmean-qstar
        else:
            condensate = qtmean-qstar
    elif flag ==2:
        qtmode =  qtmean + (sigmaqt1 - sigmaqt2) / 3.
        qtmin = max(qtmode - sigmaqt1, 0.)
        qtmax = qtmode + sigmaqt2
        if qtmax <= qstar:
            condensate = 0.
        elif (qtmode <= qstar) and (qstar < qtmax):
            constB = 2. / ((qtmax - qtmin) * (qtmax - qtmode))
            cloudf = (qtmax - qstar) * (qtmax - qstar) / ((qtmax - qtmin) * (qtmax - qtmode))
            term1 = (qstar * qstar * qstar) / 3.
            term2 = (qtmax * qstar * qstar)/2.
            term3 = (qtmax * qtmax * qtmax)/6.
            condensate = constB * (term1 - term2 + term3) - qstar * cloudf
        elif (qtmin <= qstar) and (qstar < qtmode):
            constA = 2. / ((qtmax - qtmin) * (qtmode - qtmin))
            cloudf = 1. - ((qstar - qtmin) * (qstar - qtmin) / ((qtmax - qtmin) * (qtmode - qtmin)))
            term1 = (qstar * qstar * qstar) / 3.
            term2 = (qtmin * qstar * qstar) / 2.
            term3 = (qtmin * qtmin * qtmin) / 6.
            condensate = qtmean - (constA * (term1 - term2 + term3)) - qstar * cloudf
        elif qstar <= qtmin:
            condensate = qtmean - qstar

    return condensate

@gtscript.function
def bergeron_partition(
    DTIME: Float,
    PL: Float,
    TE: Float,
    QV: Float,
    QILS: Float,
    QICN: Float,
    QLLS: Float,
    QLCN: Float,
    CF: Float,
    AF: Float,
    NL: Float,
    NI: Float,
    DQALL: Float,
    fQI: Float,
    CNV_FRC: Float,
    SRF_TYPE: Float,
    needs_preexisting: bool = False,
):
    DIFF = 0.0
    DEP  = 0.0
    QI   = QILS + QICN # neccesary because NI is for convective and large scale
    QL   = QLLS + QLCN
    QTOT = QI+QL
    FQA  = 0.0
    if QTOT > 0.0: FQA = (QICN + QILS) / QTOT
    NIX  = (1.0 - FQA) * NI

    DQALL = DQALL / DTIME
    CFALL = min(CF + AF, 1.0)
    TC    = TE - constants.t_ice
    fQI_0 = fQI

    # Completelely glaciated cloud:
    if TE >= constants.iT_ICE_MAX: # liquid cloud
        FQI = 0.
    elif TE <= constants.iT_ICE_ALL: # ice cloud
        FQI = 1.
    else: # mixed phase cloud
        fQI = 0.
        if QILS <= 0.:
            if needs_preexisting is True:
                raise NotImplementedError(
                    "Not implemented in Fortran, see Process_Library.F90"
                )
            else:
                fQI = ice_fraction(TE, CNV_FRC, SRF_TYPE)
            return
    
        QVINC =  QV
        QSLIQ = GEOS_QsatLQU( TE, PL*100.0 , DQ=DQSL )
        QSICE = GEOS_QsatICE( TE, PL*100.0 , DQ=DQSI )
        QVINC = MIN(QVINC, QSLIQ) !limit to below water saturation



def get_last(in_field: FloatField, temporary_field: FloatFieldIJ, out_field: FloatField):
    with computation(FORWARD), interval(-1, None):
        temporary_field = in_field

    with computation(PARALLEL), interval(...):
        out_field = temporary_field

def hybrid_index_2dout(
    data_field: FloatField,
    k_mask: FloatField,
    k_index_desired: FloatFieldIJ,
    out_field: FloatFieldIJ,
):
    with computation(FORWARD), interval(...):
        if k_mask == k_index_desired:
            out_field = data_field

def initial_calc(
    EIS: FloatFieldIJ,
    dw_land: Float,
    dw_ocean: Float,
    TURNRHCRIT_PARAM: Float,
    minrhcrit: FloatField,
    PLmb_at_klcl: FloatFieldIJ,
    PLmb: FloatField,
    PLEmb_top: FloatField,
    AREA: FloatFieldIJ,
    ALPHA: FloatField,
): 
    with computation(PARALLEL), interval(...):
        # Send the condensates through the pdf after convection
        facEIS = max(0.0, min(1.0, EIS/10.0))**2
        # determine combined minrhcrit in stable/unstable regimes
        minrhcrit = (1.0-dw_ocean)*(1.0-facEIS) + (1.0-dw_land)*facEIS
        # determine the turn pressure using the LCL
        if TURNRHCRIT_PARAM <= 0:
            turnrhcrit = PLmb_at_klcl -250 # implementation of for loop needs to go here
        else:
            turnrhcrit = TURNRHCRIT_PARAM

    # lower turn from maxrhcrit=1.0 # implemented in multiple "with" statements to deal with hybrid indexing
    with computation(PARALLEL), interval(0,-1):
        # Use Slingo-Ritter (1985) formulation for critical rel ative humidity
        RHCRIT = 1.0
        if PLmb <= turnrhcrit:
            RHCRIT = minrhcrit
        else:
            RHCRIT = minrhcrit + (1.0-minrhcrit)/(19.) * ((atan( (2.*(PLmb-turnrhcrit)/(PLEmb_top-turnrhcrit)-1.) *
                tan(20.*constants.PI/21.-0.5*constants.PI) )
                + 0.5*constants.PI) * 21./constants.PI - 1.)
    with computation(PARALLEL), interval(-1,None):
        # lower turn from maxrhcrit=1.0
        if PLmb <= turnrhcrit:
            RHCRIT = minrhcrit
        else:
            RHCRIT = 1.0
    with computation(PARALLEL), interval(...):
        # include grid cell area scaling and limit RHcrit to > 70%\
        ALPHA = max(0.0,min(0.30, (1.0-RHCRIT)*sqrt(sqrt(AREA/1.e10))))

def hystpdf(
    DT_MOIST: Float,
    ALPHA: FloatField,
    PDFSHAPE: Float,
    cnv_frc: FloatFieldIJ,
    srf_type: FloatFieldIJ,
    PL: FloatField,
    ZL: FloatField,
    Q: FloatField,
    QLLS: FloatField,
    QLCN: FloatField,
    QILS: FloatField,
    QICN: FloatField,
    TE: FloatField,
    CLLS: FloatField,
    CLCN: FloatField,
    NL: FloatField,
    NI: FloatField,
    WHL: FloatField,
    WQT: FloatField,
    HL2: FloatField,
    QT2: FloatField,
    HLQT: FloatField,
    W3: FloatField,
    W2: FloatField,
    MFQT3: FloatField,
    MFHL3: FloatField,
    WTHV2: FloatField,
    WQL: FloatField,
    needs_preexisting: bool = False,
    USE_BERGERON: bool = True,
):
    with computation(PARALLEL), interval(...):
        scice = 1.0 # don't ask, I don't know
        if CLCN < 1.0:
            tmpARR = 1.0 / (1.0 - CLCN)
        else:
            tmpARR = 0.0
        if CLCN > 0.0: # need to make sure this is equivilant to fortran CLCN > tiny(0.0)
            QAx = (QLCN + QICN) / CLCN
        else:
            QAx = 0.0
        CFn = CLLS * tmpARR
        QCn = (QLLS + QILS) * tmpARR
        QCi = QILS * tmpARR
        TEn = TE
        QSx, DQS = QSat_Float(TE, PL, DQSAT_trigger=True)
        QVn = (Q - QSx * CLCN) * tmpARR
        
        QT = QCn + QVn # Total LS water after microphysics

        count = 1
        while count <= 20:
            QVp = QVn
            QCp = QCn
            CFp = CFn
            TEp = TEn
            QSn, DQS = QSat_Float(TEn, PL, DQSAT_trigger=True)

            # SC_ICE option not implemented, will go here if needed
            # Fortran:
            # if(present(SC_ICE)) then
            #     PRINT *, SC_ICE
            #     scice = min(max(SC_ICE, 1.0), 1.7)
            #     qsnx= Qsn*scice !
            #     if ((QCi .ge. 0.0) .and. (Qsn .gt. Qt))  QSn=Qsnx !this way we do not evaporate preexisting ice but maintain supersat
            # end if

            if PDFSHAPE < 3: # 1 = top-hat 2 = triangulat
                sigmaqt1  = ALPHA*QSn
                sigmaqt2  = ALPHA*QSn
            elif PDFSHAPE == 4: # lognormal (sigma is dimensionless)
                sigmaqt1 =  max(ALPHA/sqrt(3.0), 0.001)
            
            if PDFSHAPE < 5:
                CFn = pdffrac(PDFSHAPE, QT, sigmaqt1, sigmaqt2, QSn, CFn)
                QCn = pdfcondensate(PDFSHAPE, QT, sigmaqt1, sigmaqt2, QSn, QCn)
            elif PDFSHAPE == 5:
                raise NotImplementedError(
                    "Function: partition_dblgss not implemented.\nReference Fortran: Process_LibraryF90"
                )
            
            if USE_BERGERON is True:
                DQCALL = QCn - QCp
                Nfac = 100. * PL * constants.R_AIR / TEn # density times conversion factor
                NLv = NL/Nfac
                NIv = NI/Nfac
                bergeron_partition
            else:
                fQi = ice_fraction(TEn, cnv_frc, srf_type)

            latent_heat_factor = (1.0 - fQi) * constants.latent_heat_vaporization / constants.cpdry + fQi*constants.latent_heat_fusion / constants.cpdry
            if PDFSHAPE == 1:
                QCn = QCp + (QCn-QCp) / (1. - (CFn * (ALPHA - 1.) - (QCn / QSn)) * DQS * latent_heat_factor)
            elif PDFSHAPE == 2:
                raise NotImplementedError(
                    "Fortran comments says code is incorrect and therefore not implemented. See Process_Library.F90"
                )
            elif PDFSHAPE == 5:
                QCn = QCp + 0.5 * (QCn - QCp)
            
            if CLCN > 0:
                QAo = QAx
            else:
                QAo = 0
            
            QVn = QVp - (QCn - QCp)
            TEn = TEp + (1.0-fQi)* constants.latent_heat_vaporization / constants.cpdry * ((QCn - QCp) * (1. - CLCN) + (QAo - QAx) * CLCN) \
                + fQi * constants.latent_heat_sublimation / constants.cpdry * ((QCn - QCp) * (1. - CLCN) + (QAo - QAx) * CLCN)
            
            PDFITERS = count
            if abs(TEn - TEp) < 0.00001: break
        
        if CLCN < 1.0:
            CLLS = CFn * (1.-CLCN)
            QCn  = QCn * (1.-CLCN)
            QAo  = QAo * CLCN
        else:
            # Special case CLCN=1, i.e., box filled with anvil.
            # - Note: no guarantee QV_box > QS_box
            CLLS = 0. #Remove any LS cloud
            QAo  = QLCN+QICN+QLLS+QILS # Add all LS condensate to anvil type
            QCn  = 0. # Remove same from new LS
            QT   = QAo + QV # Update total water
            # Now set anvil condensate to any excess of total water
            # over QSx (saturation value at top)
            QAo = max(QT - QSx, 0.)

        # Now take {\em New} condensate and partition into ice and liquid

        # large scale
        QCx = QCn - (QLLS + QILS)
        if (QCx < 0.0): # net evaporation
            dQLLS = max(QCx, -QLLS) # Water evaporates first
            dQILS = max(QCx - dQLLS, -QILS) # Then sublimation
        else:
            dQLLS = (1.0 - fQi) * QCx
            dQILS = fQi * QCx
        
        # convective
        QAx = QAo - (QLCN+QICN)
        if  QAx < 0.0: #net evaporation
            dQLCN = max(QAx        , -QLCN) # Water evaporates first
            dQICN = max(QAx - dQLCN, -QICN) # Then sublimation
        else:
            dQLCN  = (1.0 - fQi) * QAx
            dQICN  = fQi *QAx

        # Clean-up cloud if fractions are too small
        if CLCN < 1.e-5:
            dQICN = -QICN
            dQLCN = -QLCN
        if CLLS < 1.e-5:
            dQILS = -QILS
            dQLLS = -QLLS

        QICN = QICN + dQICN
        QLCN = QLCN + dQLCN
        QILS = QILS + dQILS
        QLLS = QLLS + dQLLS
        QV = QV - dQICN + dQILS + dQLCN + dQLLS
        TE = TE + constants.latent_heat_vaporization / constants.cpdry *(dQICN + dQILS + dQLCN + dQLLS) \
            + constants.latent_heat_fusion / constants.cpdry *(dQICN + dQILS)

        # We need to take care of situations where QS moves past QA
        # during QSAT iteration. This should be only when QA/AF is small
        # to begin with. Effect is to make QAo negative. So, we
        # "evaporate" offending QA's
        # We get rid of anvil fraction also, although strictly
        # speaking, PDF-wise, we should not do this.

        if QAo <= 0.:
            QV = QV + QICN + QLCN
            TE = TE - constants.latent_heat_sublimation / constants.cpdry *QICN - constants.latent_heat_vaporization / constants.cpdry *QLCN
            QICN = 0.
            QLCN = 0.
            CLCN = 0.

def meltfrz(
    dt: Float,
    cnv_frc: FloatFieldIJ,
    srf_type: FloatFieldIJ,
    T: FloatField,
    QLCN: FloatField,
    QICN: FloatField,
):
    with computation(PARALLEL), interval(...):
        if T < constants.t_ice:
            fQi  = ice_fraction( T, cnv_frc, srf_type)
            dQil = QLCN *(1.0 - exp( -dt * fQi / constants.taufrz))
            dQil = max(0., dQil)
            QICN = QICN + dQil
            QLCN = QLCN - dQil
            T = T + (constants.latent_heat_sublimation - constants.latent_heat_vaporization) * dQil / constants.cpdry
        if T > constants.t_ice:
            dQil = -QICN
            dQil = min(0., dQil)
            QICN = QICN + dQil
            QLCN = QLCN - dQil
            T = T + (constants.latent_heat_sublimation - constants.latent_heat_vaporization) * dQil / constants.cpdry

def evap(
    DT_MOIST: Float,
    CCW_EVAP_EFF: Float,
    RHCRIT: Float,
    PLmb: FloatField,
    T: FloatField,
    Q: FloatField,
    QLCN: FloatField,
    QICN: FloatField,
    CLCN: FloatField,
    NACTL: FloatField,
    NACTI: FloatField,
    QST: FloatField,
    EVAPC: FloatField,
    RADIUS: FloatField,
    QCm: FloatField,
):
    with computation(PARALLEL), interval(...):
        EVAPC = Q
        # Evaporation of cloud water. DelGenio et al formulation (Eq.s 15-17, 1996, J. Clim., 9, 270-303)
        ES = 100.* PLmb * QST  / (constants.epsilon + (1.0 - constants.epsilon) * QST)  # (100's <-^ convert from mbar to Pa)
        RHx = min(Q/QST, 1.00)
        K1 = (constants.latent_heat_vaporization**2) * constants.RHO_W / (constants.K_COND * constants.rvap * (T**2))
        K2 = constants.rvap * T * constants.RHO_W / (constants.DIFFU * (1000. / PLmb) * ES)
        # Here, DIFFU is given for 1000 mb so 1000./PLmb accounts for increased diffusivity at lower pressure
        if CLCN > 0. and QLCN > 0.:
            QCm = QLCN / CLCN
        else:
            QCm = 0.
        RADIUS = cloud_effective_radius_liquid(PLmb, T, QCm, NACTL, NACTI)
        if RHx < RHCRIT and RADIUS > 0.0:
            EVAP = CCW_EVAP_EFF * QLCN * DT_MOIST * (RHCRIT - RHx) / ((K1 + K2) * RADIUS**2)
            EVAP = min(EVAP, QLCN)
        else:
            EVAP = 0.0
        QC = QLCN + QICN
        # if QC > 0.:
        #     CLCN = CLCN * (QC - EVAP) / QC
        # Q = Q + EVAP
        # QLCN = QLCN - EVAP
        # T = T - (constants.latent_heat_vaporization/constants.cpdry) * EVAP
        # EVAPC = (Q - EVAPC) / DT_MOIST

def subl(
    DT_MOIST: Float,
    CCI_EVAP_EFF: Float,
    RHCRIT: Float,
    PLmb: FloatField,
    T: FloatField,
    Q: FloatField,
    QLCN: FloatField,
    QICN: FloatField,
    CLCN: FloatField,
    NACTL: FloatField,
    NACTI: FloatField,
    QST: FloatField,
    EVAPC: FloatField,
):
    with computation(PARALLEL), interval(...):
        # Sublimation of cloud water. DelGenio et al formulation (Eq.s 15-17, 1996, J. Clim., 9, 270-303)
        ES = 100.* PLmb * QST  / (constants.epsilon + (1.0-constants.epsilon)*QST) # (100s <-^ convert from mbar to Pa)
        RHx = min(Q/QST, 1.00)
        K1 = (constants.latent_heat_vaporization**2) * constants.RHO_I / (constants.K_COND * constants.rvap * (T**2))
        K2 = constants.rvap * T * constants.RHO_I / (constants.DIFFU * (1000./PLmb) * ES)
        # Here, DIFFU is given for 1000 mb so 1000./PLmb accounts for increased diffusivity at lower pressure
        if CLCN > 0. and QICN > 0.:
            QCm = QICN / CLCN
        else:
            QCm = 0.
        radius = cloud_effective_radius_ice(PLmb, T, QCm, NACTL, NACTI)
        if RHx < RHCRIT and radius > 0.0:
            SUBL = CCI_EVAP_EFF * QICN * DT_MOIST * (RHCRIT - RHx) / ((K1 + K2) * radius**2)
            SUBL = min(SUBL, QICN)
        else:
            SUBL = 0.0
        QC = QLCN + QICN
        if QC > 0.:
            CLCN = CLCN * (QC - SUBL) / QC
        Q = Q + SUBL
        QICN = QICN - SUBL
        T = T - (constants.latent_heat_sublimation/constants.cpdry) * SUBL