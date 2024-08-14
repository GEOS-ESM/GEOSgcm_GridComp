import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import computation, interval, PARALLEL, FORWARD, atan, sin, tan, sqrt, tanh, exp, log10
from ndsl.boilerplate import get_factories_single_tile_numpy
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField, FloatFieldIJ, Float, IntField, IntFieldIJ
from ndsl import StencilFactory, QuantityFactory, orchestrate
import numpy as np
import pyMoist.pyMoist_constants as constants

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
def dqsat(
    TL: FloatField,
    PL: FloatField,
    RAMP: FloatField = False,
    PASCALS: FloatField = False,
    QSAT: FloatField = False,
): #unfinished
    if RAMP:
        URAMP = -abs(RAMP)
    else:
        URAMP = constants.TIMX
    
    if PASCALS:
        PP = PL
    else:
        PP = PL * 100
    
    if (URAMP == constants.TIMX or URAMP == 0) and constants.UTBL == True:
        if constants.FIRST == True:
            pass #need to make new class in different file that initializes table (GEOS_Utilities.F90 subroutine ESINIT)

    if TL < constants.TMINTBL:
        TI = constants.TMINTBL
    elif TL >= (constants.TMAXTBL - 0.001):
        TI = constants.TMAXTBL - 0.001
    else:
        TI = TL
    
    TT = (TI - constants.TMINTBL) * constants.DEGSUBS + 1
    # IT = int(TT) # need to find replacement for int, perhaps Int?

    # if URAMP == constants.TMIX:
    #     DQQ = 

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
    cnv_frc: FloatFieldIJ,
    srf_type: FloatFieldIJ,
    PLmb: FloatField,
    ZL0: FloatField,
    Q: FloatField,
    QLLS: FloatField,
    QLCN: FloatField,
    QILS: FloatField,
    QICN: FloatField,
    T: FloatField,
    CLLS: FloatField,
    CLCN: FloatField,
    NACTL: FloatField,
    NACTI: FloatField,
    WSL: FloatField,
    WQT: FloatField,
    SL2: FloatField,
    QT2: FloatField,
    SLQT: FloatField,
    W3: FloatField,
    W2: FloatField,
    QT3: FloatField,
    SL3: FloatField,
    PDF_A: FloatField,
    PDFITERS: FloatField,
    WTHV2: FloatField,
    WQL: FloatField,
): #unfinished
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
        #DQS = dqsat(T, PLmb, QSn)

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

class evap_subl_pdf:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
    ):
        orchestrate(obj=self, config=stencil_factory.config.dace_config)
        self._get_last = stencil_factory.from_dims_halo(
            func=get_last,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )
        self._hybrid_index_2dout = stencil_factory.from_dims_halo(
            func=hybrid_index_2dout,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )
        self._initial_calc = stencil_factory.from_dims_halo(
            func=initial_calc,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )
        self._hystpdf = stencil_factory.from_dims_halo(
            func=hystpdf,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )
        self._meltfrz = stencil_factory.from_dims_halo(
            func=meltfrz,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )
        self._evap = stencil_factory.from_dims_halo(
            func=evap,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )
        self._subl = stencil_factory.from_dims_halo(
            func=subl,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )
        self._tmp = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        self._minrhcrit = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._PLEmb_top = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._PLmb_at_klcl = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        self._halo = stencil_factory.grid_indexing.n_halo
        self._k_mask = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._alpha = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._evapc = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        for i in range(0, self._k_mask.view[:].shape[0]):
            for j in range(0, self._k_mask.view[:].shape[1]):
                for k in range(0, self._k_mask.view[:].shape[2]):
                    self._k_mask.view[i, j, k] = k + 1
    def __call__(self,
        EIS: FloatFieldIJ,
        dw_land: Float,
        dw_ocean: Float,
        TURNRHCRIT_PARAM: Float,
        PLmb: FloatField,
        KLCL: FloatFieldIJ,
        PLEmb: FloatField,
        AREA: FloatFieldIJ,
        DT_MOIST: Float,
        CNV_FRC: FloatFieldIJ,
        SRF_TYPE: FloatFieldIJ,
        T: FloatField,
        QLCN: FloatField,
        QICN: FloatField,
        QLLS: FloatField,
        QILS: FloatField,
        CCW_EVAP_EFF: Float,
        CCI_EVAP_EFF: Float,
        Q: FloatField,
        CLCN: FloatField,
        NACTL: FloatField,
        NACTI: FloatField,
        QST: FloatField,
        RADIUS: FloatField,
        QCm: FloatField,
    ):
        # Theoretically, the following stencil and for loop should provide the same (correct) result. However, they both provide different incorrect results
        # The for loop is currently closer to being correct, with only the k level incorrect, so I am using that.
        self._get_last(PLEmb, self._tmp, self._PLEmb_top)

        # Temporary implementation of hybrid_index_2dout.py, perhaps not working as indended (backend issue), will need to be addressed at later date 
        self._hybrid_index_2dout(PLmb, self._k_mask, KLCL, self._PLmb_at_klcl)
        
        self._initial_calc(EIS=EIS, dw_land=dw_land, dw_ocean=dw_ocean, TURNRHCRIT_PARAM=TURNRHCRIT_PARAM, minrhcrit=self._minrhcrit, 
                           PLmb_at_klcl=self._PLmb_at_klcl, PLmb=PLmb, PLEmb_top=self._PLEmb_top, AREA=AREA, ALPHA=self._alpha)

        # self._meltfrz(DT_MOIST, CNV_FRC, SRF_TYPE, T, QLCN, QICN)
        # self._meltfrz(DT_MOIST, CNV_FRC, SRF_TYPE, T, QLLS, QILS)

        RHCRIT = Float(1.0)
        self._evap(DT_MOIST, CCW_EVAP_EFF, RHCRIT, PLmb, T, Q, QLCN, QICN, CLCN, NACTL, NACTI, QST, self._evapc, RADIUS, QCm)

        #self._subl(DT_MOIST, CCW_EVAP_EFF, RHCRIT, PLmb, T, Q, QLCN, QICN, CLCN, NACTL, NACTI, QST, self._evapc)