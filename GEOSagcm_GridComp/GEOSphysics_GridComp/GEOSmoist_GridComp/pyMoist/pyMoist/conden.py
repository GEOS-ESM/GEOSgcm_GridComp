import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import computation, interval, PARALLEL, sin
import pyMoist.radiation_coupling_constants as radconstants
import pyMoist.pyMoist_constants as constants
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField, Float, IntField 
from ndsl import StencilFactory, QuantityFactory 
from pyMoist.saturation.qsat import QSat, QSat_Float

@gtscript.function
def exnerfn(
    p: Float,
)-> Float:
    
    return (p / 100000.0) ** (radconstants.MAPL_RDRY / radconstants.MAPL_CPDRY)


@gtscript.function
def geos_qsat(
    temps: Float,
    ps: Float, 
)-> Float:
    
    # This is just a placeholder until GEOS_QSAT is ported!!!
    e_s0 = 6.11         # Reference vapor pressure (hPa)
    a = 17.67           # Constant for temperature in the exponential function
    b = 273.15          # Reference temperature (K)
    c = 29.65           # Temperature constant in the denominator (K)
    q_s0 = 0.62         # Ratio of molecular weights of water vapor and dry air (dimensionless)
    
    # Saturation vapor pressure calculation
    e_s = e_s0 * temps
    
    # Saturation specific humidity calculation
    q_s = q_s0 * e_s / (ps - e_s)
    
    return q_s


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
        icefrct_c = sin(0.5 * radconstants.MAPL_PI * (1.00 - (temp - constants.JaT_ICE_ALL) / (constants.JaT_ICE_MAX - constants.JaT_ICE_ALL)))
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
           icefrct_m = sin(0.5 * radconstants.MAPL_PI * (1.00 - (temp - constants.lT_ICE_ALL) / (constants.lT_ICE_MAX - constants.lT_ICE_ALL)))
        else:
            icefrct_m = 0.00
        icefrct_m = max(min(icefrct_m, 1.00), 0.00) ** constants.lICEFRPWR
    else:
        if temp <= constants.oT_ICE_ALL:
           icefrct_m = 1.000
        elif temp > constants.oT_ICE_ALL and temp <= constants.oT_ICE_MAX:
           icefrct_m = sin(0.5 * radconstants.MAPL_PI * (1.00 - (temp - constants.oT_ICE_ALL) / (constants.oT_ICE_MAX - constants.oT_ICE_ALL)))
        else:
            icefrct_m = 0.00
        icefrct_m = max(min(icefrct_m, 1.00), 0.00) ** constants.oICEFRPWR
    ice_frac = icefrct_m * (1.0-cnv_frc) + icefrct_c * cnv_frc
    return ice_frac


def condensation(
    p: FloatField, 
    thl: FloatField, 
    qt: FloatField,
    th: FloatField, 
    qv: FloatField, 
    ql: FloatField,
    qi: FloatField, 
    rvls: FloatField,
    id_check: IntField,
):

    with computation(PARALLEL), interval(...):

        tc = thl * exnerfn(p)
        
        icefrac = ice_fraction(tc, 0.0, 0.0)
        nu = icefrac
        leff = (1.0 - nu) * radconstants.MAPL_ALHL + nu * radconstants.MAPL_ALHS # Effective latent heat
        temps = tc
        ps = p
        qs = geos_qsat(temps, ps / 100.0) # Saturation specific humidity
        rvls = qs
        
        if qs >= qt:      # no condensation
            id_check = 0
            qv = qt
            qc = 0.0
            ql = 0.0
            qi = 0.0
            th = thl
        else:             # condensation
            iteration = 0
            while iteration < 10:
                temps = temps + ((tc - temps) * radconstants.MAPL_CPDRY / leff + qt - rvls) / (radconstants.MAPL_CPDRY / leff + (constants.rdry / constants.rvap) * leff * rvls / (radconstants.MAPL_RDRY * temps * temps))
                qs = geos_qsat(temps, ps / 100.0)
                rvls = qs
                iteration+=1
            qc = max(qt - qs, 0.0)
            qv = qt - qc
            ql = qc * (1.0 - nu)
            qi = nu * qc
            th = temps / exnerfn(p)
            if abs((temps - (leff / radconstants.MAPL_CPDRY) * qc) - tc) >= 1.0:
                id_check = 1
            else:
                id_check = 0
                


class Conden:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
    ):
      
        self._condensation = stencil_factory.from_dims_halo(
            func=condensation,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(self,
        pifc0_test: FloatField, 
        thl0bot_test: FloatField, 
        qt0bot_test: FloatField,
        thj_test: FloatField, 
        qvj_test: FloatField, 
        qlj_test: FloatField,
        qij_test: FloatField, 
        qse_test: FloatField,
        id_check_test: IntField,
    ):

        """
            Calculate condensation process.

            Parameters:
            p (3D in): Environmental pressure [ Pa ].
            thl (3D in): 
            qt (3D in):
            th (3D out):
            qv (3D out):
            ql (3D out):
            qi (3D out):
            rvls (3D out):
            id_check (3D out): Indicates if condensation process occurs [ N/A ].
        """
        self.qsat = QSat(
            self.stencil_factory,
            self.quantity_factory,
            formulation=formulation,
            use_table_lookup=use_table_lookup,
        )

        # These are args for QSat_Float
        self.qsat._ese
        self.qsat._esw
        self.qsat._esx

        self._condensation(
                pifc0_test, 
                thl0bot_test, 
                qt0bot_test, 
                thj_test, 
                qvj_test, 
                qlj_test, 
                qij_test, 
                qse_test, 
                id_check_test)
