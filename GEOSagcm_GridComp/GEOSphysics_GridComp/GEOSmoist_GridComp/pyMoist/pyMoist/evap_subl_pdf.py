import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import computation, interval, PARALLEL, FORWARD, atan, sin, tan, sqrt, tanh, exp
from ndsl.boilerplate import get_factories_single_tile_numpy
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField, FloatFieldIJ, Float, IntField, IntFieldIJ
from ndsl import StencilFactory, QuantityFactory, orchestrate
import numpy as np
import pyMoist.pyMoist_constants as constants



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
    ice_frac = icefrct_m * (1.0-cnv_frc) + icefrct_m * cnv_frc
    return ice_frac

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
    testminrhcrit: FloatField,
    testPLEmb: FloatField,
    testturnrhcrit: FloatField,
    testrhcritpart1: FloatField,
    testrhcritpart2: FloatField,
    testrhcrit: FloatField,
    testalpha: FloatField,
    testtancalc: FloatField,
): 
    with computation(PARALLEL), interval(...):
        # Send the condensates through the pdf after convection
        facEIS = max(0.0, min(1.0, EIS/10.0))**2
        # determine combined minrhcrit in stable/unstable regimes
        minrhcrit = (1.0-dw_ocean)*(1.0-facEIS) + (1.0-dw_land)*facEIS
        testminrhcrit = (1.0-dw_ocean)*(1.0-facEIS) + (1.0-dw_land)*facEIS
        # determine the turn pressure using the LCL
        if TURNRHCRIT_PARAM <= 0:
            turnrhcrit = PLmb_at_klcl -250 # implementation of for loop needs to go here
            testturnrhcrit = PLmb_at_klcl -250
        else:
            turnrhcrit = TURNRHCRIT_PARAM
            testturnrhcrit = TURNRHCRIT_PARAM

    # lower turn from maxrhcrit=1.0 # implemented in multiple "with" statements to deal with hybrid indexing
    with computation(PARALLEL), interval(0,-1):
        testPLEmb = PLEmb_top
        # Use Slingo-Ritter (1985) formulation for critical rel ative humidity
        RHCRIT = 1.0
        testrhcrit = 1.0
        testtancalc = 0.0
        if PLmb <= turnrhcrit:
            RHCRIT = minrhcrit
            testrhcrit = minrhcrit
        else:
            testrhcritpart1 = minrhcrit + (1.0-minrhcrit)/(19.)
            testrhcritpart2 = ((atan( (2.*(PLmb-turnrhcrit)/(PLEmb_top-turnrhcrit)-1.) *
                tan(20.*constants.PI/21.-0.5*constants.PI) )
                + 0.5*constants.PI) * 21./constants.PI - 1.)
            RHCRIT = minrhcrit + (1.0-minrhcrit)/(19.) * ((atan( (2.*(PLmb-turnrhcrit)/(PLEmb_top-turnrhcrit)-1.) *
                tan(20.*constants.PI/21.-0.5*constants.PI) )
                + 0.5*constants.PI) * 21./constants.PI - 1.)
            testrhcrit = minrhcrit + (1.0-minrhcrit)/(19.) * ((atan( (2.*(PLmb-turnrhcrit)/(PLEmb_top-turnrhcrit)-1.) *
                tan(20.*constants.PI/21.-0.5*constants.PI) )
                + 0.5*constants.PI) * 21./constants.PI - 1.)
            testtancalc = 1.0
    with computation(PARALLEL), interval(-1,None):
        testPLEmb = PLEmb_top
        testtancalc = 0.0
        # lower turn from maxrhcrit=1.0
        if PLmb <= turnrhcrit:
            RHCRIT = minrhcrit
            testrhcrit = minrhcrit
        else:
            RHCRIT = 1.0
            testrhcrit = 1.0
    with computation(PARALLEL), interval(...):
        # include grid cell area scaling and limit RHcrit to > 70%\
        ALPHA = max(0.0,min(0.30, (1.0-RHCRIT)*sqrt(sqrt(AREA/1.e10))))
        testalpha = max(0.0,min(0.30, (1.0-RHCRIT)*sqrt(sqrt(AREA/1.e10))))

def meltfrz(
    dt: Float,
    cnv_frc: FloatFieldIJ,
    srf_type: FloatFieldIJ,
    T: FloatField,
    QLCN: FloatField,
    QICN: FloatField,
    testtancalc: FloatField,
):
    with computation(PARALLEL), interval(...):
        compvar = T
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
        if T == compvar:
            testtancalc = 1
        else:
            testtancalc = 0
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
        self._meltfrz = stencil_factory.from_dims_halo(
            func=meltfrz,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )
        self._tmp = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        self._minrhcrit = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._PLEmb_top = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._PLmb_at_klcl = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        self._halo = stencil_factory.grid_indexing.n_halo
        self._k_mask = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
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
                testminrhcrit: FloatField,
                testPLEmb: FloatField,
                testturnrhcrit: FloatField,
                testrhcritpart1: FloatField,
                testrhcritpart2: FloatField,
                testrhcrit: FloatField,
                testalpha: FloatField,
                testtancalc: FloatField,
                 ):
        # Theoretically, the following stencil and for loop should provide the same (correct) result. However, they both provide different incorrect results
        # The for loop is currently closer to being correct, with only the k level incorrect, so I am using that.
        self._get_last(PLEmb, self._tmp, self._PLEmb_top)

        # Temporary implementation of hybrid_index_2dout.py, perhaps not working as indended (backend issue), will need to be addressed at later date 
        self._hybrid_index_2dout(PLmb, self._k_mask, KLCL, self._PLmb_at_klcl)
        
        self._initial_calc(EIS=EIS, dw_land=dw_land, dw_ocean=dw_ocean, TURNRHCRIT_PARAM=TURNRHCRIT_PARAM, minrhcrit=self._minrhcrit, 
                           PLmb_at_klcl=self._PLmb_at_klcl, PLmb=PLmb, PLEmb_top=self._PLEmb_top, AREA=AREA, testminrhcrit=testminrhcrit, testPLEmb=testPLEmb, 
                           testturnrhcrit=testturnrhcrit, testrhcritpart1=testrhcritpart1, testrhcritpart2=testrhcritpart2,
                           testrhcrit=testrhcrit, testalpha=testalpha, testtancalc=testtancalc)

        self._meltfrz(DT_MOIST, CNV_FRC, SRF_TYPE, T, QLCN, QICN, testtancalc)

"""

if __name__ == "__main__": #Should only run when testing. DOE NOT WORK IN CURRENT STATE ---> need to change PLEmb[...] to PLEmb.view[...] to handle a Quantity instead of np array, currently line 150

    domain = (24, 24, 72)
    nhalo = 0
    stcil_fctry, ijk_qty_fctry = get_factories_single_tile_numpy(
        domain[0], domain[1], domain[2], nhalo
    )

    #Initialize dummy variables of same name and arbitrary domain
    EIS = ijk_qty_fctry.zeros([X_DIM, Y_DIM], "n/a")
    EIS.view[:, :] = 42
    dw_land = 25.
    dw_ocean = 60.
    TURNRHCRIT_PARAM = -9999.
    PLmb = ijk_qty_fctry.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
    PLmb.view[:,:,:] = 850
    KLCL = ijk_qty_fctry.zeros([X_DIM, Y_DIM], "n/a")
    KLCL.view[:,:] = 2
    PLEmb = ijk_qty_fctry.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
    PLEmb.data[:,:,-2] = 300
    AREA = ijk_qty_fctry.zeros([X_DIM, Y_DIM], "n/a")
    AREA.view[:,:] = domain[0] * domain[1]

    testminrhcrit = ijk_qty_fctry.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
    testPLEmb = ijk_qty_fctry.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
    testturnrhcrit = ijk_qty_fctry.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
    testrhcritpart1 = ijk_qty_fctry.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
    testrhcritpart2 = ijk_qty_fctry.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
    testrhcrit = ijk_qty_fctry.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
    testalpha = ijk_qty_fctry.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
    testtancalc = ijk_qty_fctry.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

    code = evap_subl_pdf(stcil_fctry, ijk_qty_fctry)
    code(EIS, dw_land, dw_ocean, TURNRHCRIT_PARAM, PLmb, KLCL, PLEmb, AREA, testminrhcrit, testPLEmb, testturnrhcrit, testrhcritpart1, testrhcritpart2, testrhcrit, testalpha, testtancalc)

    print(f"Output:\n{testalpha}\n")

"""