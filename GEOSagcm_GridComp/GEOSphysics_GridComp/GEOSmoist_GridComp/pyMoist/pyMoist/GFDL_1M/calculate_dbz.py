from gt4py.cartesian.gtscript import f64

from ndsl import StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.gt4py import PARALLEL, computation, exp, interval, log10, tanh
from ndsl.dsl.typing import FloatField


def calculate_dbz(
    p_mb: FloatField,
    t: FloatField,
    vapor: FloatField,
    rain: FloatField,
    snow: FloatField,
    graupel: FloatField,
    model_has_snow: bool,
    ivarint: bool,
    iliqskin: bool,
    dbz: FloatField,
):
    from __externals__ import (
        ALPHA,
        CELKEL,
        FACTOR_G,
        FACTOR_R,
        FACTOR_S,
        GAMMA_SEVEN,
        GON,
        PI,
        R1,
        RD,
        RHO_G,
        RHO_R,
        RHO_S,
        RHOWAT,
        RN0_G,
        RN0_R,
        RN0_S,
        RON,
        RON2,
        RON_CONST1R,
        RON_CONST2R,
        RON_DELQR0,
        RON_MIN,
        RON_QR0,
        SON,
    )

    with computation(PARALLEL), interval(...):
        p_pascals = 100 * p_mb
        # Force all Q arrays to be 0.0 or greater.
        if vapor < 0:
            vapor = 0
        if rain < 0:
            rain = 0
        if snow < 0:
            snow = 0
        if graupel < 0:
            graupel = 0

        if model_has_snow is False:
            # COmpute snow based on rain and temperature
            if t < CELKEL:
                snow = rain
                rain = f64(0.0)  # is swapping rain to 64 bit here going to mess up things later?

        VIRTUAL_T = t * (f64(0.622) + vapor) / (f64(0.622) * (f64(1.0) + vapor))
        RHOAIR = p_pascals / (RD * VIRTUAL_T)

        # Adjust factor for brightband, where snow or graupel particle
        # scatters like liquid water (alpha=1.0) because it is assumed to
        # have a liquid skin.
        if iliqskin is True and t > CELKEL:
            FACTORB_S = FACTOR_S / ALPHA
            FACTORB_G = FACTOR_G / ALPHA
        else:
            FACTORB_S = FACTOR_S
            FACTORB_G = FACTOR_G

        # Calculate variable intercept parameters
        if ivarint is True:
            temp_c = min(f64(-0.001), t - CELKEL)
            sonv = min(f64(2.0e8), f64(2.0e6) * exp(f64(-0.12e0) * temp_c))
            gonv = GON
            if graupel > R1:
                gonv = f64(2.38) * (PI * RHO_G / (RHOAIR * graupel)) ** f64(0.92)
                gonv = max(f64(1.0e4), min(gonv, GON))
            ronv = RON2
            if rain > R1:
                ronv = RON_CONST1R * tanh((RON_QR0 - rain) / RON_DELQR0) + RON_CONST2R
        else:
            RONV = RN0_R
            SONV = RN0_S
            GONV = RN0_G

        # Total equivalent reflectivity factor (z_e, in mm^6 m^-3) is
        # the sum of z_e for each hydrometeor species:

        z_e = (
            FACTOR_R * (RHOAIR * rain) ** f64(1.75) / RONV ** f64(0.75)
            + FACTORB_S * (RHOAIR * snow) ** f64(1.75) / SONV ** f64(0.75)
            + FACTORB_G * (RHOAIR * graupel) ** f64(1.75) / GONV ** f64(0.75)
        )

        # Adjust small values of Z_e so that dBZ is no lower than -30
        z_e = max(z_e, f64(0.001))

        # Convert to dBZ
        dbz = f64(10.0) * log10(z_e)


class CalculateDBZ:
    """
    This routine computes equivalent reflectivity factor (in dBZ) at
    each model grid point.  In calculating Ze, the RIP algorithm makes
    assumptions consistent with those made in an early version
    (ca. 1996) of the bulk mixed-phase microphysical scheme in the MM5
    model (i.e., the scheme known as "Resiner-2").  For each species:

    1. Particles are assumed to be spheres of constant density.  The
    densities of rain drops, snow particles, and graupel particles are
    taken to be rho_r = rho_l = 1000 kg m^-3, rho_s = 100 kg m^-3, and
    rho_g = 400 kg m^-3, respectively. (l refers to the density of
    liquid water.)

    2. The size distribution (in terms of the actual diameter of the
    particles, rather than the melted diameter or the equivalent solid
    ice sphere diameter) is assumed to follow an exponential
    distribution of the form N(D) = N_0 * exp( lambda*D ).

    3. If ivarint=0, the intercept parameters are assumed constant
    (as in early Reisner-2), with values of 8x10^6, 2x10^7,
    and 4x10^6 m^-4, for rain, snow, and graupel, respectively.
    If ivarint=1, variable intercept parameters are used, as
    calculated in Thompson, Rasmussen, and Manning (2004, Monthly
    Weather Review, Vol. 132, No. 2, pp. 519-542.)

    4. If iliqskin=1, frozen particles that are at a temperature above
    freezing are assumed to scatter as a liquid particle.

    More information on the derivation of simulated reflectivity in
    RIP can be found in Stoelinga (2005, unpublished write-up).
    Contact Mark Stoelinga (stoeling@atmos.washington.edu) for a copy.

    Arguments:
        p (in): pressure (millibars)
        t (in): temperature (Kelvin)
        vapor (inout): water vapor mixing ratio (kg/kg)
        rain (inout): rain mass fraction (unitless)
        snow (inout): snow mass fraction (unitless)
        graupel (inout): graupel mass fraction (unitless)
        rain_is_snow (ino): control parameter, turns all rain into snow if true
        ivarint (in): control parameter (see above)
        iliqskin (in): control paramter (see above)
        dbz (inout): radar reflectivity (dBZ)
    """

    def __init__(
        self,
        stencil_factory: StencilFactory,
    ):

        # Constants used to calculate variable intercepts
        R1 = f64(1.0e-150)
        RON = f64(8.0e6)
        RON2 = f64(1.0e10)
        SON = f64(2.0e7)
        GON = f64(5.0e7)
        RON_MIN = f64(8.0e6)
        RON_QR0 = f64(0.00010)
        RON_DELQR0 = f64(0.25) * RON_QR0
        RON_CONST1R = (RON2 - RON_MIN) * f64(0.5)
        RON_CONST2R = (RON2 + RON_MIN) * f64(0.5)

        # Constant intercepts
        RN0_R = f64(8.0e6)
        RN0_S = f64(2.0e7)
        RN0_G = f64(4.0e6)

        # Other constants
        GAMMA_SEVEN = f64(720.0)
        RHOWAT = f64(1000.0)
        RHO_R = f64(RHOWAT)
        RHO_S = f64(100.0)
        RHO_G = f64(400.0)
        ALPHA = f64(0.224)
        CELKEL = f64(273.15)
        PI = f64(3.141592653589793)
        RD = f64(287.04)

        FACTOR_R = GAMMA_SEVEN * f64(1.0e18) * (f64(1.0) / (PI * RHO_R)) ** f64(1.75)
        FACTOR_S = (
            GAMMA_SEVEN * f64(1.0e18) * (f64(1.0) / (PI * RHO_S)) ** f64(1.75) * (RHO_S / RHOWAT) ** 2 * ALPHA
        )
        FACTOR_G = (
            GAMMA_SEVEN * f64(1.0e18) * (f64(1.0) / (PI * RHO_G)) ** f64(1.75) * (RHO_G / RHOWAT) ** 2 * ALPHA
        )

        # Construct stencil
        self.calculate_dbz = stencil_factory.from_dims_halo(
            func=calculate_dbz,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "R1": R1,
                "RON": RON,
                "RON2": RON2,
                "SON": SON,
                "GON": GON,
                "RON_MIN": RON_MIN,
                "RON_QR0": RON_QR0,
                "RON_DELQR0": RON_DELQR0,
                "RON_CONST1R": RON_CONST1R,
                "RON_CONST2R": RON_CONST2R,
                "RN0_R": RN0_R,
                "RN0_S": RN0_S,
                "RN0_G": RN0_G,
                "GAMMA_SEVEN": GAMMA_SEVEN,
                "RHOWAT": RHOWAT,
                "RHO_R": RHO_R,
                "RHO_S": RHO_S,
                "RHO_G": RHO_G,
                "ALPHA": ALPHA,
                "CELKEL": CELKEL,
                "PI": PI,
                "RD": RD,
                "FACTOR_R": FACTOR_R,
                "FACTOR_S": FACTOR_S,
                "FACTOR_G": FACTOR_G,
            },
        )

    def __call__(
        self,
        p: FloatField,
        t: FloatField,
        vapor: FloatField,
        rain: FloatField,
        snow: FloatField,
        graupel: FloatField,
        model_has_snow: bool,
        ivarint: bool,
        iliqskin: bool,
        dbz: FloatField,
    ):
        self.calculate_dbz(
            p=p,
            t=t,
            vapor=vapor,
            rain=rain,
            snow=snow,
            graupel=graupel,
            rain_is_snow=model_has_snow,
            ivarint=ivarint,
            iliqskin=iliqskin,
            dbz=dbz,
        )
