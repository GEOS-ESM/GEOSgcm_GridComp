import copy

from gt4py.cartesian.gtscript import (
    PARALLEL,
    computation,
    exp,
    f64,
    interval,
    log,
    sqrt,
)

import pyMoist.constants as constants
from ndsl import QuantityFactory, StencilFactory, orchestrate
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ, Int
from pyMoist.field_types import FloatField_NModes
from pyMoist.numerical_recipes import Erf


# 64 bit
ZERO_PAR = f64(1.0e-6)  # small non-zero value
AI = f64(0.0000594)
BI = f64(3.33)
CI = f64(0.0264)
DI = f64(0.0033)

BETAAI = f64(-2.262e3)
GAMAI = f64(5.113e6)
DELTAAI = f64(2.809e3)
DENSIC = f64(917.0)  # Ice crystal density in kgm-3

# Default precision
NN_MIN = 100.0e6
NN_MAX = 1000.0e6

# ACTFRAC_Mat constants - all 64 bit
PI = f64(3.141592653589793e00)
TWOPI = f64(2.0e00) * PI
SQRT2 = f64(1.414213562e00)
THREESQRT2BY2 = f64(1.5e00) * SQRT2

AVGNUM = f64(6.0221367e23)  # [1/mol]
RGASJMOL = f64(8.31451e00)  # [j/mol/k]
WMOLMASS = f64(18.01528e-03)  # molar mass of h2o [kg/mol]
AMOLMASS = f64(28.966e-03)  # molar mass of air     [kg/mol]
ASMOLMSS = f64(132.1406e-03)  # molar mass of nh42so4 [kg/mol]
DENH2O = f64(1.00e03)  # density of water [kg/m^3]
DENAMSUL = f64(1.77e03)  # density of pure ammonium sulfate [kg/m^3]
XNUAMSUL = f64(3.00e00)  # number of ions formed when the salt is dissolved in water [1]
PHIAMSUL = f64(1.000e00)  # osmotic coefficient value in a-r 1998. [1]
GRAVITY = f64(9.81e00)  # grav. accel. at the earth's surface [m/s/s]
HEATVAP = f64(40.66e03) / WMOLMASS  # latent heat of vap. for water and tnbp [j/kg]
CPAIR = f64(1006.0e00)  # heat capacity of air [j/kg/k]
T0DIJ = f64(273.15e00)  # reference temp. for dv [k]
P0DIJ = f64(101325.0e00)  # reference pressure for dv [pa]
DIJH2O0 = f64(0.211e-04)  # reference value of dv [m^2/s] (p&k,2nd ed., p.503)
DELTAV = f64(1.096e-07)  # vapor jump length [m]
DELTAT = f64(2.160e-07)  # thermal jump length [m]
ALPHAC = f64(1.000e00)  # condensation mass accommodation coefficient [1]
ALPHAT = f64(0.960e00)  # thermal accommodation coefficient [1]


def aer_activation_stencil(
    aero_dgn: FloatField_NModes,
    aero_num: FloatField_NModes,
    aero_sigma: FloatField_NModes,
    aero_hygroscopicity: FloatField_NModes,
    frland: FloatFieldIJ,
    nn_ocean: Float,
    nn_land: Float,
    t: FloatField,
    plo: FloatField,
    qicn: FloatField,
    qils: FloatField,
    qlcn: FloatField,
    qlls: FloatField,
    vvel: FloatField,
    tke: FloatField,
    nwfa: FloatField,
    nacti: FloatField,
    nactl: FloatField,
    nact: FloatField_NModes,
    ni: FloatField_NModes,
    rg: FloatField_NModes,
    sig0: FloatField_NModes,
    bibar: FloatField_NModes,
):
    """
    Compute the aerosol activation stencil.

    Parameters:
    aero_dgn (4D in): AeroProps aerosol geometric mean diameter.
    aero_num (4D in): AeroProps aerosol number concentration.
    nacti (3D out): Activated aerosol number concentration.
    t (3D in): Temperature field.
    plo (3D in): Pressure field.
    qicn (3D in): Ice cloud number concentration.
    qils (3D in): Ice liquid water content.
    qlcn (3D in): Cloud number concentration.
    qlls (3D in): Liquid water content.
    nn_land (1D in): Land-based aerosol activation number.
    frland (2D in): Fraction of land.
    nn_ocean (1D in): Ocean-based aerosol activation number.
    aero_hygroscopicity (4D in): Aerosol hygroscopicity parameter.
    nwfa (3D out): Number of activated aerosols.
    nactl (3D out): Number of activated aerosols in liquid clouds.
    vvel (3D in): Vertical velocity field.
    tke (3D in): Turbulent kinetic energy field.
    aero_sigma (4D in): AeroProps aerosol geometric standard deviation.
    nact (4D temporary in): Activated aerosol number concentration field.
    ni (4D temporary in): AeroProp ice crystal number concentration field.
    rg (4D temporary in): AeroProp geometric mean radius of aerosols.
    sig0 (4D temporary in): AeroProp aerosol geometric standard deviation field.
    bibar (4D temporary in): AeroProp Hygroscopicity parameter field.

    Returns:
    None
    """
    with computation(PARALLEL), interval(...):
        # Compute nwfa
        # Fortran AeroProps aero_kap is aero_hygroscopicity
        nfaux = 0.0
        n = 0
        while n < constants.N_MODES:
            if aero_hygroscopicity[0, 0, 0][n] > 0.4:
                nfaux += aero_num[0, 0, 0][n]
            n += 1
        nwfa = nfaux

        # Determine aerosol number concentration at cloud base
        nactl = nn_land * frland + nn_ocean * (1.0 - frland)
        nacti = nn_land * frland + nn_ocean * (1.0 - frland)

        tk = t  # [K]
        press = plo  # [Pa]
        air_den = press * 28.8e-3 / 8.31 / tk  # [kg/m3]
        qi = (qicn + qils) * 1.0e3  # [g/kg]
        ql = (qlcn + qlls) * 1.0e3  # [g/kg]
        wupdraft = vvel + sqrt(tke)

        # Liquid Clouds Calculation
        if (
            (tk >= (constants.MAPL_TICE - 40.0))
            and (plo > 10000.0)
            and (0.1 < wupdraft)
            and (wupdraft < 100.0)
        ):
            n = 0
            while n < constants.N_MODES:
                ni[0, 0, 0][n] = max(
                    aero_num[0, 0, 0][n] * air_den, ZERO_PAR
                )  # unit: [m-3]
                rg[0, 0, 0][n] = max(
                    aero_dgn[0, 0, 0][n] * 0.5 * 1.0e6, ZERO_PAR
                )  # unit: [um]
                sig0[0, 0, 0][n] = aero_sigma[0, 0, 0][n]  # unit: [um]
                bibar[0, 0, 0][n] = max(aero_hygroscopicity[0, 0, 0][n], ZERO_PAR)
                n += 1

            """
            12-12-06, DLW: Routine to calculate the activated fraction of the number
            and mass concentrations, as well as the number and mass
            concentrations activated for each of nmodes modes. The
            minimum dry radius for activation for each mode is also returned.

            Each mode is assumed to potentially contains 5 chemical species:
                (1) sulfate
                (2) BC
                (3) OC
                (4) mineral dust
                (5) sea salt

            The aerosol activation parameterizations are described in

                1. Abdul-Razzak et al.   1998, JGR, vol.103, p.6123-6131.
                2. Abdul-Razzak and Ghan 2000, JGR, vol.105, p.6837-6844.

            and values for many of the required parameters were taken from

                3. Ghan et al. 2001, JGR vol 106, p.5295-5316.

            With the density of sea salt set to the value used in ref. 3 (1900 kg/m^3),
            this routine yields values for the hygroscopicity parameters Bi in
            agreement with ref. 3.

            This routine is for the multiple-aerosol type parameterization.

            Original Fortran: ACTFRAC_Mat
            """
            # rdrp is the radius value used in eqs.(17) & (18) and was adjusted to
            # yield eta and zeta values close to those given in
            # a-z et al. 1998 figure 5.

            # tuned to approximate the results in figures 1-5 in a-z et al. 1998.
            rdrp = f64(0.105e-06)  # [m]

            # These variables are common to all modes and need only be computed once.
            dv = (
                DIJH2O0 * (P0DIJ / plo) * (tk / T0DIJ) ** f64(1.94e00)
            )  # [m^2/s] (p&k,2nd ed., p.503)
            surten = f64(76.10e-3) - f64(0.155e-3) * (tk - f64(273.15e00))  # [j/m^2]
            wpe = exp(
                f64(77.34491296)
                - f64(7235.424651) / tk
                - f64(8.2) * log(tk)
                + tk * f64(5.7113e-3)
            )  # [pa]
            dumw = sqrt(TWOPI * WMOLMASS / RGASJMOL / tk)  # [s/m]
            dvprime = dv / (
                (rdrp / (rdrp + DELTAV)) + (dv * dumw / (rdrp * ALPHAC))
            )  # [m^2/s] - eq. (17)
            xka = (f64(5.69) + f64(0.017) * (tk - f64(273.15))) * f64(
                418.4e-5
            )  # [j/m/s/k] (0.0238 j/m/s/k at 273.15 k)
            duma = sqrt(TWOPI * AMOLMASS / RGASJMOL / tk)  # [s/m]
            xkaprime = xka / (
                (rdrp / (rdrp + DELTAT))
                + (xka * duma / (rdrp * ALPHAT * DENH2O * CPAIR))
            )  # [j/m/s/k]
            g = f64(1.0) / (
                (DENH2O * RGASJMOL * tk) / (wpe * dvprime * WMOLMASS)
                + ((HEATVAP * DENH2O) / (xkaprime * tk))
                * ((HEATVAP * WMOLMASS) / (RGASJMOL * tk) - f64(1.0))
            )  # [m^2/s]
            a = (f64(2.0) * surten * WMOLMASS) / (DENH2O * RGASJMOL * tk)  # [m]
            alpha = (GRAVITY / (RGASJMOL * tk)) * (
                (WMOLMASS * HEATVAP) / (CPAIR * tk) - AMOLMASS
            )  # [1/m]
            gamma = (RGASJMOL * tk) / (wpe * WMOLMASS) + (
                WMOLMASS * HEATVAP * HEATVAP
            ) / (
                CPAIR * plo * AMOLMASS * tk
            )  # [m^3/kg]
            dum = sqrt(alpha * wupdraft / g)  # [1/m]
            zeta = f64(2.0) * a * dum / f64(3.0)  # [1]

            # These variables must be computed for each mode
            n = 0
            while n < constants.N_MODES:
                xlogsigm = log(sig0[0, 0, 0][n])  # [1]
                smax = f64(0.0)  # [1]
                sm = (f64(2.0) / sqrt(bibar[0, 0, 0][n])) * (
                    a / (3.0 * rg[0, 0, 0][n])
                ) ** f64(
                    1.5
                )  # [1]
                eta = dum ** 3 / (TWOPI * DENH2O * gamma * ni[0, 0, 0][n])  # [1]
                f1 = f64(0.5) * exp(2.50 * xlogsigm ** 2)  # [1]
                f2 = f64(1.0) + 0.25 * xlogsigm  # [1]
                smax = (
                    smax
                    + (
                        f1 * (zeta / eta) ** f64(1.5)
                        + f2 * (sm ** 2 / (eta + f64(3.0) * zeta)) ** f64(0.75)
                    )
                    / sm ** 2
                )  # [1] - eq. (6)
                n += 1

            smax = f64(1.0e00) / sqrt(smax)  # [1]
            n = 0
            u: f64 = f64(0.0)
            while n < constants.N_MODES:
                sm = (f64(2.0) / sqrt(bibar[0, 0, 0][n])) * (
                    a / (3.0 * rg[0, 0, 0][n])
                ) ** f64(
                    1.5
                )  # [1]
                xlogsigm = log(sig0[0, 0, 0][n])  # [1]
                ac = rg[0, 0, 0][n] * (sm / smax) ** f64(0.66666666666666667)  # [um]
                u = log(ac / rg[0, 0, 0][n]) / (SQRT2 * xlogsigm)  # [1]
                fracactn = f64(0.5) * (f64(1.0) - Erf(u))  # [1]
                nact[0, 0, 0][n] = fracactn * ni[0, 0, 0][n]  # [#/m^3]
                n += 1

            numbinit = 0.0
            nactl = 0.0
            n = 0
            while n < constants.N_MODES:
                numbinit += aero_num[0, 0, 0][n] * air_den
                nactl += nact[0, 0, 0][n]
                n += 1
            nactl = min(nactl, 0.99 * numbinit)

        # Ice Clouds Calculation
        if (tk <= constants.MAPL_TICE) and (
            qi > constants.FLOAT_TINY or ql > constants.FLOAT_TINY
        ):
            numbinit = 0.0
            n = 0
            while n < constants.N_MODES:
                # diameters > 0.5 microns
                if aero_dgn[0, 0, 0][n] >= 0.5e-6:
                    numbinit += aero_num[0, 0, 0][n]
                n += 1
            numbinit *= air_den  # [#/m3]
            # Number of activated IN following deMott (2010) [#/m3]
            nacti = (
                AI
                * ((constants.MAPL_TICE - tk) ** BI)
                * (numbinit ** (CI * (constants.MAPL_TICE - tk) + DI))
            )

        # apply limits for NACTL/NACTI
        if nactl < NN_MIN:
            nactl = NN_MIN
        if nactl > NN_MAX:
            nactl = NN_MAX
        if nacti < NN_MIN:
            nacti = NN_MIN
        if nacti > NN_MAX:
            nacti = NN_MAX


class AerActivation:
    """
    Class for aerosol activation computation.

    Attributes:
    stencil_factory (StencilFactory): Factory for creating stencil computations.
    quantity_factory (QuantityFactory): Factory for creating quantities.
    n_modes (Int): Number of aerosol modes.
    USE_AERSOL_NN (bool): Flag indicating whether to use neural network for aerosol.
    """

    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        n_modes: Int,
        USE_AERSOL_NN: bool,
    ) -> None:
        """
        Initialize the AerActivation class.

        Parameters:
        stencil_factory (StencilFactory): Factory for creating stencil computations.
        quantity_factory (QuantityFactory): Factory for creating quantities.
        n_modes (Int): Number of aerosol modes.
        USE_AERSOL_NN (bool): Flag indicating whether to use neural network for aerosol.

        Raises:
        NotImplementedError: If the number of modes is not equal to the expected number.
        NotImplementedError: If the neural network for aerosol is not implemented.
        """
        orchestrate(obj=self, config=stencil_factory.config.dace_config)

        if constants.N_MODES != n_modes:
            raise NotImplementedError(
                f"Coding limitation: {constants.N_MODES} modes are expected, "
                f"getting {n_modes}"
            )

        if not USE_AERSOL_NN:
            raise NotImplementedError("Non NN Aerosol not implemented")

        # Temporary buffers
        nmodes_quantity_factory = AerActivation.make_nmodes_quantity_factory(
            quantity_factory
        )
        self._nact = nmodes_quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_DIM, "n_modes"],
            units="n/a",
            dtype=Float,
        )
        self._ni = nmodes_quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_DIM, "n_modes"],
            units="n/a",
            dtype=Float,
        )
        self._rg = nmodes_quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_DIM, "n_modes"],
            units="n/a",
            dtype=Float,
        )
        self._sig0 = nmodes_quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_DIM, "n_modes"],
            units="n/a",
            dtype=Float,
        )
        self._bibar = nmodes_quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_DIM, "n_modes"],
            units="n/a",
            dtype=Float,
        )

        # Stencil
        self.aer_activation = stencil_factory.from_dims_halo(
            func=aer_activation_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    @staticmethod
    def make_nmodes_quantity_factory(ijk_quantity_factory: QuantityFactory):
        nmodes_quantity_factory = copy.deepcopy(ijk_quantity_factory)
        nmodes_quantity_factory.set_extra_dim_lengths(
            **{
                "n_modes": constants.N_MODES,
            }
        )
        return nmodes_quantity_factory

    def __call__(
        self,
        aero_dgn: FloatField_NModes,
        aero_num: FloatField_NModes,
        aero_hygroscopicity: FloatField_NModes,
        aero_sigma: FloatField_NModes,
        frland: FloatFieldIJ,
        nn_ocean: Float,
        nn_land: Float,
        t: FloatField,
        plo: FloatField,
        qicn: FloatField,
        qils: FloatField,
        qlcn: FloatField,
        qlls: FloatField,
        vvel: FloatField,
        tke: FloatField,
        nwfa: FloatField,
        nacti: FloatField,
        nactl: FloatField,
    ) -> None:
        """
        Compute aerosol activation by calling the stencil function.

        Parameters:
        aero_dgn (4D in): AeroProps aerosol geometric mean diameter.
        aero_num (4D in): AeroProps aerosol number concentration.
        aero_hygroscopicity (4D in): Aerosol hygroscopicity parameter.
        aero_sigma (4D in): AeroProps aerosol geometric standard deviation.
        frland (2D in): Fraction of land.
        nn_ocean (1D in): Ocean-based aerosol activation number.
        nn_land (1D in): Land-based aerosol activation number.
        t (3D in): Temperature field.
        plo (3D in): Pressure field.
        qicn (3D in): Ice cloud number concentration.
        qils (3D in): Ice liquid water content.
        qlcn (3D in): Cloud number concentration.
        qlls (3D in): Liquid water content.
        vvel (3D in): Vertical velocity field.
        tke (3D in): Turbulent kinetic energy field.
        nwfa (3D out): Number of activated aerosols.
        nactl (3D out): Number of activated aerosols in liquid clouds.
        nacti (3D out): Activated aerosol number concentration.

        Returns:
        None
        """
        self.aer_activation(
            aero_dgn,
            aero_num,
            aero_sigma,
            aero_hygroscopicity,
            frland,
            nn_ocean,
            nn_land,
            t,
            plo,
            qicn,
            qils,
            qlcn,
            qlls,
            vvel,
            tke,
            nwfa,
            nacti,
            nactl,
            self._nact,
            self._ni,
            self._rg,
            self._sig0,
            self._bibar,
        )
