import copy

from gt4py.cartesian.gtscript import PARALLEL, computation, exp, interval, log, sqrt

import pyMoist.aer_activation_constants as constants
from ndsl import QuantityFactory, StencilFactory, orchestrate
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ, Int
from pyMoist.numerical_recipes import Erf
from pyMoist.types import FloatField_NModes


def aer_activation_stencil(
    aero_dgn: FloatField_NModes,
    aero_num: FloatField_NModes,
    nacti: FloatField,
    t: FloatField,
    plo: FloatField,
    qicn: FloatField,
    qils: FloatField,
    qlcn: FloatField,
    qlls: FloatField,
    nn_land: Float,
    frland: FloatFieldIJ,
    nn_ocean: Float,
    aero_hygroscopicity: FloatField_NModes,
    nwfa: FloatField,
    nactl: FloatField,
    vvel: FloatField,
    tke: FloatField,
    aero_sigma: FloatField_NModes,
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
        while n < constants.n_modes:
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
            while n < constants.n_modes:
                ni[0, 0, 0][n] = max(
                    aero_num[0, 0, 0][n] * air_den, constants.ZERO_PAR
                )  # unit: [m-3]
                rg[0, 0, 0][n] = max(
                    aero_dgn[0, 0, 0][n] * 0.5 * 1.0e6, constants.ZERO_PAR
                )  # unit: [um]
                sig0[0, 0, 0][n] = aero_sigma[0, 0, 0][n]  # unit: [um]
                bibar[0, 0, 0][n] = max(
                    aero_hygroscopicity[0, 0, 0][n], constants.ZERO_PAR
                )
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
            """
            # rdrp is the radius value used in eqs.(17) & (18) and was adjusted to
            # yield eta and zeta values close to those given in
            # a-z et al. 1998 figure 5.

            # tuned to approximate the results in figures 1-5 in a-z et al. 1998.
            rdrp = 0.105e-06  # [m]

            # These variables are common to all modes and need only be computed once.
            dv = (
                constants.DIJH2O0
                * (constants.P0DIJ / plo)
                * (tk / constants.T0DIJ) ** 1.94e00
            )  # [m^2/s] (p&k,2nd ed., p.503)
            surten = 76.10e-3 - 0.155e-3 * (tk - 273.15e00)  # [j/m^2]
            wpe = exp(
                77.34491296 - 7235.424651 / tk - 8.2 * log(tk) + tk * 5.7113e-3
            )  # [pa]
            dumw = sqrt(
                constants.TWOPI * constants.WMOLMASS / constants.RGASJMOL / tk
            )  # [s/m]
            dvprime = dv / (
                (rdrp / (rdrp + constants.DELTAV))
                + (dv * dumw / (rdrp * constants.ALPHAC))
            )  # [m^2/s] - eq. (17)
            xka = (
                5.69 + 0.017 * (tk - 273.15)
            ) * 418.4e-5  # [j/m/s/k] (0.0238 j/m/s/k at 273.15 k)
            duma = sqrt(
                constants.TWOPI * constants.AMOLMASS / constants.RGASJMOL / tk
            )  # [s/m]
            xkaprime = xka / (
                (rdrp / (rdrp + constants.DELTAT))
                + (
                    xka
                    * duma
                    / (rdrp * constants.ALPHAT * constants.DENH2O * constants.CPAIR)
                )
            )  # [j/m/s/k]
            g = 1.0 / (
                (constants.DENH2O * constants.RGASJMOL * tk)
                / (wpe * dvprime * constants.WMOLMASS)
                + ((constants.HEATVAP * constants.DENH2O) / (xkaprime * tk))
                * (
                    (constants.HEATVAP * constants.WMOLMASS) / (constants.RGASJMOL * tk)
                    - 1.0
                )
            )  # [m^2/s]
            a = (2.0 * surten * constants.WMOLMASS) / (
                constants.DENH2O * constants.RGASJMOL * tk
            )  # [m]
            alpha = (constants.GRAVITY / (constants.RGASJMOL * tk)) * (
                (constants.WMOLMASS * constants.HEATVAP) / (constants.CPAIR * tk)
                - constants.AMOLMASS
            )  # [1/m]
            gamma = (constants.RGASJMOL * tk) / (wpe * constants.WMOLMASS) + (
                constants.WMOLMASS * constants.HEATVAP * constants.HEATVAP
            ) / (constants.CPAIR * plo * constants.AMOLMASS * tk)  # [m^3/kg]
            dum = sqrt(alpha * wupdraft / g)  # [1/m]
            zeta = 2.0 * a * dum / 3.0  # [1]

            # These variables must be computed for each mode
            n = 0
            while n < constants.n_modes:
                xlogsigm = log(sig0[0, 0, 0][n])  # [1]
                smax = 0.0  # [1]
                sm = (2.0 / sqrt(bibar[0, 0, 0][n])) * (
                    a / (3.0 * rg[0, 0, 0][n])
                ) ** 1.5  # [1]
                eta = dum**3 / (
                    constants.TWOPI * constants.DENH2O * gamma * ni[0, 0, 0][n]
                )  # [1]
                f1 = 0.5 * exp(2.50 * xlogsigm**2)  # [1]
                f2 = 1.0 + 0.25 * xlogsigm  # [1]
                smax = (
                    smax
                    + (
                        f1 * (zeta / eta) ** 1.5
                        + f2 * (sm**2 / (eta + 3.0 * zeta)) ** 0.75
                    )
                    / sm**2
                )  # [1] - eq. (6)
                n += 1

            smax = 1.0e00 / sqrt(smax)  # [1]
            n = 0
            u = 0.0
            while n < constants.n_modes:
                sm = (2.0 / sqrt(bibar[0, 0, 0][n])) * (
                    a / (3.0 * rg[0, 0, 0][n])
                ) ** 1.5  # [1]
                xlogsigm = log(sig0[0, 0, 0][n])  # [1]
                ac = rg[0, 0, 0][n] * (sm / smax) ** 0.66666666666666667  # [um]
                u = log(ac / rg[0, 0, 0][n]) / (constants.SQRT2 * xlogsigm)  # [1]
                fracactn = 0.5 * (1.0 - Erf(u))  # [1]
                nact[0, 0, 0][n] = fracactn * ni[0, 0, 0][n]  # [#/m^3]
                n += 1

            numbinit = 0.0
            nactl = 0.0
            n = 0
            while n < constants.n_modes:
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
            while n < constants.n_modes:
                # diameters > 0.5 microns
                if aero_dgn[0, 0, 0][n] >= 0.5e-6:
                    numbinit += aero_num[0, 0, 0][n]
                n += 1
            numbinit *= air_den  # [#/m3]
            # Number of activated IN following deMott (2010) [#/m3]
            nacti = (
                constants.AI
                * ((constants.MAPL_TICE - tk) ** constants.BI)
                * (
                    numbinit
                    ** (constants.CI * (constants.MAPL_TICE - tk) + constants.DI)
                )
            )

        # apply limits for NACTL/NACTI
        if nactl < constants.NN_MIN:
            nactl = constants.NN_MIN
        if nactl > constants.NN_MAX:
            nactl = constants.NN_MAX
        if nacti < constants.NN_MIN:
            nacti = constants.NN_MIN
        if nacti > constants.NN_MAX:
            nacti = constants.NN_MAX


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

        if constants.n_modes != n_modes:
            raise NotImplementedError(
                f"Coding limitation: 14 modes are expected, getting {n_modes}"
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
                "n_modes": constants.n_modes,
            }
        )
        return nmodes_quantity_factory

    def __call__(
        self,
        aero_dgn: FloatField_NModes,
        aero_num: FloatField_NModes,
        nacti: FloatField,
        t: FloatField,
        plo: FloatField,
        qicn: FloatField,
        qils: FloatField,
        qlcn: FloatField,
        qlls: FloatField,
        nn_land: Float,
        frland: FloatFieldIJ,
        nn_ocean: Float,
        aero_hygroscopicity: FloatField_NModes,
        nwfa: FloatField,
        nactl: FloatField,
        vvel: FloatField,
        tke: FloatField,
        aero_sigma: FloatField_NModes,
    ) -> None:
        """
        Compute aerosol activation by calling the stencil function.

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

        Returns:
        None
        """
        self.aer_activation(
            aero_dgn,
            aero_num,
            nacti,
            t,
            plo,
            qicn,
            qils,
            qlcn,
            qlls,
            nn_land,
            frland,
            nn_ocean,
            aero_hygroscopicity,
            nwfa,
            nactl,
            vvel,
            tke,
            aero_sigma,
            self._nact,
            self._ni,
            self._rg,
            self._sig0,
            self._bibar,
        )
