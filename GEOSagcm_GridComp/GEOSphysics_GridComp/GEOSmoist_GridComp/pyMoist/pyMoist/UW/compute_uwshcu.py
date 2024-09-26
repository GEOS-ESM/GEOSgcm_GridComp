import copy

import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import computation, interval, PARALLEL, FORWARD, sin
import pyMoist.pyMoist_constants as constants
import pyMoist.radiation_coupling_constants as radconstants
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField, Float, Int, IntField
from ndsl import StencilFactory, QuantityFactory
from pyMoist.types import FloatField_NTracers
from pyMoist.saturation.qsat import QSat, QSat_Float, FloatField_Extra_Dim
from pyMoist.saturation.formulation import SaturationFormulation


@gtscript.function
def slope(
    kmask: Float,
    field: Float,
    field_above: Float,
    field_below: Float,
    p0: Float,
    p0_above: Float,
    p0_below: Float,
):
    """
    Calculates slope of a given field.

    Parameters:
    kmask (Float): K-level mask (e.g., 1 for k=0, 2 for k=1,71, 3 for k=72).
    field (Float): Field of interest.
    field_above (Float): 1 k-level above field (e.g., field[0,0,1]).
    field_below (Float): 1 k-level below field (e.g., field[0,0,-1]).
    p0 (Float): Pressure
    p0_above (Float): 1 k-level above p0 (e.g., p0[0,0,1]).
    p0_below (Float): 1 k-level below p0 (e.g., p0[0,0,-1]).

    Returns:
    Slope: Slope of the field of interest.
    """
    if kmask == 1:
        value = (field_above - field) / (p0_above - p0)
        if value > 0.0:
            slope = max(0.0, value)
        else:
            slope = min(0.0, value)
    elif kmask == 2:
        above_value = (field_above - field) / (p0_above - p0)
        below_value = (field - field_below) / (p0 - p0_below)
        if above_value > 0.0:
            slope = max(0.0, min(above_value, below_value))
        else:
            slope = min(0.0, max(above_value, below_value))
    else:
        above_value = (field_above - field) / (p0_above - p0)
        below_value = (field - field_below) / (p0 - p0_below)
        if above_value > 0.0:
            slope = max(0.0, min(above_value, below_value))
        else:
            slope = min(0.0, max(above_value, below_value))

    return slope


@gtscript.function
def exnerfn(
    p: Float,
) -> Float:

    return (p / 100000.0) ** (radconstants.MAPL_RGAS / constants.cp)


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
        icefrct_c = sin(
            0.5
            * constants.PI
            * (
                1.00
                - (temp - constants.JaT_ICE_ALL)
                / (constants.JaT_ICE_MAX - constants.JaT_ICE_ALL)
            )
        )
    else:
        icefrct_c = 0.00
    icefrct_c = max(min(icefrct_c, 1.00), 0.00) ** constants.aICEFRPWR
    # Sigmoidal functions like figure 6b/6c of Hu et al 2010, doi:10.1029/2009JD012384
    if srf_type == 2.0:
        if temp <= constants.JiT_ICE_ALL:
            icefrct_m = 1.000
        elif temp > constants.JiT_ICE_ALL and temp <= constants.JiT_ICE_MAX:
            icefrct_m = 1.00 - (temp - constants.JiT_ICE_ALL) / (
                constants.JiT_ICE_MAX - constants.JiT_ICE_ALL
            )
        else:
            icefrct_m = 0.00
        icefrct_m = max(min(icefrct_m, 1.00), 0.00) ** constants.iICEFRPWR
    elif srf_type > 1.0:
        if temp <= constants.lT_ICE_ALL:
            icefrct_m = 1.000
        elif temp > constants.lT_ICE_ALL and temp <= constants.lT_ICE_MAX:
            icefrct_m = sin(
                0.5
                * constants.PI
                * (
                    1.00
                    - (temp - constants.lT_ICE_ALL)
                    / (constants.lT_ICE_MAX - constants.lT_ICE_ALL)
                )
            )
        else:
            icefrct_m = 0.00
        icefrct_m = max(min(icefrct_m, 1.00), 0.00) ** constants.lICEFRPWR
    else:
        if temp <= constants.oT_ICE_ALL:
            icefrct_m = 1.000
        elif temp > constants.oT_ICE_ALL and temp <= constants.oT_ICE_MAX:
            icefrct_m = sin(
                0.5
                * constants.PI
                * (
                    1.00
                    - (temp - constants.oT_ICE_ALL)
                    / (constants.oT_ICE_MAX - constants.oT_ICE_ALL)
                )
            )
        else:
            icefrct_m = 0.00
        icefrct_m = max(min(icefrct_m, 1.00), 0.00) ** constants.oICEFRPWR
    ice_frac = icefrct_m * (1.0 - cnv_frc) + icefrct_c * cnv_frc
    return ice_frac


@gtscript.function
def conden(
    p: Float,
    thl: Float,
    qt: Float,
    ese: FloatField_Extra_Dim,
    esx: FloatField_Extra_Dim,
):

    tc = thl * exnerfn(p)

    nu = ice_fraction(tc, 0.0, 0.0)
    leff = (1.0 - nu) * constants.xlv + (nu * constants.xls)  # Effective latent heat
    temps = tc
    ps = p
    qs, _ = QSat_Float(ese, esx, temps, ps / 100.0)  # Saturation specific humidity
    rvls = qs

    if qs >= qt:  # no condensation
        id_check = 0
        qv = qt
        qc = 0.0
        ql = 0.0
        qi = 0.0
        th = thl
    else:  # condensation
        iteration = 0
        while iteration < 10:
            temps = temps + ((tc - temps) * constants.cp / leff + qt - rvls) / (
                constants.cp / leff
                + constants.ep2 * leff * rvls / (constants.r * temps * temps)
            )
            qs, _ = QSat_Float(ese, esx, temps, ps / 100.0)
            rvls = qs
            iteration += 1
        qc = max(qt - qs, 0.0)
        qv = qt - qc
        ql = qc * (1.0 - nu)
        qi = nu * qc
        th = temps / exnerfn(p)
        if abs((temps - (leff / constants.cp) * qc) - tc) >= 1.0:
            id_check = 1
        else:
            id_check = 0

    return th, qv, ql, qi, rvls, id_check


def compute_uwshcu(
    dotransport: Float,
    pifc0_in: FloatField,
    pmid0_in: FloatField,
    zmid0_in: FloatField,
    exnmid0_in: FloatField,
    u0_in: FloatField,
    v0_in: FloatField,
    qv0_in: FloatField,
    ql0_in: FloatField,
    qi0_in: FloatField,
    th0_in: FloatField,
    tr0_inout: FloatField_NTracers,
    tr0: FloatField_NTracers,
    ssthl0: FloatField,
    ssqt0: FloatField,
    ssu0: FloatField,
    ssv0: FloatField,
    sstr0: FloatField_NTracers,
    thj: FloatField,
    qvj: FloatField,
    qlj: FloatField,
    qij: FloatField,
    qse: FloatField,
    id_check: IntField,
    kmask: FloatField,
    ese: FloatField_Extra_Dim,
    esx: FloatField_Extra_Dim,
):
    """
    University of Washington Shallow Convection Scheme

    Described in Park and Bretherton. 2008. J. Climate :

    'The University of Washington shallow convection and
    moist turbulent schemes and their impact on climate
    simulations with the Community Atmosphere Model'

    Coded in CESM by Sungsu Park. Oct.2005. May.2008.

    Coded in GEOS by Nathan Arnold. July 2016.

    For general questions, email sungsup@ucar.edu or
    sungsu@atmos.washington.edu

    For GEOS-specific questions, email nathan.arnold@nasa.gov
    """

    """
    Add description of variables
    """

    # Start Main Calculation
    with computation(PARALLEL), interval(0, 1):
        pmid0 = pmid0_in
        pmid0_above = pmid0_in[0, 0, 1]
        u0 = u0_in
        u0_above = u0_in[0, 0, 1]
        v0 = v0_in
        v0_above = v0_in[0, 0, 1]
        qv0 = qv0_in
        qv0_above = qv0_in[0, 0, 1]
        ql0 = ql0_in
        ql0_above = ql0_in[0, 0, 1]
        qi0 = qi0_in
        qi0_above = qi0_in[0, 0, 1]
        qt0 = qv0 + ql0 + qi0
        qt0_above = qv0_above + ql0_above + qi0_above
        exnmid0 = exnmid0_in
        exnmid0_above = exnmid0_in[0, 0, 1]
        t0 = th0_in * exnmid0
        t0_above = th0_in[0, 0, 1] * exnmid0_above
        thl0 = (
            t0
            - ((constants.xlv * ql0) / constants.cp)
            - ((constants.xls * qi0) / constants.cp)
        ) / exnmid0
        thl0_above = (
            t0_above
            - ((constants.xlv * ql0_above) / constants.cp)
            - ((constants.xls * qi0_above) / constants.cp)
        ) / exnmid0_above

        if dotransport == 1.0:
            n = 0
            # Loop over tracers
            while n < constants.ncnst:
                tr0[0, 0, 0][n] = tr0_inout[0, 0, 0][n]
                n += 1

        ssthl0 = slope(
            kmask, thl0, thl0_above, thl0_above, pmid0, pmid0_above, pmid0_above
        )
        ssqt0 = slope(kmask, qt0, qt0_above, qt0_above, pmid0, pmid0_above, pmid0_above)
        ssu0 = slope(kmask, u0, u0_above, u0_above, pmid0, pmid0_above, pmid0_above)
        ssv0 = slope(kmask, v0, v0_above, v0_above, pmid0, pmid0_above, pmid0_above)

        if dotransport == 1.0:
            n = 0
            while n < constants.ncnst:
                if True:
                    sstr0[0, 0, 0][n] = slope(
                        kmask,
                        tr0[0, 0, 0][n],
                        tr0[0, 0, 1][n],
                        tr0[0, 0, 1][n],
                        pmid0,
                        pmid0_above,
                        pmid0_above,
                    )
                n += 1

    with computation(PARALLEL), interval(1, -1):
        pmid0 = pmid0_in
        pmid0_above = pmid0_in[0, 0, 1]
        pmid0_below = pmid0_in[0, 0, -1]
        zmid0 = zmid0_in
        u0 = u0_in
        u0_above = u0_in[0, 0, 1]
        u0_below = u0_in[0, 0, -1]
        v0 = v0_in
        v0_above = v0_in[0, 0, 1]
        v0_below = v0_in[0, 0, -1]
        qv0 = qv0_in
        qv0_above = qv0_in[0, 0, 1]
        qv0_below = qv0_in[0, 0, -1]
        ql0 = ql0_in
        ql0_above = ql0_in[0, 0, 1]
        ql0_below = ql0_in[0, 0, -1]
        qi0 = qi0_in
        qi0_above = qi0_in[0, 0, 1]
        qi0_below = qi0_in[0, 0, -1]

        if dotransport == 1.0:
            n = 0
            # Loop over tracers
            while n < constants.ncnst:
                tr0[0, 0, 0][n] = tr0_inout[0, 0, 0][n]
                n += 1

        exnmid0 = exnmid0_in
        exnmid0_above = exnmid0_in[0, 0, 1]
        exnmid0_below = exnmid0_in[0, 0, -1]
        t0 = th0_in * exnmid0
        t0_above = th0_in[0, 0, 1] * exnmid0_above
        t0_below = th0_in[0, 0, -1] * exnmid0_below
        qt0 = qv0 + ql0 + qi0
        qt0_above = qv0_above + ql0_above + qi0_above
        qt0_below = qv0_below + ql0_below + qi0_below
        thl0 = (
            t0
            - ((constants.xlv * ql0) / constants.cp)
            - ((constants.xls * qi0) / constants.cp)
        ) / exnmid0
        thl0_above = (
            t0_above
            - ((constants.xlv * ql0_above) / constants.cp)
            - ((constants.xls * qi0_above) / constants.cp)
        ) / exnmid0_above
        thl0_below = (
            t0_below
            - ((constants.xlv * ql0_below) / constants.cp)
            - ((constants.xls * qi0_below) / constants.cp)
        ) / exnmid0_below

        ssthl0 = slope(
            kmask, thl0, thl0_above, thl0_below, pmid0, pmid0_above, pmid0_below
        )
        ssqt0 = slope(kmask, qt0, qt0_above, qt0_below, pmid0, pmid0_above, pmid0_below)
        ssu0 = slope(kmask, u0, u0_above, u0_below, pmid0, pmid0_above, pmid0_below)
        ssv0 = slope(kmask, v0, v0_above, v0_below, pmid0, pmid0_above, pmid0_below)

        if dotransport == 1.0:
            n = 0
            while n < constants.ncnst:
                if True:
                    sstr0[0, 0, 0][n] = slope(
                        kmask,
                        tr0[0, 0, 0][n],
                        tr0[0, 0, 1][n],
                        tr0[0, 0, -1][n],
                        pmid0,
                        pmid0_above,
                        pmid0_below,
                    )
                n += 1

    with computation(PARALLEL), interval(-1, None):
        pmid0 = pmid0_in[0, 0, -1]
        pmid0_above = pmid0_in
        pmid0_below = pmid0_in[0, 0, -2]
        u0 = u0_in[0, 0, -1]
        u0_above = u0_in
        u0_below = u0_in[0, 0, -2]
        v0 = v0_in[0, 0, -1]
        v0_above = v0_in
        v0_below = v0_in[0, 0, -2]
        qv0 = qv0_in[0, 0, -1]
        qv0_above = qv0_in
        qv0_below = qv0_in[0, 0, -2]
        ql0 = ql0_in[0, 0, -1]
        ql0_above = ql0_in
        ql0_below = ql0_in[0, 0, -2]
        qi0 = qi0_in[0, 0, -1]
        qi0_above = qi0_in
        qi0_below = qi0_in[0, 0, -2]

        exnmid0 = exnmid0_in[0, 0, -1]
        exnmid0_above = exnmid0_in
        exnmid0_below = exnmid0_in[0, 0, -2]
        t0 = th0_in[0, 0, -1] * exnmid0
        t0_above = th0_in * exnmid0_above
        t0_below = th0_in[0, 0, -2] * exnmid0_below
        qt0 = qv0 + ql0 + qi0
        qt0_above = qv0_above + ql0_above + qi0_above
        qt0_below = qv0_below + ql0_below + qi0_below
        thl0 = (
            t0
            - ((constants.xlv * ql0) / constants.cp)
            - ((constants.xls * qi0) / constants.cp)
        ) / exnmid0
        thl0_above = (
            t0_above
            - ((constants.xlv * ql0_above) / constants.cp)
            - ((constants.xls * qi0_above) / constants.cp)
        ) / exnmid0_above
        thl0_below = (
            t0_below
            - ((constants.xlv * ql0_below) / constants.cp)
            - ((constants.xls * qi0_below) / constants.cp)
        ) / exnmid0_below

        if dotransport == 1.0:
            n = 0
            # Loop over tracers
            while n < constants.ncnst:
                tr0[0, 0, 0][n] = tr0_inout[0, 0, 0][n]
                n += 1

        ssthl0 = slope(
            kmask, thl0, thl0_above, thl0_below, pmid0, pmid0_above, pmid0_below
        )
        ssqt0 = slope(kmask, qt0, qt0_above, qt0_below, pmid0, pmid0_above, pmid0_below)
        ssu0 = slope(kmask, u0, u0_above, u0_below, pmid0, pmid0_above, pmid0_below)
        ssv0 = slope(kmask, v0, v0_above, v0_below, pmid0, pmid0_above, pmid0_below)

        # Calculate slope for each tracer by hand
        if dotransport == 1.0:
            n = 0
            while n < constants.ncnst:
                if True:
                    sstr0[0, 0, 0][n] = slope(
                        kmask,
                        tr0[0, 0, -1][n],
                        tr0[0, 0, 0][n],
                        tr0[0, 0, -2][n],
                        pmid0,
                        pmid0_above,
                        pmid0_below,
                    )
                n += 1

    with computation(PARALLEL), interval(...):
        pmid0 = pmid0_in
        exnmid0 = exnmid0_in
        qv0 = qv0_in
        ql0 = ql0_in
        qi0 = qi0_in
        qt0 = qv0 + ql0 + qi0
        t0 = th0_in * exnmid0
        thl0 = (
            t0 - constants.xlv * ql0 / constants.cp - constants.xls * qi0 / constants.cp
        ) / exnmid0
        zvir = 0.609  # r_H2O/r_air-1
        pifc0 = pifc0_in
        thl0bot = thl0 + ssthl0 * (pifc0 - pmid0)
        qt0bot = qt0 + ssqt0 * (pifc0 - pmid0)

        thj, qvj, qlj, qij, qse, id_check = conden(pifc0, thl0bot, qt0bot, ese, esx)
        # Raise an error if id_check = 1
        thv0bot = thj * (1.0 + zvir * qvj - qlj - qij)
        thvl0bot = thl0bot * (1.0 + zvir * qt0bot)

        thl0top = thl0 + ssthl0 * (pifc0_in[0, 0, 1] - pmid0)
        qt0top = qt0 + ssqt0 * (pifc0_in[0, 0, 1] - pmid0)

    with computation(PARALLEL), interval(0, -1):
        thj, qvj, qlj, qij, qse, id_check = conden(
            pifc0_in[0, 0, 1], thl0top, qt0top, ese, esx
        )
        # Raise an error if id_check = 1
        thv0top = thj * (1.0 + zvir * qvj - qlj - qij)
        thvl0top = thl0top * (1.0 + zvir * qt0top)

    with computation(PARALLEL), interval(-1, None):
        thv0top = thv0bot
        thvl0top = thvl0bot


class ComputeUwshcu:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        ncnst: Int,
        formulation: SaturationFormulation = SaturationFormulation.Staars,
    ) -> None:
        """
        Initialize the ComputeUwshcu class.

        Parameters:
        stencil_factory (StencilFactory): Factory for creating stencil computations.
        quantity_factory (QuantityFactory): Factory for creating quantities.
        ncnst (Int): Number of tracers.
        """

        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory
        grid_indexing = stencil_factory.grid_indexing
        self._compute_uwshcu = self.stencil_factory.from_dims_halo(
            func=compute_uwshcu,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._k_mask = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        for i in range(0, self._k_mask.view[:].shape[0]):
            for j in range(0, self._k_mask.view[:].shape[1]):
                for k in range(0, self._k_mask.view[:].shape[2]):
                    if k == 0:
                        self._k_mask.view[i, j, k] = 1
                    elif k == 71:
                        self._k_mask.view[i, j, k] = 3
                    else:
                        self._k_mask.view[i, j, k] = 2

    @staticmethod
    def make_ntracers_quantity_factory(ijk_quantity_factory: QuantityFactory):
        ntracers_quantity_factory = copy.deepcopy(ijk_quantity_factory)
        ntracers_quantity_factory.set_extra_dim_lengths(
            **{
                "ntracers": constants.ncnst,
            }
        )
        return ntracers_quantity_factory

    def __call__(
        self,
        dotransport: Float,
        pifc0_in: FloatField,
        pmid0_in: FloatField,
        zmid0_in: FloatField,
        exnmid0_in: FloatField,
        u0_in: FloatField,
        v0_in: FloatField,
        qv0_in: FloatField,
        ql0_in: FloatField,
        qi0_in: FloatField,
        th0_in: FloatField,
        tr0_inout: FloatField_NTracers,
        tr0_test: FloatField_NTracers,
        ssthl0_test: FloatField,
        ssqt0_test: FloatField,
        ssu0_test: FloatField,
        ssv0_test: FloatField,
        sstr0_test: FloatField_NTracers,
        thj_test: FloatField,
        qvj_test: FloatField,
        qlj_test: FloatField,
        qij_test: FloatField,
        qse_test: FloatField,
        id_check_test: IntField,
        formulation: SaturationFormulation = SaturationFormulation.Staars,
    ):

        self.qsat = QSat(
            self.stencil_factory,
            self.quantity_factory,
            formulation=formulation,
        )

        self._compute_uwshcu(
            dotransport=dotransport,
            pifc0_in=pifc0_in,
            pmid0_in=pmid0_in,
            zmid0_in=zmid0_in,
            exnmid0_in=exnmid0_in,
            u0_in=u0_in,
            v0_in=v0_in,
            qv0_in=qv0_in,
            ql0_in=ql0_in,
            qi0_in=qi0_in,
            th0_in=th0_in,
            tr0_inout=tr0_inout,
            tr0=tr0_test,
            ssthl0=ssthl0_test,
            ssqt0=ssqt0_test,
            ssu0=ssu0_test,
            ssv0=ssv0_test,
            sstr0=sstr0_test,
            thj=thj_test,
            qvj=qvj_test,
            qlj=qlj_test,
            qij=qij_test,
            qse=qse_test,
            id_check=id_check_test,
            kmask=self._k_mask,
            ese=self.qsat.ese,
            esx=self.qsat.esx,
        )
