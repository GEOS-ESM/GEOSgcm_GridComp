import copy
from typing import Optional

from gt4py.cartesian.gtscript import i32

from ndsl import QuantityFactory, StencilFactory, orchestrate
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.gt4py import PARALLEL, GlobalTable, computation, floor, function, interval
from ndsl.dsl.typing import Float, FloatField, Int
from pyMoist.field_types import GlobalTable_saturaion_tables
from pyMoist.saturation_tables.constants import (
    DEGSUBS,
    ERFAC,
    ESFAC,
    MAPL_TICE,
    MAX_MIXING_RATIO,
    TABLESIZE,
    TMAXTBL,
    TMINLQU,
    TMINTBL,
    TMIX,
)
from pyMoist.saturation_tables.formulation import SaturationFormulation
from pyMoist.saturation_tables.tables.main import get_table


@function
def saturation_specific_humidity_frozen_surface(
    ese: GlobalTable_saturaion_tables,  # type: ignore
    frz: Float,
    t: Float,
    p: Float,
    pressure_correction: bool,
    compute_dq: bool,
):
    """
    Computes saturation specific humitidy over liquid water surface, using
    data from saturation pressure tables.

    Arguments:
        ese (in): saturation pressure table in Pascals, specifics unknown
        frz (in): saturation pressure at reference temperature (273.16 K)
        t (in): temperature in Kelvin
        p (in): pressure in Pascals
        pressure_correction (in): trigger for pressure correction
        compute_dq (in): trigger for computing derivative saturation specific humidity

    Returns:
        qsat (out): saturation specific humidity
        dqsat (out): derivative saturation specific humidity with respect to temperature
    """
    dq = 0.0
    if t <= TMINTBL:
        qs = ese.A[0]  # type: ignore
        if compute_dq:
            ddq = 0.0
    elif t >= MAPL_TICE:
        qs = frz
        if compute_dq:
            ddq = 0.0
    else:
        t = (t - TMINTBL) * DEGSUBS + 1
        t_integer = int(floor(t))
        ddq = ese.A[t_integer] - ese.A[t_integer - 1]  # type: ignore
        qs = (t - t_integer) * ddq + ese.A[t_integer - 1]  # type: ignore

    if pressure_correction == True:
        if p > qs:
            dd = ESFAC / (p - (1.0 - ESFAC) * qs)
            qs = qs * dd
            if compute_dq:
                dq = ddq * ERFAC * p * dd * dd
        else:
            qs = MAX_MIXING_RATIO
            if compute_dq:
                dq = 0.0
    else:
        if compute_dq:
            dq = ddq

    return qs, dq


@function
def saturation_specific_humidity_liquid_surface(
    esw: GlobalTable_saturaion_tables,  # type: ignore
    lqu: Float,
    t: Float,
    p: Float,
    pressure_correction: bool = False,
    compute_dq: bool = False,
):
    """
    Computes saturation specific humitidy over liquid water surface, using
    data from saturation pressure tables.

    Arguments:
        esw (in): saturation pressure table in Pascals, specifics unknown
        lqu (in): saturation pressure at reference temperature (233.16 K)
        t (in): temperature in Kelvin
        p (in): pressure in Pascals
        pressure_correction (in): trigger for pressure correction
        compute_dq (in): trigger for computing derivative saturation specific humidity

    Returns:
        qsat (out): saturation specific humidity
        dqsat (out): derivative saturation specific humidity with respect to temperature
    """
    dqsat = 0.0
    if t <= TMINLQU:
        qsat = lqu
        if compute_dq == True:
            ddq = 0.0
    elif t >= TMAXTBL:
        TABLESIZE_MINUS_1: i32 = TABLESIZE - 1
        qsat = esw.A[TABLESIZE_MINUS_1]  # type: ignore
        if compute_dq == True:
            ddq = 0.0
    else:
        t = (t - TMINTBL) * DEGSUBS + 1
        t_integer = int(floor(t))
        ddq = esw.A[t_integer] - esw.A[t_integer - 1]  # type: ignore
        qsat = (t - t_integer) * ddq + esw.A[t_integer - 1]  # type: ignore

    if pressure_correction == True:
        if p > qsat:
            dd = ESFAC / (p - (1.0 - ESFAC) * qsat)
            qsat = qsat * dd
            if compute_dq == True:
                dqsat = ddq * ERFAC * p * dd * dd
        else:
            qsat = MAX_MIXING_RATIO
            if compute_dq == True:
                dqsat = 0.0
    else:
        if compute_dq == True:
            dqsat = ddq

    return qsat, dqsat


# Function version of GEOS_Qsat
@function
def saturation_specific_humidity(
    t: Float,
    p: Float,
    ese: GlobalTable_saturaion_tables,  # type: ignore
    esx: GlobalTable_saturaion_tables,  # type: ignore
    use_ramp: bool = False,
    ramp: Float = -999.0,
):
    """
    Compute saturation specific humidity and derivative saturation specific humidity
    with respect to temperature from saturation pressure tables.

    Tables must be initalized before use.

    Arguments:
        t (in): temperature in Kelvin
        p (in): pressure in Pascals
        ese (in): saturation pressure table in Pascals, specifics unknown
        esx (in): saturation presure table in Pascals, specifics unknown
        use_ramp (in): trigger for "ramp" option. details unknown
        ramp (in): parameter used for "ramp" option. details unknown

    Returns:
        qsat (out): saturation specific humidity
        dqsat (out): derivative saturation specific humidity with respect to temperature
    """
    if use_ramp == True:
        uramp = -abs(ramp)
    else:
        uramp = TMIX

    if t <= TMINTBL:
        t = TMINTBL
    elif t >= TMAXTBL - 0.001:
        t = TMAXTBL - 0.001

    t = (t - TMINTBL) * DEGSUBS + 1
    t_integer = i32(floor(t))
    IT_MINUS_1 = t_integer - 1

    if uramp == TMIX:
        dq = esx.A[t_integer] - esx.A[IT_MINUS_1]  # type: ignore
        qsat = (t - t_integer) * dq + esx.A[IT_MINUS_1]  # type: ignore
    else:
        dq = ese.A[t_integer] - ese.A[IT_MINUS_1]  # type: ignore
        qsat = (t - t_integer) * dq + ese.A[IT_MINUS_1]  # type: ignore

    if p <= qsat:
        qsat = MAX_MIXING_RATIO
        dqsat = 0.0
    else:
        dd = 1.0 / (p - (1.0 - ESFAC) * qsat)
        qsat = ESFAC * qsat * dd
        dqsat = ESFAC * dq * DEGSUBS * p * (dd * dd)

    return qsat, dqsat


# Stencils implement GEOS_Qsat subroutine from GEOS_Utilities.F90
def QSat_FloatField(
    ese: GlobalTable_saturaion_tables,  # type: ignore
    esx: GlobalTable_saturaion_tables,  # type: ignore
    T: FloatField,
    PL: FloatField,
    QSAT: FloatField,
    DQSAT: FloatField,
    RAMP: FloatField,
):
    from __externals__ import FILL_DQSAT, USE_PASCALS, USE_RAMP

    with computation(PARALLEL), interval(...):
        if USE_RAMP:
            URAMP = -abs(RAMP)
        else:
            URAMP = TMIX

        if USE_PASCALS:
            PP = PL
        else:
            PP = PL * 100.0

        if T <= TMINTBL:
            TI = TMINTBL
        elif T >= TMAXTBL - 0.001:
            TI = TMAXTBL - 0.001
        else:
            TI = T

        TI = (TI - TMINTBL) * DEGSUBS + 1
        IT = i32(floor(TI))
        IT_MINUS_1: i32 = (
            IT - 1
        )  # dace backend does not allow for [IT - 1] indexing because of cast to int

        if URAMP == TMIX:
            DQ = esx[0][IT] - esx[0][IT_MINUS_1]  # type: ignore
            QSAT = (TI - IT) * DQ + esx[0][IT_MINUS_1]  # type: ignore
        else:
            DQ = ese[0][IT] - ese[0][IT_MINUS_1]  # type: ignore
            QSAT = (TI - IT) * DQ + ese[0][IT_MINUS_1]  # type: ignore
        if PP <= QSAT:
            QSAT = MAX_MIXING_RATIO
            if FILL_DQSAT:
                DQSAT = 0.0
        else:
            DD = 1.0 / (PP - (1.0 - ESFAC) * QSAT)
            QSAT = ESFAC * QSAT * DD
            if FILL_DQSAT:
                DQSAT = ESFAC * DQ * DEGSUBS * PP * (DD * DD)


class QSat:
    """Traditional satuation specific humidity.
    In Fortran: GEOS_Utilities:QSat

    Uses various formulations of the saturation vapor pressure to compute
    the saturation specific humidity for temperature T and pressure PL.

    The function get_table called in __init__ creates tables using exact
    calculations (equivalent to UTBL=False in Fortran).
    `get_table` can be called again for further exact computations.
    QSat_FloatField computes QSat for an entire field using table lookups
    (equivalent to UTBL=True in Fortran).

    For temperatures <= TMIX (-20C)
    the calculation is done over ice; for temperatures >= ZEROC (0C)
    the calculation is done over liquid water; and in between these values,
    it interpolates linearly between the two.

    The optional RAMP is the width of this
    ice/water ramp (i.e., TMIX = ZEROC-RAMP); its default is 20.

    If PASCALS is true, PL is
    assumed to be in Pa; if false or not present, it is assumed to be in mb.
    """

    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        formulation: SaturationFormulation = SaturationFormulation.Staars,
        use_ramp: bool = False,
        use_pascals: bool = False,
        fill_dqsat: bool = False,
    ) -> None:
        self.extra_dim_quantity_factory = self.make_extra_dim_quantity_factory(
            quantity_factory
        )

        self.ese = self.extra_dim_quantity_factory.zeros([Z_DIM, "table_axis"], "n/a")
        self.esw = self.extra_dim_quantity_factory.zeros([Z_DIM, "table_axis"], "n/a")
        self.esx = self.extra_dim_quantity_factory.zeros([Z_DIM, "table_axis"], "n/a")

        self.table = get_table(formulation)
        self.ese.view[:] = self.table.ese
        self.esw.view[:] = self.table.esw
        self.esx.view[:] = self.table.esx

        self._RAMP = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._DQSAT = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        self.QSat = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        orchestrate(obj=self, config=stencil_factory.config.dace_config)
        self._QSat_FloatField = stencil_factory.from_dims_halo(
            func=QSat_FloatField,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "USE_RAMP": use_ramp,
                "USE_PASCALS": use_pascals,
                "FILL_DQSAT": fill_dqsat,
            },
        )

    @staticmethod
    def make_extra_dim_quantity_factory(ijk_quantity_factory: QuantityFactory):
        extra_dim_quantity_factory = copy.deepcopy(ijk_quantity_factory)
        extra_dim_quantity_factory.set_extra_dim_lengths(
            **{
                "table_axis": TABLESIZE,
            }
        )
        return extra_dim_quantity_factory

    def __call__(
        self,
        T: FloatField,
        PL: FloatField,
        DQSat: Optional[FloatField] = None,
        use_table_lookup: bool = True,
    ):
        if DQSat:
            dqsat = DQSat
        else:
            dqsat = self._DQSAT

        if use_table_lookup:
            self._QSat_FloatField(
                self.ese,
                self.esx,
                T,
                PL,
                self.QSat,
                dqsat,
                self._RAMP,
            )

        if not use_table_lookup:
            raise NotImplementedError(
                "Saturation calculation: exact formulation not available,"
                " only table look up"
            )
