import copy
from typing import Optional

import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import (
    PARALLEL,
    computation,
    floor,
    interval,
    i32,
)

from ndsl import QuantityFactory, StencilFactory, orchestrate
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float, FloatField
from pyMoist.saturation.constants import (
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
from pyMoist.saturation.formulation import SaturationFormulation
from pyMoist.saturation.table import get_table


# FloatField with extra dimension initialized to handle table data
# This is a temporary solution. This solution creates a KxTABLESIZE array
# At some point a proper solution will be implemented, enabling a 1xTABLESIZE array
# Allocation is currently fixed to TABLESIZE constant. Fortran has some examples of
# flexible table sixes (larger than TABLESIZE, with increased granulatiry).
# Current implementation does not allow for flexible table sizes
# so if needed this will have to be implemented in another way.
# DEV NOTE: we have to "type ignore" the indexation due to a misread of mypy
FloatField_Extra_Dim = gtscript.Field[gtscript.K, (Float, (int(TABLESIZE)))]


# The following two functions, QSat_Float_Liquid and QSat_Float_Ice are,
# together, a non-equivalent replacement for QSat_Float.
# QSat_Float interpolates a smooth transition between liquid water and ice
# (based on the parameter TMIX), while the Liquic/Ice versions have a
# hard transition at zero C.
@gtscript.function
def QSat_Float_Liquid(
    esw: FloatField_Extra_Dim,  # type: ignore  # type: ignore
    estlqu: Float,
    TL: Float,
    PL: Float = -999.0,
    DQ_trigger: bool = False,
):
    # qsatlqu.code with UTBL = True
    if TL <= TMINLQU:
        QS = estlqu
        if DQ_trigger:
            DDQ = 0.0
    elif TL >= TMAXTBL:
        QS = esw[0][TABLESIZE - 1]  # type: ignore
        if DQ_trigger:
            DDQ = 0.0
    else:
        TT = (TL - TMINTBL) * DEGSUBS + 1
        IT = i32(TT)
        IT_MINUS_1 = (
            IT - 1
        )  # dace backend does not allow for [IT - 1] indexing because of cast to int
        DDQ = esw[0][IT] - esw[0][IT_MINUS_1]  # type: ignore
        QS = (TT - IT) * DDQ + esw[0][IT_MINUS_1]  # type: ignore

    if PL != -999.0:
        if PL > QS:
            DD = ESFAC / (PL - (1.0 - ESFAC) * QS)
            QS = QS * DD
            if DQ_trigger:
                DQ = DDQ * ERFAC * PL * DD * DD
        else:
            QS = MAX_MIXING_RATIO
            if DQ_trigger:
                DQ = 0.0
    else:
        if DQ_trigger:
            DQ = DDQ

    return QS, DQ


@gtscript.function
def QSat_Float_Ice(
    ese: FloatField_Extra_Dim,  # type: ignore
    estfrz: Float,
    TL: Float,
    PL: Float = -999.0,
    DQ_trigger: bool = False,
):
    # qsatice.code with UTBL = True
    if TL <= TMINTBL:
        QS = ese[0][0]  # type: ignore
        if DQ_trigger:
            DDQ = 0.0
    elif TL >= MAPL_TICE:
        QS = estfrz
        if DQ_trigger:
            DDQ = 0.0
    else:
        TT = (TL - TMINTBL) * DEGSUBS + 1
        IT = i32(floor(TT))
        IT_MINUS_1 = (
            IT - 1
        )  # dace backend does not allow for [IT - 1] indexing because of cast to int
        DDQ = ese[0][IT] - ese[0][IT_MINUS_1]  # type: ignore
        QS = (TT - IT) * DDQ + ese[0][IT_MINUS_1]  # type: ignore

    if PL != -999.0:
        if PL > QS:
            DD = ESFAC / (PL - (1.0 - ESFAC) * QS)
            QS = QS * DD
            if DQ_trigger:
                DQ = DDQ * ERFAC * PL * DD * DD
        else:
            QS = MAX_MIXING_RATIO
            if DQ_trigger:
                DQ = 0.0
    else:
        if DQ_trigger:
            DQ = DDQ

    return QS, DQ


# Function version of QSat_table
@gtscript.function
def QSat_Float(
    ese: FloatField_Extra_Dim,  # type: ignore
    esx: FloatField_Extra_Dim,  # type: ignore
    T: Float,
    PL: Float,
    RAMP: Float = -999.0,
    PASCALS_trigger: bool = False,
    RAMP_trigger: bool = False,
    DQSAT_trigger: bool = False,
):
    if RAMP_trigger:
        URAMP = -abs(RAMP)
    else:
        URAMP = TMIX

    if PASCALS_trigger:
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
    IT_MINUS_1 = (
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
        if DQSAT_trigger:
            DQSAT = 0.0
        else:
            DQSAT = -999.0
    else:
        DD = 1.0 / (PP - (1.0 - ESFAC) * QSAT)
        QSAT = ESFAC * QSAT * DD
        if DQSAT_trigger:
            DQSAT = ESFAC * DQ * DEGSUBS * PP * (DD * DD)
        else:
            DQSAT = -999.0

    return QSAT, DQSAT


# Stencils implement QSAT0 function from GEOS_Utilities.F90
def QSat_FloatField(
    ese: FloatField_Extra_Dim,  # type: ignore
    esx: FloatField_Extra_Dim,  # type: ignore
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
        IT_MINUS_1 = (
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
