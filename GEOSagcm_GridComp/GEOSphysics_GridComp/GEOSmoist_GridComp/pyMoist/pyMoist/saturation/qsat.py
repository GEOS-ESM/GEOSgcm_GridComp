import os
from ndsl import StencilFactory, QuantityFactory, orchestrate
from ndsl.dsl.typing import Float, FloatField, Int, IntField
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from gt4py.cartesian.gtscript import computation, interval, PARALLEL, floor
import gt4py.cartesian.gtscript as gtscript
import copy
from typing import Optional
import xarray as xr
from pyMoist.saturation.formulation import SaturationFormulation
from pyMoist.saturation.table import get_table
from pyMoist.saturation.constants import (
    TMIX,
    TMINTBL,
    TMAXTBL,
    TMINLQU,
    MAPL_TICE,
    DEGSUBS,
    MAX_MIXING_RATIO,
    ESFAC,
    ERFAC,
    TABLESIZE,
)

# FloatField with extra dimension initialized to handle table data
# This is a temporary solution. This solution creates a KxTABLESIZE array
# At some point a proper solution will be implemented, enabling a 1xTABLESIZE array
# Allocation is currently fixed to TABLESIZE constant. Fortran has some examples of
# flexible table sixes (larger than TABLESIZE, with increased granulatiry).
# Current implementation does not allow for flexible table sizes
# so if needed this will have to be implemented in another way.
FloatField_Extra_Dim = gtscript.Field[gtscript.K, (Float, (int(TABLESIZE)))]


# The following two functions, QSat_Float_Liquid and QSat_Float_Ice are,
# together, a non-equivalent replacement for QSat_Float.
# QSat_Float interpolates a smooth transition between liquid water and ice
# (based on the parameter TMIX), while the Liquic/Ice versions have a
# hard transition at zero C.
@gtscript.function
def QSat_Float_Liquid(
    esw: FloatField_Extra_Dim,
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
        QS = esw[0][TABLESIZE - 1]
        if DQ_trigger:
            DDQ = 0.0
    else:
        TT = (TL - TMINTBL) * DEGSUBS + 1
        IT = int(TT)
        IT_PLUS_1 = (
            IT + 1
        )  # dace backend does not allow for [IT + 1] indexing because of cast to int
        DDQ = esw[0][IT_PLUS_1] - esw[0][IT]
        QS = (TT - IT) * DDQ + esw[0][IT]

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
    ese: FloatField_Extra_Dim,
    estfrz: Float,
    TL: Float,
    PL: Float = -999.0,
    DQ_trigger: bool = False,
):
    # qsatice.code with UTBL = True
    if TL <= TMINTBL:
        QS = ese[0][0]
        if DQ_trigger:
            DDQ = 0.0
    elif TL >= MAPL_TICE:
        QS = estfrz
        if DQ_trigger:
            DDQ = 0.0
    else:
        TT = (TL - TMINTBL) * DEGSUBS + 1
        IT = int(TT)
        IT_PLUS_1 = (
            IT + 1
        )  # dace backend does not allow for [IT + 1] indexing because of cast to int
        DDQ = ese[0][IT_PLUS_1] - ese[0][IT]
        QS = (TT - IT) * DDQ + ese[0][IT]

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
    ese: FloatField_Extra_Dim,
    esx: FloatField_Extra_Dim,
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
    IT = int(TI)
    IT_PLUS_1 = (
        IT + 1
    )  # dace backend does not allow for [IT + 1] indexing because of cast to int

    if URAMP == TMIX:
        DQ = esx[0][IT_PLUS_1] - esx[0][IT]  # should be esx[0][IT + 1] - esx[0][IT]
        QSAT = (TI - IT) * DQ + esx[0][IT]
    else:
        DQ = ese[0][IT_PLUS_1] - ese[0][IT]  # should be ese[0][IT + 1] - ese[0][IT]
        QSAT = (TI - IT) * DQ + ese[0][IT]

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
    ese: FloatField_Extra_Dim,
    esx: FloatField_Extra_Dim,
    T: FloatField,
    PL: FloatField,
    QSAT: FloatField,
    DQSAT: FloatField,
    RAMP: FloatField,
    PASCALS_trigger: bool,
    RAMP_trigger: bool,
    DQSAT_trigger: bool,
    IT_Float: FloatField,
    IFELSE: FloatField,
    QSAT_HALFWAY: FloatField,
    TI: FloatField,
):
    with computation(PARALLEL), interval(...):
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
        IT_Float = int(floor(TI))
        IT = int(floor(TI))
        IT_MINUS_1 = (
            IT - 1
        )  # dace backend does not allow for [IT + 1] indexing because of cast to int

        if URAMP == TMIX:
            DQ = esx[0][IT] - esx[0][IT_MINUS_1]
            QSAT = (TI - IT) * DQ + esx[0][IT_MINUS_1]
        else:
            DQ = ese[0][IT] - ese[0][IT_MINUS_1]
            QSAT = (TI - IT) * DQ + ese[0][IT_MINUS_1]
        QSAT_HALFWAY = QSAT
        if PP <= QSAT:
            IFELSE = 1
            QSAT = MAX_MIXING_RATIO
            if DQSAT_trigger:
                DQSAT = 0.0
        else:
            IFELSE = 0
            DD = 1.0 / (PP - (1.0 - ESFAC) * QSAT)
            QSAT = ESFAC * QSAT * DD
            if DQSAT_trigger:
                DQSAT = ESFAC * DQ * DEGSUBS * PP * (DD * DD)


class QSat:
    """Traditional satuation specific humidity.
    In Fortran: GEOS_Utilities:QSat

    Uses various formulations of the saturation vapor pressure to compute the saturation specific
    humidity for temperature T and pressure PL.

    The function get_table called in __init__ creates tables using exact calculations (equivalent
    to UTBL=False in Fortran). get_table can be called again for further exact computations.
    QSat_FloatField computes QSat for an entire field using table lookups
    (equivalent to UTBL=True in Fortran).

    For temperatures <= TMIX (-20C)
    the calculation is done over ice; for temperatures >= ZEROC (0C) the calculation
    is done over liquid water; and in between these values,
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
        table_method: Float = 0,
    ) -> None:
        self.extra_dim_quantity_factory = self.make_extra_dim_quantity_factory(
            quantity_factory
        )

        self.ese = self.extra_dim_quantity_factory.zeros([Z_DIM, "table_axis"], "n/a")
        self.esw = self.extra_dim_quantity_factory.zeros([Z_DIM, "table_axis"], "n/a")
        self.esx = self.extra_dim_quantity_factory.zeros([Z_DIM, "table_axis"], "n/a")

        if table_method == 0:
            with xr.open_dataset(
                os.path.join(os.path.dirname(__file__), "netCDFs", "QSat_Tables.nc")
            ) as ds:
                ese_array = ds.data_vars["ese"].values[0, 0, :]
                esw_array = ds.data_vars["esw"].values[0, 0, :]
                esx_array = ds.data_vars["esx"].values[0, 0, :]

                self.ese.view[:] = ese_array
                self.esw.view[:] = esw_array
                self.esx.view[:] = esx_array

        if table_method == 1:
            raise NotImplementedError(
                """Calculated tables are incorrect due to differences between C and Fortran calculations.\n
                For now the tables are being read in from a netCDF file.\n
                See https://github.com/GEOS-ESM/SMT-Nebulae/issues/88 for more information."""
            )
            self.table = get_table(formulation)

            self.ese.view[:] = self.table.ese
            self.esw.view[:] = self.table.esw
            self.esx.view[:] = self.table.esx

        self._RAMP = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._DQSAT = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        self.QSat = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        # For testing
        self._IT = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._IFELSE = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._QSAT_HALFWAY = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._TI = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        orchestrate(obj=self, config=stencil_factory.config.dace_config)
        self._QSat_FloatField = stencil_factory.from_dims_halo(
            func=QSat_FloatField,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
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
        RAMP: Optional[FloatField] = None,
        PASCALS: bool = False,
        DQSAT: bool = False,
        use_table_lookup: bool = True,
        ese: Optional[FloatField_Extra_Dim] = None,
        esw: Optional[FloatField_Extra_Dim] = None,
        esx: Optional[FloatField_Extra_Dim] = None,
    ):
        if RAMP is None:
            RAMP = self._RAMP
            RAMP_trigger = False
        else:
            RAMP_trigger = True

        if ese is not None:
            self.ese = ese
        if esw is not None:
            self.esw = esw
        if esx is not None:
            self.esx = esx

        if use_table_lookup:
            self._QSat_FloatField(
                self.ese,
                self.esx,
                T,
                PL,
                self.QSat,
                self._DQSAT,
                RAMP,
                PASCALS,
                RAMP_trigger,
                DQSAT,
                self._IT,
                self._IFELSE,
                self._QSAT_HALFWAY,
                self._TI,
            )

        if not use_table_lookup:
            raise NotImplementedError(
                "Saturation calculation: exact formulation not available, only table look up"
            )
