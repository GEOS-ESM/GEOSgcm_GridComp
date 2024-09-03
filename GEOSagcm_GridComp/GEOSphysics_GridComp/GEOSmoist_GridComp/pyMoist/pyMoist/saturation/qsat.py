from ndsl import StencilFactory, QuantityFactory, orchestrate
from ndsl.dsl.typing import Float, FloatField, Int, IntField
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from gt4py.cartesian.gtscript import computation, interval, PARALLEL
import gt4py.cartesian.gtscript as gtscript
import copy
from typing import Optional
from pyMoist.saturation.formulation import SaturationFormulation
from pyMoist.saturation.table import get_table
from pyMoist.saturation.constants import (
    TMIX,
    TMINTBL,
    TMAXTBL,
    DEGSUBS,
    MAX_MIXING_RATIO,
    ESFAC,
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

# Function version of QSat_table
@gtscript.function
def QSat_Float(
    ese: FloatField_Extra_Dim,
    esw: FloatField_Extra_Dim,
    esx: FloatField_Extra_Dim,
    T: Float,
    PL: Float,
    RAMP: Float = -999.,
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
        PP = PL*100.

    
    if T <= TMINTBL:
        TI = TMINTBL
    elif T >= TMAXTBL-.001:
        TI = TMAXTBL-.001
    else:
        TI = T
    
    TI = (TI - TMINTBL)*DEGSUBS+1
    IT = int(TI)

    if URAMP==TMIX:
        DQ = esx[0][IT] - esx[0][IT]
        QSAT = (TI-IT)*DQ + esx[0][IT]
    else:
        DQ    = ese[0][IT] - ese[0][IT]
        QSAT  = (TI-IT)*DQ + ese[0][IT]

    if PP <= QSAT:
        QSAT = MAX_MIXING_RATIO
        if DQSAT_trigger: DQSAT = 0.0
        else: DQSAT = -999.
    else:
        DD = 1.0/(PP - (1.0-ESFAC)*QSAT)
        QSAT = ESFAC*QSAT*DD
        if DQSAT_trigger: DQSAT = ESFAC*DQ*DEGSUBS*PP*(DD*DD)
        else: DQSAT = -999.

    return QSAT, DQSAT


# Stencils implement QSAT0 function from GEOS_Utilities.F90
def QSat_FloatField(
    ese: FloatField_Extra_Dim,
    esw: FloatField_Extra_Dim,
    esx: FloatField_Extra_Dim,
    T: FloatField,
    PL: FloatField,
    QSAT: FloatField,
    DQSAT: FloatField,
    RAMP: FloatField,
    PASCALS_trigger: bool,
    RAMP_trigger: bool,
    DQSAT_trigger: bool,
):
    with computation(PARALLEL), interval(...):
        if RAMP_trigger:
            URAMP = -abs(RAMP)
        else:
            URAMP = TMIX
        
        if PASCALS_trigger:
            PP = PL
        else:
            PP = PL*100.

        
        if T <= TMINTBL:
            TI = TMINTBL
        elif T >= TMAXTBL-.001:
            TI = TMAXTBL-.001
        else:
            TI = T
        
        TI = (TI - TMINTBL)*DEGSUBS+1
        IT = int(TI)

        if URAMP==TMIX:
            DQ = esx[0][IT] - esx[0][IT]
            QSAT = (TI-IT)*DQ + esx[0][IT]
        else:
            DQ    = ese[0][IT] - ese[0][IT]
            QSAT  = (TI-IT)*DQ + ese[0][IT]

        if PP <= QSAT:
            QSAT = MAX_MIXING_RATIO
            if DQSAT_trigger: DQSAT = 0.0
        else:
            DD = 1.0/(PP - (1.0-ESFAC)*QSAT)
            QSAT = ESFAC*QSAT*DD
            if DQSAT_trigger: DQSAT = ESFAC*DQ*DEGSUBS*PP*(DD*DD)


class QSat:
    """Traditional satuation specific humidity.
    In Fortran: GEOS_Utilities:QSat

    Uses various formulations of the saturation vapor pressure to compute the saturation specific
    humidity for temperature T and pressure PL.

    For temperatures <= TMIX (-20C)
    the calculation is done over ice; for temperatures >= ZEROC (0C) the calculation
    is done over liquid water; and in between these values,
    it interpolates linearly between the two.

    The optional RAMP is the width of this
    ice/water ramp (i.e., TMIX = ZEROC-RAMP); its default is 20.

    If PASCALS is true, PL is
    assumed to be in Pa; if false or not present, it is assumed to be in mb.

    Another compile time choice is whether to use the exact formulation
    or a table look-up.
    If UTBL is true, tabled values of the saturation vapor pressures
    are used. These tables are automatically generated at a 0.1K resolution
    for whatever vapor pressure formulation is being used.
    """

    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        formulation: SaturationFormulation = SaturationFormulation.Staars,
        use_table_lookup: bool = True,
    ) -> None:

        self.table = get_table(formulation)

        self.extra_dim_quantity_factory = self.make_extra_dim_quantity_factory(
            quantity_factory
        )

        self._ese = self.extra_dim_quantity_factory.zeros([Z_DIM, "table_axis"], "n/a")
        self._esw = self.extra_dim_quantity_factory.zeros([Z_DIM, "table_axis"], "n/a")
        self._esx = self.extra_dim_quantity_factory.zeros([Z_DIM, "table_axis"], "n/a")
        self._ese.view[:] = self.table.ese
        self._esw.view[:] = self.table.esw
        self._esx.view[:] = self.table.esx

        self._RAMP = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._DQSAT = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        self.QSat = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

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
    ):

        if RAMP is None:
            RAMP = self._RAMP
            RAMP_trigger = False
        else:
            RAMP_trigger = True

        if use_table_lookup:
            self._QSat_FloatField(self._ese, self._esw, self._esx, T, PL, self.QSat, self._DQSAT,
                             RAMP, PASCALS, RAMP_trigger, DQSAT)

        if not use_table_lookup:
            raise NotImplementedError(
                "Saturation calculation: exact formulation not available, only table look up"
            )