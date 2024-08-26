from ndsl import StencilFactory, QuantityFactory, orchestrate
from ndsl.dsl.typing import Float, FloatField, Int, IntField
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from gt4py.cartesian.gtscript import computation, interval, PARALLEL
import gt4py.cartesian.gtscript as gtscript
import copy
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

# Stencils implement QSAT0 function from GEOS_Utilities.F90

def QSat_table(
    BLE: FloatField_Extra_Dim,
    BLW: FloatField_Extra_Dim,
    BLX: FloatField_Extra_Dim,
    T: FloatField,
    PL: FloatField,
    QSAT: FloatField,
    RAMP: FloatField,
    DQSAT: FloatField,
    IT: FloatField,
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
        IT_USE = int(TI)

        if URAMP==TMIX:
            DQ = BLX[0][IT_USE] - BLX[0][IT_USE]
            QSAT = (TI-IT_USE)*DQ + BLX[0][IT_USE]
        else:
            DQ    = BLE[0][IT_USE] - BLE[0][IT_USE]
            QSAT  = (TI-IT_USE)*DQ + BLE[0][IT_USE]

        if DQSAT_trigger == True:
            DQSAT = DQ*DEGSUBS

        if PP <= QSAT:
            QSAT = MAX_MIXING_RATIO
            if DQSAT_trigger: DQSAT = 0.0
        else:
            DD = 1.0/(PP - (1.0-ESFAC)*QSAT)
            QSAT = ESFAC*QSAT*DD
            if DQSAT_trigger: DQSAT = ESFAC*DQSAT*PP*(DD*DD)


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

        self._ble = self.extra_dim_quantity_factory.zeros([Z_DIM, "table_axis"], "n/a")
        self._blw = self.extra_dim_quantity_factory.zeros([Z_DIM, "table_axis"], "n/a")
        self._blx = self.extra_dim_quantity_factory.zeros([Z_DIM, "table_axis"], "n/a")
        self._ble.view[:] = self.table.ble
        self._blw.view[:] = self.table.blw
        self._blx.view[:] = self.table.blx

        self.RAMP = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.DQSAT = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        self._IT = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._TI = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        orchestrate(obj=self, config=stencil_factory.config.dace_config)
        self._QSat_table = stencil_factory.from_dims_halo(
            func=QSat_table,
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
            QSAT: FloatField,
            IT: FloatField,
            ESTBLE_TEST: FloatField_Extra_Dim,
            ESTBLW_TEST: FloatField_Extra_Dim,
            ESTBLX_TEST: FloatField_Extra_Dim,
            TEMP_IN_QSAT: FloatField_Extra_Dim,
            use_table_lookup: bool = True,
            PASCALS_trigger: bool = False,
            RAMP_trigger: bool = False,
            DQSAT_trigger: bool = False,
    ):

        ESTBLE_TEST.view[0] = self.table.ble
        ESTBLW_TEST.view[0] = self.table.blw
        ESTBLX_TEST.view[0] = self.table.blx
        TEMP_IN_QSAT.view[0] = self.table._TI[:]
        print(self.table._TI)

        if use_table_lookup:
            self._QSat_table(self._ble, self._blw, self._blx, T, PL, QSAT, self.RAMP, self.DQSAT,
                             IT, PASCALS_trigger, RAMP_trigger, DQSAT_trigger)

        if not use_table_lookup:
            raise NotImplementedError(
                "Saturation calculation: exact formulation not available, only table look up"
            )