from ndsl import StencilFactory, QuantityFactory, orchestrate
from ndsl.dsl.typing import Float, FloatField, Int, FloatFieldI, IntField
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.boilerplate import get_factories_single_tile_numpy
from gt4py.cartesian.gtscript import computation, interval, PARALLEL, floor
import gt4py.cartesian.gtscript as gtscript
from pyMoist.saturation.formulation import SaturationFormulation
from pyMoist.saturation.table import get_table
from pyMoist.saturation.constants import (
    TMIX,
    TMINTBL,
    TMAXTBL,
    DEGSUBS,
    MAX_MIXING_RATIO,
    ESFAC,
)

# Stencils implement QSAT0 function from GEOS_Utilities.F90

def _QSat_table_part1(
    TL: FloatField,
    PL: FloatField,
    TI: FloatField,
    IT: FloatField,
    PASCALS: bool,
    RAMP: Float,
):
    with computation(PARALLEL), interval(...):
        if RAMP != -999:
            URAMP = -abs(RAMP)
        else:
            URAMP = TMIX
        
        if PASCALS == True:
            PP = PL
        else:
            PP = PL*100.

        
        if TL <= TMINTBL:
            TI = TMINTBL
        elif TL >= TMAXTBL-.001:
            TI = TMAXTBL-.001
        else:
            TI = TL
        
        TI = (TI - TMINTBL)*DEGSUBS+1
        IT = floor(TI)


def _QSat_table_part2(
    TL: FloatField,
    PL: FloatField,
    TI: FloatField,
    IT: IntField,
    ITP1: IntField,
    BLE: FloatFieldI,
    BLW: FloatFieldI,
    BLX: FloatFieldI,
    QSAT: FloatField,
    PASCALS: bool,
    RAMP: Float,
):
    with computation(PARALLEL), interval(...):
        if RAMP != -999:
            URAMP = -abs(RAMP)
        else:
            URAMP = TMIX
        
        if PASCALS == True:
            PP = PL
        else:
            PP = PL*100.
        
        if URAMP==TMIX:
            DQ = BLX.A[ITP1] - BLX.A[IT]
            QSAT = (TI-IT)*DQ + BLX.A[IT]
        else:
            DQ    = BLE.A[ITP1] - BLE.A[IT]
            QSAT  = (TI-IT)*DQ + BLE.A[IT]

        # if DQSAT != -999:
        #     DQSAT = DQ*DEGSUBS

        if PP <= QSAT:
            QSAT = MAX_MIXING_RATIO
            # if DQSAT != -999:
            #     DQSAT = 0.0
        else:
            DD = 1.0/(PP - (1.0-ESFAC)*QSAT)
            QSAT = ESFAC*QSAT*DD
            # if DQSAT != -999:
            #     DQSAT = ESFAC*DQSAT*PP*(DD*DD)


class QSat:
    """Traditional satuation specific humidity.
    In Fortran: GEOS_Utilities:QSat

    Uses various formulations of the saturation vapor pressure to compute the saturation specific
    humidity for temperature TL and pressure PL.

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
        
        orchestrate(obj=self, config=stencil_factory.config.dace_config)
        self._QSat_table_part1 = stencil_factory.from_dims_halo(
            func=_QSat_table_part1,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )
        self._QSat_table_part2 = stencil_factory.from_dims_halo(
            func=_QSat_table_part2,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )
        
        self._TI = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._IT = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._IT_plus1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        self._table = get_table(formulation)

    def __call__(
            self,
            TL: FloatField,
            PL: FloatField,
            QSAT: FloatField,
            use_table_lookup: bool = True,
            PASCALS: bool = False,
            RAMP: Float = -999,
    ) -> Float:
        
        #Enable only for testing
        use_table_lookup: bool = True
        self._PASCALS = False
        self._RAMP = -999
        self._DQSAT = -999

        if use_table_lookup:
            self._QSat_table_part1(TL, PL, self._TI, self._IT, self._PASCALS, self._RAMP)
            #
            for i in range(0, self._IT.view[:].shape[0]):
                for j in range(0, self._IT.view[:].shape[1]):
                    for k in range(0, self._IT.view[:].shape[2]):
                        self._IT[i,j,k] = int(self._IT[i,j,k])
                        self._IT_plus1[i,j,k] = int(self._IT[i,j,k])
                        print(self._IT[i,j,k])
            self._QSat_table_part2(TL, PL, self._TI, self._IT, self._IT_plus1, self._table.ble, self._table.blw,
                                   self._table.blx, QSAT, self._PASCALS, self._RAMP)

        if not use_table_lookup:
            raise NotImplementedError(
                "Saturation calculation: exact formulation not available, only table look up"
            )
