"""This function identifies the index of the lifted
condensation level (LCL) within a column. LCL is 
the point at which a lifted parcel becomes saturated."""

from gt4py.cartesian.gtscript import (
    FORWARD,
    PARALLEL,
    computation,
    interval,
    log,
)
import pyMoist.pyMoist_constants as constants
from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField, FloatFieldIJ
from pyMoist.field_types import FloatField_Extra_Dim
from pyMoist.saturation.formulation import SaturationFormulation
from pyMoist.saturation.qsat import QSat, QSat_Float


def _get_last(
    in_field: FloatField, temporary_field: FloatFieldIJ, out_field: FloatField
):
    with computation(FORWARD), interval(-1, None):
        temporary_field = in_field

    with computation(PARALLEL), interval(...):
        out_field = temporary_field


def _find_klcl(
    T: FloatField,
    P: FloatField,
    Q: FloatField,
    T_top: FloatFieldIJ,
    P_top: FloatFieldIJ,
    Q_top: FloatFieldIJ,
    ese: FloatField_Extra_Dim,
    esw: FloatField_Extra_Dim,
    esx: FloatField_Extra_Dim,
    PLCL: FloatFieldIJ,
):
    with computation(FORWARD), interval(-1, None):
        T_top = T
        P_top = P
        Q_top = Q

    with computation(FORWARD), interval(-1, None):
        QSat_function, _ = QSat_Float(ese, esw, esx, T_top, P_top)
        RHSFC = 100.0 * Q_top / QSat_function
        TLCL = (
            1 / ((1.0 / (T_top - 55.0)) - (log(max(0.1, RHSFC) / 100.0) / 2840.0))
        ) + 55.0
        Rm = (1.0 - Q_top) * constants.rdry + Q_top * constants.rvap
        Cpm = (1.0 - Q_top) * constants.cpdry + Q_top * constants.cpvap
        PLCL = P_top * ((TLCL / T_top) ** (Cpm / Rm))


class FindKLCL:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        formulation: SaturationFormulation = SaturationFormulation.Staars,
        use_table_lookup: bool = True,
    ) -> None:

        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory
        self._T_top = self.quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        self._P_top = self.quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        self._Q_top = self.quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        self.KLCL = self.quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        self.PLCL = self.quantity_factory.zeros([X_DIM, Y_DIM], "n/a")

        self.get_last = self.stencil_factory.from_dims_halo(
            func=_get_last,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )
        self.find_klcl = self.stencil_factory.from_dims_halo(
            func=_find_klcl,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        T: FloatField,
        P: FloatField,
        Q: FloatField,
        formulation: SaturationFormulation = SaturationFormulation.Staars,
        use_table_lookup: bool = True,
    ):

        self.qsat = QSat(
            self.stencil_factory,
            self.quantity_factory,
            formulation=formulation,
            use_table_lookup=use_table_lookup,
        )

        self.find_klcl(
            T,
            P,
            Q,
            self._T_top,
            self._P_top,
            self._Q_top,
            self.qsat.ese,
            self.qsat.esw,
            self.qsat.esx,
            self.PLCL,
        )

        # TODO: Workaround for hybrid indexing (mask)
        for i in range(0, P.view[:].shape[0]):
            for j in range(0, P.view[:].shape[1]):
                for k in reversed(range(0, P.view[:].shape[2])):
                    if P.view[i, j, k] > self.PLCL.view[i, j]:
                        self.KLCL.view[:][i, j] = k
