from typing import Any

from ndsl import Namelist, StencilFactory
from ndsl.boilerplate import get_factories_single_tile
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.gt4py import FORWARD, PARALLEL, computation, interval
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.field_types import GlobalTable_saturaion_tables
from pyMoist.saturation_tables.qsat_functions import (
    saturation_specific_humidity,
    saturation_specific_humidity_frozen_surface,
    saturation_specific_humidity_liquid_surface,
)
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable


def _stencil(
    t: FloatField,
    p: FloatField,
    out_qsatice: FloatField,
    out_qsatlqu: FloatField,
    out_dqsi: FloatField,
    out_dqsl: FloatField,
    out_qs: FloatField,
    out_dq: FloatField,
    ese: GlobalTable_saturaion_tables,
    esw: GlobalTable_saturaion_tables,
    esx: GlobalTable_saturaion_tables,
    frz: Float,
    lqu: Float,
):
    with computation(PARALLEL), interval(...):
        out_qsatice, out_dqsi = saturation_specific_humidity_frozen_surface(ese, frz, t, p * 100, True, True)
        out_qsatlqu, out_dqsl = saturation_specific_humidity_liquid_surface(esw, lqu, t, p * 100, True, True)

        out_qs, out_dq = saturation_specific_humidity(t=t, p=p * 100, ese=ese, esx=esx)


def _stencil_2d(
    t_2d: FloatFieldIJ,
    p_2d: FloatFieldIJ,
    out_qsatice_2d: FloatFieldIJ,
    out_qsatlqu_2d: FloatFieldIJ,
    out_dqsi_2d: FloatFieldIJ,
    out_dqsl_2d: FloatFieldIJ,
    out_qs_2d: FloatFieldIJ,
    out_dq_2d: FloatFieldIJ,
    ese: GlobalTable_saturaion_tables,
    esw: GlobalTable_saturaion_tables,
    esx: GlobalTable_saturaion_tables,
    frz: Float,
    lqu: Float,
):
    with computation(FORWARD), interval(0, 1):
        out_qsatice_2d, out_dqsi_2d = saturation_specific_humidity_frozen_surface(
            ese, frz, t_2d, p_2d * 100, True, True
        )
        out_qsatlqu_2d, out_dqsl_2d = saturation_specific_humidity_liquid_surface(
            esw, lqu, t_2d, p_2d * 100, True, True
        )

        out_qs_2d, out_dq_2d = saturation_specific_humidity(
            t=t_2d,
            p=p_2d * 100,
            ese=ese,
            esx=esx,
        )


class Translateqsat_functions(TranslateFortranData2Py):
    def __init__(self, grid: Grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory
        self._grid = grid
        self.max_error = 1e-9

        self.in_vars["data_vars"] = {
            "t": grid.compute_dict() | {"serialname": "T"},
            "p": grid.compute_dict() | {"serialname": "PLmb"},
            "t_2d": grid.compute_dict() | {"serialname": "TEMP_ARRAY"},
            "p_2d": grid.compute_dict() | {"serialname": "PRES_ARRAY"},
        }

        # Set Up Outputs
        self.out_vars: dict[str, Any] = {
            "SER_QSATLQU": {},
            "SER_QSATICE": {},
            "SER_DQSL": {},
            "SER_DQSI": {},
            "SER_QSATLQU_2D": {},
            "SER_QSATICE_2D": {},
            "SER_DQSL_2D": {},
            "SER_DQSI_2D": {},
            "SER_QSAT0_QS_2D": {},
            "SER_QSAT0_DQ_2D": {},
            "SER_DQSAT0_QS_2D": {},
            "SER_DQSAT0_DQ_2D": {},
            "SER_QSAT_QS": {},
            "SER_QSAT_DQ": {},
            "SER_DQSAT_QS": {},
            "SER_DQSAT_DQ": {},
        }

    def compute(self, inputs):
        # Initalize tables
        tables = SaturationVaporPressureTable(self.stencil_factory.backend)

        # Get input data
        t_2d = inputs.pop("TEMP_ARRAY")
        p_2d = inputs.pop("PRES_ARRAY")

        domain_2d = (t_2d.shape[0], t_2d.shape[1], 1)
        stencil_factory_2d, quantity_factory_2d = get_factories_single_tile(
            domain_2d[0],
            domain_2d[1],
            domain_2d[2],
            0,
            self.stencil_factory.backend,
        )

        t = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        t.field[:] = inputs.pop("T")
        p = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        p.field[:] = inputs.pop("PLmb")
        out_qsatlqu = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        out_qsatice = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        out_dqsl = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        out_dqsi = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        out_qsatlqu_2d = quantity_factory_2d.zeros([X_DIM, Y_DIM], "n/a")
        out_qsatice_2d = quantity_factory_2d.zeros([X_DIM, Y_DIM], "n/a")
        out_dqsl_2d = quantity_factory_2d.zeros([X_DIM, Y_DIM], "n/a")
        out_dqsi_2d = quantity_factory_2d.zeros([X_DIM, Y_DIM], "n/a")

        out_qsat_qs_2d = quantity_factory_2d.zeros([X_DIM, Y_DIM], "n/a")
        out_qsat_dq_2d = quantity_factory_2d.zeros([X_DIM, Y_DIM], "n/a")

        out_qs = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        out_dq = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        stencil = self.stencil_factory.from_dims_halo(
            func=_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        origin, domain = self.stencil_factory.grid_indexing.get_origin_domain(dims=[X_DIM, Y_DIM, Z_DIM])

        stencil_2d = stencil_factory_2d.from_origin_domain(
            func=_stencil_2d,
            origin=(0, 0, 0),
            domain=domain_2d,
        )

        stencil(
            t=t,
            p=p,
            out_qsatice=out_qsatice,
            out_qsatlqu=out_qsatlqu,
            out_dqsi=out_dqsi,
            out_dqsl=out_dqsl,
            out_qs=out_qs,
            out_dq=out_dq,
            ese=tables.ese,
            esw=tables.esw,
            esx=tables.esx,
            frz=tables.frz,
            lqu=tables.lqu,
        )

        stencil_2d(
            t_2d=t_2d,
            p_2d=p_2d,
            out_qsatice_2d=out_qsatice_2d,
            out_qsatlqu_2d=out_qsatlqu_2d,
            out_dqsi_2d=out_dqsi_2d,
            out_dqsl_2d=out_dqsl_2d,
            out_qs_2d=out_qsat_qs_2d,
            out_dq_2d=out_qsat_dq_2d,
            ese=tables.ese,
            esw=tables.esw,
            esx=tables.esx,
            frz=tables.frz,
            lqu=tables.lqu,
        )

        inputs.update(
            {
                "SER_QSATICE": out_qsatice.field,
                "SER_QSATLQU": out_qsatlqu.field,
                "SER_DQSI": out_dqsi.field,
                "SER_DQSL": out_dqsl.field,
                "SER_QSATICE_2D": out_qsatice_2d.field,
                "SER_QSATLQU_2D": out_qsatlqu_2d.field,
                "SER_DQSI_2D": out_dqsi_2d.field,
                "SER_DQSL_2D": out_dqsl_2d.field,
                "SER_QSAT0_QS_2D": out_qsat_qs_2d.field,
                "SER_QSAT0_DQ_2D": out_qsat_dq_2d.field,
                "SER_DQSAT0_QS_2D": out_qsat_qs_2d.field,
                "SER_DQSAT0_DQ_2D": out_qsat_dq_2d.field,
                "SER_QSAT_QS": out_qs.field,
                "SER_QSAT_DQ": out_dq.field,
                "SER_DQSAT_QS": out_qs.field,
                "SER_DQSAT_DQ": out_dq.field,
            }
        )

        return inputs
