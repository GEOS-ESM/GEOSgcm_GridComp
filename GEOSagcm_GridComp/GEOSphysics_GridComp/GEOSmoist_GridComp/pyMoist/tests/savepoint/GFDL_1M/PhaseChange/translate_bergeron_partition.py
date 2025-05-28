from ndsl import Namelist, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.gt4py import PARALLEL, computation, interval
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.field_types import GlobalTable_saturaion_tables
from pyMoist.GFDL_1M.PhaseChange.hydrostatic_pdf import bergeron_partition
from pyMoist.saturation_tables.formulation import SaturationFormulation
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable


def _stencil(
    DT_MOIST: Float,
    PL: FloatField,
    TEn: FloatField,
    QT: FloatField,
    QILS: FloatField,
    QICN: FloatField,
    QLLS: FloatField,
    QLCN: FloatField,
    NIv: FloatField,
    DQCALL: FloatField,
    CNVFRC: FloatFieldIJ,
    SRF_TYPE: FloatFieldIJ,
    ese: GlobalTable_saturaion_tables,
    esw: GlobalTable_saturaion_tables,
    estfrz: Float,
    estlqu: Float,
    fQi: FloatField,
):
    with computation(PARALLEL), interval(...):
        fQi, DQCALL = bergeron_partition(
            DT_MOIST,
            PL,
            TEn,
            QT,
            QILS,
            QICN,
            QLLS,
            QLCN,
            NIv,
            DQCALL,
            CNVFRC,
            SRF_TYPE,
            ese,
            esw,
            estfrz,
            estlqu,
        )


class Translatebergeron_partition(TranslateFortranData2Py):
    def __init__(self, grid: Grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "PL": grid.compute_dict(),
            "TEn": grid.compute_dict(),
            "QT": grid.compute_dict(),
            "QILS": grid.compute_dict(),
            "QICN": grid.compute_dict(),
            "QLLS": grid.compute_dict(),
            "QLCN": grid.compute_dict(),
            "NIv": grid.compute_dict(),
            "DQCALL": grid.compute_dict(),
            "CNVFRC": grid.compute_dict(),
            "SRF_TYPE": grid.compute_dict(),
            "fQi": grid.compute_dict(),
        }

        self.in_vars["parameters"] = [
            "DT",
        ]

        self.out_vars = self.in_vars["data_vars"].copy()
        # self.out_vars = {
        #     "DQCALL": {} | {"serialname": "DQCALL"},
        #     "fQi": {},
        # }

    def compute_from_storage(self, inputs):
        stencil = self.stencil_factory.from_dims_halo(
            func=_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        tables = SaturationVaporPressureTable(self.stencil_factory.backend)

        f_qi = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dq_all_out = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        stencil(
            DT_MOIST=inputs.pop("DT"),
            ese=tables.ese,
            esw=tables.esw,
            estfrz=tables.frz,
            estlqu=tables.lqu,
            **inputs
        )

        # inputs.update(
        #     {
        #         "DQCALL": dq_all_out.field,
        #         "fQi": f_qi.field,
        #     }
        # )

        return inputs
