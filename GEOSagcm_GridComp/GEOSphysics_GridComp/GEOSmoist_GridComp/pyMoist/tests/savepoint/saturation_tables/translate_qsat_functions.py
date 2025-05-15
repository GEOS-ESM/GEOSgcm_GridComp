from ndsl import Namelist, StencilFactory
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.saturation_tables.qsat_functions import QSat_Float_Liquid, QSat_Float_Ice
from ndsl.dsl.gt4py import GlobalTable, PARALLEL, computation, interval
from ndsl.dsl.typing import Float, FloatField
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from pyMoist.saturation_tables.constants import TABLESIZE


GlobalTable_saturaion_tables = GlobalTable[(Float, (int(TABLESIZE)))]


def _stencil(
    t: FloatField,
    p: FloatField,
    out_qsatice: FloatField,
    out_qsatlqu: FloatField,
    out_dqsi: FloatField,
    out_dqsl: FloatField,
    ese: GlobalTable_saturaion_tables,
    esw: GlobalTable_saturaion_tables,
    esx: GlobalTable_saturaion_tables,
    frz: Float,
    lqu: Float,
):
    with computation(PARALLEL), interval(...):
        out_qsatice, out_dqsi = QSat_Float_Ice(ese, frz, t, p, True, True)
        out_qsatlqu, out_dqsl = QSat_Float_Liquid(esw, lqu, t, p, True, True)


def _stencil_2d(
    t: FloatField,
    p: FloatField,
    out_qsatice: FloatField,
    out_qsatlqu: FloatField,
    out_dqsi: FloatField,
    out_dqsl: FloatField,
    ese: GlobalTable_saturaion_tables,
    esw: GlobalTable_saturaion_tables,
    esx: GlobalTable_saturaion_tables,
    frz: Float,
    lqu: Float,
):
    with computation(PARALLEL), interval(0, 1):
        out_qsatice, out_dqsi = QSat_Float_Ice(ese, frz, t, p, True, True)
        out_qsatlqu, out_dqsl = QSat_Float_Liquid(esw, lqu, t, p, True, True)


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
        self.out_vars = {
            "SER_QSATLQU": {},
            "SER_QSATICE": {},
            "SER_DQSL": {},
            "SER_DQSI": {},
            "SER_QSATLQU_2D": {},
            "SER_QSATICE_2D": {},
            "SER_DQSL_2D": {},
            "SER_DQSI_2D": {},
        }

    def compute(self, inputs):
        # Initalize tables
        tables = SaturationVaporPressureTable(self.stencil_factory.backend)

        # Get input data
        t = inputs.pop("T")
        p = inputs.pop("PLmb")
        t_2d = inputs.pop("TEMP_ARRAY")
        p_2d = inputs.pop("PRES_ARRAY")

        out_qsatlqu = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        out_qsatice = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        out_dqsl = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        out_dqsi = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        domain_2d = (t_2d.shape[0], t_2d.shape[1])
        out_qsatlqu_2d = self.quantity_factory.zeros(
            [domain_2d[0], domain_2d[1]], "n/a"
        )
        out_qsatice_2d = self.quantity_factory.zeros(
            [domain_2d[0], domain_2d[1]], "n/a"
        )
        out_dqsl_2d = self.quantity_factory.zeros([domain_2d[0], domain_2d[1]], "n/a")
        out_dqsi_2d = self.quantity_factory.zeros([domain_2d[0], domain_2d[1]], "n/a")

        origin, domain = self.stencil_factory.grid_indexing.get_origin_domain(
            dims=[X_DIM, Y_DIM, Z_DIM]
        )

        stencil = self.stencil_factory.from_origin_domain(
            func=_stencil,
            origin=(0, 0, 0),
            domain=domain,
        )

        stencil_2d = self.stencil_factory.from_origin_domain(
            func=_stencil_2d,
            origin=(0, 0),
            domain=domain_2d,
        )

        # stencil(
        #     t,
        #     p,
        #     out_qsatice,
        #     out_qsatlqu,
        #     out_dqsi,
        #     out_dqsl,
        #     tables.ese,
        #     tables.esw,
        #     tables.esx,
        #     tables.frz,
        #     tables.lqu,
        # )

        # stencil_2d(
        #     t_2d,
        #     p_2d,
        #     out_qsatlqu_2d,
        #     out_qsatice_2d,
        #     out_dqsl_2d,
        #     out_dqsi_2d,
        #     tables.ese,
        #     tables.esw,
        #     tables.esx,
        #     tables.frz,
        #     tables.lqu,
        # )

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
            }
        )

        return inputs
