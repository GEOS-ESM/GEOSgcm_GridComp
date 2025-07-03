from ndsl import Namelist, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.constants import FLOAT_TINY
from pyMoist.GFDL_1M.PhaseChange.hydrostatic_pdf import hydrostatic_pdf
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable


class Translatehydrostatic_pdf(TranslateFortranData2Py):
    def __init__(self, grid: Grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "cnv_frc": grid.compute_dict() | {"serialname": "CNV_FRC"},
            "srf_type": grid.compute_dict() | {"serialname": "SRF_TYPE"},
            "alpha": grid.compute_dict() | {"serialname": "ALPHA"},
            "p_mb": grid.compute_dict() | {"serialname": "PLmb"},
            "q": grid.compute_dict() | {"serialname": "Q"},
            "qlls": grid.compute_dict() | {"serialname": "QLLS"},
            "qlcn": grid.compute_dict() | {"serialname": "QLCN"},
            "qils": grid.compute_dict() | {"serialname": "QILS"},
            "qicn": grid.compute_dict() | {"serialname": "QICN"},
            "t": grid.compute_dict() | {"serialname": "T"},
            "clls": grid.compute_dict() | {"serialname": "CLLS"},
            "clcn": grid.compute_dict() | {"serialname": "CLCN"},
            "nacti": grid.compute_dict() | {"serialname": "NACTI"},
            "rhx": grid.compute_dict() | {"serialname": "RHX"},
        }

        self.in_vars["parameters"] = [
            "DT_MOIST",
            "PDFSHAPE",
            "USE_BERGERON",
        ]

        self.out_vars = self.in_vars["data_vars"].copy()
        del (
            self.out_vars["cnv_frc"],
            self.out_vars["srf_type"],
            self.out_vars["alpha"],
            self.out_vars["p_mb"],
            self.out_vars["nacti"],
        )

    def compute_from_storage(self, inputs):
        _hydrostatic_pdf = self.stencil_factory.from_dims_halo(
            func=hydrostatic_pdf,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "DT_MOIST": inputs.pop("DT_MOIST"),
                "PDF_SHAPE": inputs.pop("PDFSHAPE"),
                "USE_BERGERON": bool(inputs.pop("USE_BERGERON")),
                "FLOAT_TINY": FLOAT_TINY,
            },
        )

        tables = SaturationVaporPressureTable(self.stencil_factory.backend)

        _hydrostatic_pdf(
            ese=tables.ese,
            esw=tables.esw,
            esx=tables.esx,
            estfrz=tables.frz,
            estlqu=tables.lqu,
            **inputs,
        )

        return inputs
