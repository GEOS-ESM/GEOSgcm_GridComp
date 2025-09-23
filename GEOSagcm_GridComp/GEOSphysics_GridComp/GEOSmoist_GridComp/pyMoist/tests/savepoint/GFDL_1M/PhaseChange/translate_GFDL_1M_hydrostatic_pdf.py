from ndsl import Namelist, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.constants import FLOAT_TINY
from pyMoist.GFDL_1M.PhaseChange.hydrostatic_pdf import hydrostatic_pdf
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable


class TranslateGFDL_1M_hydrostatic_pdf(TranslateFortranData2Py):
    def __init__(self, grid: Grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "convection_fraction": grid.compute_dict() | {"serialname": "CNV_FRC"},
            "surface_type": grid.compute_dict() | {"serialname": "SRF_TYPE"},
            "alpha": grid.compute_dict() | {"serialname": "ALPHA"},
            "p_mb": grid.compute_dict() | {"serialname": "PLmb"},
            "vapor": grid.compute_dict() | {"serialname": "Q"},
            "large_scale_liquid": grid.compute_dict() | {"serialname": "QLLS"},
            "convective_liquid": grid.compute_dict() | {"serialname": "QLCN"},
            "large_scale_ice": grid.compute_dict() | {"serialname": "QILS"},
            "convective_ice": grid.compute_dict() | {"serialname": "QICN"},
            "t": grid.compute_dict() | {"serialname": "T"},
            "large_scale_cloud_fraction": grid.compute_dict() | {"serialname": "CLLS"},
            "convective_cloud_fraction": grid.compute_dict() | {"serialname": "CLCN"},
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
            self.out_vars["convection_fraction"],
            self.out_vars["surface_type"],
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

        # NOTE errors in this test trace back to the saturation_specific_humidity call in the while loop,
        # which returns errors < 5 ULP for dqs. These exist because the GEOS_DQSAT and GEOS_QSAT
        # fortran functions were combined into a single saturation_specific_humidity NDSL function.
        # The two fortran functions were nearly identical, but had slightly different order of operations,
        # producing results which occasionally differ by a few ULP. These errors are considered acceptable
        # because the absolute and relative errors are both extremely small
        # (10-7-10-13 and 10-5, respectively)

        _hydrostatic_pdf(
            ese=tables.ese,
            esw=tables.esw,
            esx=tables.esx,
            estfrz=tables.frz,
            estlqu=tables.lqu,
            **inputs,
        )

        return inputs
