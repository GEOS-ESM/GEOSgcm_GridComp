from f90nml import Namelist

from ndsl import StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.constants import FLOAT_TINY
from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.GFDL_1M.locals import GFDL1MLocals
from pyMoist.GFDL_1M.PhaseChange.hydrostatic_pdf import hydrostatic_pdf
from pyMoist.GFDL_1M.state import GFDL1MState
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable


class TranslateGFDL_1M_HydrostaticPDF(TranslateFortranData2Py):
    def __init__(self, grid: Grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "convection_fraction": {},
            "surface_type": {},
            "local_alpha": {},
            "local_p_mb": {},
            "mixing_ratio_vapor": {},
            "mixing_ratio_large_scale_liquid": {},
            "mixing_ratio_convective_liquid": {},
            "mixing_ratio_large_scale_ice": {},
            "mixing_ratio_convective_ice": {},
            "t": {},
            "cloud_fraction_large_scale": {},
            "cloud_fraction_convective": {},
            "concentration_liquid": {},
            "relative_humidity_after_pdf": {},
            "hydrostatic_pdf_iterations": {},
        }

        self.out_vars = self.in_vars["data_vars"].copy()
        del (
            self.out_vars["convection_fraction"],
            self.out_vars["surface_type"],
            self.out_vars["local_alpha"],
            self.out_vars["local_p_mb"],
            self.out_vars["concentration_liquid"],
        )

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GFDL_1M-constants")

    def compute(self, inputs):
        # initialize constants
        config = GFDL1MConfig(**self.constants)

        # initialize dataclasses
        state = GFDL1MState.zeros(self.quantity_factory)
        locals_ = GFDL1MLocals.make_as_state(self.quantity_factory)

        # Internal from wrapper class needed for this test
        alpha = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        # Initialize saturation tables
        saturation_tables = SaturationVaporPressureTable(self.stencil_factory.backend)

        # fill relevant parts of dataclasses
        state.convection_fraction.field[:] = inputs["convection_fraction"]
        state.surface_type.field[:] = inputs["surface_type"]
        alpha.field[:] = inputs["local_alpha"]
        locals_.p_mb.field[:] = inputs["local_p_mb"]
        state.mixing_ratio.vapor.field[:] = inputs["mixing_ratio_vapor"]
        state.mixing_ratio.large_scale_liquid.field[:] = inputs["mixing_ratio_large_scale_liquid"]
        state.mixing_ratio.convective_liquid.field[:] = inputs["mixing_ratio_convective_liquid"]
        state.mixing_ratio.large_scale_ice.field[:] = inputs["mixing_ratio_large_scale_ice"]
        state.mixing_ratio.convective_ice.field[:] = inputs["mixing_ratio_convective_ice"]
        state.t.field[:] = inputs["t"]
        state.cloud_fraction.large_scale.field[:] = inputs["cloud_fraction_large_scale"]
        state.cloud_fraction.convective.field[:] = inputs["cloud_fraction_convective"]
        state.concentration.liquid.field[:] = inputs["concentration_liquid"]
        state.relative_humidity_after_pdf.field[:] = inputs["relative_humidity_after_pdf"]
        state.hydrostatic_pdf_iterations.field[:] = inputs["hydrostatic_pdf_iterations"]

        # construct test stencil
        code = self.stencil_factory.from_dims_halo(
            func=hydrostatic_pdf,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "DT_MOIST": config.DT_MOIST,
                "PDF_SHAPE": config.PDFSHAPE,
                "USE_BERGERON": config.USE_BERGERON,
                "FLOAT_TINY": FLOAT_TINY,
            },
        )
        code(
            alpha=alpha,
            convection_fraction=state.convection_fraction,
            surface_type=state.surface_type,
            p_mb=locals_.p_mb,
            vapor=state.mixing_ratio.vapor,
            large_scale_liquid=state.mixing_ratio.large_scale_liquid,
            convective_liquid=state.mixing_ratio.convective_liquid,
            large_scale_ice=state.mixing_ratio.large_scale_ice,
            convective_ice=state.mixing_ratio.convective_ice,
            t=state.t,
            large_scale_cloud_fraction=state.cloud_fraction.large_scale,
            convective_cloud_fraction=state.cloud_fraction.convective,
            ice_concentration=state.concentration.ice,
            relative_humidity=state.relative_humidity_after_pdf,
            pdf_iters=state.hydrostatic_pdf_iterations,
            ese=saturation_tables.ese,
            esw=saturation_tables.esw,
            esx=saturation_tables.esx,
            estfrz=saturation_tables.frz,
            estlqu=saturation_tables.lqu,
        )

        return {
            "mixing_ratio_vapor": state.mixing_ratio.vapor.field[:],
            "mixing_ratio_large_scale_liquid": state.mixing_ratio.large_scale_liquid.field[:],
            "mixing_ratio_convective_liquid": state.mixing_ratio.convective_liquid.field[:],
            "mixing_ratio_large_scale_ice": state.mixing_ratio.large_scale_ice.field[:],
            "mixing_ratio_convective_ice": state.mixing_ratio.convective_ice.field[:],
            "t": state.t.field[:],
            "cloud_fraction_large_scale": state.cloud_fraction.large_scale.field[:],
            "cloud_fraction_convective": state.cloud_fraction.convective.field[:],
            "relative_humidity_after_pdf": state.relative_humidity_after_pdf.field[:],
            "hydrostatic_pdf_iterations": state.hydrostatic_pdf_iterations.field[:],
        }
