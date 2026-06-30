import copy
import time

from f90nml import Namelist
from ndsl import StencilFactory, ndsl_log
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py

from pyMoist.microphysics.GFDL_1M.config import GFDL1MConfig
from pyMoist.microphysics.GFDL_1M.locals import GFDL1MLocals
from pyMoist.microphysics.GFDL_1M.PhaseChange.phase_change import PhaseChange
from pyMoist.microphysics.GFDL_1M.state import GFDL1MState
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable


class TranslateGFDL_1M_PhaseChange(TranslateFortranData2Py):
    def __init__(self, grid: Grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.quantity_factory = grid.quantity_factory

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "estimated_inversion_strength": {},
            "convection_fraction": {},
            "surface_type": {},
            "t": {},
            "mixing_ratio_vapor": {},
            "concentration_liquid": {},
            "concentration_ice": {},
            "mixing_ratio_convective_liquid": {},
            "mixing_ratio_convective_ice": {},
            "mixing_ratio_large_scale_ice": {},
            "mixing_ratio_large_scale_liquid": {},
            "cloud_fraction_large_scale": {},
            "cloud_fraction_convective": {},
            "local_lcl_level": {},
            "local_p_mb": {},
            "local_p_interface_mb": {},
            "area": {},
            "hydrostatic_pdf_iterations": {},
            "local_saturation_specific_humidity": {},
            "cloud_liquid_evaporation": {},
            "cloud_ice_sublimation": {},
            "relative_humidity_after_pdf": {},
            "critical_relative_humidity_for_pdf": {},
        }

        self.out_vars = self.in_vars["data_vars"].copy()

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GFDL_1M-constants")

    def compute(self, inputs):
        # initialize constants
        config = GFDL1MConfig(**self.constants)

        # initialize dataclasses
        state = GFDL1MState.zeros(self.quantity_factory)
        locals_ = GFDL1MLocals.make_as_state(self.quantity_factory)

        # Initialize saturation tables
        saturation_tables = SaturationVaporPressureTable(self.stencil_factory.backend)

        state.estimated_inversion_strength.field[:] = inputs["estimated_inversion_strength"]
        state.convection_fraction.field[:] = inputs["convection_fraction"]
        state.surface_type.field[:] = inputs["surface_type"]
        state.t.field[:] = inputs["t"]
        state.mixing_ratio.vapor.field[:] = inputs["mixing_ratio_vapor"]
        state.concentration.liquid.field[:] = inputs["concentration_liquid"]
        state.concentration.ice.field[:] = inputs["concentration_ice"]
        state.mixing_ratio.convective_liquid.field[:] = inputs["mixing_ratio_convective_liquid"]
        state.mixing_ratio.convective_ice.field[:] = inputs["mixing_ratio_convective_ice"]
        state.mixing_ratio.large_scale_ice.field[:] = inputs["mixing_ratio_large_scale_ice"]
        state.mixing_ratio.large_scale_liquid.field[:] = inputs["mixing_ratio_large_scale_liquid"]
        state.cloud_fraction.large_scale.field[:] = inputs["cloud_fraction_large_scale"]
        state.cloud_fraction.convective.field[:] = inputs["cloud_fraction_convective"]
        locals_.lcl_level.field[:] = inputs["local_lcl_level"] - 1
        locals_.p_mb.field[:] = inputs["local_p_mb"]
        locals_.p_interface_mb.field[:] = inputs["local_p_interface_mb"]
        state.area.field[:] = inputs["area"]
        state.hydrostatic_pdf_iterations.field[:] = inputs["hydrostatic_pdf_iterations"]
        locals_.saturation_specific_humidity.field[:] = inputs["local_saturation_specific_humidity"]
        state.cloud_liquid_evaporation.field[:] = inputs["cloud_liquid_evaporation"]
        state.cloud_ice_sublimation.field[:] = inputs["cloud_ice_sublimation"]
        state.relative_humidity_after_pdf.field[:] = inputs["relative_humidity_after_pdf"]
        state.critical_relative_humidity_for_pdf.field[:] = inputs["critical_relative_humidity_for_pdf"]

        # construct test stencil
        code = PhaseChange(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            saturation_tables=saturation_tables,
        )

        # run test code
        code(
            t=state.t,
            mixing_ratio_vapor=state.mixing_ratio.vapor,
            mixing_ratio_large_scale_liquid=state.mixing_ratio.large_scale_liquid,
            mixing_ratio_convective_liquid=state.mixing_ratio.convective_liquid,
            mixing_ratio_large_scale_ice=state.mixing_ratio.large_scale_ice,
            mixing_ratio_convective_ice=state.mixing_ratio.convective_ice,
            cloud_fraction_large_scale=state.cloud_fraction.large_scale,
            cloud_fraction_convective=state.cloud_fraction.convective,
            concentration_ice=state.concentration.ice,
            concentration_liquid=state.concentration.liquid,
            relative_humidity_after_pdf=state.relative_humidity_after_pdf,
            estimated_inversion_strength=state.estimated_inversion_strength,
            area=state.area,
            critical_relative_humidity_for_pdf=state.critical_relative_humidity_for_pdf,
            pdf_iters=state.hydrostatic_pdf_iterations,
            cloud_liquid_evaporation=state.cloud_liquid_evaporation,
            cloud_ice_sublimation=state.cloud_ice_sublimation,
            convection_fraction=state.convection_fraction,
            surface_type=state.surface_type,
            local_lcl_level=locals_.lcl_level,
            local_p_mb=locals_.p_mb,
            local_p_interface_mb=locals_.p_interface_mb,
            local_saturation_specific_humidity=locals_.saturation_specific_humidity,
        )

        outputs = {
            "estimated_inversion_strength": copy.copy(state.estimated_inversion_strength.field[:]),
            "convection_fraction": copy.copy(state.convection_fraction.field[:]),
            "surface_type": copy.copy(state.surface_type.field[:]),
            "t": copy.copy(state.t.field[:]),
            "mixing_ratio_vapor": copy.copy(state.mixing_ratio.vapor.field[:]),
            "concentration_liquid": copy.copy(state.concentration.liquid.field[:]),
            "concentration_ice": copy.copy(state.concentration.ice.field[:]),
            "mixing_ratio_convective_liquid": copy.copy(state.mixing_ratio.convective_liquid.field[:]),
            "mixing_ratio_convective_ice": copy.copy(state.mixing_ratio.convective_ice.field[:]),
            "mixing_ratio_large_scale_ice": copy.copy(state.mixing_ratio.large_scale_ice.field[:]),
            "mixing_ratio_large_scale_liquid": copy.copy(state.mixing_ratio.large_scale_liquid.field[:]),
            "cloud_fraction_large_scale": copy.copy(state.cloud_fraction.large_scale.field[:]),
            "cloud_fraction_convective": copy.copy(state.cloud_fraction.convective.field[:]),
            "local_lcl_level": copy.copy(locals_.lcl_level.field[:]) + 1,
            "local_p_mb": copy.copy(locals_.p_mb.field[:]),
            "local_p_interface_mb": copy.copy(locals_.p_interface_mb.field[:]),
            "area": copy.copy(state.area.field[:]),
            "hydrostatic_pdf_iterations": copy.copy(state.hydrostatic_pdf_iterations.field[:]),
            "local_saturation_specific_humidity": copy.copy(locals_.saturation_specific_humidity.field[:]),
            "cloud_liquid_evaporation": copy.copy(state.cloud_liquid_evaporation.field[:]),
            "cloud_ice_sublimation": copy.copy(state.cloud_ice_sublimation.field[:]),
            "relative_humidity_after_pdf": copy.copy(state.relative_humidity_after_pdf.field[:]),
            "critical_relative_humidity_for_pdf": copy.copy(state.critical_relative_humidity_for_pdf.field[:]),
        }

        # Micro-bench
        ts = 0
        if ts > 0:
            s = time.perf_counter()
            for _ in range(ts):
                code(
                    t=state.t,
                    mixing_ratio_vapor=state.mixing_ratio.vapor,
                    mixing_ratio_large_scale_liquid=state.mixing_ratio.large_scale_liquid,
                    mixing_ratio_convective_liquid=state.mixing_ratio.convective_liquid,
                    mixing_ratio_large_scale_ice=state.mixing_ratio.large_scale_ice,
                    mixing_ratio_convective_ice=state.mixing_ratio.convective_ice,
                    cloud_fraction_large_scale=state.cloud_fraction.large_scale,
                    cloud_fraction_convective=state.cloud_fraction.convective,
                    concentration_ice=state.concentration.ice,
                    concentration_liquid=state.concentration.liquid,
                    relative_humidity_after_pdf=state.relative_humidity_after_pdf,
                    estimated_inversion_strength=state.estimated_inversion_strength,
                    area=state.area,
                    critical_relative_humidity_for_pdf=state.critical_relative_humidity_for_pdf,
                    pdf_iters=state.hydrostatic_pdf_iterations,
                    cloud_liquid_evaporation=state.cloud_liquid_evaporation,
                    cloud_ice_sublimation=state.cloud_ice_sublimation,
                    convection_fraction=state.convection_fraction,
                    surface_type=state.surface_type,
                    local_lcl_level=locals_.lcl_level,
                    local_p_mb=locals_.p_mb,
                    local_p_interface_mb=locals_.p_interface_mb,
                    local_saturation_specific_humidity=locals_.saturation_specific_humidity,
                )
            e = time.perf_counter()
            ndsl_log.info(f"GFDL1M Phase Change micro bench: {(e - s) / ts:.4f}s")

        return outputs
