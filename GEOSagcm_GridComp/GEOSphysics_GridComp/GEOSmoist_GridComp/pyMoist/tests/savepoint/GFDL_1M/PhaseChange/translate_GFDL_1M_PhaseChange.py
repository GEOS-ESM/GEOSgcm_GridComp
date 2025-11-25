from f90nml import Namelist

from ndsl import StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.GFDL_1M.state import GFDL1MState
from pyMoist.GFDL_1M.locals import GFDL1MLocals
from pyMoist.GFDL_1M.config import GFDL1MConfig
from ndsl.stencils.testing.savepoint import DataLoader
from pyMoist.GFDL_1M.PhaseChange.sublimate import sublimate


class TranslateGFDL_1M_PhaseChange(TranslateFortranData2Py):
    def __init__(self, grid: Grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
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
            "p_interface_mb": {},
            "area": {},
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

        # initalize constants
        config = GFDL1MConfig(**self.constants)

        # initalize dataclasses
        state = GFDL1MState.zeros(self.quantity_factory)
        locals = GFDL1MLocals.zeros(self.quantity_factory)

        # Initalize saturation tables
        self.saturation_tables = SaturationVaporPressureTable(self.stencil_factory.backend)

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
        locals.lcl_level.field[:] = inputs["local_lcl_level"]
        locals.p_mb.field[:] = inputs["local_p_mb"]
        locals.p_interface_mb.field[:] = inputs["p_interface_mb"]
        state.area.field[:] = inputs["area"]
        locals.saturation_specific_humidity.field[:] = inputs["local_saturation_specific_humidity"]
        state.cloud_liquid_evaporation.field[:] = inputs["cloud_liquid_evaporation"]
        state.cloud_ice_sublimation.field[:] = inputs["cloud_ice_sublimation"]
        state.relative_humidity_after_pdf.field[:] = inputs["relative_humidity_after_pdf"]
        state.critical_relative_humidity_for_pdf.field[:] = inputs["critical_relative_humidity_for_pdf"]

        # construct test stencil
        code = self.stencil_factory.from_dims_halo(
            func=sublimate,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"DT_MOIST": config.DT_MOIST, "CCI_EVAP_EFF": config.CCI_EVAP_EFF},
        )
        code(
            p_mb=locals.p_mb,
            t=state.t,
            vapor=state.mixing_ratio.vapor,
            convective_liquid=state.mixing_ratio.convective_liquid,
            convective_ice=state.mixing_ratio.convective_ice,
            convective_cloud_fraction=state.cloud_fraction.convective,
            liquid_concentration=state.concentration.liquid,
            ice_concentration=state.concentration.ice,
            saturation_specific_humidity=locals.saturation_specific_humidity,
            sublimation=state.cloud_ice_sublimation,
        )

        return {
            "estimated_inversion_strength": state.estimated_inversion_strength.field[:],
            "convection_fraction": state.convection_fraction.field[:],
            "surface_type": state.surface_type.field[:],
            "t": state.t.field[:],
            "mixing_ratio_vapor": state.mixing_ratio.vapor.field[:],
            "concentration_liquid": state.concentration.liquid.field[:],
            "concentration_ice": state.concentration.ice.field[:],
            "mixing_ratio_convective_liquid": state.mixing_ratio.convective_liquid.field[:],
            "mixing_ratio_convective_ice": state.mixing_ratio.convective_ice.field[:],
            "mixing_ratio_large_scale_ice": state.mixing_ratio.large_scale_ice.field[:],
            "mixing_ratio_large_scale_liquid": state.mixing_ratio.large_scale_liquid.field[:],
            "cloud_fraction_large_scale": state.cloud_fraction.large_scale.field[:],
            "cloud_fraction_convective": state.cloud_fraction.convective.field[:],
            "local_lcl_level": locals.lcl_level.field[:],
            "local_p_mb": locals.p_mb.field[:],
            "p_interface_mb": locals.p_interface_mb.field[:],
            "area": state.area.field[:],
            "local_saturation_specific_humidity": locals.saturation_specific_humidity.field[:],
            "cloud_liquid_evaporation": state.cloud_liquid_evaporation.field[:],
            "cloud_ice_sublimation": state.cloud_ice_sublimation.field[:],
            "relative_humidity_after_pdf": state.relative_humidity_after_pdf.field[:],
            "critical_relative_humidity_for_pdf": state.critical_relative_humidity_for_pdf.field[:],
        }
