from f90nml import Namelist

from ndsl import StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.GFDL_1M.locals import GFDL1MLocals
from pyMoist.GFDL_1M.PhaseChange.sublimate import sublimate
from pyMoist.GFDL_1M.state import GFDL1MState
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable


class TranslateGFDL_1M_Sublimate(TranslateFortranData2Py):
    def __init__(self, grid: Grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "local_p_mb": {},
            "t": {},
            "mixing_ratio_vapor": {},
            "mixing_ratio_convective_liquid": {},
            "mixing_ratio_convective_ice": {},
            "cloud_fraction_convective": {},
            "concentration_liquid": {},
            "concentration_ice": {},
            "local_saturation_specific_humidity": {},
            "cloud_ice_sublimation": {},
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

        locals.p_mb.field[:] = inputs["local_p_mb"]
        state.t.field[:] = inputs["t"]
        state.mixing_ratio.vapor.field[:] = inputs["mixing_ratio_vapor"]
        state.mixing_ratio.convective_liquid.field[:] = inputs["mixing_ratio_convective_liquid"]
        state.mixing_ratio.convective_ice.field[:] = inputs["mixing_ratio_convective_ice"]
        state.cloud_fraction.convective.field[:] = inputs["cloud_fraction_convective"]
        state.concentration.liquid.field[:] = inputs["concentration_liquid"]
        state.concentration.ice.field[:] = inputs["concentration_ice"]
        locals.saturation_specific_humidity.field[:] = inputs["local_saturation_specific_humidity"]
        state.cloud_ice_sublimation.field[:] = inputs["cloud_ice_sublimation"]

        # construct test stencil
        code = self.stencil_factory.from_dims_halo(
            func=sublimate,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "DT_MOIST": config.DT_MOIST,
                "CCI_EVAP_EFF": config.CCI_EVAP_EFF,
            },
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
            "local_p_mb": locals.p_mb.field[:],
            "t": state.t.field[:],
            "mixing_ratio_vapor": state.mixing_ratio.vapor.field[:],
            "mixing_ratio_convective_liquid": state.mixing_ratio.convective_liquid.field[:],
            "mixing_ratio_convective_ice": state.mixing_ratio.convective_ice.field[:],
            "cloud_fraction_convective": state.cloud_fraction.convective.field[:],
            "concentration_liquid": state.concentration.liquid.field[:],
            "concentration_ice": state.concentration.ice.field[:],
            "local_saturation_specific_humidity": locals.saturation_specific_humidity.field[:],
            "cloud_ice_sublimation": state.cloud_ice_sublimation.field[:],
        }
