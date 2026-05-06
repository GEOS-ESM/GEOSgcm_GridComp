from f90nml import Namelist
from ndsl import StencilFactory
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py

from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.constants import MAXENS1, MAXENS2, MAXENS3, NUMBER_OF_PLUMES
from pyMoist.convection.GF_2020.cumulus_parameterization.diurnal_cycle import DiurnalCycle
from pyMoist.convection.GF_2020.cumulus_parameterization.locals import GF2020CumulusParameterizationLocals
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import GF2020PlumeDependentConstants
from pyMoist.convection.GF_2020.cumulus_parameterization.setup.set_constants import set_constants
from pyMoist.convection.GF_2020.cumulus_parameterization.state import GF2020CumulusParameterizationState


class TestCore:
    def __init__(
        self,
        grid: Grid,
        stencil_factory: StencilFactory,
        in_vars: dict,
        out_vars: dict,
    ):
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        in_vars["data_vars"] = {
            "error_code": {},
            "updraft_lfc_level": {},
            "cloud_top_level": {},
            "pbl_level": {},
            "grid_length": {},
            "ocean_fraction": {},
            "local_geopotential_height_cloud_levels_forced": {},
            "topography_height_no_negative": {},
            "t_old": {},
            "local_t_new": {},
            "local_t_cloud_levels_forced": {},
            "vapor_old": {},
            "local_vapor_forced": {},
            "u": {},
            "v": {},
            "local_vertical_velocity_2d": {},
            "local_cape_removal_time_scale": {},
            "cape_removal_time_scale": {},
            "local_pbl_time_scale": {},
            "pbl_time_scale": {},
            "local_cloud_work_function_1_pbl": {},
            "local_cloud_work_function_1_fa": {},
        }

        out_vars.update(in_vars["data_vars"])

    def __call__(self, constants: dict, cu_param_constants: dict, plume: str, **inputs):
        # initialize constants
        config = GF2020Config(**constants)
        cumulus_parameterization_config = GF2020CumulusParameterizationConfig(**cu_param_constants)
        plume_dependent_constants = GF2020PlumeDependentConstants()
        plume_dependent_constants = set_constants(cumulus_parameterization_config, plume_dependent_constants, plume)

        # initialize dataclasses
        state = GF2020CumulusParameterizationState.zeros(
            self.quantity_factory,
            data_dimensions={
                "plumes": NUMBER_OF_PLUMES,
                "convection_tracers": config.NUMBER_OF_TRACERS,
            },
        )

        locals = GF2020CumulusParameterizationLocals.zeros(
            self.quantity_factory,
            data_dimensions={
                "ensemble_1": MAXENS1,
                "ensemble_2": MAXENS2,
                "ensemble_3": MAXENS3,
                "ensemble_members": MAXENS1 * MAXENS2 * MAXENS3,
                "convection_tracers": config.NUMBER_OF_TRACERS,
            },
        )

        # fill relevant parts of dataclasses
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["error_code"]
        state.output.updraft_lfc_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["updraft_lfc_level"] - 1
        state.output.cloud_top_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["cloud_top_level"] - 1
        state.input_output.pbl_level.data[:] = inputs["pbl_level"] - 1
        state.input_output.grid_length.data[:] = inputs["grid_length"]
        state.input.ocean_fraction.data[:] = inputs["ocean_fraction"]
        locals.geopotential_height_cloud_levels_forced.data[:] = inputs["local_geopotential_height_cloud_levels_forced"]
        state.input_output.topography_height_no_negative.data[:] = inputs["topography_height_no_negative"]
        state.input_output.t_old.data[:] = inputs["t_old"]
        locals.t_new.data[:] = inputs["local_t_new"]
        locals.t_cloud_levels_forced.data[:] = inputs["local_t_cloud_levels_forced"]
        state.input_output.vapor_old.data[:] = inputs["vapor_old"]
        locals.vapor_forced.data[:] = inputs["local_vapor_forced"]
        state.input_output.u.data[:] = inputs["u"]
        state.input_output.v.data[:] = inputs["v"]
        locals.vertical_velocity_2d.data[:] = inputs["local_vertical_velocity_2d"]
        locals.cape_removal_time_scale.data[:] = inputs["local_cape_removal_time_scale"]
        state.output.cape_removal_time_scale.data[:] = inputs["cape_removal_time_scale"]
        locals.pbl_time_scale.data[:] = inputs["local_pbl_time_scale"]
        state.output.pbl_time_scale.data[:] = inputs["pbl_time_scale"]
        locals.cloud_workfunction_1_pbl.data[:] = inputs["local_cloud_work_function_1_pbl"]
        locals.cloud_workfunction_1_fa.data[:] = inputs["local_cloud_work_function_1_fa"]

        # initialize test code
        code = DiurnalCycle(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        # call test code
        if plume_dependent_constants.ENABLE_PLUME == 1:
            code(
                error_code=state.output.error_code,
                updraft_lfc_level=state.output.updraft_lfc_level,
                cloud_top_level=state.output.cloud_top_level,
                pbl_level=state.input_output.pbl_level,
                grid_length=state.input_output.grid_length,
                ocean_fraction=state.input.ocean_fraction,
                topography_height_no_negative=state.input_output.topography_height_no_negative,
                geopotential_height_cloud_levels_forced=locals.geopotential_height_cloud_levels_forced,
                t_old=state.input_output.t_old,
                t_new=locals.t_new,
                t_cloud_levels_forced=locals.t_cloud_levels_forced,
                vapor_old=state.input_output.vapor_old,
                vapor_forced=locals.vapor_forced,
                u=state.input_output.u,
                v=state.input_output.v,
                vertical_velocity_2d=locals.vertical_velocity_2d,
                cape_removal_time_scale=locals.cape_removal_time_scale,
                cape_removal_time_scale_from_state=state.output.cape_removal_time_scale,
                pbl_time_scale=locals.pbl_time_scale,
                pbl_time_scale_from_state=state.output.pbl_time_scale,
                cloud_work_function_1_pbl=locals.cloud_workfunction_1_pbl,
                cloud_work_function_1_fa=locals.cloud_workfunction_1_fa,
                plume_dependent_constants=plume_dependent_constants,
            )

        # write output
        outputs = {
            "error_code": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "updraft_lfc_level": state.output.updraft_lfc_level.field[:, :, plume_dependent_constants.PLUME_INDEX] + 1,
            "cloud_top_level": state.output.cloud_top_level.field[:, :, plume_dependent_constants.PLUME_INDEX] + 1,
            "pbl_level": state.input_output.pbl_level.field[:] + 1,
            "grid_length": state.input_output.grid_length.field[:],
            "ocean_fraction": state.input.ocean_fraction.field[:],
            "local_geopotential_height_cloud_levels_forced": locals.geopotential_height_cloud_levels_forced.field[:],
            "topography_height_no_negative": state.input_output.topography_height_no_negative.field[:],
            "t_old": state.input_output.t_old.field[:],
            "local_t_new": locals.t_new.field[:],
            "local_t_cloud_levels_forced": locals.t_cloud_levels_forced.field[:],
            "vapor_old": state.input_output.vapor_old.field[:],
            "local_vapor_forced": locals.vapor_forced.field[:],
            "u": state.input_output.u.field[:],
            "v": state.input_output.v.field[:],
            "local_vertical_velocity_2d": locals.vertical_velocity_2d.field[:],
            "local_cape_removal_time_scale": locals.cape_removal_time_scale.field[:],
            "cape_removal_time_scale": state.output.cape_removal_time_scale.field[:],
            "local_pbl_time_scale": locals.pbl_time_scale.field[:],
            "pbl_time_scale": state.output.pbl_time_scale.field[:],
            "local_cloud_work_function_1_pbl": locals.cloud_workfunction_1_pbl.field[:],
            "local_cloud_work_function_1_fa": locals.cloud_workfunction_1_fa.field[:],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_DiurnalCycle_shallow(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        _namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)

        self.test_core = TestCore(grid, stencil_factory, self.in_vars, self.out_vars)

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF2020-constants")
        self.cu_param_constants = data_loader.load("GF2020_CumulusParameterization-constants")

    def compute_func(self, **inputs):
        outputs = self.test_core(self.constants, self.cu_param_constants, "shallow", **inputs)

        return outputs


class TranslateGF2020_CumulusParameterization_DiurnalCycle_mid(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        _namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)

        self.test_core = TestCore(grid, stencil_factory, self.in_vars, self.out_vars)

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF2020-constants")
        self.cu_param_constants = data_loader.load("GF2020_CumulusParameterization-constants")

    def compute_func(self, **inputs):
        outputs = self.test_core(self.constants, self.cu_param_constants, "mid", **inputs)

        return outputs


class TranslateGF2020_CumulusParameterization_DiurnalCycle_deep(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        _namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)

        self.test_core = TestCore(grid, stencil_factory, self.in_vars, self.out_vars)

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF2020-constants")
        self.cu_param_constants = data_loader.load("GF2020_CumulusParameterization-constants")

    def compute_func(self, **inputs):
        outputs = self.test_core(self.constants, self.cu_param_constants, "deep", **inputs)

        return outputs
