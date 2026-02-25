from f90nml import Namelist

from ndsl import StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.constants import (
    MAXENS1,
    MAXENS2,
    MAXENS3,
    NUMBER_OF_PLUMES,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.get_levels import find_lcl
from pyMoist.convection.GF_2020.cumulus_parameterization.locals import GF2020CumulusParameterizationLocals
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.setup.set_constants import set_constants
from pyMoist.convection.GF_2020.cumulus_parameterization.state import GF2020CumulusParameterizationState


class TestCore:
    def __init__(
        self,
        grid: Grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
        in_vars: dict,
        out_vars: dict,
    ):
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        in_vars["data_vars"] = {
            "p_forced": {},
            "local_p_cloud_levels": {},
            "local_t_excess": {},
            "local_t_cloud_levels": {},
            "t_perturbation": {},
            "local_vapor_excess": {},
            "local_vapor_cloud_levels_forced": {},
            "omega": {},
            "air_density": {},
            "local_geopotential_height_cloud_levels": {},
            "topography_height_no_negative": {},
            "ocean_fraction": {},
            "updraft_origin_level": {},
            "grid_length": {},
            "lcl_level": {},
            "error_code": {},
        }

        out_vars.update(in_vars["data_vars"])

    def __call__(self, constants: dict, cu_param_constants: dict, plume: str, **inputs):
        # initialize constants
        config = GF2020Config(SINGLE_COLUMN_MODE=False, **constants)
        cumulus_parameterization_config = GF2020CumulusParameterizationConfig(**cu_param_constants)
        plume_dependent_constants = GF2020PlumeDependentConstants()
        plume_dependent_constants = set_constants(
            cumulus_parameterization_config, plume_dependent_constants, plume
        )

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
        state.input_output.p_forced.data[:] = inputs["p_forced"]
        locals.p_cloud_levels.data[:] = inputs["local_p_cloud_levels"]
        locals.t_excess.data[:] = inputs["local_t_excess"]
        locals.t_cloud_levels.data[:] = inputs["local_t_cloud_levels"]
        state.output.t_perturbation.data[:] = inputs["t_perturbation"]
        locals.vapor_excess.data[:] = inputs["local_vapor_excess"]
        locals.vapor_cloud_levels.data[:] = inputs["local_vapor_cloud_levels_forced"]
        state.input_output.omega.data[:] = inputs["omega"]
        state.input_output.air_density.data[:] = inputs["air_density"]
        locals.geopotential_height_cloud_levels.data[:] = inputs["local_geopotential_height_cloud_levels"]
        state.input_output.topography_height_no_negative.data[:] = inputs["topography_height_no_negative"]
        state.input.ocean_fraction.data[:] = inputs["ocean_fraction"]
        state.output.updraft_origin_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["updraft_origin_level"] - 1
        )
        state.input_output.grid_length.data[:] = inputs["grid_length"]
        state.output.lcl_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["lcl_level"] - 1
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["error_code"]

        code = self.stencil_factory.from_dims_halo(
            func=find_lcl,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "BOUNDARY_CONDITION_METHOD": cumulus_parameterization_config.BOUNDARY_CONDITION_METHOD,
                "ADV_TRIGGER": config.ADV_TRIGGER,
            },
        )

        if plume_dependent_constants.ENABLE_PLUME == 1:
            code(
                p=state.input_output.p_forced,
                p_cloud_levels=locals.p_cloud_levels,
                t_excess=locals.t_excess,
                t_cloud_levels_forced=locals.t_cloud_levels,
                t_perturbation=state.output.t_perturbation,
                vapor_excess=locals.vapor_excess,
                vapor_cloud_levels_forced=locals.vapor_cloud_levels,
                omega=state.input_output.omega,
                air_density=state.input_output.air_density,
                geopotential_height_cloud_levels=locals.geopotential_height_cloud_levels,
                topography_height_no_negative=state.input_output.topography_height_no_negative,
                ocean_fraction=state.input.ocean_fraction,
                updraft_origin_level=state.output.updraft_origin_level,
                grid_length=state.input_output.grid_length,
                lcl_level=state.output.lcl_level,
                error_code=state.output.error_code,
                AVERAGE_LAYER_DEPTH=plume_dependent_constants.AVERAGE_LAYER_DEPTH,
                plume=plume_dependent_constants.PLUME_INDEX,
            )

        outputs = {
            "p_forced": state.input_output.p_forced.field[:],
            "local_p_cloud_levels": locals.p_cloud_levels.field[:],
            "local_t_excess": locals.t_excess.field[:],
            "local_t_cloud_levels": locals.t_cloud_levels.field[:],
            "t_perturbation": state.output.t_perturbation.field[:],
            "local_vapor_excess": locals.vapor_excess.field[:],
            "local_vapor_cloud_levels_forced": locals.vapor_cloud_levels.field[:],
            "omega": state.input_output.omega.field[:],
            "air_density": state.input_output.air_density.field[:],
            "local_geopotential_height_cloud_levels": locals.geopotential_height_cloud_levels.field[:],
            "topography_height_no_negative": state.input_output.topography_height_no_negative.field[:],
            "ocean_fraction": state.input.ocean_fraction.field[:],
            "updraft_origin_level": state.output.updraft_origin_level.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ]
            + 1,
            "grid_length": state.input_output.grid_length.field[:],
            "lcl_level": state.output.lcl_level.field[:, :, plume_dependent_constants.PLUME_INDEX] + 1,
            "error_code": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_GetLCL_shallow(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        self.test_core = TestCore(grid, namelist, stencil_factory, self.in_vars, self.out_vars)

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF2020-constants")
        self.cu_param_constants = data_loader.load("GF2020_CumulusParameterization-constants")

    def compute_func(self, **inputs):
        outputs = self.test_core(self.constants, self.cu_param_constants, "shallow", **inputs)

        return outputs


class TranslateGF2020_CumulusParameterization_GetLCL_mid(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        self.test_core = TestCore(grid, namelist, stencil_factory, self.in_vars, self.out_vars)

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF2020-constants")
        self.cu_param_constants = data_loader.load("GF2020_CumulusParameterization-constants")

    def compute_func(self, **inputs):
        outputs = self.test_core(self.constants, self.cu_param_constants, "mid", **inputs)

        return outputs


class TranslateGF2020_CumulusParameterization_GetLCL_deep(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        self.test_core = TestCore(grid, namelist, stencil_factory, self.in_vars, self.out_vars)

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF2020-constants")
        self.cu_param_constants = data_loader.load("GF2020_CumulusParameterization-constants")

    def compute_func(self, **inputs):
        outputs = self.test_core(self.constants, self.cu_param_constants, "deep", **inputs)

        return outputs
