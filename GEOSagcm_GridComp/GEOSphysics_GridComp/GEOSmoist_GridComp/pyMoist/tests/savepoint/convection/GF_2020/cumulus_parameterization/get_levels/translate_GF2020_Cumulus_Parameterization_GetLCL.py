from f90nml import Namelist
from ndsl import StencilFactory
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.convection.GF_2020.cumulus_parameterization.state import GF2020CumulusParameterizationState
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.locals import GF2020CumulusParameterizationLocals
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.constants import MAXENS1, MAXENS2, MAXENS3
from pyMoist.convection.GF_2020.cumulus_parameterization.get_levels.get_levels import (
    GetLCL,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.setup.set_constants import set_constants


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

        self.in_vars["data_vars"] = {
            "p_forced_getlcl": {},
            "local_p_cloud_levels_getlcl": {},
            "local_t_excess_getlcl": {},
            "local_t_cloud_levels_getlcl": {},
            "t_perturbation_getlcl": {},
            "local_vapor_excess_getlcl": {},
            "local_vapor_cloud_levels_forced_getlcl": {},
            "omega_getlcl": {},
            "air_density_getlcl": {},
            "local_geopotential_height_cloud_levels_getlcl": {},
            "topography_height_no_negative_getlcl": {},
            "ocean_fraction_getlcl": {},
            "local_updraft_origin_level_getlcl": {},
            "grid_length_getlcl": {},
            "lcl_level_getlcl": {},
            "error_code_getlcl": {},
        }

        self.out_vars = self.in_vars["data_vars"].copy()

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF2020-constants")
        self.cu_param_constants = data_loader.load("GF2020_CumulusParameterization-constants")

    def compute_func(self, **inputs):
        # initalize constants
        config = GF2020Config(SINGLE_COLUMN_MODE=False, **self.constants)
        cumulus_parameterization_config = GF2020CumulusParameterizationConfig(**self.cu_param_constants)
        plume_dependent_constants = GF2020PlumeDependentConstants()
        plume_dependent_constants = set_constants(
            cumulus_parameterization_config, plume_dependent_constants, "shallow"
        )

        # initalize dataclasses
        state = GF2020CumulusParameterizationState.zeros(
            self.quantity_factory,
            data_dimensions={
                "plumes": 3,
            },
        )

        locals = GF2020CumulusParameterizationLocals.zeros(
            self.quantity_factory,
            data_dimensions={
                "ensemble_members": MAXENS1 * MAXENS2 * MAXENS3,
            },
        )

        # fill relevant parts of dataclasses
        state.input_output.p_forced.data[:] = inputs["p_forced_getlcl"]
        locals.p_cloud_levels.data[:] = inputs["local_p_cloud_levels_getlcl"]
        locals.t_excess.data[:] = inputs["local_t_excess_getlcl"]
        locals.t_cloud_levels.data[:] = inputs["local_t_cloud_levels_getlcl"]
        state.output.t_perturbation.data[:] = inputs["t_perturbation_getlcl"]
        locals.vapor_excess.data[:] = inputs["local_vapor_excess_getlcl"]
        locals.vapor_cloud_levels.data[:] = inputs["local_vapor_cloud_levels_forced_getlcl"]
        state.input_output.omega.data[:] = inputs["omega_getlcl"]
        state.input_output.air_density.data[:] = inputs["air_density_getlcl"]
        locals.geopotential_height_cloud_levels.data[:] = inputs[
            "local_geopotential_height_cloud_levels_getlcl"
        ]
        state.input_output.topography_height_no_negative.data[:] = inputs[
            "topography_height_no_negative_getlcl"
        ]
        state.input.ocean_fraction.data[:] = inputs["ocean_fraction_getlcl"]
        locals.updraft_origin_level.data[:] = inputs["local_updraft_origin_level_getlcl"] - 1
        state.input_output.grid_length.data[:] = inputs["grid_length_getlcl"]
        state.output.lcl_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["lcl_level_getlcl"]
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "error_code_getlcl"
        ]

        # initalize test code
        code = GetLCL(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        # call test code
        if plume_dependent_constants.ENABLE_PLUME == 1:
            code(
                state=state,
                locals=locals,
                plume_dependent_constants=plume_dependent_constants,
            )

            locals.updraft_origin_level.field[:] += 1
            state.output.lcl_level.field[:, :, plume_dependent_constants.PLUME_INDEX] += 1

        # write output
        outputs = {
            "p_forced_getlcl": state.input_output.p_forced.field[:],
            "local_p_cloud_levels_getlcl": locals.p_cloud_levels.field[:],
            "local_t_excess_getlcl": locals.t_excess.field[:],
            "local_t_cloud_levels_getlcl": locals.t_cloud_levels.field[:],
            "t_perturbation_getlcl": state.output.t_perturbation.field[:],
            "local_vapor_excess_getlcl": locals.vapor_excess.field[:],
            "local_vapor_cloud_levels_forced_getlcl": locals.vapor_cloud_levels.field[:],
            "omega_getlcl": state.input_output.omega.field[:],
            "air_density_getlcl": state.input_output.air_density.field[:],
            "local_geopotential_height_cloud_levels_getlcl": locals.geopotential_height_cloud_levels.field[:],
            "topography_height_no_negative_getlcl": state.input_output.topography_height_no_negative.field[:],
            "ocean_fraction_getlcl": state.input.ocean_fraction.field[:],
            "local_updraft_origin_level_getlcl": locals.updraft_origin_level.field[:],
            "grid_length_getlcl": state.input_output.grid_length.field[:],
            "lcl_level_getlcl": state.output.lcl_level.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "error_code_getlcl": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
        }

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

        self.in_vars["data_vars"] = {
            "p_forced_getlcl": {},
            "local_p_cloud_levels_getlcl": {},
            "local_t_excess_getlcl": {},
            "local_t_cloud_levels_getlcl": {},
            "t_perturbation_getlcl": {},
            "local_vapor_excess_getlcl": {},
            "local_vapor_cloud_levels_forced_getlcl": {},
            "omega_getlcl": {},
            "air_density_getlcl": {},
            "local_geopotential_height_cloud_levels_getlcl": {},
            "topography_height_no_negative_getlcl": {},
            "ocean_fraction_getlcl": {},
            "local_updraft_origin_level_getlcl": {},
            "grid_length_getlcl": {},
            "lcl_level_getlcl": {},
            "error_code_getlcl": {},
        }

        self.out_vars = self.in_vars["data_vars"].copy()

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF2020-constants")
        self.cu_param_constants = data_loader.load("GF2020_CumulusParameterization-constants")

    def compute_func(self, **inputs):
        # initalize constants
        config = GF2020Config(SINGLE_COLUMN_MODE=False, **self.constants)
        cumulus_parameterization_config = GF2020CumulusParameterizationConfig(**self.cu_param_constants)
        plume_dependent_constants = GF2020PlumeDependentConstants()
        plume_dependent_constants = set_constants(
            cumulus_parameterization_config, plume_dependent_constants, "mid"
        )

        # initalize dataclasses
        state = GF2020CumulusParameterizationState.zeros(
            self.quantity_factory,
            data_dimensions={
                "plumes": 3,
            },
        )

        locals = GF2020CumulusParameterizationLocals.zeros(
            self.quantity_factory,
            data_dimensions={
                "ensemble_members": MAXENS1 * MAXENS2 * MAXENS3,
            },
        )

        # fill relevant parts of dataclasses
        state.input_output.p_forced.data[:] = inputs["p_forced_getlcl"]
        locals.p_cloud_levels.data[:] = inputs["local_p_cloud_levels_getlcl"]
        locals.t_excess.data[:] = inputs["local_t_excess_getlcl"]
        locals.t_cloud_levels.data[:] = inputs["local_t_cloud_levels_getlcl"]
        state.output.t_perturbation.data[:] = inputs["t_perturbation_getlcl"]
        locals.vapor_excess.data[:] = inputs["local_vapor_excess_getlcl"]
        locals.vapor_cloud_levels.data[:] = inputs["local_vapor_cloud_levels_forced_getlcl"]
        state.input_output.omega.data[:] = inputs["omega_getlcl"]
        state.input_output.air_density.data[:] = inputs["air_density_getlcl"]
        locals.geopotential_height_cloud_levels.data[:] = inputs[
            "local_geopotential_height_cloud_levels_getlcl"
        ]
        state.input_output.topography_height_no_negative.data[:] = inputs[
            "topography_height_no_negative_getlcl"
        ]
        state.input.ocean_fraction.data[:] = inputs["ocean_fraction_getlcl"]
        locals.updraft_origin_level.data[:] = inputs["local_updraft_origin_level_getlcl"] - 1
        state.input_output.grid_length.data[:] = inputs["grid_length_getlcl"]
        state.output.lcl_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["lcl_level_getlcl"]
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "error_code_getlcl"
        ]

        # initalize test code
        code = GetLCL(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        # call test code
        if plume_dependent_constants.ENABLE_PLUME == 1:
            code(
                state=state,
                locals=locals,
                plume_dependent_constants=plume_dependent_constants,
            )

            locals.updraft_origin_level.field[:] += 1
            state.output.lcl_level.field[:, :, plume_dependent_constants.PLUME_INDEX] += 1

        # write output
        outputs = {
            "p_forced_getlcl": state.input_output.p_forced.field[:],
            "local_p_cloud_levels_getlcl": locals.p_cloud_levels.field[:],
            "local_t_excess_getlcl": locals.t_excess.field[:],
            "local_t_cloud_levels_getlcl": locals.t_cloud_levels.field[:],
            "t_perturbation_getlcl": state.output.t_perturbation.field[:],
            "local_vapor_excess_getlcl": locals.vapor_excess.field[:],
            "local_vapor_cloud_levels_forced_getlcl": locals.vapor_cloud_levels.field[:],
            "omega_getlcl": state.input_output.omega.field[:],
            "air_density_getlcl": state.input_output.air_density.field[:],
            "local_geopotential_height_cloud_levels_getlcl": locals.geopotential_height_cloud_levels.field[:],
            "topography_height_no_negative_getlcl": state.input_output.topography_height_no_negative.field[:],
            "ocean_fraction_getlcl": state.input.ocean_fraction.field[:],
            "local_updraft_origin_level_getlcl": locals.updraft_origin_level.field[:],
            "grid_length_getlcl": state.input_output.grid_length.field[:],
            "lcl_level_getlcl": state.output.lcl_level.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "error_code_getlcl": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
        }

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

        self.in_vars["data_vars"] = {
            "p_forced_getlcl": {},
            "local_p_cloud_levels_getlcl": {},
            "local_t_excess_getlcl": {},
            "local_t_cloud_levels_getlcl": {},
            "t_perturbation_getlcl": {},
            "local_vapor_excess_getlcl": {},
            "local_vapor_cloud_levels_forced_getlcl": {},
            "omega_getlcl": {},
            "air_density_getlcl": {},
            "local_geopotential_height_cloud_levels_getlcl": {},
            "topography_height_no_negative_getlcl": {},
            "ocean_fraction_getlcl": {},
            "local_updraft_origin_level_getlcl": {},
            "grid_length_getlcl": {},
            "lcl_level_getlcl": {},
            "error_code_getlcl": {},
        }

        self.out_vars = self.in_vars["data_vars"].copy()

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF2020-constants")
        self.cu_param_constants = data_loader.load("GF2020_CumulusParameterization-constants")

    def compute_func(self, **inputs):
        # initalize constants
        config = GF2020Config(SINGLE_COLUMN_MODE=False, **self.constants)
        cumulus_parameterization_config = GF2020CumulusParameterizationConfig(**self.cu_param_constants)
        plume_dependent_constants = GF2020PlumeDependentConstants()
        plume_dependent_constants = set_constants(
            cumulus_parameterization_config, plume_dependent_constants, "deep"
        )

        # initalize dataclasses
        state = GF2020CumulusParameterizationState.zeros(
            self.quantity_factory,
            data_dimensions={
                "plumes": 3,
            },
        )

        locals = GF2020CumulusParameterizationLocals.zeros(
            self.quantity_factory,
            data_dimensions={
                "ensemble_members": MAXENS1 * MAXENS2 * MAXENS3,
            },
        )

        # fill relevant parts of dataclasses
        state.input_output.p_forced.data[:] = inputs["p_forced_getlcl"]
        locals.p_cloud_levels.data[:] = inputs["local_p_cloud_levels_getlcl"]
        locals.t_excess.data[:] = inputs["local_t_excess_getlcl"]
        locals.t_cloud_levels.data[:] = inputs["local_t_cloud_levels_getlcl"]
        state.output.t_perturbation.data[:] = inputs["t_perturbation_getlcl"]
        locals.vapor_excess.data[:] = inputs["local_vapor_excess_getlcl"]
        locals.vapor_cloud_levels.data[:] = inputs["local_vapor_cloud_levels_forced_getlcl"]
        state.input_output.omega.data[:] = inputs["omega_getlcl"]
        state.input_output.air_density.data[:] = inputs["air_density_getlcl"]
        locals.geopotential_height_cloud_levels.data[:] = inputs[
            "local_geopotential_height_cloud_levels_getlcl"
        ]
        state.input_output.topography_height_no_negative.data[:] = inputs[
            "topography_height_no_negative_getlcl"
        ]
        state.input.ocean_fraction.data[:] = inputs["ocean_fraction_getlcl"]
        locals.updraft_origin_level.data[:] = inputs["local_updraft_origin_level_getlcl"] - 1
        state.input_output.grid_length.data[:] = inputs["grid_length_getlcl"]
        state.output.lcl_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["lcl_level_getlcl"]
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "error_code_getlcl"
        ]

        # initalize test code
        code = GetLCL(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        # call test code
        if plume_dependent_constants.ENABLE_PLUME == 1:
            code(
                state=state,
                locals=locals,
                plume_dependent_constants=plume_dependent_constants,
            )

            locals.updraft_origin_level.field[:] += 1
            state.output.lcl_level.field[:, :, plume_dependent_constants.PLUME_INDEX] += 1

        # write output
        outputs = {
            "p_forced_getlcl": state.input_output.p_forced.field[:],
            "local_p_cloud_levels_getlcl": locals.p_cloud_levels.field[:],
            "local_t_excess_getlcl": locals.t_excess.field[:],
            "local_t_cloud_levels_getlcl": locals.t_cloud_levels.field[:],
            "t_perturbation_getlcl": state.output.t_perturbation.field[:],
            "local_vapor_excess_getlcl": locals.vapor_excess.field[:],
            "local_vapor_cloud_levels_forced_getlcl": locals.vapor_cloud_levels.field[:],
            "omega_getlcl": state.input_output.omega.field[:],
            "air_density_getlcl": state.input_output.air_density.field[:],
            "local_geopotential_height_cloud_levels_getlcl": locals.geopotential_height_cloud_levels.field[:],
            "topography_height_no_negative_getlcl": state.input_output.topography_height_no_negative.field[:],
            "ocean_fraction_getlcl": state.input.ocean_fraction.field[:],
            "local_updraft_origin_level_getlcl": locals.updraft_origin_level.field[:],
            "grid_length_getlcl": state.input_output.grid_length.field[:],
            "lcl_level_getlcl": state.output.lcl_level.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "error_code_getlcl": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
        }

        return outputs
