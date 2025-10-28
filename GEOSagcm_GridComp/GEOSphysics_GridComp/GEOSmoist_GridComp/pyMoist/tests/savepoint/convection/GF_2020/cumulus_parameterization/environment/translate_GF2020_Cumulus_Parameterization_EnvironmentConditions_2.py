from ndsl import Namelist, StencilFactory
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
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.convection.GF_2020.cumulus_parameterization.constants import MAXENS1, MAXENS2, MAXENS3
from pyMoist.convection.GF_2020.cumulus_parameterization.environment.environment import EnvironmentConditions
from ndsl.dsl.typing import Int
from pyMoist.convection.GF_2020.cumulus_parameterization.setup.set_constants import set_constants


class TranslateGF2020_CumulusParameterization_EnvironmentConditions_2_shallow(TranslateFortranData2Py):
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
            "geopotential_height": {"serialname": "geopotential_height_env_cond"},
            "saturation_mixing_ratio_forced": {
                "serialname": "local_env_saturation_mixing_ratio_forced_env_cond"
            },
            "moist_static_energy_forced": {"serialname": "local_env_moist_static_energy_forced_env_cond"},
            "saturation_moist_static_energy_forced": {
                "serialname": "local_env_saturation_moist_static_energy_forced_env_cond"
            },
            "t_new": {"serialname": "local_t_new_env_cond"},
            "vapor_new": {"serialname": "local_vapor_new_env_cond"},
            "p": {"serialname": "p_env_cond"},
            "topography_height_no_negative": {"serialname": "topography_height_no_negative_env_cond"},
            "p_surface": {"serialname": "p_surface_env_cond"},
            "error_code": {"serialname": "error_code_env_cond"},
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

        # initalize saturation tables
        saturation_tables = SaturationVaporPressureTable(self.stencil_factory.backend)

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
        state.input_output.geopotential_height.data[:] = inputs["geopotential_height"]
        locals.environment_saturation_mixing_ratio_forced.data[:] = inputs["saturation_mixing_ratio_forced"]
        locals.environment_moist_static_energy_forced.data[:] = inputs["moist_static_energy_forced"]
        locals.environment_saturation_moist_static_energy_forced.data[:] = inputs[
            "saturation_moist_static_energy_forced"
        ]

        locals.t_new.data[:] = inputs["t_new"]
        locals.vapor_new.data[:] = inputs["vapor_new"]
        state.input_output.p.data[:] = inputs["p"]
        state.input_output.topography_height_no_negative.data[:] = inputs["topography_height_no_negative"]
        state.input_output.p_surface.data[:] = inputs["p_surface"]
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["error_code"]

        environment_conditions = EnvironmentConditions(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        if plume_dependent_constants.ENABLE_PLUME == 1:
            environment_conditions(
                state=state,
                locals=locals,
                saturation_tables=saturation_tables,
                plume_dependent_constants=plume_dependent_constants,
                data_type=1,
            )

        outputs = {
            # state fields
            "geopotential_height": state.input_output.geopotential_height.field[:],
            "saturation_mixing_ratio_forced": locals.environment_saturation_mixing_ratio_forced.field[:],
            "moist_static_energy_forced": locals.environment_moist_static_energy_forced.field[:],
            "saturation_moist_static_energy_forced": locals.environment_saturation_moist_static_energy_forced.field[
                :
            ],
            "t_old": locals.t_new.field[:],
            "vapor_old": locals.vapor_new.field[:],
            "p": state.input_output.p.field[:],
            "topography_height_no_negative": state.input_output.topography_height_no_negative.field[:],
            "p_surface": state.input_output.p_surface.field[:],
            "error_code": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_EnvironmentConditions_2_mid(TranslateFortranData2Py):
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
            "geopotential_height": {"serialname": "geopotential_height_env_cond"},
            "saturation_mixing_ratio_forced": {
                "serialname": "local_env_saturation_mixing_ratio_forced_env_cond"
            },
            "moist_static_energy_forced": {"serialname": "local_env_moist_static_energy_forced_env_cond"},
            "saturation_moist_static_energy_forced": {
                "serialname": "local_env_saturation_moist_static_energy_forced_env_cond"
            },
            "t_new": {"serialname": "local_t_new_env_cond"},
            "vapor_new": {"serialname": "local_vapor_new_env_cond"},
            "p": {"serialname": "p_env_cond"},
            "topography_height_no_negative": {"serialname": "topography_height_no_negative_env_cond"},
            "p_surface": {"serialname": "p_surface_env_cond"},
            "error_code": {"serialname": "error_code_env_cond"},
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

        # initalize saturation tables
        saturation_tables = SaturationVaporPressureTable(self.stencil_factory.backend)

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
        state.input_output.geopotential_height.data[:] = inputs["geopotential_height"]
        locals.environment_saturation_mixing_ratio_forced.data[:] = inputs["saturation_mixing_ratio_forced"]
        locals.environment_moist_static_energy_forced.data[:] = inputs["moist_static_energy_forced"]
        locals.environment_saturation_moist_static_energy_forced.data[:] = inputs[
            "saturation_moist_static_energy_forced"
        ]

        locals.t_new.data[:] = inputs["t_new"]
        locals.vapor_new.data[:] = inputs["vapor_new"]
        state.input_output.p.data[:] = inputs["p"]
        state.input_output.topography_height_no_negative.data[:] = inputs["topography_height_no_negative"]
        state.input_output.p_surface.data[:] = inputs["p_surface"]
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["error_code"]

        environment_conditions = EnvironmentConditions(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        if plume_dependent_constants.ENABLE_PLUME == 1:
            environment_conditions(
                state=state,
                locals=locals,
                saturation_tables=saturation_tables,
                plume_dependent_constants=plume_dependent_constants,
                data_type=1,
            )

        outputs = {
            # state fields
            "geopotential_height": state.input_output.geopotential_height.field[:],
            "saturation_mixing_ratio_forced": locals.environment_saturation_mixing_ratio_forced.field[:],
            "moist_static_energy_forced": locals.environment_moist_static_energy_forced.field[:],
            "saturation_moist_static_energy_forced": locals.environment_saturation_moist_static_energy_forced.field[
                :
            ],
            "t_old": locals.t_new.field[:],
            "vapor_old": locals.vapor_new.field[:],
            "p": state.input_output.p.field[:],
            "topography_height_no_negative": state.input_output.topography_height_no_negative.field[:],
            "p_surface": state.input_output.p_surface.field[:],
            "error_code": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_EnvironmentConditions_2_deep(TranslateFortranData2Py):
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
            "geopotential_height": {"serialname": "geopotential_height_env_cond"},
            "saturation_mixing_ratio_forced": {
                "serialname": "local_env_saturation_mixing_ratio_forced_env_cond"
            },
            "moist_static_energy_forced": {"serialname": "local_env_moist_static_energy_forced_env_cond"},
            "saturation_moist_static_energy_forced": {
                "serialname": "local_env_saturation_moist_static_energy_forced_env_cond"
            },
            "t_new": {"serialname": "local_t_new_env_cond"},
            "vapor_new": {"serialname": "local_vapor_new_env_cond"},
            "p": {"serialname": "p_env_cond"},
            "topography_height_no_negative": {"serialname": "topography_height_no_negative_env_cond"},
            "p_surface": {"serialname": "p_surface_env_cond"},
            "error_code": {"serialname": "error_code_env_cond"},
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

        # initalize saturation tables
        saturation_tables = SaturationVaporPressureTable(self.stencil_factory.backend)

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
        state.input_output.geopotential_height.data[:] = inputs["geopotential_height"]
        locals.environment_saturation_mixing_ratio_forced.data[:] = inputs["saturation_mixing_ratio_forced"]
        locals.environment_moist_static_energy_forced.data[:] = inputs["moist_static_energy_forced"]
        locals.environment_saturation_moist_static_energy_forced.data[:] = inputs[
            "saturation_moist_static_energy_forced"
        ]

        locals.t_new.data[:] = inputs["t_new"]
        locals.vapor_new.data[:] = inputs["vapor_new"]
        state.input_output.p.data[:] = inputs["p"]
        state.input_output.topography_height_no_negative.data[:] = inputs["topography_height_no_negative"]
        state.input_output.p_surface.data[:] = inputs["p_surface"]
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["error_code"]

        environment_conditions = EnvironmentConditions(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        if plume_dependent_constants.ENABLE_PLUME == 1:
            environment_conditions(
                state=state,
                locals=locals,
                saturation_tables=saturation_tables,
                plume_dependent_constants=plume_dependent_constants,
                data_type=1,
            )

        outputs = {
            # state fields
            "geopotential_height": state.input_output.geopotential_height.field[:],
            "saturation_mixing_ratio_forced": locals.environment_saturation_mixing_ratio_forced.field[:],
            "moist_static_energy_forced": locals.environment_moist_static_energy_forced.field[:],
            "saturation_moist_static_energy_forced": locals.environment_saturation_moist_static_energy_forced.field[
                :
            ],
            "t_old": locals.t_new.field[:],
            "vapor_old": locals.vapor_new.field[:],
            "p": state.input_output.p.field[:],
            "topography_height_no_negative": state.input_output.topography_height_no_negative.field[:],
            "p_surface": state.input_output.p_surface.field[:],
            "error_code": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
        }

        return outputs
