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
from pyMoist.convection.GF_2020.cumulus_parameterization.environment.environment import EnvironmentCloudLevels
from ndsl.dsl.typing import Int
from pyMoist.convection.GF_2020.cumulus_parameterization.setup.set_constants import set_constants


class TranslateGF2020_CumulusParameterization_EnvironmentCloudLevels_2_shallow(TranslateFortranData2Py):
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
            "local_t_new": {"serialname": "local_t_new_env_clev"},
            "local_env_saturation_mixing_ratio_forced": {
                "serialname": "local_env_saturation_mixing_ratio_forced_env_clev"
            },
            "local_vapor_new": {"serialname": "local_vapor_new_env_cond"},
            "local_env_moist_static_energy_forced": {
                "serialname": "local_env_moist_static_energy_forced_env_clev"
            },
            "local_env_saturation_moist_static_energy_forced": {
                "serialname": "local_env_saturation_moist_static_energy_forced_env_clev"
            },
            "geopotential_height": {"serialname": "geopotential_height_env_clev"},
            "p": {"serialname": "p_env_clev"},
            "local_env_saturation_mixing_ratio_cloud_levels_forced": {
                "serialname": "local_env_saturation_mixing_ratio_cloud_levels_forced_env_clev"
            },
            "local_vapor_cloud_levels_forced": {"serialname": "local_vapor_cloud_levels_forced_env_clev"},
            "local_env_moist_static_energy_cloud_levels_forced": {
                "serialname": "local_env_moist_static_energy_cloud_levels_forced_env_clev"
            },
            "u": {"serialname": "u_env_clev"},
            "v": {"serialname": "v_env_clev"},
            "local_u_cloud_levels": {"serialname": "local_u_cloud_levels_env_clev"},
            "local_v_cloud_levels": {"serialname": "local_v_cloud_levels_env_clev"},
            "local_env_saturation_moist_static_energy_cloud_levels_forced": {
                "serialname": "local_env_saturation_moist_static_energy_cloud_levels_forced_env_clev"
            },
            "local_geopotential_height_cloud_levels_forced": {
                "serialname": "local_geopotential_height_cloud_levels_forced_env_clev"
            },
            "p_output": {"serialname": "p_output_env_clev"},
            "local_gamma_cloud_levels_forced": {"serialname": "local_gamma_cloud_levels_forced_env_clev"},
            "local_t_cloud_levels_forced": {"serialname": "local_t_cloud_levels_forced_env_clev"},
            "p_surface": {"serialname": "p_surface_env_clev"},
            "t_surface": {"serialname": "t_surface_env_clev"},
            "error_code": {"serialname": "error_code_env_clev"},
            "topography_height_no_negative": {"serialname": "topography_height_no_negative_env_clev"},
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
        locals.t_new.data[:] = inputs["local_t_new"]
        locals.environment_saturation_mixing_ratio_forced.data[:] = inputs[
            "local_env_saturation_mixing_ratio_forced"
        ]
        locals.vapor_new.data[:] = inputs["local_vapor_new"]
        locals.environment_moist_static_energy_forced.data[:] = inputs["local_env_moist_static_energy_forced"]
        locals.environment_saturation_moist_static_energy_forced.data[:] = inputs[
            "local_env_saturation_moist_static_energy_forced"
        ]
        state.input_output.geopotential_height.data[:] = inputs["geopotential_height"]
        state.input_output.p.data[:] = inputs["p"]
        locals.environment_saturation_mixing_ratio_cloud_levels.data[:] = inputs[
            "local_env_saturation_mixing_ratio_cloud_levels_forced"
        ]
        locals.vapor_cloud_levels_forced.data[:] = inputs["local_vapor_cloud_levels_forced"]
        locals.environment_moist_static_energy_cloud_levels.data[:] = inputs[
            "local_env_moist_static_energy_cloud_levels_forced"
        ]
        state.input_output.u.data[:] = inputs["u"]
        state.input_output.v.data[:] = inputs["v"]
        locals.u_cloud_levels.data[:] = inputs["local_u_cloud_levels"]
        locals.v_cloud_levels.data[:] = inputs["local_v_cloud_levels"]
        locals.environment_saturation_moist_static_energy_cloud_levels_forced.data[:] = inputs[
            "local_env_saturation_moist_static_energy_cloud_levels_forced"
        ]
        locals.geopotential_height_cloud_levels_forced.data[:] = inputs[
            "local_geopotential_height_cloud_levels_forced"
        ]
        state.output.p.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs["p_output"]
        locals.gamma_cloud_levels_forced.data[:] = inputs["local_gamma_cloud_levels_forced"]
        locals.t_cloud_levels_forced.data[:] = inputs["local_t_cloud_levels_forced"]
        state.input_output.p_surface.data[:] = inputs["p_surface"]
        state.input_output.t_surface.data[:] = inputs["t_surface"]
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["error_code"]
        state.input_output.topography_height_no_negative.data[:] = inputs["topography_height_no_negative"]

        environment_cloud_levels = EnvironmentCloudLevels(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        if plume_dependent_constants.ENABLE_PLUME == 1:
            environment_cloud_levels(
                state=state,
                locals=locals,
                saturation_tables=saturation_tables,
                plume_dependent_constants=plume_dependent_constants,
                data_type=1,
            )

        outputs = {
            # state fields
            "local_t_new": locals.t_new.field[:],
            "local_env_saturation_mixing_ratio_forced": locals.environment_saturation_mixing_ratio_forced.field[
                :
            ],
            "local_vapor_new": locals.vapor_new.field[:],
            "local_env_moist_static_energy_forced": locals.environment_moist_static_energy_forced.field[:],
            "local_env_saturation_moist_static_energy_forced": locals.environment_saturation_moist_static_energy_forced.field[
                :
            ],
            "geopotential_height": state.input_output.geopotential_height.field[:],
            "p": state.input_output.p.field[:],
            "local_env_saturation_mixing_ratio_cloud_levels_forced": locals.environment_saturation_mixing_ratio_cloud_levels.field[
                :
            ],
            "local_vapor_cloud_levels_forced": locals.vapor_cloud_levels_forced.field[:],
            "local_env_moist_static_energy_cloud_levels_forced": locals.environment_moist_static_energy_cloud_levels.field[
                :
            ],
            "u": state.input_output.u.field[:],
            "v": state.input_output.v.field[:],
            "local_u_cloud_levels": locals.u_cloud_levels.field[:],
            "local_v_cloud_levels": locals.v_cloud_levels.field[:],
            "local_env_saturation_moist_static_energy_cloud_levels_forced": locals.environment_saturation_moist_static_energy_cloud_levels_forced.field[
                :
            ],
            "local_geopotential_height_cloud_levels_forced": locals.geopotential_height_cloud_levels_forced.field[
                :
            ],
            "p_output": state.output.p.field[:, :, :, plume_dependent_constants.PLUME_INDEX],
            "local_gamma_cloud_levels_forced": locals.gamma_cloud_levels_forced.field[:],
            "local_t_cloud_levels_forced": locals.t_cloud_levels_forced.field[:],
            "p_surface": state.input_output.p_surface.field[:],
            "t_surface": state.input_output.t_surface.field[:],
            "error_code": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "topography_height_no_negative": state.input_output.topography_height_no_negative.field[:],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_EnvironmentCloudLevels_2_mid(TranslateFortranData2Py):
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
            "local_t_new": {"serialname": "local_t_new_env_clev"},
            "local_env_saturation_mixing_ratio_forced": {
                "serialname": "local_env_saturation_mixing_ratio_forced_env_clev"
            },
            "local_vapor_new": {"serialname": "local_vapor_new_env_cond"},
            "local_env_moist_static_energy_forced": {
                "serialname": "local_env_moist_static_energy_forced_env_clev"
            },
            "local_env_saturation_moist_static_energy_forced": {
                "serialname": "local_env_saturation_moist_static_energy_forced_env_clev"
            },
            "geopotential_height": {"serialname": "geopotential_height_env_clev"},
            "p": {"serialname": "p_env_clev"},
            "local_env_saturation_mixing_ratio_cloud_levels_forced": {
                "serialname": "local_env_saturation_mixing_ratio_cloud_levels_forced_env_clev"
            },
            "local_vapor_cloud_levels_forced": {"serialname": "local_vapor_cloud_levels_forced_env_clev"},
            "local_env_moist_static_energy_cloud_levels_forced": {
                "serialname": "local_env_moist_static_energy_cloud_levels_forced_env_clev"
            },
            "u": {"serialname": "u_env_clev"},
            "v": {"serialname": "v_env_clev"},
            "local_u_cloud_levels": {"serialname": "local_u_cloud_levels_env_clev"},
            "local_v_cloud_levels": {"serialname": "local_v_cloud_levels_env_clev"},
            "local_env_saturation_moist_static_energy_cloud_levels_forced": {
                "serialname": "local_env_saturation_moist_static_energy_cloud_levels_forced_env_clev"
            },
            "local_geopotential_height_cloud_levels_forced": {
                "serialname": "local_geopotential_height_cloud_levels_forced_env_clev"
            },
            "p_output": {"serialname": "p_output_env_clev"},
            "local_gamma_cloud_levels_forced": {"serialname": "local_gamma_cloud_levels_forced_env_clev"},
            "local_t_cloud_levels_forced": {"serialname": "local_t_cloud_levels_forced_env_clev"},
            "p_surface": {"serialname": "p_surface_env_clev"},
            "t_surface": {"serialname": "t_surface_env_clev"},
            "error_code": {"serialname": "error_code_env_clev"},
            "topography_height_no_negative": {"serialname": "topography_height_no_negative_env_clev"},
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
        locals.t_new.data[:] = inputs["local_t_new"]
        locals.environment_saturation_mixing_ratio_forced.data[:] = inputs[
            "local_env_saturation_mixing_ratio_forced"
        ]
        locals.vapor_new.data[:] = inputs["local_vapor_new"]
        locals.environment_moist_static_energy_forced.data[:] = inputs["local_env_moist_static_energy_forced"]
        locals.environment_saturation_moist_static_energy_forced.data[:] = inputs[
            "local_env_saturation_moist_static_energy_forced"
        ]
        state.input_output.geopotential_height.data[:] = inputs["geopotential_height"]
        state.input_output.p.data[:] = inputs["p"]
        locals.environment_saturation_mixing_ratio_cloud_levels.data[:] = inputs[
            "local_env_saturation_mixing_ratio_cloud_levels_forced"
        ]
        locals.vapor_cloud_levels_forced.data[:] = inputs["local_vapor_cloud_levels_forced"]
        locals.environment_moist_static_energy_cloud_levels.data[:] = inputs[
            "local_env_moist_static_energy_cloud_levels_forced"
        ]
        state.input_output.u.data[:] = inputs["u"]
        state.input_output.v.data[:] = inputs["v"]
        locals.u_cloud_levels.data[:] = inputs["local_u_cloud_levels"]
        locals.v_cloud_levels.data[:] = inputs["local_v_cloud_levels"]
        locals.environment_saturation_moist_static_energy_cloud_levels_forced.data[:] = inputs[
            "local_env_saturation_moist_static_energy_cloud_levels_forced"
        ]
        locals.geopotential_height_cloud_levels_forced.data[:] = inputs[
            "local_geopotential_height_cloud_levels_forced"
        ]
        state.output.p.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs["p_output"]
        locals.gamma_cloud_levels_forced.data[:] = inputs["local_gamma_cloud_levels_forced"]
        locals.t_cloud_levels_forced.data[:] = inputs["local_t_cloud_levels_forced"]
        state.input_output.p_surface.data[:] = inputs["p_surface"]
        state.input_output.t_surface.data[:] = inputs["t_surface"]
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["error_code"]
        state.input_output.topography_height_no_negative.data[:] = inputs["topography_height_no_negative"]

        environment_cloud_levels = EnvironmentCloudLevels(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        if plume_dependent_constants.ENABLE_PLUME == 1:
            environment_cloud_levels(
                state=state,
                locals=locals,
                saturation_tables=saturation_tables,
                plume_dependent_constants=plume_dependent_constants,
                data_type=1,
            )

        outputs = {
            # state fields
            "local_t_new": locals.t_new.field[:],
            "local_env_saturation_mixing_ratio_forced": locals.environment_saturation_mixing_ratio_forced.field[
                :
            ],
            "local_vapor_new": locals.vapor_new.field[:],
            "local_env_moist_static_energy_forced": locals.environment_moist_static_energy_forced.field[:],
            "local_env_saturation_moist_static_energy_forced": locals.environment_saturation_moist_static_energy_forced.field[
                :
            ],
            "geopotential_height": state.input_output.geopotential_height.field[:],
            "p": state.input_output.p.field[:],
            "local_env_saturation_mixing_ratio_cloud_levels_forced": locals.environment_saturation_mixing_ratio_cloud_levels.field[
                :
            ],
            "local_vapor_cloud_levels_forced": locals.vapor_cloud_levels_forced.field[:],
            "local_env_moist_static_energy_cloud_levels_forced": locals.environment_moist_static_energy_cloud_levels.field[
                :
            ],
            "u": state.input_output.u.field[:],
            "v": state.input_output.v.field[:],
            "local_u_cloud_levels": locals.u_cloud_levels.field[:],
            "local_v_cloud_levels": locals.v_cloud_levels.field[:],
            "local_env_saturation_moist_static_energy_cloud_levels_forced": locals.environment_saturation_moist_static_energy_cloud_levels_forced.field[
                :
            ],
            "local_geopotential_height_cloud_levels_forced": locals.geopotential_height_cloud_levels_forced.field[
                :
            ],
            "p_output": state.output.p.field[:, :, :, plume_dependent_constants.PLUME_INDEX],
            "local_gamma_cloud_levels_forced": locals.gamma_cloud_levels_forced.field[:],
            "local_t_cloud_levels_forced": locals.t_cloud_levels_forced.field[:],
            "p_surface": state.input_output.p_surface.field[:],
            "t_surface": state.input_output.t_surface.field[:],
            "error_code": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "topography_height_no_negative": state.input_output.topography_height_no_negative.field[:],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_EnvironmentCloudLevels_2_deep(TranslateFortranData2Py):
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
            "local_t_new": {"serialname": "local_t_new_env_clev"},
            "local_env_saturation_mixing_ratio_forced": {
                "serialname": "local_env_saturation_mixing_ratio_forced_env_clev"
            },
            "local_vapor_new": {"serialname": "local_vapor_new_env_cond"},
            "local_env_moist_static_energy_forced": {
                "serialname": "local_env_moist_static_energy_forced_env_clev"
            },
            "local_env_saturation_moist_static_energy_forced": {
                "serialname": "local_env_saturation_moist_static_energy_forced_env_clev"
            },
            "geopotential_height": {"serialname": "geopotential_height_env_clev"},
            "p": {"serialname": "p_env_clev"},
            "local_env_saturation_mixing_ratio_cloud_levels_forced": {
                "serialname": "local_env_saturation_mixing_ratio_cloud_levels_forced_env_clev"
            },
            "local_vapor_cloud_levels_forced": {"serialname": "local_vapor_cloud_levels_forced_env_clev"},
            "local_env_moist_static_energy_cloud_levels_forced": {
                "serialname": "local_env_moist_static_energy_cloud_levels_forced_env_clev"
            },
            "u": {"serialname": "u_env_clev"},
            "v": {"serialname": "v_env_clev"},
            "local_u_cloud_levels": {"serialname": "local_u_cloud_levels_env_clev"},
            "local_v_cloud_levels": {"serialname": "local_v_cloud_levels_env_clev"},
            "local_env_saturation_moist_static_energy_cloud_levels_forced": {
                "serialname": "local_env_saturation_moist_static_energy_cloud_levels_forced_env_clev"
            },
            "local_geopotential_height_cloud_levels_forced": {
                "serialname": "local_geopotential_height_cloud_levels_forced_env_clev"
            },
            "p_output": {"serialname": "p_output_env_clev"},
            "local_gamma_cloud_levels_forced": {"serialname": "local_gamma_cloud_levels_forced_env_clev"},
            "local_t_cloud_levels_forced": {"serialname": "local_t_cloud_levels_forced_env_clev"},
            "p_surface": {"serialname": "p_surface_env_clev"},
            "t_surface": {"serialname": "t_surface_env_clev"},
            "error_code": {"serialname": "error_code_env_clev"},
            "topography_height_no_negative": {"serialname": "topography_height_no_negative_env_clev"},
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
        locals.t_new.data[:] = inputs["local_t_new"]
        locals.environment_saturation_mixing_ratio_forced.data[:] = inputs[
            "local_env_saturation_mixing_ratio_forced"
        ]
        locals.vapor_new.data[:] = inputs["local_vapor_new"]
        locals.environment_moist_static_energy_forced.data[:] = inputs["local_env_moist_static_energy_forced"]
        locals.environment_saturation_moist_static_energy_forced.data[:] = inputs[
            "local_env_saturation_moist_static_energy_forced"
        ]
        state.input_output.geopotential_height.data[:] = inputs["geopotential_height"]
        state.input_output.p.data[:] = inputs["p"]
        locals.environment_saturation_mixing_ratio_cloud_levels.data[:] = inputs[
            "local_env_saturation_mixing_ratio_cloud_levels_forced"
        ]
        locals.vapor_cloud_levels_forced.data[:] = inputs["local_vapor_cloud_levels_forced"]
        locals.environment_moist_static_energy_cloud_levels.data[:] = inputs[
            "local_env_moist_static_energy_cloud_levels_forced"
        ]
        state.input_output.u.data[:] = inputs["u"]
        state.input_output.v.data[:] = inputs["v"]
        locals.u_cloud_levels.data[:] = inputs["local_u_cloud_levels"]
        locals.v_cloud_levels.data[:] = inputs["local_v_cloud_levels"]
        locals.environment_saturation_moist_static_energy_cloud_levels_forced.data[:] = inputs[
            "local_env_saturation_moist_static_energy_cloud_levels_forced"
        ]
        locals.geopotential_height_cloud_levels_forced.data[:] = inputs[
            "local_geopotential_height_cloud_levels_forced"
        ]
        state.output.p.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs["p_output"]
        locals.gamma_cloud_levels_forced.data[:] = inputs["local_gamma_cloud_levels_forced"]
        locals.t_cloud_levels_forced.data[:] = inputs["local_t_cloud_levels_forced"]
        state.input_output.p_surface.data[:] = inputs["p_surface"]
        state.input_output.t_surface.data[:] = inputs["t_surface"]
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["error_code"]
        state.input_output.topography_height_no_negative.data[:] = inputs["topography_height_no_negative"]

        environment_cloud_levels = EnvironmentCloudLevels(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        if plume_dependent_constants.ENABLE_PLUME == 1:
            environment_cloud_levels(
                state=state,
                locals=locals,
                saturation_tables=saturation_tables,
                plume_dependent_constants=plume_dependent_constants,
                data_type=1,
            )

        outputs = {
            # state fields
            "local_t_new": locals.t_new.field[:],
            "local_env_saturation_mixing_ratio_forced": locals.environment_saturation_mixing_ratio_forced.field[
                :
            ],
            "local_vapor_new": locals.vapor_new.field[:],
            "local_env_moist_static_energy_forced": locals.environment_moist_static_energy_forced.field[:],
            "local_env_saturation_moist_static_energy_forced": locals.environment_saturation_moist_static_energy_forced.field[
                :
            ],
            "geopotential_height": state.input_output.geopotential_height.field[:],
            "p": state.input_output.p.field[:],
            "local_env_saturation_mixing_ratio_cloud_levels_forced": locals.environment_saturation_mixing_ratio_cloud_levels.field[
                :
            ],
            "local_vapor_cloud_levels_forced": locals.vapor_cloud_levels_forced.field[:],
            "local_env_moist_static_energy_cloud_levels_forced": locals.environment_moist_static_energy_cloud_levels.field[
                :
            ],
            "u": state.input_output.u.field[:],
            "v": state.input_output.v.field[:],
            "local_u_cloud_levels": locals.u_cloud_levels.field[:],
            "local_v_cloud_levels": locals.v_cloud_levels.field[:],
            "local_env_saturation_moist_static_energy_cloud_levels_forced": locals.environment_saturation_moist_static_energy_cloud_levels_forced.field[
                :
            ],
            "local_geopotential_height_cloud_levels_forced": locals.geopotential_height_cloud_levels_forced.field[
                :
            ],
            "p_output": state.output.p.field[:, :, :, plume_dependent_constants.PLUME_INDEX],
            "local_gamma_cloud_levels_forced": locals.gamma_cloud_levels_forced.field[:],
            "local_t_cloud_levels_forced": locals.t_cloud_levels_forced.field[:],
            "p_surface": state.input_output.p_surface.field[:],
            "t_surface": state.input_output.t_surface.field[:],
            "error_code": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "topography_height_no_negative": state.input_output.topography_height_no_negative.field[:],
        }

        return outputs
