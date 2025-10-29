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


class TranslateGF2020_CumulusParameterization_EnvironmentCloudLevels_1_shallow(TranslateFortranData2Py):
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
            "p": {"serialname": "p_env_clev"},
            "p_surface": {"serialname": "p_surface_env_clev"},
            "p_cloud_levels": {"serialname": "local_p_cloud_levels_env_clev"},
            "topography_height_no_negative": {"serialname": "topography_height_no_negative_env_clev"},
            "geopotential_height": {"serialname": "local_geopotential_height_env_clev"},
            "geopotential_height_cloud_levels": {
                "serialname": "local_geopotential_height_cloud_levels_env_clev"
            },
            "t_old": {"serialname": "t_old_env_clev"},
            "t_surface": {"serialname": "t_surface_env_clev"},
            "t_cloud_levels": {"serialname": "local_t_cloud_levels_env_clev"},
            "vapor_old": {"serialname": "vapor_old_env_clev"},
            "vapor_cloud_levels": {"serialname": "local_vapor_cloud_levels_env_clev"},
            "u": {"serialname": "u_env_clev"},
            "v": {"serialname": "v_env_clev"},
            "u_cloud_levels": {"serialname": "local_u_cloud_levels_env_clev"},
            "v_cloud_levels": {"serialname": "local_v_cloud_levels_env_clev"},
            "environment_saturation_mixing_ratio": {
                "serialname": "local_env_saturation_mixing_ratio_env_clev"
            },
            "environment_saturation_mixing_ratio_cloud_levels": {
                "serialname": "local_env_saturation_mixing_ratio_cloud_levels_env_clev"
            },
            "environment_moist_static_energy": {"serialname": "local_env_moist_static_energy_env_clev"},
            "environment_moist_static_energy_cloud_levels": {
                "serialname": "local_env_moist_static_energy_cloud_levels_env_clev"
            },
            "environment_saturation_moist_static_energy": {
                "serialname": "local_env_saturation_moist_static_energy_env_clev"
            },
            "environment_saturation_moist_static_energy_cloud_levels": {
                "serialname": "local_env_saturation_moist_static_energy_cloud_levels_env_clev"
            },
            "gamma_cloud_levels": {"serialname": "local_gamma_cloud_levels_env_clev"},
            "error_code": {"serialname": "error_code_env_clev"},
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
        state.input_output.p.data[:] = inputs["p"]
        state.input_output.p_surface.data[:] = inputs["p_surface"]
        locals.p_cloud_levels.data[:] = inputs["p_cloud_levels"]
        state.input_output.topography_height_no_negative.data[:] = inputs["topography_height_no_negative"]
        locals.geopotential_height.data[:] = inputs["geopotential_height"]
        locals.geopotential_height_cloud_levels.data[:] = inputs["geopotential_height_cloud_levels"]
        state.input_output.t_old.data[:] = inputs["t_old"]
        state.input_output.t_surface.data[:] = inputs["t_surface"]
        locals.t_cloud_levels.data[:] = inputs["t_cloud_levels"]
        state.input_output.vapor_old.data[:] = inputs["vapor_old"]
        locals.vapor_cloud_levels.data[:] = inputs["vapor_cloud_levels"]
        state.input_output.u.data[:] = inputs["u"]
        state.input_output.v.data[:] = inputs["v"]
        locals.u_cloud_levels.data[:] = inputs["u_cloud_levels"]
        locals.v_cloud_levels.data[:] = inputs["v_cloud_levels"]
        locals.environment_saturation_mixing_ratio.data[:] = inputs["environment_saturation_mixing_ratio"]
        locals.environment_saturation_mixing_ratio_cloud_levels.data[:] = inputs[
            "environment_saturation_mixing_ratio_cloud_levels"
        ]
        locals.environment_moist_static_energy.data[:] = inputs["environment_moist_static_energy"]
        locals.environment_moist_static_energy_cloud_levels.data[:] = inputs[
            "environment_moist_static_energy_cloud_levels"
        ]
        locals.environment_saturation_moist_static_energy.data[:] = inputs[
            "environment_saturation_moist_static_energy"
        ]
        locals.environment_saturation_moist_static_energy_cloud_levels.data[:] = inputs[
            "environment_saturation_moist_static_energy_cloud_levels"
        ]
        locals.gamma_cloud_levels.data[:] = inputs["gamma_cloud_levels"]
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["error_code"]

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
                data_type=0,
            )

        outputs = {
            # state fields
            "p": state.input_output.p.field[:],
            "p_surface": state.input_output.p_surface.field[:],
            "p_cloud_levels": locals.p_cloud_levels.field[:],
            "topography_height_no_negative": state.input_output.topography_height_no_negative.field[:],
            "geopotential_height": locals.geopotential_height.field[:],
            "geopotential_height_cloud_levels": locals.geopotential_height_cloud_levels.field[:],
            "t_old": state.input_output.t_old.field[:],
            "t_surface": state.input_output.t_surface.field[:],
            "t_cloud_levels": locals.t_cloud_levels.field[:],
            "vapor_old": state.input_output.vapor_old.field[:],
            "vapor_cloud_levels": locals.vapor_cloud_levels.field[:],
            "u": state.input_output.u.field[:],
            "v": state.input_output.v.field[:],
            "u_cloud_levels": locals.u_cloud_levels.field[:],
            "v_cloud_levels": locals.v_cloud_levels.field[:],
            "environment_saturation_mixing_ratio": locals.environment_saturation_mixing_ratio.field[:],
            "environment_saturation_mixing_ratio_cloud_levels": locals.environment_saturation_mixing_ratio_cloud_levels.field[
                :
            ],
            "environment_moist_static_energy": locals.environment_moist_static_energy.field[:],
            "environment_moist_static_energy_cloud_levels": locals.environment_moist_static_energy_cloud_levels.field[
                :
            ],
            "environment_saturation_moist_static_energy": locals.environment_saturation_moist_static_energy.field[
                :
            ],
            "environment_saturation_moist_static_energy_cloud_levels": locals.environment_saturation_moist_static_energy_cloud_levels.field[
                :
            ],
            "gamma_cloud_levels": locals.gamma_cloud_levels.field[:],
            "error_code": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_EnvironmentCloudLevels_1_mid(TranslateFortranData2Py):
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
            "p": {"serialname": "p_env_clev"},
            "p_surface": {"serialname": "p_surface_env_clev"},
            "p_cloud_levels": {"serialname": "local_p_cloud_levels_env_clev"},
            "topography_height_no_negative": {"serialname": "topography_height_no_negative_env_clev"},
            "geopotential_height": {"serialname": "local_geopotential_height_env_clev"},
            "geopotential_height_cloud_levels": {
                "serialname": "local_geopotential_height_cloud_levels_env_clev"
            },
            "t_old": {"serialname": "t_old_env_clev"},
            "t_surface": {"serialname": "t_surface_env_clev"},
            "t_cloud_levels": {"serialname": "local_t_cloud_levels_env_clev"},
            "vapor_old": {"serialname": "vapor_old_env_clev"},
            "vapor_cloud_levels": {"serialname": "local_vapor_cloud_levels_env_clev"},
            "u": {"serialname": "u_env_clev"},
            "v": {"serialname": "v_env_clev"},
            "u_cloud_levels": {"serialname": "local_u_cloud_levels_env_clev"},
            "v_cloud_levels": {"serialname": "local_v_cloud_levels_env_clev"},
            "environment_saturation_mixing_ratio": {
                "serialname": "local_env_saturation_mixing_ratio_env_clev"
            },
            "environment_saturation_mixing_ratio_cloud_levels": {
                "serialname": "local_env_saturation_mixing_ratio_cloud_levels_env_clev"
            },
            "environment_moist_static_energy": {"serialname": "local_env_moist_static_energy_env_clev"},
            "environment_moist_static_energy_cloud_levels": {
                "serialname": "local_env_moist_static_energy_cloud_levels_env_clev"
            },
            "environment_saturation_moist_static_energy": {
                "serialname": "local_env_saturation_moist_static_energy_env_clev"
            },
            "environment_saturation_moist_static_energy_cloud_levels": {
                "serialname": "local_env_saturation_moist_static_energy_cloud_levels_env_clev"
            },
            "gamma_cloud_levels": {"serialname": "local_gamma_cloud_levels_env_clev"},
            "error_code": {"serialname": "error_code_env_clev"},
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
        state.input_output.p.data[:] = inputs["p"]
        state.input_output.p_surface.data[:] = inputs["p_surface"]
        locals.p_cloud_levels.data[:] = inputs["p_cloud_levels"]
        state.input_output.topography_height_no_negative.data[:] = inputs["topography_height_no_negative"]
        locals.geopotential_height.data[:] = inputs["geopotential_height"]
        locals.geopotential_height_cloud_levels.data[:] = inputs["geopotential_height_cloud_levels"]
        state.input_output.t_old.data[:] = inputs["t_old"]
        state.input_output.t_surface.data[:] = inputs["t_surface"]
        locals.t_cloud_levels.data[:] = inputs["t_cloud_levels"]
        state.input_output.vapor_old.data[:] = inputs["vapor_old"]
        locals.vapor_cloud_levels.data[:] = inputs["vapor_cloud_levels"]
        state.input_output.u.data[:] = inputs["u"]
        state.input_output.v.data[:] = inputs["v"]
        locals.u_cloud_levels.data[:] = inputs["u_cloud_levels"]
        locals.v_cloud_levels.data[:] = inputs["v_cloud_levels"]
        locals.environment_saturation_mixing_ratio.data[:] = inputs["environment_saturation_mixing_ratio"]
        locals.environment_saturation_mixing_ratio_cloud_levels.data[:] = inputs[
            "environment_saturation_mixing_ratio_cloud_levels"
        ]
        locals.environment_moist_static_energy.data[:] = inputs["environment_moist_static_energy"]
        locals.environment_moist_static_energy_cloud_levels.data[:] = inputs[
            "environment_moist_static_energy_cloud_levels"
        ]
        locals.environment_saturation_moist_static_energy.data[:] = inputs[
            "environment_saturation_moist_static_energy"
        ]
        locals.environment_saturation_moist_static_energy_cloud_levels.data[:] = inputs[
            "environment_saturation_moist_static_energy_cloud_levels"
        ]
        locals.gamma_cloud_levels.data[:] = inputs["gamma_cloud_levels"]
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["error_code"]

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
                data_type=0,
            )

        outputs = {
            # state fields
            "p": state.input_output.p.field[:],
            "p_surface": state.input_output.p_surface.field[:],
            "p_cloud_levels": locals.p_cloud_levels.field[:],
            "topography_height_no_negative": state.input_output.topography_height_no_negative.field[:],
            "geopotential_height": locals.geopotential_height.field[:],
            "geopotential_height_cloud_levels": locals.geopotential_height_cloud_levels.field[:],
            "t_old": state.input_output.t_old.field[:],
            "t_surface": state.input_output.t_surface.field[:],
            "t_cloud_levels": locals.t_cloud_levels.field[:],
            "vapor_old": state.input_output.vapor_old.field[:],
            "vapor_cloud_levels": locals.vapor_cloud_levels.field[:],
            "u": state.input_output.u.field[:],
            "v": state.input_output.v.field[:],
            "u_cloud_levels": locals.u_cloud_levels.field[:],
            "v_cloud_levels": locals.v_cloud_levels.field[:],
            "environment_saturation_mixing_ratio": locals.environment_saturation_mixing_ratio.field[:],
            "environment_saturation_mixing_ratio_cloud_levels": locals.environment_saturation_mixing_ratio_cloud_levels.field[
                :
            ],
            "environment_moist_static_energy": locals.environment_moist_static_energy.field[:],
            "environment_moist_static_energy_cloud_levels": locals.environment_moist_static_energy_cloud_levels.field[
                :
            ],
            "environment_saturation_moist_static_energy": locals.environment_saturation_moist_static_energy.field[
                :
            ],
            "environment_saturation_moist_static_energy_cloud_levels": locals.environment_saturation_moist_static_energy_cloud_levels.field[
                :
            ],
            "gamma_cloud_levels": locals.gamma_cloud_levels.field[:],
            "error_code": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_EnvironmentCloudLevels_1_deep(TranslateFortranData2Py):
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
            "p": {"serialname": "p_env_clev"},
            "p_surface": {"serialname": "p_surface_env_clev"},
            "p_cloud_levels": {"serialname": "local_p_cloud_levels_env_clev"},
            "topography_height_no_negative": {"serialname": "topography_height_no_negative_env_clev"},
            "geopotential_height": {"serialname": "local_geopotential_height_env_clev"},
            "geopotential_height_cloud_levels": {
                "serialname": "local_geopotential_height_cloud_levels_env_clev"
            },
            "t_old": {"serialname": "t_old_env_clev"},
            "t_surface": {"serialname": "t_surface_env_clev"},
            "t_cloud_levels": {"serialname": "local_t_cloud_levels_env_clev"},
            "vapor_old": {"serialname": "vapor_old_env_clev"},
            "vapor_cloud_levels": {"serialname": "local_vapor_cloud_levels_env_clev"},
            "u": {"serialname": "u_env_clev"},
            "v": {"serialname": "v_env_clev"},
            "u_cloud_levels": {"serialname": "local_u_cloud_levels_env_clev"},
            "v_cloud_levels": {"serialname": "local_v_cloud_levels_env_clev"},
            "environment_saturation_mixing_ratio": {
                "serialname": "local_env_saturation_mixing_ratio_env_clev"
            },
            "environment_saturation_mixing_ratio_cloud_levels": {
                "serialname": "local_env_saturation_mixing_ratio_cloud_levels_env_clev"
            },
            "environment_moist_static_energy": {"serialname": "local_env_moist_static_energy_env_clev"},
            "environment_moist_static_energy_cloud_levels": {
                "serialname": "local_env_moist_static_energy_cloud_levels_env_clev"
            },
            "environment_saturation_moist_static_energy": {
                "serialname": "local_env_saturation_moist_static_energy_env_clev"
            },
            "environment_saturation_moist_static_energy_cloud_levels": {
                "serialname": "local_env_saturation_moist_static_energy_cloud_levels_env_clev"
            },
            "gamma_cloud_levels": {"serialname": "local_gamma_cloud_levels_env_clev"},
            "error_code": {"serialname": "error_code_env_clev"},
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
        state.input_output.p.data[:] = inputs["p"]
        state.input_output.p_surface.data[:] = inputs["p_surface"]
        locals.p_cloud_levels.data[:] = inputs["p_cloud_levels"]
        state.input_output.topography_height_no_negative.data[:] = inputs["topography_height_no_negative"]
        locals.geopotential_height.data[:] = inputs["geopotential_height"]
        locals.geopotential_height_cloud_levels.data[:] = inputs["geopotential_height_cloud_levels"]
        state.input_output.t_old.data[:] = inputs["t_old"]
        state.input_output.t_surface.data[:] = inputs["t_surface"]
        locals.t_cloud_levels.data[:] = inputs["t_cloud_levels"]
        state.input_output.vapor_old.data[:] = inputs["vapor_old"]
        locals.vapor_cloud_levels.data[:] = inputs["vapor_cloud_levels"]
        state.input_output.u.data[:] = inputs["u"]
        state.input_output.v.data[:] = inputs["v"]
        locals.u_cloud_levels.data[:] = inputs["u_cloud_levels"]
        locals.v_cloud_levels.data[:] = inputs["v_cloud_levels"]
        locals.environment_saturation_mixing_ratio.data[:] = inputs["environment_saturation_mixing_ratio"]
        locals.environment_saturation_mixing_ratio_cloud_levels.data[:] = inputs[
            "environment_saturation_mixing_ratio_cloud_levels"
        ]
        locals.environment_moist_static_energy.data[:] = inputs["environment_moist_static_energy"]
        locals.environment_moist_static_energy_cloud_levels.data[:] = inputs[
            "environment_moist_static_energy_cloud_levels"
        ]
        locals.environment_saturation_moist_static_energy.data[:] = inputs[
            "environment_saturation_moist_static_energy"
        ]
        locals.environment_saturation_moist_static_energy_cloud_levels.data[:] = inputs[
            "environment_saturation_moist_static_energy_cloud_levels"
        ]
        locals.gamma_cloud_levels.data[:] = inputs["gamma_cloud_levels"]
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["error_code"]

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
                data_type=0,
            )

        outputs = {
            # state fields
            "p": state.input_output.p.field[:],
            "p_surface": state.input_output.p_surface.field[:],
            "p_cloud_levels": locals.p_cloud_levels.field[:],
            "topography_height_no_negative": state.input_output.topography_height_no_negative.field[:],
            "geopotential_height": locals.geopotential_height.field[:],
            "geopotential_height_cloud_levels": locals.geopotential_height_cloud_levels.field[:],
            "t_old": state.input_output.t_old.field[:],
            "t_surface": state.input_output.t_surface.field[:],
            "t_cloud_levels": locals.t_cloud_levels.field[:],
            "vapor_old": state.input_output.vapor_old.field[:],
            "vapor_cloud_levels": locals.vapor_cloud_levels.field[:],
            "u": state.input_output.u.field[:],
            "v": state.input_output.v.field[:],
            "u_cloud_levels": locals.u_cloud_levels.field[:],
            "v_cloud_levels": locals.v_cloud_levels.field[:],
            "environment_saturation_mixing_ratio": locals.environment_saturation_mixing_ratio.field[:],
            "environment_saturation_mixing_ratio_cloud_levels": locals.environment_saturation_mixing_ratio_cloud_levels.field[
                :
            ],
            "environment_moist_static_energy": locals.environment_moist_static_energy.field[:],
            "environment_moist_static_energy_cloud_levels": locals.environment_moist_static_energy_cloud_levels.field[
                :
            ],
            "environment_saturation_moist_static_energy": locals.environment_saturation_moist_static_energy.field[
                :
            ],
            "environment_saturation_moist_static_energy_cloud_levels": locals.environment_saturation_moist_static_energy_cloud_levels.field[
                :
            ],
            "gamma_cloud_levels": locals.gamma_cloud_levels.field[:],
            "error_code": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
        }

        return outputs
