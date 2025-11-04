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
from pyMoist.convection.GF_2020.cumulus_parameterization.environment.environment import (
    EnvironmentConditions,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.setup.set_constants import set_constants


class TranslateGF2020_CumulusParameterization_EnvironmentConditions_1_shallow(TranslateFortranData2Py):
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
            "local_geopotential_height_env_cond": {},
            "local_env_saturation_mixing_ratio_env_cond": {},
            "local_env_moist_static_energy_env_cond": {},
            "local_env_saturation_moist_static_energy_env_cond": {},
            "t_old_env_cond": {},
            "vapor_old_env_cond": {},
            "p_forced_env_cond": {},
            "topography_height_no_negative_env_cond": {},
            "p_surface_env_cond": {},
            "error_code_env_cond": {},
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
        locals.geopotential_height.data[:] = inputs["local_geopotential_height_env_cond"]
        locals.environment_saturation_mixing_ratio.data[:] = inputs[
            "local_env_saturation_mixing_ratio_env_cond"
        ]
        locals.environment_moist_static_energy.data[:] = inputs["local_env_moist_static_energy_env_cond"]
        locals.environment_saturation_moist_static_energy.data[:] = inputs[
            "local_env_saturation_moist_static_energy_env_cond"
        ]
        state.input_output.t_old.data[:] = inputs["t_old_env_cond"]
        state.input_output.vapor_old.data[:] = inputs["vapor_old_env_cond"]
        state.input_output.p_forced.data[:] = inputs["p_forced_env_cond"]
        state.input_output.topography_height_no_negative.data[:] = inputs[
            "topography_height_no_negative_env_cond"
        ]
        state.input_output.p_surface.data[:] = inputs["p_surface_env_cond"]
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "error_code_env_cond"
        ]

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
                plume_dependent_constants=plume_dependent_constants,
                data_type=0,
            )

        outputs = {
            # state fields
            "local_geopotential_height_env_cond": locals.geopotential_height.field[:],
            "local_env_saturation_mixing_ratio_env_cond": locals.environment_saturation_mixing_ratio.field[:],
            "local_env_moist_static_energy_env_cond": locals.environment_moist_static_energy.field[:],
            "local_env_saturation_moist_static_energy_env_cond": locals.environment_saturation_moist_static_energy.field[
                :
            ],
            "t_old_env_cond": state.input_output.t_old.field[:],
            "vapor_old_env_cond": state.input_output.vapor_old.field[:],
            "p_forced_env_cond": state.input_output.p_forced.field[:],
            "topography_height_no_negative_env_cond": state.input_output.topography_height_no_negative.field[
                :
            ],
            "p_surface_env_cond": state.input_output.p_surface.field[:],
            "error_code_env_cond": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_EnvironmentConditions_1_mid(TranslateFortranData2Py):
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
            "local_geopotential_height_env_cond": {},
            "local_env_saturation_mixing_ratio_env_cond": {},
            "local_env_moist_static_energy_env_cond": {},
            "local_env_saturation_moist_static_energy_env_cond": {},
            "t_old_env_cond": {},
            "vapor_old_env_cond": {},
            "p_forced_env_cond": {},
            "topography_height_no_negative_env_cond": {},
            "p_surface_env_cond": {},
            "error_code_env_cond": {},
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
        locals.geopotential_height.data[:] = inputs["local_geopotential_height_env_cond"]
        locals.environment_saturation_mixing_ratio.data[:] = inputs[
            "local_env_saturation_mixing_ratio_env_cond"
        ]
        locals.environment_moist_static_energy.data[:] = inputs["local_env_moist_static_energy_env_cond"]
        locals.environment_saturation_moist_static_energy.data[:] = inputs[
            "local_env_saturation_moist_static_energy_env_cond"
        ]
        state.input_output.t_old.data[:] = inputs["t_old_env_cond"]
        state.input_output.vapor_old.data[:] = inputs["vapor_old_env_cond"]
        state.input_output.p_forced.data[:] = inputs["p_forced_env_cond"]
        state.input_output.topography_height_no_negative.data[:] = inputs[
            "topography_height_no_negative_env_cond"
        ]
        state.input_output.p_surface.data[:] = inputs["p_surface_env_cond"]
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "error_code_env_cond"
        ]

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
                plume_dependent_constants=plume_dependent_constants,
                data_type=0,
            )

        outputs = {
            # state fields
            "local_geopotential_height_env_cond": locals.geopotential_height.field[:],
            "local_env_saturation_mixing_ratio_env_cond": locals.environment_saturation_mixing_ratio.field[:],
            "local_env_moist_static_energy_env_cond": locals.environment_moist_static_energy.field[:],
            "local_env_saturation_moist_static_energy_env_cond": locals.environment_saturation_moist_static_energy.field[
                :
            ],
            "t_old_env_cond": state.input_output.t_old.field[:],
            "vapor_old_env_cond": state.input_output.vapor_old.field[:],
            "p_forced_env_cond": state.input_output.p_forced.field[:],
            "topography_height_no_negative_env_cond": state.input_output.topography_height_no_negative.field[
                :
            ],
            "p_surface_env_cond": state.input_output.p_surface.field[:],
            "error_code_env_cond": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_EnvironmentConditions_1_deep(TranslateFortranData2Py):
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
            "local_geopotential_height_env_cond": {},
            "local_env_saturation_mixing_ratio_env_cond": {},
            "local_env_moist_static_energy_env_cond": {},
            "local_env_saturation_moist_static_energy_env_cond": {},
            "t_old_env_cond": {},
            "vapor_old_env_cond": {},
            "p_forced_env_cond": {},
            "topography_height_no_negative_env_cond": {},
            "p_surface_env_cond": {},
            "error_code_env_cond": {},
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
        locals.geopotential_height.data[:] = inputs["local_geopotential_height_env_cond"]
        locals.environment_saturation_mixing_ratio.data[:] = inputs[
            "local_env_saturation_mixing_ratio_env_cond"
        ]
        locals.environment_moist_static_energy.data[:] = inputs["local_env_moist_static_energy_env_cond"]
        locals.environment_saturation_moist_static_energy.data[:] = inputs[
            "local_env_saturation_moist_static_energy_env_cond"
        ]
        state.input_output.t_old.data[:] = inputs["t_old_env_cond"]
        state.input_output.vapor_old.data[:] = inputs["vapor_old_env_cond"]
        state.input_output.p_forced.data[:] = inputs["p_forced_env_cond"]
        state.input_output.topography_height_no_negative.data[:] = inputs[
            "topography_height_no_negative_env_cond"
        ]
        state.input_output.p_surface.data[:] = inputs["p_surface_env_cond"]
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "error_code_env_cond"
        ]

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
                plume_dependent_constants=plume_dependent_constants,
                data_type=0,
            )

        outputs = {
            # state fields
            "local_geopotential_height_env_cond": locals.geopotential_height.field[:],
            "local_env_saturation_mixing_ratio_env_cond": locals.environment_saturation_mixing_ratio.field[:],
            "local_env_moist_static_energy_env_cond": locals.environment_moist_static_energy.field[:],
            "local_env_saturation_moist_static_energy_env_cond": locals.environment_saturation_moist_static_energy.field[
                :
            ],
            "t_old_env_cond": state.input_output.t_old.field[:],
            "vapor_old_env_cond": state.input_output.vapor_old.field[:],
            "p_forced_env_cond": state.input_output.p_forced.field[:],
            "topography_height_no_negative_env_cond": state.input_output.topography_height_no_negative.field[
                :
            ],
            "p_surface_env_cond": state.input_output.p_surface.field[:],
            "error_code_env_cond": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
        }

        return outputs
