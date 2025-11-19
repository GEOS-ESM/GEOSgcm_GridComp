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
from pyMoist.convection.GF_2020.cumulus_parameterization.environment.environment import EnvironmentCloudLevels
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
            "local_t_new_env_clev": {},
            "local_env_saturation_mixing_ratio_forced_env_clev": {},
            "local_vapor_forced_env_cond": {},
            "local_env_moist_static_energy_forced_env_clev": {},
            "local_env_saturation_moist_static_energy_forced_env_clev": {},
            "geopotential_height_forced_env_clev": {},
            "p_forced_env_clev": {},
            "local_env_saturation_mixing_ratio_cloud_levels_forced_env_clev": {},
            "local_vapor_cloud_levels_forced_env_clev": {},
            "local_env_moist_static_energy_cloud_levels_forced_env_clev": {},
            "u_env_clev": {},
            "v_env_clev": {},
            "local_u_cloud_levels_env_clev": {},
            "local_v_cloud_levels_env_clev": {},
            "local_env_saturation_moist_static_energy_cloud_levels_forced_env_clev": {},
            "local_geopotential_height_cloud_levels_forced_env_clev": {},
            "p_cloud_levels_forced_env_clev": {},
            "local_gamma_cloud_levels_forced_env_clev": {},
            "local_t_cloud_levels_forced_env_clev": {},
            "p_surface_env_clev": {},
            "t_surface_env_clev": {},
            "error_code_env_clev": {},
            "topography_height_no_negative_env_clev": {},
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
        locals.t_new.data[:] = inputs["local_t_new_env_clev"]
        locals.environment_saturation_mixing_ratio_forced.data[:] = inputs[
            "local_env_saturation_mixing_ratio_forced_env_clev"
        ]
        locals.vapor_forced.data[:] = inputs["local_vapor_forced_env_cond"]
        locals.environment_moist_static_energy_forced.data[:] = inputs[
            "local_env_moist_static_energy_forced_env_clev"
        ]
        locals.environment_saturation_moist_static_energy_forced.data[:] = inputs[
            "local_env_saturation_moist_static_energy_forced_env_clev"
        ]
        state.input_output.geopotential_height_forced.data[:] = inputs["geopotential_height_forced_env_clev"]
        state.input_output.p_forced.data[:] = inputs["p_forced_env_clev"]
        locals.environment_saturation_mixing_ratio_cloud_levels_forced.data[:] = inputs[
            "local_env_saturation_mixing_ratio_cloud_levels_forced_env_clev"
        ]
        locals.vapor_cloud_levels_forced.data[:] = inputs["local_vapor_cloud_levels_forced_env_clev"]
        locals.environment_moist_static_energy_cloud_levels_forced.data[:] = inputs[
            "local_env_moist_static_energy_cloud_levels_forced_env_clev"
        ]
        state.input_output.u.data[:] = inputs["u_env_clev"]
        state.input_output.v.data[:] = inputs["v_env_clev"]
        locals.u_cloud_levels.data[:] = inputs["local_u_cloud_levels_env_clev"]
        locals.v_cloud_levels.data[:] = inputs["local_v_cloud_levels_env_clev"]
        locals.environment_saturation_moist_static_energy_cloud_levels_forced.data[:] = inputs[
            "local_env_saturation_moist_static_energy_cloud_levels_forced_env_clev"
        ]
        locals.geopotential_height_cloud_levels_forced.data[:] = inputs[
            "local_geopotential_height_cloud_levels_forced_env_clev"
        ]
        state.output.p_cloud_levels_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "p_cloud_levels_forced_env_clev"
        ]
        locals.gamma_cloud_levels_forced.data[:] = inputs["local_gamma_cloud_levels_forced_env_clev"]
        locals.t_cloud_levels_forced.data[:] = inputs["local_t_cloud_levels_forced_env_clev"]
        state.input_output.p_surface.data[:] = inputs["p_surface_env_clev"]
        state.input_output.t_surface.data[:] = inputs["t_surface_env_clev"]
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "error_code_env_clev"
        ]
        state.input_output.topography_height_no_negative.data[:] = inputs[
            "topography_height_no_negative_env_clev"
        ]

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
                plume_dependent_constants=plume_dependent_constants,
                data_type=1,
            )

        outputs = {
            "local_t_new_env_clev": locals.t_new.field[:],
            "local_env_saturation_mixing_ratio_forced_env_clev": locals.environment_saturation_mixing_ratio_forced.field[
                :
            ],
            "local_vapor_forced_env_cond": locals.vapor_forced.field[:],
            "local_env_moist_static_energy_forced_env_clev": locals.environment_moist_static_energy_forced.field[
                :
            ],
            "local_env_saturation_moist_static_energy_forced_env_clev": locals.environment_saturation_moist_static_energy_forced.field[
                :
            ],
            "geopotential_height_forced_env_clev": state.input_output.geopotential_height_forced.field[:],
            "p_forced_env_clev": state.input_output.p_forced.field[:],
            "local_env_saturation_mixing_ratio_cloud_levels_forced_env_clev": locals.environment_saturation_mixing_ratio_cloud_levels_forced.field[
                :
            ],
            "local_vapor_cloud_levels_forced_env_clev": locals.vapor_cloud_levels_forced.field[:],
            "local_env_moist_static_energy_cloud_levels_forced_env_clev": locals.environment_moist_static_energy_cloud_levels_forced.field[
                :
            ],
            "u_env_clev": state.input_output.u.field[:],
            "v_env_clev": state.input_output.v.field[:],
            "local_u_cloud_levels_env_clev": locals.u_cloud_levels.field[:],
            "local_v_cloud_levels_env_clev": locals.v_cloud_levels.field[:],
            "local_env_saturation_moist_static_energy_cloud_levels_forced_env_clev": locals.environment_saturation_moist_static_energy_cloud_levels_forced.field[
                :
            ],
            "local_geopotential_height_cloud_levels_forced_env_clev": locals.geopotential_height_cloud_levels_forced.field[
                :
            ],
            "p_cloud_levels_forced_env_clev": state.output.p_cloud_levels_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_gamma_cloud_levels_forced_env_clev": locals.gamma_cloud_levels_forced.field[:],
            "local_t_cloud_levels_forced_env_clev": locals.t_cloud_levels_forced.field[:],
            "p_surface_env_clev": state.input_output.p_surface.field[:],
            "t_surface_env_clev": state.input_output.t_surface.field[:],
            "error_code_env_clev": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "topography_height_no_negative_env_clev": state.input_output.topography_height_no_negative.field[
                :
            ],
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
            "local_t_new_env_clev": {},
            "local_env_saturation_mixing_ratio_forced_env_clev": {},
            "local_vapor_forced_env_cond": {},
            "local_env_moist_static_energy_forced_env_clev": {},
            "local_env_saturation_moist_static_energy_forced_env_clev": {},
            "geopotential_height_forced_env_clev": {},
            "p_forced_env_clev": {},
            "local_env_saturation_mixing_ratio_cloud_levels_forced_env_clev": {},
            "local_vapor_cloud_levels_forced_env_clev": {},
            "local_env_moist_static_energy_cloud_levels_forced_env_clev": {},
            "u_env_clev": {},
            "v_env_clev": {},
            "local_u_cloud_levels_env_clev": {},
            "local_v_cloud_levels_env_clev": {},
            "local_env_saturation_moist_static_energy_cloud_levels_forced_env_clev": {},
            "local_geopotential_height_cloud_levels_forced_env_clev": {},
            "p_cloud_levels_forced_env_clev": {},
            "local_gamma_cloud_levels_forced_env_clev": {},
            "local_t_cloud_levels_forced_env_clev": {},
            "p_surface_env_clev": {},
            "t_surface_env_clev": {},
            "error_code_env_clev": {},
            "topography_height_no_negative_env_clev": {},
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
        locals.t_new.data[:] = inputs["local_t_new_env_clev"]
        locals.environment_saturation_mixing_ratio_forced.data[:] = inputs[
            "local_env_saturation_mixing_ratio_forced_env_clev"
        ]
        locals.vapor_forced.data[:] = inputs["local_vapor_forced_env_cond"]
        locals.environment_moist_static_energy_forced.data[:] = inputs[
            "local_env_moist_static_energy_forced_env_clev"
        ]
        locals.environment_saturation_moist_static_energy_forced.data[:] = inputs[
            "local_env_saturation_moist_static_energy_forced_env_clev"
        ]
        state.input_output.geopotential_height_forced.data[:] = inputs["geopotential_height_forced_env_clev"]
        state.input_output.p_forced.data[:] = inputs["p_forced_env_clev"]
        locals.environment_saturation_mixing_ratio_cloud_levels_forced.data[:] = inputs[
            "local_env_saturation_mixing_ratio_cloud_levels_forced_env_clev"
        ]
        locals.vapor_cloud_levels_forced.data[:] = inputs["local_vapor_cloud_levels_forced_env_clev"]
        locals.environment_moist_static_energy_cloud_levels_forced.data[:] = inputs[
            "local_env_moist_static_energy_cloud_levels_forced_env_clev"
        ]
        state.input_output.u.data[:] = inputs["u_env_clev"]
        state.input_output.v.data[:] = inputs["v_env_clev"]
        locals.u_cloud_levels.data[:] = inputs["local_u_cloud_levels_env_clev"]
        locals.v_cloud_levels.data[:] = inputs["local_v_cloud_levels_env_clev"]
        locals.environment_saturation_moist_static_energy_cloud_levels_forced.data[:] = inputs[
            "local_env_saturation_moist_static_energy_cloud_levels_forced_env_clev"
        ]
        locals.geopotential_height_cloud_levels_forced.data[:] = inputs[
            "local_geopotential_height_cloud_levels_forced_env_clev"
        ]
        state.output.p_cloud_levels_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "p_cloud_levels_forced_env_clev"
        ]
        locals.gamma_cloud_levels_forced.data[:] = inputs["local_gamma_cloud_levels_forced_env_clev"]
        locals.t_cloud_levels_forced.data[:] = inputs["local_t_cloud_levels_forced_env_clev"]
        state.input_output.p_surface.data[:] = inputs["p_surface_env_clev"]
        state.input_output.t_surface.data[:] = inputs["t_surface_env_clev"]
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "error_code_env_clev"
        ]
        state.input_output.topography_height_no_negative.data[:] = inputs[
            "topography_height_no_negative_env_clev"
        ]

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
                plume_dependent_constants=plume_dependent_constants,
                data_type=1,
            )

        outputs = {
            "local_t_new_env_clev": locals.t_new.field[:],
            "local_env_saturation_mixing_ratio_forced_env_clev": locals.environment_saturation_mixing_ratio_forced.field[
                :
            ],
            "local_vapor_forced_env_cond": locals.vapor_forced.field[:],
            "local_env_moist_static_energy_forced_env_clev": locals.environment_moist_static_energy_forced.field[
                :
            ],
            "local_env_saturation_moist_static_energy_forced_env_clev": locals.environment_saturation_moist_static_energy_forced.field[
                :
            ],
            "geopotential_height_forced_env_clev": state.input_output.geopotential_height_forced.field[:],
            "p_forced_env_clev": state.input_output.p_forced.field[:],
            "local_env_saturation_mixing_ratio_cloud_levels_forced_env_clev": locals.environment_saturation_mixing_ratio_cloud_levels_forced.field[
                :
            ],
            "local_vapor_cloud_levels_forced_env_clev": locals.vapor_cloud_levels_forced.field[:],
            "local_env_moist_static_energy_cloud_levels_forced_env_clev": locals.environment_moist_static_energy_cloud_levels_forced.field[
                :
            ],
            "u_env_clev": state.input_output.u.field[:],
            "v_env_clev": state.input_output.v.field[:],
            "local_u_cloud_levels_env_clev": locals.u_cloud_levels.field[:],
            "local_v_cloud_levels_env_clev": locals.v_cloud_levels.field[:],
            "local_env_saturation_moist_static_energy_cloud_levels_forced_env_clev": locals.environment_saturation_moist_static_energy_cloud_levels_forced.field[
                :
            ],
            "local_geopotential_height_cloud_levels_forced_env_clev": locals.geopotential_height_cloud_levels_forced.field[
                :
            ],
            "p_cloud_levels_forced_env_clev": state.output.p_cloud_levels_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_gamma_cloud_levels_forced_env_clev": locals.gamma_cloud_levels_forced.field[:],
            "local_t_cloud_levels_forced_env_clev": locals.t_cloud_levels_forced.field[:],
            "p_surface_env_clev": state.input_output.p_surface.field[:],
            "t_surface_env_clev": state.input_output.t_surface.field[:],
            "error_code_env_clev": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "topography_height_no_negative_env_clev": state.input_output.topography_height_no_negative.field[
                :
            ],
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
            "local_t_new_env_clev": {},
            "local_env_saturation_mixing_ratio_forced_env_clev": {},
            "local_vapor_forced_env_cond": {},
            "local_env_moist_static_energy_forced_env_clev": {},
            "local_env_saturation_moist_static_energy_forced_env_clev": {},
            "geopotential_height_forced_env_clev": {},
            "p_forced_env_clev": {},
            "local_env_saturation_mixing_ratio_cloud_levels_forced_env_clev": {},
            "local_vapor_cloud_levels_forced_env_clev": {},
            "local_env_moist_static_energy_cloud_levels_forced_env_clev": {},
            "u_env_clev": {},
            "v_env_clev": {},
            "local_u_cloud_levels_env_clev": {},
            "local_v_cloud_levels_env_clev": {},
            "local_env_saturation_moist_static_energy_cloud_levels_forced_env_clev": {},
            "local_geopotential_height_cloud_levels_forced_env_clev": {},
            "p_cloud_levels_forced_env_clev": {},
            "local_gamma_cloud_levels_forced_env_clev": {},
            "local_t_cloud_levels_forced_env_clev": {},
            "p_surface_env_clev": {},
            "t_surface_env_clev": {},
            "error_code_env_clev": {},
            "topography_height_no_negative_env_clev": {},
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
        locals.t_new.data[:] = inputs["local_t_new_env_clev"]
        locals.environment_saturation_mixing_ratio_forced.data[:] = inputs[
            "local_env_saturation_mixing_ratio_forced_env_clev"
        ]
        locals.vapor_forced.data[:] = inputs["local_vapor_forced_env_cond"]
        locals.environment_moist_static_energy_forced.data[:] = inputs[
            "local_env_moist_static_energy_forced_env_clev"
        ]
        locals.environment_saturation_moist_static_energy_forced.data[:] = inputs[
            "local_env_saturation_moist_static_energy_forced_env_clev"
        ]
        state.input_output.geopotential_height_forced.data[:] = inputs["geopotential_height_forced_env_clev"]
        state.input_output.p_forced.data[:] = inputs["p_forced_env_clev"]
        locals.environment_saturation_mixing_ratio_cloud_levels_forced.data[:] = inputs[
            "local_env_saturation_mixing_ratio_cloud_levels_forced_env_clev"
        ]
        locals.vapor_cloud_levels_forced.data[:] = inputs["local_vapor_cloud_levels_forced_env_clev"]
        locals.environment_moist_static_energy_cloud_levels_forced.data[:] = inputs[
            "local_env_moist_static_energy_cloud_levels_forced_env_clev"
        ]
        state.input_output.u.data[:] = inputs["u_env_clev"]
        state.input_output.v.data[:] = inputs["v_env_clev"]
        locals.u_cloud_levels.data[:] = inputs["local_u_cloud_levels_env_clev"]
        locals.v_cloud_levels.data[:] = inputs["local_v_cloud_levels_env_clev"]
        locals.environment_saturation_moist_static_energy_cloud_levels_forced.data[:] = inputs[
            "local_env_saturation_moist_static_energy_cloud_levels_forced_env_clev"
        ]
        locals.geopotential_height_cloud_levels_forced.data[:] = inputs[
            "local_geopotential_height_cloud_levels_forced_env_clev"
        ]
        state.output.p_cloud_levels_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "p_cloud_levels_forced_env_clev"
        ]
        locals.gamma_cloud_levels_forced.data[:] = inputs["local_gamma_cloud_levels_forced_env_clev"]
        locals.t_cloud_levels_forced.data[:] = inputs["local_t_cloud_levels_forced_env_clev"]
        state.input_output.p_surface.data[:] = inputs["p_surface_env_clev"]
        state.input_output.t_surface.data[:] = inputs["t_surface_env_clev"]
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "error_code_env_clev"
        ]
        state.input_output.topography_height_no_negative.data[:] = inputs[
            "topography_height_no_negative_env_clev"
        ]

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
                plume_dependent_constants=plume_dependent_constants,
                data_type=1,
            )

        outputs = {
            "local_t_new_env_clev": locals.t_new.field[:],
            "local_env_saturation_mixing_ratio_forced_env_clev": locals.environment_saturation_mixing_ratio_forced.field[
                :
            ],
            "local_vapor_forced_env_cond": locals.vapor_forced.field[:],
            "local_env_moist_static_energy_forced_env_clev": locals.environment_moist_static_energy_forced.field[
                :
            ],
            "local_env_saturation_moist_static_energy_forced_env_clev": locals.environment_saturation_moist_static_energy_forced.field[
                :
            ],
            "geopotential_height_forced_env_clev": state.input_output.geopotential_height_forced.field[:],
            "p_forced_env_clev": state.input_output.p_forced.field[:],
            "local_env_saturation_mixing_ratio_cloud_levels_forced_env_clev": locals.environment_saturation_mixing_ratio_cloud_levels_forced.field[
                :
            ],
            "local_vapor_cloud_levels_forced_env_clev": locals.vapor_cloud_levels_forced.field[:],
            "local_env_moist_static_energy_cloud_levels_forced_env_clev": locals.environment_moist_static_energy_cloud_levels_forced.field[
                :
            ],
            "u_env_clev": state.input_output.u.field[:],
            "v_env_clev": state.input_output.v.field[:],
            "local_u_cloud_levels_env_clev": locals.u_cloud_levels.field[:],
            "local_v_cloud_levels_env_clev": locals.v_cloud_levels.field[:],
            "local_env_saturation_moist_static_energy_cloud_levels_forced_env_clev": locals.environment_saturation_moist_static_energy_cloud_levels_forced.field[
                :
            ],
            "local_geopotential_height_cloud_levels_forced_env_clev": locals.geopotential_height_cloud_levels_forced.field[
                :
            ],
            "p_cloud_levels_forced_env_clev": state.output.p_cloud_levels_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_gamma_cloud_levels_forced_env_clev": locals.gamma_cloud_levels_forced.field[:],
            "local_t_cloud_levels_forced_env_clev": locals.t_cloud_levels_forced.field[:],
            "p_surface_env_clev": state.input_output.p_surface.field[:],
            "t_surface_env_clev": state.input_output.t_surface.field[:],
            "error_code_env_clev": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "topography_height_no_negative_env_clev": state.input_output.topography_height_no_negative.field[
                :
            ],
        }

        return outputs
