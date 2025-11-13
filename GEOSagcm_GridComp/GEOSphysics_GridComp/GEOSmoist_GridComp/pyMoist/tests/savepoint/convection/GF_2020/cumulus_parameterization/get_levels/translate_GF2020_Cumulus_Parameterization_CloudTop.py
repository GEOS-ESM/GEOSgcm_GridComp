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
    CloudTop,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.setup.set_constants import set_constants


class TranslateGF2020_CumulusParameterization_CloudTop_shallow(TranslateFortranData2Py):
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
            "error_code_cloudtop": {},
            "cloud_top_cloudtop": {},
            "entrainment_rate_cloudtop": {},
            "local_moist_static_energy_origin_level_forced_cloudtop": {},
            "local_env_moist_static_energy_forced_cloudtop": {},
            "local_env_saturation_moist_static_energy_cloud_levels_forced_cloudtop": {},
            "local_geopotential_height_cloud_levels_forced_cloudtop": {},
            "updraft_lfc_level_cloudtop": {},
            "local_cloud_moist_static_energy_forced_transported_cloudtop": {},
            "p_cloud_levels_forced_cloudtop": {},
            "local_geopotential_height_cloud_levels_cloudtop": {},
            "last_error_code_cloudtop": {},
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
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "error_code_cloudtop"
        ]
        state.output.cloud_top.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "cloud_top_cloudtop"
        ]
        state.output.entrainment_rate.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "entrainment_rate_cloudtop"
        ]
        locals.moist_static_energy_origin_level_forced.data[:] = inputs[
            "local_moist_static_energy_origin_level_forced_cloudtop"
        ]
        locals.environment_moist_static_energy_forced.data[:] = inputs[
            "local_env_moist_static_energy_forced_cloudtop"
        ]
        locals.environment_saturation_mixing_ratio_cloud_levels_forced.data[:] = inputs[
            "local_env_saturation_moist_static_energy_cloud_levels_forced_cloudtop"
        ]
        locals.geopotential_height_cloud_levels_forced.data[:] = inputs[
            "local_geopotential_height_cloud_levels_forced_cloudtop"
        ]
        state.output.updraft_lfc_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["updraft_lfc_level_cloudtop"] - 1
        )
        locals.cloud_moist_static_energy_forced_transported.data[:] = inputs[
            "local_cloud_moist_static_energy_forced_transported_cloudtop"
        ]
        state.output.p_cloud_levels_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "p_cloud_levels_forced_cloudtop"
        ]
        locals.geopotential_height_cloud_levels.data[:] = inputs[
            "local_geopotential_height_cloud_levels_cloudtop"
        ]
        state.input.last_error_code.data[:] = inputs["last_error_code_cloudtop"]

        # initalize test code
        code = CloudTop(
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

            state.output.cloud_top.field[:, :, plume_dependent_constants.PLUME_INDEX] += 1

        # write output
        outputs = {
            "error_code_cloudtop": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "cloud_top_cloudtop": state.output.cloud_top.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "entrainment_rate_cloudtop": state.output.entrainment_rate.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_moist_static_energy_origin_level_forced_cloudtop": locals.moist_static_energy_origin_level_forced.field[
                :
            ],
            "local_env_moist_static_energy_forced_cloudtop": locals.environment_moist_static_energy_forced.field[
                :
            ],
            "local_env_saturation_moist_static_energy_cloud_levels_forced_cloudtop": locals.environment_saturation_mixing_ratio_cloud_levels_forced.field[
                :
            ],
            "local_geopotential_height_cloud_levels_forced_cloudtop": locals.geopotential_height_cloud_levels_forced.field[
                :
            ],
            "updraft_lfc_level_cloudtop": state.output.updraft_lfc_level.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ]
            + 1,
            "local_cloud_moist_static_energy_forced_transported_cloudtop": locals.cloud_moist_static_energy_forced_transported.field[
                :
            ],
            "p_cloud_levels_forced_cloudtop": state.output.p_cloud_levels_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_geopotential_height_cloud_levels_cloudtop": locals.geopotential_height_cloud_levels.field[
                :
            ],
            "last_error_code_cloudtop": state.input.last_error_code.field[:],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_CloudTop_mid(TranslateFortranData2Py):
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
            "error_code_cloudtop": {},
            "cloud_top_cloudtop": {},
            "entrainment_rate_cloudtop": {},
            "local_moist_static_energy_origin_level_forced_cloudtop": {},
            "local_env_moist_static_energy_forced_cloudtop": {},
            "local_env_saturation_moist_static_energy_cloud_levels_forced_cloudtop": {},
            "local_geopotential_height_cloud_levels_forced_cloudtop": {},
            "updraft_lfc_level_cloudtop": {},
            "local_cloud_moist_static_energy_forced_transported_cloudtop": {},
            "p_cloud_levels_forced_cloudtop": {},
            "local_geopotential_height_cloud_levels_cloudtop": {},
            "last_error_code_cloudtop": {},
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
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "error_code_cloudtop"
        ]
        state.output.cloud_top.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "cloud_top_cloudtop"
        ]
        state.output.entrainment_rate.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "entrainment_rate_cloudtop"
        ]
        locals.moist_static_energy_origin_level_forced.data[:] = inputs[
            "local_moist_static_energy_origin_level_forced_cloudtop"
        ]
        locals.environment_moist_static_energy_forced.data[:] = inputs[
            "local_env_moist_static_energy_forced_cloudtop"
        ]
        locals.environment_saturation_mixing_ratio_cloud_levels_forced.data[:] = inputs[
            "local_env_saturation_moist_static_energy_cloud_levels_forced_cloudtop"
        ]
        locals.geopotential_height_cloud_levels_forced.data[:] = inputs[
            "local_geopotential_height_cloud_levels_forced_cloudtop"
        ]
        state.output.updraft_lfc_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["updraft_lfc_level_cloudtop"] - 1
        )
        locals.cloud_moist_static_energy_forced_transported.data[:] = inputs[
            "local_cloud_moist_static_energy_forced_transported_cloudtop"
        ]
        state.output.p_cloud_levels_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "p_cloud_levels_forced_cloudtop"
        ]
        locals.geopotential_height_cloud_levels.data[:] = inputs[
            "local_geopotential_height_cloud_levels_cloudtop"
        ]
        state.input.last_error_code.data[:] = inputs["last_error_code_cloudtop"]

        # initalize test code
        code = CloudTop(
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

            state.output.cloud_top.field[:, :, plume_dependent_constants.PLUME_INDEX] += 1

        # write output
        outputs = {
            "error_code_cloudtop": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "cloud_top_cloudtop": state.output.cloud_top.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "entrainment_rate_cloudtop": state.output.entrainment_rate.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_moist_static_energy_origin_level_forced_cloudtop": locals.moist_static_energy_origin_level_forced.field[
                :
            ],
            "local_env_moist_static_energy_forced_cloudtop": locals.environment_moist_static_energy_forced.field[
                :
            ],
            "local_env_saturation_moist_static_energy_cloud_levels_forced_cloudtop": locals.environment_saturation_mixing_ratio_cloud_levels_forced.field[
                :
            ],
            "local_geopotential_height_cloud_levels_forced_cloudtop": locals.geopotential_height_cloud_levels_forced.field[
                :
            ],
            "updraft_lfc_level_cloudtop": state.output.updraft_lfc_level.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ]
            + 1,
            "local_cloud_moist_static_energy_forced_transported_cloudtop": locals.cloud_moist_static_energy_forced_transported.field[
                :
            ],
            "p_cloud_levels_forced_cloudtop": state.output.p_cloud_levels_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_geopotential_height_cloud_levels_cloudtop": locals.geopotential_height_cloud_levels.field[
                :
            ],
            "last_error_code_cloudtop": state.input.last_error_code.field[:],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_CloudTop_deep(TranslateFortranData2Py):
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
            "error_code_cloudtop": {},
            "cloud_top_cloudtop": {},
            "entrainment_rate_cloudtop": {},
            "local_moist_static_energy_origin_level_forced_cloudtop": {},
            "local_env_moist_static_energy_forced_cloudtop": {},
            "local_env_saturation_moist_static_energy_cloud_levels_forced_cloudtop": {},
            "local_geopotential_height_cloud_levels_forced_cloudtop": {},
            "updraft_lfc_level_cloudtop": {},
            "local_cloud_moist_static_energy_forced_transported_cloudtop": {},
            "p_cloud_levels_forced_cloudtop": {},
            "local_geopotential_height_cloud_levels_cloudtop": {},
            "last_error_code_cloudtop": {},
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
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "error_code_cloudtop"
        ]
        state.output.cloud_top.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "cloud_top_cloudtop"
        ]
        state.output.entrainment_rate.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "entrainment_rate_cloudtop"
        ]
        locals.moist_static_energy_origin_level_forced.data[:] = inputs[
            "local_moist_static_energy_origin_level_forced_cloudtop"
        ]
        locals.environment_moist_static_energy_forced.data[:] = inputs[
            "local_env_moist_static_energy_forced_cloudtop"
        ]
        locals.environment_saturation_mixing_ratio_cloud_levels_forced.data[:] = inputs[
            "local_env_saturation_moist_static_energy_cloud_levels_forced_cloudtop"
        ]
        locals.geopotential_height_cloud_levels_forced.data[:] = inputs[
            "local_geopotential_height_cloud_levels_forced_cloudtop"
        ]
        state.output.updraft_lfc_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["updraft_lfc_level_cloudtop"] - 1
        )
        locals.cloud_moist_static_energy_forced_transported.data[:] = inputs[
            "local_cloud_moist_static_energy_forced_transported_cloudtop"
        ]
        state.output.p_cloud_levels_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "p_cloud_levels_forced_cloudtop"
        ]
        locals.geopotential_height_cloud_levels.data[:] = inputs[
            "local_geopotential_height_cloud_levels_cloudtop"
        ]
        state.input.last_error_code.data[:] = inputs["last_error_code_cloudtop"]

        # initalize test code
        code = CloudTop(
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

            state.output.cloud_top.field[:, :, plume_dependent_constants.PLUME_INDEX] += 1

        # write output
        outputs = {
            "error_code_cloudtop": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "cloud_top_cloudtop": state.output.cloud_top.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "entrainment_rate_cloudtop": state.output.entrainment_rate.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_moist_static_energy_origin_level_forced_cloudtop": locals.moist_static_energy_origin_level_forced.field[
                :
            ],
            "local_env_moist_static_energy_forced_cloudtop": locals.environment_moist_static_energy_forced.field[
                :
            ],
            "local_env_saturation_moist_static_energy_cloud_levels_forced_cloudtop": locals.environment_saturation_mixing_ratio_cloud_levels_forced.field[
                :
            ],
            "local_geopotential_height_cloud_levels_forced_cloudtop": locals.geopotential_height_cloud_levels_forced.field[
                :
            ],
            "updraft_lfc_level_cloudtop": state.output.updraft_lfc_level.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ]
            + 1,
            "local_cloud_moist_static_energy_forced_transported_cloudtop": locals.cloud_moist_static_energy_forced_transported.field[
                :
            ],
            "p_cloud_levels_forced_cloudtop": state.output.p_cloud_levels_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_geopotential_height_cloud_levels_cloudtop": locals.geopotential_height_cloud_levels.field[
                :
            ],
            "last_error_code_cloudtop": state.input.last_error_code.field[:],
        }

        return outputs
