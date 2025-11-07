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
from pyMoist.convection.GF_2020.cumulus_parameterization.entrainment.entrainment import (
    EntrainmentRates,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.setup.set_constants import set_constants


class TranslateGF2020_CumulusParameterization_EntrainmentRates_shallow(TranslateFortranData2Py):
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
            "error_code_entrainmentrates": {},
            "local_vapor_cloud_levels_forced_ecv": {},
            "local_env_saturation_mixing_ratio_cloud_levels_forced_entrainmentrates": {},
            "lcl_level_entrainmentrates": {},
            "entrainment_rate_entrainmentrates": {},
            "local_updraft_detrainment_function_entrainmentrates": {},
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
            "error_code_entrainmentrates"
        ]
        locals.vapor_cloud_levels_forced.data[:] = inputs["local_vapor_cloud_levels_forced_ecv"]
        locals.environment_saturation_mixing_ratio_cloud_levels_forced.data[:] = inputs[
            "local_env_saturation_mixing_ratio_cloud_levels_forced_entrainmentrates"
        ]
        state.output.lcl_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["lcl_level_entrainmentrates"] - 1
        )
        state.output.entrainment_rate.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "entrainment_rate_entrainmentrates"
        ]
        locals.detrainment_function_updraft.data[:] = inputs[
            "local_updraft_detrainment_function_entrainmentrates"
        ]

        # initalize test code
        code = EntrainmentRates(
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

        # write output
        outputs = {
            "error_code_entrainmentrates": state.output.error_code.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_vapor_cloud_levels_forced_ecv": locals.vapor_cloud_levels_forced.field[:],
            "local_env_saturation_mixing_ratio_cloud_levels_forced_entrainmentrates": locals.environment_saturation_mixing_ratio_cloud_levels_forced.field[
                :
            ],
            "lcl_level_entrainmentrates": state.output.lcl_level.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ]
            + 1,
            "entrainment_rate_entrainmentrates": state.output.entrainment_rate.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_updraft_detrainment_function_entrainmentrates": locals.detrainment_function_updraft.field[
                :
            ],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_EntrainmentRates_mid(TranslateFortranData2Py):
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
            "error_code_entrainmentrates": {},
            "local_vapor_cloud_levels_forced_ecv": {},
            "local_env_saturation_mixing_ratio_cloud_levels_forced_entrainmentrates": {},
            "lcl_level_entrainmentrates": {},
            "entrainment_rate_entrainmentrates": {},
            "local_updraft_detrainment_function_entrainmentrates": {},
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
            "error_code_entrainmentrates"
        ]
        locals.vapor_cloud_levels_forced.data[:] = inputs["local_vapor_cloud_levels_forced_ecv"]
        locals.environment_saturation_mixing_ratio_cloud_levels_forced.data[:] = inputs[
            "local_env_saturation_mixing_ratio_cloud_levels_forced_entrainmentrates"
        ]
        state.output.lcl_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["lcl_level_entrainmentrates"] - 1
        )
        state.output.entrainment_rate.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "entrainment_rate_entrainmentrates"
        ]
        locals.detrainment_function_updraft.data[:] = inputs[
            "local_updraft_detrainment_function_entrainmentrates"
        ]

        # initalize test code
        code = EntrainmentRates(
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

        # write output
        outputs = {
            "error_code_entrainmentrates": state.output.error_code.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_vapor_cloud_levels_forced_ecv": locals.vapor_cloud_levels_forced.field[:],
            "local_env_saturation_mixing_ratio_cloud_levels_forced_entrainmentrates": locals.environment_saturation_mixing_ratio_cloud_levels_forced.field[
                :
            ],
            "lcl_level_entrainmentrates": state.output.lcl_level.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ]
            + 1,
            "entrainment_rate_entrainmentrates": state.output.entrainment_rate.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_updraft_detrainment_function_entrainmentrates": locals.detrainment_function_updraft.field[
                :
            ],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_EntrainmentRates_deep(TranslateFortranData2Py):
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
            "error_code_entrainmentrates": {},
            "local_vapor_cloud_levels_forced_ecv": {},
            "local_env_saturation_mixing_ratio_cloud_levels_forced_entrainmentrates": {},
            "lcl_level_entrainmentrates": {},
            "entrainment_rate_entrainmentrates": {},
            "local_updraft_detrainment_function_entrainmentrates": {},
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
            "error_code_entrainmentrates"
        ]
        locals.vapor_cloud_levels_forced.data[:] = inputs["local_vapor_cloud_levels_forced_ecv"]
        locals.environment_saturation_mixing_ratio_cloud_levels_forced.data[:] = inputs[
            "local_env_saturation_mixing_ratio_cloud_levels_forced_entrainmentrates"
        ]
        state.output.lcl_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["lcl_level_entrainmentrates"] - 1
        )
        state.output.entrainment_rate.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "entrainment_rate_entrainmentrates"
        ]
        locals.detrainment_function_updraft.data[:] = inputs[
            "local_updraft_detrainment_function_entrainmentrates"
        ]

        # initalize test code
        code = EntrainmentRates(
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

        # write output
        outputs = {
            "error_code_entrainmentrates": state.output.error_code.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_vapor_cloud_levels_forced_ecv": locals.vapor_cloud_levels_forced.field[:],
            "local_env_saturation_mixing_ratio_cloud_levels_forced_entrainmentrates": locals.environment_saturation_mixing_ratio_cloud_levels_forced.field[
                :
            ],
            "lcl_level_entrainmentrates": state.output.lcl_level.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ]
            + 1,
            "entrainment_rate_entrainmentrates": state.output.entrainment_rate.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_updraft_detrainment_function_entrainmentrates": locals.detrainment_function_updraft.field[
                :
            ],
        }

        return outputs
