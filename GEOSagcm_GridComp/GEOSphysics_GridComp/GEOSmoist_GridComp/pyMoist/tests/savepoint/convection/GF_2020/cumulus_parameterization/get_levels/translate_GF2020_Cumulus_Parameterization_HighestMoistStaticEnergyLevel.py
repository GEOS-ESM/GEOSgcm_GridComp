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
from pyMoist.convection.GF_2020.cumulus_parameterization.constants import MAXENS1, MAXENS2, MAXENS3
from pyMoist.convection.GF_2020.cumulus_parameterization.get_levels.get_levels import (
    HighestMoistStaticEnergyLevel,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.setup.set_constants import set_constants


class TranslateGF2020_CumulusParameterization_HighestMoistStaticEnergyLevel_shallow(TranslateFortranData2Py):
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
            "local_updraft_origin_level": {},
            "error_code": {},
            "local_env_moist_static_energy_cloud_levels_forced": {},
            "local_maximum_updraft_origin_level": {},
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
        locals.updraft_origin_level.data[:] = inputs["local_updraft_origin_level"]
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["error_code"]
        locals.environment_moist_static_energy_cloud_levels_forced.data[:] = inputs[
            "local_env_moist_static_energy_cloud_levels_forced"
        ]
        locals.maximum_updraft_origin_level.data[:] = inputs["local_maximum_updraft_origin_level"] - 1

        # initalize test code
        highest_moist_static_energy_level = HighestMoistStaticEnergyLevel(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        # call test code
        if plume_dependent_constants.ENABLE_PLUME == 1:
            highest_moist_static_energy_level(
                state=state,
                locals=locals,
                plume_dependent_constants=plume_dependent_constants,
            )

            locals.updraft_origin_level.field[:] = locals.updraft_origin_level.field[:] + 1

        # write output
        outputs = {
            "local_updraft_origin_level": locals.updraft_origin_level.field[:],
            "error_code": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "local_env_moist_static_energy_cloud_levels_forced": locals.environment_moist_static_energy_cloud_levels_forced.field[
                :
            ],
            "local_maximum_updraft_origin_level": locals.maximum_updraft_origin_level.field[:] + 1,
        }

        return outputs


class TranslateGF2020_CumulusParameterization_HighestMoistStaticEnergyLevel_mid(TranslateFortranData2Py):
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
            "local_updraft_origin_level": {},
            "error_code": {},
            "local_env_moist_static_energy_cloud_levels_forced": {},
            "local_maximum_updraft_origin_level": {},
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
        locals.updraft_origin_level.data[:] = inputs["local_updraft_origin_level"]
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["error_code"]
        locals.environment_moist_static_energy_cloud_levels_forced.data[:] = inputs[
            "local_env_moist_static_energy_cloud_levels_forced"
        ]
        locals.maximum_updraft_origin_level.data[:] = inputs["local_maximum_updraft_origin_level"] - 1

        # initalize test code
        highest_moist_static_energy_level = HighestMoistStaticEnergyLevel(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        # call test code
        if plume_dependent_constants.ENABLE_PLUME == 1:
            highest_moist_static_energy_level(
                state=state,
                locals=locals,
                plume_dependent_constants=plume_dependent_constants,
            )

            locals.updraft_origin_level.field[:] = locals.updraft_origin_level.field[:] + 1

        # write output
        outputs = {
            "local_updraft_origin_level": locals.updraft_origin_level.field[:],
            "error_code": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "local_env_moist_static_energy_cloud_levels_forced": locals.environment_moist_static_energy_cloud_levels_forced.field[
                :
            ],
            "local_maximum_updraft_origin_level": locals.maximum_updraft_origin_level.field[:] + 1,
        }

        return outputs


class TranslateGF2020_CumulusParameterization_HighestMoistStaticEnergyLevel_deep(TranslateFortranData2Py):
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
            "local_updraft_origin_level": {},
            "error_code": {},
            "local_env_moist_static_energy_cloud_levels_forced": {},
            "local_maximum_updraft_origin_level": {},
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
        locals.updraft_origin_level.data[:] = inputs["local_updraft_origin_level"]
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["error_code"]
        locals.environment_moist_static_energy_cloud_levels_forced.data[:] = inputs[
            "local_env_moist_static_energy_cloud_levels_forced"
        ]
        locals.maximum_updraft_origin_level.data[:] = inputs["local_maximum_updraft_origin_level"] - 1

        # initalize test code
        highest_moist_static_energy_level = HighestMoistStaticEnergyLevel(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        # call test code
        if plume_dependent_constants.ENABLE_PLUME == 1:
            highest_moist_static_energy_level(
                state=state,
                locals=locals,
                plume_dependent_constants=plume_dependent_constants,
            )

            locals.updraft_origin_level.field[:] = locals.updraft_origin_level.field[:] + 1

        # write output
        outputs = {
            "local_updraft_origin_level": locals.updraft_origin_level.field[:],
            "error_code": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "local_env_moist_static_energy_cloud_levels_forced": locals.environment_moist_static_energy_cloud_levels_forced.field[
                :
            ],
            "local_maximum_updraft_origin_level": locals.maximum_updraft_origin_level.field[:] + 1,
        }

        return outputs
