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
from pyMoist.convection.GF_2020.cumulus_parameterization.moist_static_energy.moist_static_energy import (
    FirstGuessMoistStaticEnergy,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.setup.set_constants import set_constants


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
            "error_code_firstmse": {},
            "local_start_level_firstmse": {},
            "cloud_top_firstmse": {},
            "mass_detrainment_updraft_forced_firstmse": {},
            "mass_entrainment_updraft_forced_firstmse": {},
            "local_noramlized_massflux_updraft_firstmse": {},
            "local_cloud_moist_static_energy_forced_firstmse": {},
            "noramlized_massflux_updraft_forced_firstmse": {},
            "local_env_moist_static_energy_forced_firstmse": {},
            "local_vapor_excess_firstmse": {},
            "local_t_excess_firstmse": {},
            "local_add_buoyancy_firstmse": {},
            "local_env_saturation_moist_static_energy_cloud_levels_forced_firstmse": {},
        }

        out_vars.update(in_vars["data_vars"])

    def __call__(self, constants: dict, cu_param_constants: dict, plume: str, **inputs):
        # initalize constants
        config = GF2020Config(SINGLE_COLUMN_MODE=False, **constants)
        cumulus_parameterization_config = GF2020CumulusParameterizationConfig(**cu_param_constants)
        plume_dependent_constants = GF2020PlumeDependentConstants()
        plume_dependent_constants = set_constants(
            cumulus_parameterization_config, plume_dependent_constants, plume
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
            "error_code_firstmse"
        ]
        locals.start_level.data[:] = inputs["local_start_level_firstmse"] - 1
        state.output.cloud_top.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["cloud_top_firstmse"] - 1
        )
        state.output.mass_detrainment_updraft_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["mass_detrainment_updraft_forced_firstmse"]
        )
        state.output.mass_entrainment_updraft_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["mass_entrainment_updraft_forced_firstmse"]
        )
        locals.normalized_massflux_updraft.data[:] = inputs["local_noramlized_massflux_updraft_firstmse"]
        locals.cloud_moist_static_energy_forced.data[:] = inputs[
            "local_cloud_moist_static_energy_forced_firstmse"
        ]
        state.output.normalized_massflux_updraft_forced.data[
            :, :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["noramlized_massflux_updraft_forced_firstmse"]
        locals.environment_moist_static_energy_forced.data[:] = inputs[
            "local_env_moist_static_energy_forced_firstmse"
        ]
        locals.vapor_excess.data[:] = inputs["local_vapor_excess_firstmse"]
        locals.t_excess.data[:] = inputs["local_t_excess_firstmse"]
        locals.add_buoyancy.data[:] = inputs["local_add_buoyancy_firstmse"]
        locals.environment_saturation_moist_static_energy_cloud_levels_forced.data[:] = inputs[
            "local_env_saturation_moist_static_energy_cloud_levels_forced_firstmse"
        ]

        # initalize test code
        code = FirstGuessMoistStaticEnergy(
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
            "error_code_firstmse": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "local_start_level_firstmse": locals.start_level.field[:] + 1,
            "cloud_top_firstmse": state.output.cloud_top.field[:, :, plume_dependent_constants.PLUME_INDEX]
            + 1,
            "mass_detrainment_updraft_forced_firstmse": state.output.mass_detrainment_updraft_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "mass_entrainment_updraft_forced_firstmse": state.output.mass_entrainment_updraft_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_noramlized_massflux_updraft_firstmse": locals.normalized_massflux_updraft.field[:],
            "local_cloud_moist_static_energy_forced_firstmse": locals.cloud_moist_static_energy_forced,
            "noramlized_massflux_updraft_forced_firstmse": state.output.normalized_massflux_updraft_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_env_moist_static_energy_forced_firstmse": locals.environment_moist_static_energy_forced.field[
                :
            ],
            "local_vapor_excess_firstmse": locals.vapor_excess.field[:],
            "local_t_excess_firstmse": locals.t_excess.field[:],
            "local_add_buoyancy_firstmse": locals.add_buoyancy.field[:],
            "local_env_saturation_moist_static_energy_cloud_levels_forced_firstmse": locals.environment_saturation_moist_static_energy_cloud_levels_forced.field[
                :
            ],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_FirstGuessMoistStaticEnergy_shallow(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_FirstGuessMoistStaticEnergy_mid(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_FirstGuessMoistStaticEnergy_deep(TranslateFortranData2Py):
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
