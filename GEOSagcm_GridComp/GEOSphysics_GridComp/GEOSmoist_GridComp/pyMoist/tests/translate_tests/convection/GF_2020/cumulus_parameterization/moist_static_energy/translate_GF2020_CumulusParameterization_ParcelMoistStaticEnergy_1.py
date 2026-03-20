from f90nml import Namelist

from ndsl import StencilFactory
from ndsl.constants import I_DIM, J_DIM, K_DIM
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.constants import (
    MAXENS1,
    MAXENS2,
    MAXENS3,
    NUMBER_OF_PLUMES,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.locals import GF2020CumulusParameterizationLocals
from pyMoist.convection.GF_2020.cumulus_parameterization.moist_static_energy import parcel_moist_static_energy
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.setup.set_constants import set_constants
from pyMoist.convection.GF_2020.cumulus_parameterization.state import GF2020CumulusParameterizationState


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
            "error_code": {},
            "ocean_fraction": {},
            "local_vapor_excess": {},
            "local_t_excess": {},
            "local_add_buoy": {},
            "p_forced": {},
            "local_env_moist_static_energy_cloud_levels": {},
            "local_env_moist_static_energy_cloud_levels_forced": {},
            "local_moist_static_energy_origin_level": {},
            "local_moist_static_energy_origin_level_forced": {},
            "updraft_origin_level": {},
            "t_perturbation": {},
        }

        out_vars.update(in_vars["data_vars"])

    def __call__(self, constants: dict, cu_param_constants: dict, plume: str, **inputs):
        # initialize constants
        config = GF2020Config(SINGLE_COLUMN_MODE=False, **constants)
        cumulus_parameterization_config = GF2020CumulusParameterizationConfig(**cu_param_constants)
        plume_dependent_constants = GF2020PlumeDependentConstants()
        plume_dependent_constants = set_constants(
            cumulus_parameterization_config, plume_dependent_constants, plume
        )

        # initialize dataclasses
        state = GF2020CumulusParameterizationState.zeros(
            self.quantity_factory,
            data_dimensions={
                "plumes": NUMBER_OF_PLUMES,
                "convection_tracers": config.NUMBER_OF_TRACERS,
            },
        )

        locals = GF2020CumulusParameterizationLocals.zeros(
            self.quantity_factory,
            data_dimensions={
                "ensemble_1": MAXENS1,
                "ensemble_2": MAXENS2,
                "ensemble_3": MAXENS3,
                "ensemble_members": MAXENS1 * MAXENS2 * MAXENS3,
                "convection_tracers": config.NUMBER_OF_TRACERS,
            },
        )

        # fill relevant parts of dataclasses
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["error_code"]
        state.input.ocean_fraction.data[:] = inputs["ocean_fraction"]
        locals.vapor_excess.data[:] = inputs["local_vapor_excess"]
        locals.t_excess.data[:] = inputs["local_t_excess"]
        locals.add_buoyancy.data[:] = inputs["local_add_buoy"]
        state.input_output.p_forced.data[:] = inputs["p_forced"]
        locals.environment_moist_static_energy_cloud_levels.data[:] = inputs[
            "local_env_moist_static_energy_cloud_levels"
        ]
        locals.environment_moist_static_energy_cloud_levels_forced.data[:] = inputs[
            "local_env_moist_static_energy_cloud_levels_forced"
        ]
        locals.moist_static_energy_origin_level.data[:] = inputs["local_moist_static_energy_origin_level"]
        locals.moist_static_energy_origin_level_forced.data[:] = inputs[
            "local_moist_static_energy_origin_level_forced"
        ]
        state.output.updraft_origin_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["updraft_origin_level"] - 1
        )
        state.output.t_perturbation.data[:] = inputs["t_perturbation"]

        code = self.stencil_factory.from_dims_halo(
            func=parcel_moist_static_energy,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "BOUNDARY_CONDITION_METHOD": cumulus_parameterization_config.BOUNDARY_CONDITION_METHOD,
            },
        )

        if plume_dependent_constants.ENABLE_PLUME == 1:
            code(
                error_code=state.output.error_code,
                t_excess=locals.t_excess,
                vapor_excess=locals.vapor_excess,
                add_buoyancy=locals.add_buoyancy,
                ocean_fraction=state.input.ocean_fraction,
                updraft_origin_level=state.output.updraft_origin_level,
                p=state.input_output.p_forced,
                environmenet_moist_static_energy=locals.environment_moist_static_energy_cloud_levels,
                environmenet_moist_static_energy_forced=locals.environment_moist_static_energy_cloud_levels_forced,
                t_perturbation=state.output.t_perturbation,
                moist_static_energy_origin_level=locals.moist_static_energy_origin_level,
                moist_static_energy_origin_level_forced=locals.moist_static_energy_origin_level_forced,
                AVERAGE_LAYER_DEPTH=plume_dependent_constants.AVERAGE_LAYER_DEPTH,
                plume=plume_dependent_constants.PLUME_INDEX,
            )

        outputs = {
            "error_code": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "ocean_fraction": state.input.ocean_fraction.field[:],
            "local_vapor_excess": locals.vapor_excess.field[:],
            "local_t_excess": locals.t_excess.field[:],
            "local_add_buoy": locals.add_buoyancy.field[:],
            "p_forced": state.input_output.p_forced.field[:],
            "local_env_moist_static_energy_cloud_levels": locals.environment_moist_static_energy_cloud_levels.field[
                :
            ],
            "local_env_moist_static_energy_cloud_levels_forced": locals.environment_moist_static_energy_cloud_levels_forced.field[
                :
            ],
            "local_moist_static_energy_origin_level": locals.moist_static_energy_origin_level.field[:],
            "local_moist_static_energy_origin_level_forced": locals.moist_static_energy_origin_level_forced.field[
                :
            ],
            "updraft_origin_level": state.output.updraft_origin_level.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ]
            + 1,
            "t_perturbation": state.output.t_perturbation.field[:],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_ParcelMoistStaticEnergy_1_shallow(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_ParcelMoistStaticEnergy_1_mid(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_ParcelMoistStaticEnergy_1_deep(TranslateFortranData2Py):
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
