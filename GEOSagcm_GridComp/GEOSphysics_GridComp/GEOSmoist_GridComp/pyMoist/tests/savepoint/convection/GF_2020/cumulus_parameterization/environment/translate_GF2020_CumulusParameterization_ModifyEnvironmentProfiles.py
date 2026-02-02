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
from pyMoist.convection.GF_2020.cumulus_parameterization.constants import (
    MAXENS1,
    MAXENS2,
    MAXENS3,
    NUMBER_OF_PLUMES,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.environment import (
    modify_environment_profiles,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.setup.set_constants import set_constants
from ndsl.constants import X_DIM, Y_DIM, Z_DIM


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
            "cloud_top_level": {},
            "updraft_origin_level": {},
            "ocean_fraction": {},
            "p_forced": {},
            "local_t_new": {},
            "local_t_modified": {},
            "local_vapor_forced": {},
            "local_vapor_modified": {},
            "local_env_moist_static_energy_forced": {},
            "local_env_moist_static_energy_modified": {},
            "local_moist_static_energy_origin_level_forced": {},
            "local_moist_static_energy_origin_level_modified": {},
            "local_partition_liquid_ice": {},
            "local_del_moist_static_energy_cloud_ensemble": {},
            "local_del_t_cloud_ensemble": {},
            "local_del_vapor_cloud_ensemble": {},
            "local_del_cloud_liquid_cloud_ensemble": {},
            "local_del_u_cloud_ensemble": {},
            "local_del_v_cloud_ensemble": {},
            "local_moist_static_energy_tendency_from_environmental_subsidence": {},
            "local_vapor_tendency_from_environmental_subsidence": {},
            "local_t_tendency_from_environmental_subsidence": {},
            "local_arbitrary_numerical_parameter": {},
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
        state.output.cloud_top_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["cloud_top_level"] - 1
        )
        state.output.updraft_origin_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["updraft_origin_level"] - 1
        )
        state.input.ocean_fraction.data[:] = inputs["ocean_fraction"]
        state.input_output.p_forced.data[:] = inputs["p_forced"]
        locals.t_new.data[:] = inputs["local_t_new"]
        locals.t_modified.data[:] = inputs["local_t_modified"]
        locals.vapor_forced.data[:] = inputs["local_vapor_forced"]
        locals.vapor_modified.data[:] = inputs["local_vapor_modified"]
        locals.environment_moist_static_energy_forced.data[:] = inputs["local_env_moist_static_energy_forced"]
        locals.environment_moist_static_energy_modified.data[:] = inputs[
            "local_env_moist_static_energy_modified"
        ]
        locals.moist_static_energy_origin_level_forced.data[:] = inputs[
            "local_moist_static_energy_origin_level_forced"
        ]
        locals.moist_static_energy_origin_level_modified.data[:] = inputs[
            "local_moist_static_energy_origin_level_modified"
        ]
        locals.partition_liquid_ice.data[:] = inputs["local_partition_liquid_ice"]
        locals.del_moist_static_energy_cloud_ensemble.data[:] = inputs[
            "local_del_moist_static_energy_cloud_ensemble"
        ]
        locals.del_t_cloud_ensemble.data[:] = inputs["local_del_t_cloud_ensemble"]
        locals.del_vapor_cloud_ensemble.data[:] = inputs["local_del_vapor_cloud_ensemble"]
        locals.del_cloud_liquid_cloud_ensemble.data[:] = inputs["local_del_cloud_liquid_cloud_ensemble"]
        locals.del_u_cloud_ensemble.data[:] = inputs["local_del_u_cloud_ensemble"]
        locals.del_v_cloud_ensemble.data[:] = inputs["local_del_v_cloud_ensemble"]
        locals.moist_static_energy_tendency_from_environmental_subsidence.data[:] = inputs[
            "local_moist_static_energy_tendency_from_environmental_subsidence"
        ]
        locals.vapor_tendency_from_environmental_subsidence.data[:] = inputs[
            "local_vapor_tendency_from_environmental_subsidence"
        ]
        locals.t_tendency_from_environmental_subsidence.data[:] = inputs[
            "local_t_tendency_from_environmental_subsidence"
        ]
        locals.arbitrary_numerical_parameter.data[:] = inputs["local_arbitrary_numerical_parameter"]

        code = self.stencil_factory.from_dims_halo(
            func=modify_environment_profiles,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "COUPLE_MICROPHYSICS": cumulus_parameterization_config.COUPLE_MICROPHYSICS,
                "BOUNDARY_CONDITION_METHOD": cumulus_parameterization_config.BOUNDARY_CONDITION_METHOD,
            },
        )

        if plume_dependent_constants.ENABLE_PLUME == 1:
            code(
                error_code=state.output.error_code,
                cloud_top_level=state.output.cloud_top_level,
                updraft_origin_level=state.output.updraft_origin_level,
                ocean_fraction=state.input.ocean_fraction,
                p_forced=state.input_output.p_forced,
                t_new=locals.t_new,
                t_modified=locals.t_modified,
                vapor_forced=locals.vapor_forced,
                vapor_modified=locals.vapor_modified,
                environment_moist_static_energy_forced=locals.environment_moist_static_energy_forced,
                environment_moist_static_energy_modified=locals.environment_moist_static_energy_modified,
                moist_static_energy_origin_level_forced=locals.moist_static_energy_origin_level_forced,
                moist_static_energy_origin_level_modified=locals.moist_static_energy_origin_level_modified,
                partition_liquid_ice=locals.partition_liquid_ice,
                del_moist_static_energy_cloud_ensemble=locals.del_moist_static_energy_cloud_ensemble,
                del_t_cloud_ensemble=locals.del_t_cloud_ensemble,
                del_vapor_cloud_ensemble=locals.del_vapor_cloud_ensemble,
                del_cloud_liquid_cloud_ensemble=locals.del_cloud_liquid_cloud_ensemble,
                del_u_cloud_ensemble=locals.del_u_cloud_ensemble,
                del_v_cloud_ensemble=locals.del_v_cloud_ensemble,
                moist_static_energy_tendency_from_environmental_subsidence=locals.moist_static_energy_tendency_from_environmental_subsidence,
                vapor_tendency_from_environmental_subsidence=locals.vapor_tendency_from_environmental_subsidence,
                t_tendency_from_environmental_subsidence=locals.t_tendency_from_environmental_subsidence,
                arbitrary_numerical_parameter=locals.arbitrary_numerical_parameter,
                AVERAGE_LAYER_DEPTH=plume_dependent_constants.AVERAGE_LAYER_DEPTH,
                plume=plume_dependent_constants.PLUME_INDEX,
            )

        outputs = {
            "error_code": state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX],
            "cloud_top_level": state.output.cloud_top_level.data[:, :, plume_dependent_constants.PLUME_INDEX]
            + 1,
            "updraft_origin_level": state.output.updraft_origin_level.data[
                :, :, plume_dependent_constants.PLUME_INDEX
            ]
            + 1,
            "ocean_fraction": state.input.ocean_fraction.data[:],
            "p_forced": state.input_output.p_forced.data[:],
            "local_t_new": locals.t_new.data[:],
            "local_t_modified": locals.t_modified.data[:],
            "local_vapor_forced": locals.vapor_forced.data[:],
            "local_vapor_modified": locals.vapor_modified.data[:],
            "local_env_moist_static_energy_forced": locals.environment_moist_static_energy_forced.data[:],
            "local_env_moist_static_energy_modified": locals.environment_moist_static_energy_modified.data[:],
            "local_moist_static_energy_origin_level_forced": locals.moist_static_energy_origin_level_forced.data[
                :
            ],
            "local_moist_static_energy_origin_level_modified": locals.moist_static_energy_origin_level_modified.data[
                :
            ],
            "local_partition_liquid_ice": locals.partition_liquid_ice.data[:],
            "local_del_moist_static_energy_cloud_ensemble": locals.del_moist_static_energy_cloud_ensemble.data[
                :
            ],
            "local_del_t_cloud_ensemble": locals.del_t_cloud_ensemble.data[:],
            "local_del_vapor_cloud_ensemble": locals.del_vapor_cloud_ensemble.data[:],
            "local_del_cloud_liquid_cloud_ensemble": locals.del_cloud_liquid_cloud_ensemble.data[:],
            "local_del_u_cloud_ensemble": locals.del_u_cloud_ensemble.data[:],
            "local_del_v_cloud_ensemble": locals.del_v_cloud_ensemble.data[:],
            "local_moist_static_energy_tendency_from_environmental_subsidence": locals.moist_static_energy_tendency_from_environmental_subsidence.data[
                :
            ],
            "local_vapor_tendency_from_environmental_subsidence": locals.vapor_tendency_from_environmental_subsidence.data[
                :
            ],
            "local_t_tendency_from_environmental_subsidence": locals.t_tendency_from_environmental_subsidence.data[
                :
            ],
            "local_arbitrary_numerical_parameter": locals.arbitrary_numerical_parameter.data[:],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_ModifyEnvironmentProfiles_shallow(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_ModifyEnvironmentProfiles_mid(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_ModifyEnvironmentProfiles_deep(TranslateFortranData2Py):
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
