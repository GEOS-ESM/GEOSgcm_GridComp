from f90nml import Namelist
from ndsl import StencilFactory
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.convection.GF_2020.cumulus_parameterization.state import (
    GF2020CumulusParameterizationState,
)
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import (
    GF2020CumulusParameterizationConfig,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.locals import (
    GF2020CumulusParameterizationLocals,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.constants import (
    MAXENS1,
    MAXENS2,
    MAXENS3,
    NUMBER_OF_PLUMES,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.vertical_discretization import VerticalDiscretization
from pyMoist.convection.GF_2020.cumulus_parameterization.setup.set_constants import (
    set_constants,
)


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
            "error_code_vertdisc": {},
            "cloud_top_level_vertdisc": {},
            "p_cloud_levels_forced_vertdisc": {},
            "normalized_massflux_updraft_forced_vertdisc": {},
            "normalized_massflux_downdraft_forced_vertdisc": {},
            "local_environment_massflux_vertdisc": {},
            "u_vertdisc": {},
            "v_vertdisc": {},
            "local_u_cloud_levels_vertdisc": {},
            "local_v_cloud_levels_vertdisc": {},
            "local_u_c_vertdisc": {},
            "local_v_c_vertdisc": {},
            "local_u_c_downdraft_vertdisc": {},
            "local_v_c_downdraft_vertdisc": {},
            "local_env_moist_static_energy_cloud_levels_forced_vertdisc": {},
            "epsilon_forced_vertdisc": {},
            "local_del_u_cloud_ensemble_vertdisc": {},
            "local_del_v_cloud_ensemble_vertdisc": {},
            "local_del_moist_static_energy_cloud_ensemble_vertdisc": {},
            "local_del_t_cloud_ensemble_vertdisc": {},
            "local_del_vapor_cloud_ensemble_vertdisc": {},
            "local_del_cloud_liquid_cloud_ensemble_vertdisc": {},
            "local_del_buoyancy_cloud_ensemble_vertdisc": {},
            "local_t_tendency_from_environmental_subsidence_vertdisc": {},
            "local_moist_static_energy_tendency_from_environmental_subsidence_vertdisc": {},
            "local_vapor_tendency_from_environmental_subsidence_vertdisc": {},
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
            },
        )

        locals = GF2020CumulusParameterizationLocals.zeros(
            self.quantity_factory,
            data_dimensions={
                "ensemble_1": MAXENS1,
                "ensemble_2": MAXENS2,
                "ensemble_3": MAXENS3,
                "ensemble_members": MAXENS1 * MAXENS2 * MAXENS3,
            },
        )

        # fill relevant parts of dataclasses
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "error_code_vertdisc"
        ]
        state.output.cloud_top_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["cloud_top_level_vertdisc"] - 1
        )
        state.output.p_cloud_levels_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "p_cloud_levels_forced_vertdisc"
        ]
        state.output.normalized_massflux_updraft_forced.data[
            :, :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["normalized_massflux_updraft_forced_vertdisc"]
        state.output.normalized_massflux_downdraft_forced.data[
            :, :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["normalized_massflux_downdraft_forced_vertdisc"]
        locals.environment_massflux.data[:] = inputs["local_environment_massflux_vertdisc"]
        state.input_output.u.data[:] = inputs["u_vertdisc"]
        state.input_output.v.data[:] = inputs["v_vertdisc"]
        locals.u_cloud_levels.data[:] = inputs["local_u_cloud_levels_vertdisc"]
        locals.v_cloud_levels.data[:] = inputs["local_v_cloud_levels_vertdisc"]
        locals.u_c.data[:] = inputs["local_u_c_vertdisc"]
        locals.v_c.data[:] = inputs["local_v_c_vertdisc"]
        locals.u_c_downdraft.data[:] = inputs["local_u_c_downdraft_vertdisc"]
        locals.v_c_downdraft.data[:] = inputs["local_v_c_downdraft_vertdisc"]
        locals.environment_moist_static_energy_cloud_levels_forced.data[:] = inputs[
            "local_env_moist_static_energy_cloud_levels_forced_vertdisc"
        ]
        state.output.epsilon_forced.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "epsilon_forced_vertdisc"
        ]
        locals.del_u_cloud_ensemble.data[:] = inputs["local_del_u_cloud_ensemble_vertdisc"]
        locals.del_v_cloud_ensemble.data[:] = inputs["local_del_v_cloud_ensemble_vertdisc"]
        locals.del_moist_static_energy_cloud_ensemble.data[:] = inputs[
            "local_del_moist_static_energy_cloud_ensemble_vertdisc"
        ]
        locals.del_t_cloud_ensemble.data[:] = inputs["local_del_t_cloud_ensemble_vertdisc"]
        locals.del_vapor_cloud_ensemble.data[:] = inputs["local_del_vapor_cloud_ensemble_vertdisc"]
        locals.del_cloud_liquid_cloud_ensemble.data[:] = inputs[
            "local_del_cloud_liquid_cloud_ensemble_vertdisc"
        ]
        locals.del_buoyancy_cloud_ensemble.data[:] = inputs["local_del_buoyancy_cloud_ensemble_vertdisc"]
        locals.t_tendency_from_environmental_subsidence.data[:] = inputs[
            "local_t_tendency_from_environmental_subsidence_vertdisc"
        ]
        locals.moist_static_energy_tendency_from_environmental_subsidence.data[:] = inputs[
            "local_moist_static_energy_tendency_from_environmental_subsidence_vertdisc"
        ]
        locals.vapor_tendency_from_environmental_subsidence.data[:] = inputs[
            "local_vapor_tendency_from_environmental_subsidence_vertdisc"
        ]

        # initalize test code
        code = VerticalDiscretization(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        # call test code
        if plume_dependent_constants.ENABLE_PLUME == 1:
            code(
                error_code=state.output.error_code,
                cloud_top_level=state.output.cloud_top_level,
                p_cloud_levels_forced=state.output.p_cloud_levels_forced,
                normalized_massflux_updraft_forced=state.output.normalized_massflux_updraft_forced,
                normalized_massflux_downdraft_forced=state.output.normalized_massflux_downdraft_forced,
                environment_massflux=locals.environment_massflux,
                u=state.input_output.u,
                v=state.input_output.v,
                u_cloud_levels=locals.u_cloud_levels,
                v_cloud_levels=locals.v_cloud_levels,
                u_c=locals.u_c,
                v_c=locals.v_c,
                u_c_downdraft=locals.u_c_downdraft,
                v_c_downdraft=locals.v_c_downdraft,
                del_u_cloud_ensemble=locals.del_u_cloud_ensemble,
                del_v_cloud_ensemble=locals.del_v_cloud_ensemble,
                del_moist_static_energy_cloud_ensemble=locals.del_moist_static_energy_cloud_ensemble,
                del_t_cloud_ensemble=locals.del_t_cloud_ensemble,
                del_vapor_cloud_ensemble=locals.del_vapor_cloud_ensemble,
                del_cloud_liquid_cloud_ensemble=locals.del_cloud_liquid_cloud_ensemble,
                del_buoyancy_cloud_ensemble=locals.del_buoyancy_cloud_ensemble,
                t_tendency_from_environmental_subsidence=locals.t_tendency_from_environmental_subsidence,
                moist_static_energy_tendency_from_environmental_subsidence=locals.moist_static_energy_tendency_from_environmental_subsidence,
                vapor_tendency_from_environmental_subsidence=locals.vapor_tendency_from_environmental_subsidence,
                epsilon_forced=state.output.epsilon_forced,
                plume_dependent_constants=plume_dependent_constants,
            )

        # write output
        outputs = {
            "error_code_vertdisc": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "cloud_top_level_vertdisc": state.output.cloud_top_level.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ]
            + 1,
            "p_cloud_levels_forced_vertdisc": state.output.p_cloud_levels_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "normalized_massflux_updraft_forced_vertdisc": state.output.normalized_massflux_updraft_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "normalized_massflux_downdraft_forced_vertdisc": state.output.normalized_massflux_downdraft_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_environment_massflux_vertdisc": locals.environment_massflux.field[:],
            "u_vertdisc": state.input_output.u.field[:],
            "v_vertdisc": state.input_output.v.field[:],
            "local_u_cloud_levels_vertdisc": locals.u_cloud_levels.field[:],
            "local_v_cloud_levels_vertdisc": locals.v_cloud_levels.field[:],
            "local_u_c_vertdisc": locals.u_c.field[:],
            "local_v_c_vertdisc": locals.v_c.field[:],
            "local_u_c_downdraft_vertdisc": locals.u_c_downdraft.field[:],
            "local_v_c_downdraft_vertdisc": locals.v_c_downdraft.field[:],
            "local_env_moist_static_energy_cloud_levels_forced_vertdisc": locals.environment_moist_static_energy_cloud_levels_forced.field[
                :
            ],
            "epsilon_forced_vertdisc": state.output.epsilon_forced.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_del_u_cloud_ensemble_vertdisc": locals.del_u_cloud_ensemble.field[:],
            "local_del_v_cloud_ensemble_vertdisc": locals.del_v_cloud_ensemble.field[:],
            "local_del_moist_static_energy_cloud_ensemble_vertdisc": locals.del_moist_static_energy_cloud_ensemble.field[
                :
            ],
            "local_del_t_cloud_ensemble_vertdisc": locals.del_t_cloud_ensemble.field[:],
            "local_del_vapor_cloud_ensemble_vertdisc": locals.del_vapor_cloud_ensemble.field[:],
            "local_del_cloud_liquid_cloud_ensemble_vertdisc": locals.del_cloud_liquid_cloud_ensemble.field[:],
            "local_del_buoyancy_cloud_ensemble_vertdisc": locals.del_buoyancy_cloud_ensemble.field[:],
            "local_t_tendency_from_environmental_subsidence_vertdisc": locals.t_tendency_from_environmental_subsidence.field[
                :
            ],
            "local_moist_static_energy_tendency_from_environmental_subsidence_vertdisc": locals.moist_static_energy_tendency_from_environmental_subsidence.field[
                :
            ],
            "local_vapor_tendency_from_environmental_subsidence_vertdisc": locals.vapor_tendency_from_environmental_subsidence.field[
                :
            ],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_VerticalDiscretization_shallow(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_VerticalDiscretization_mid(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_VerticalDiscretization_deep(TranslateFortranData2Py):
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
