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
)
from pyMoist.convection.GF_2020.cumulus_parameterization.updraft import (
    updraft_moist_static_energy_and_momentum_budget,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.setup.set_constants import (
    set_constants,
)
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
            "error_code_upmsemomentumo": {},
            "local_start_level_upmsemomentumo": {},
            "cloud_top_level_upmsemomentumo": {},
            "p_forced_upmsemomentumo": {},
            "local_env_moist_static_energy_upmsemomentumo": {},
            "local_env_moist_static_energy_forced_upmsemomentumo": {},
            "local_env_moist_static_energy_cloud_levels_upmsemomentumo": {},
            "local_env_moist_static_energy_cloud_levels_forced_upmsemomentumo": {},
            "local_env_saturation_moist_static_energy_cloud_levels_upmsemomentumo": {},
            "local_env_saturation_moist_static_energy_cloud_levels_forced_upmsemomentumo": {},
            "local_cloud_moist_static_energy_upmsemomentumo": {},
            "local_cloud_moist_static_energy_forced_upmsemomentumo": {},
            "local_normalized_massflux_updraft_upmsemomentumo": {},
            "local_normalized_massflux_updraft_forced_upmsemomentumo": {},
            "local_mass_entrainment_updraft_upmsemomentumo": {},
            "local_mass_detrainment_updraft_upmsemomentumo": {},
            "local_mass_entrainment_u_updraft_upmsemomentumo": {},
            "local_mass_detrainment_u_updraft_upmsemomentumo": {},
            "mass_detrainment_updraft_forced_upmsemomentumo": {},
            "mass_entrainment_updraft_forced_upmsemomentumo": {},
            "u_upmsemomentumo": {},
            "v_upmsemomentumo": {},
            "local_u_c_upmsemomentumo": {},
            "local_v_c_upmsemomentumo": {},
            "local_u_cloud_levels_upmsemomentumo": {},
            "local_v_cloud_levels_upmsemomentumo": {},
            "local_partition_liquid_ice_upmsemomentumo": {},
            "cloud_liquid_after_rain_forced_upmsemomentumo": {},
            "local_vapor_excess_upmsemomentumo": {},
            "local_t_excess_upmsemomentumo": {},
            "local_add_buoyancy_upmsemomentumo": {},
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
            "error_code_upmsemomentumo"
        ]
        locals.start_level.data[:] = inputs["local_start_level_upmsemomentumo"] - 1
        state.output.cloud_top_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["cloud_top_level_upmsemomentumo"] - 1
        )
        state.input_output.p_forced.data[:] = inputs["p_forced_upmsemomentumo"]
        locals.environment_moist_static_energy.data[:] = inputs[
            "local_env_moist_static_energy_upmsemomentumo"
        ]
        locals.environment_moist_static_energy_forced.data[:] = inputs[
            "local_env_moist_static_energy_forced_upmsemomentumo"
        ]
        locals.environment_moist_static_energy_cloud_levels.data[:] = inputs[
            "local_env_moist_static_energy_cloud_levels_upmsemomentumo"
        ]
        locals.environment_moist_static_energy_cloud_levels_forced.data[:] = inputs[
            "local_env_moist_static_energy_cloud_levels_forced_upmsemomentumo"
        ]
        locals.environment_saturation_moist_static_energy_cloud_levels.data[:] = inputs[
            "local_env_saturation_moist_static_energy_cloud_levels_upmsemomentumo"
        ]
        locals.environment_saturation_moist_static_energy_cloud_levels_forced.data[:] = inputs[
            "local_env_saturation_moist_static_energy_cloud_levels_forced_upmsemomentumo"
        ]
        locals.cloud_moist_static_energy.data[:] = inputs["local_cloud_moist_static_energy_upmsemomentumo"]
        locals.cloud_moist_static_energy_forced.data[:] = inputs[
            "local_cloud_moist_static_energy_forced_upmsemomentumo"
        ]
        locals.normalized_massflux_updraft.data[:] = inputs[
            "local_normalized_massflux_updraft_upmsemomentumo"
        ]
        locals.normalized_massflux_updraft_forced.data[:] = inputs[
            "local_normalized_massflux_updraft_forced_upmsemomentumo"
        ]
        locals.mass_entrainment_updraft.data[:] = inputs["local_mass_entrainment_updraft_upmsemomentumo"]
        locals.mass_detrainment_updraft.data[:] = inputs["local_mass_detrainment_updraft_upmsemomentumo"]
        locals.mass_entrainment_u_updraft.data[:] = inputs["local_mass_entrainment_u_updraft_upmsemomentumo"]
        locals.mass_detrainment_u_updraft.data[:] = inputs["local_mass_detrainment_u_updraft_upmsemomentumo"]
        state.output.mass_detrainment_updraft_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["mass_detrainment_updraft_forced_upmsemomentumo"]
        )
        state.output.mass_entrainment_updraft_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["mass_entrainment_updraft_forced_upmsemomentumo"]
        )
        state.input_output.u.data[:] = inputs["u_upmsemomentumo"]
        state.input_output.v.data[:] = inputs["v_upmsemomentumo"]
        locals.u_c.data[:] = inputs["local_u_c_upmsemomentumo"]
        locals.v_c.data[:] = inputs["local_v_c_upmsemomentumo"]
        locals.u_cloud_levels.data[:] = inputs["local_u_cloud_levels_upmsemomentumo"]
        locals.v_cloud_levels.data[:] = inputs["local_v_cloud_levels_upmsemomentumo"]
        locals.partition_liquid_ice.data[:] = inputs["local_partition_liquid_ice_upmsemomentumo"]
        state.output.cloud_liquid_after_rain_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["cloud_liquid_after_rain_forced_upmsemomentumo"]
        )
        locals.vapor_excess.data[:] = inputs["local_vapor_excess_upmsemomentumo"]
        locals.t_excess.data[:] = inputs["local_t_excess_upmsemomentumo"]
        locals.add_buoyancy.data[:] = inputs["local_add_buoyancy_upmsemomentumo"]

        # initalize test code
        code = self.stencil_factory.from_dims_halo(
            func=updraft_moist_static_energy_and_momentum_budget,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "USE_LINEAR_SUBCLOUD_MOISTURE_FLUXES": cumulus_parameterization_config.USE_LINEAR_SUBCLOUD_MOISTURE_FLUXES,
                "PRESSURE_GRADIENT_CONSTANT": cumulus_parameterization_config.PRESSURE_GRADIENT_CONSTANT,
            },
        )

        # call test code
        if plume_dependent_constants.ENABLE_PLUME == 1:
            if cumulus_parameterization_config.FIRST_GUESS_W == 0:
                code(
                    error_code=state.output.error_code,
                    start_level=locals.start_level,
                    cloud_top_level=state.output.cloud_top_level,
                    p_forced=state.input_output.p_forced,
                    environment_moist_static_energy=locals.environment_moist_static_energy,
                    environment_moist_static_energy_forced=locals.environment_moist_static_energy_forced,
                    environment_moist_static_energy_cloud_levels=locals.environment_moist_static_energy_cloud_levels,
                    environment_moist_static_energy_cloud_levels_forced=locals.environment_moist_static_energy_cloud_levels_forced,
                    environment_saturation_moist_static_energy_cloud_levels=locals.environment_saturation_moist_static_energy_cloud_levels,
                    environment_saturation_moist_static_energy_cloud_levels_forced=locals.environment_saturation_moist_static_energy_cloud_levels_forced,
                    cloud_moist_static_energy=locals.cloud_moist_static_energy,
                    cloud_moist_static_energy_forced=locals.cloud_moist_static_energy_forced,
                    normalized_massflux_updraft=locals.normalized_massflux_updraft,
                    normalized_massflux_updraft_forced=locals.normalized_massflux_updraft_forced,
                    mass_entrainment_updraft=locals.mass_entrainment_updraft,
                    mass_detrainment_updraft=locals.mass_detrainment_updraft,
                    mass_entrainment_u_updraft=locals.mass_entrainment_u_updraft,
                    mass_detrainment_u_updraft=locals.mass_detrainment_u_updraft,
                    mass_detrainment_updraft_forced=state.output.mass_detrainment_updraft_forced,
                    mass_entrainment_updraft_forced=state.output.mass_entrainment_updraft_forced,
                    u=state.input_output.u,
                    v=state.input_output.v,
                    u_c=locals.u_c,
                    v_c=locals.v_c,
                    u_cloud_levels=locals.u_cloud_levels,
                    v_cloud_levels=locals.v_cloud_levels,
                    partition_liquid_ice=locals.partition_liquid_ice,
                    cloud_liquid_after_rain_forced=state.output.cloud_liquid_after_rain_forced,
                    vapor_excess=locals.vapor_excess,
                    t_excess=locals.t_excess,
                    add_buoyancy=locals.add_buoyancy,
                    plume=plume_dependent_constants.PLUME_INDEX,
                )

        # write output
        outputs = {
            "error_code_upmsemomentumo": state.output.error_code.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_start_level_upmsemomentumo": locals.start_level.field[:] + 1,
            "cloud_top_level_upmsemomentumo": state.output.cloud_top_level.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ]
            + 1,
            "p_forced_upmsemomentumo": state.input_output.p_forced.field[:],
            "local_env_moist_static_energy_upmsemomentumo": locals.environment_moist_static_energy.field[:],
            "local_env_moist_static_energy_forced_upmsemomentumo": locals.environment_moist_static_energy_forced.field[
                :
            ],
            "local_env_moist_static_energy_cloud_levels_upmsemomentumo": locals.environment_moist_static_energy_cloud_levels.field[
                :
            ],
            "local_env_moist_static_energy_cloud_levels_forced_upmsemomentumo": locals.environment_moist_static_energy_cloud_levels_forced.field[
                :
            ],
            "local_env_saturation_moist_static_energy_cloud_levels_upmsemomentumo": locals.environment_saturation_moist_static_energy_cloud_levels.field[
                :
            ],
            "local_env_saturation_moist_static_energy_cloud_levels_forced_upmsemomentumo": locals.environment_saturation_moist_static_energy_cloud_levels_forced.field[
                :
            ],
            "local_cloud_moist_static_energy_upmsemomentumo": locals.cloud_moist_static_energy.field[:],
            "local_cloud_moist_static_energy_forced_upmsemomentumo": locals.cloud_moist_static_energy_forced.field[
                :
            ],
            "local_normalized_massflux_updraft_upmsemomentumo": locals.normalized_massflux_updraft.field[:],
            "local_normalized_massflux_updraft_forced_upmsemomentumo": locals.normalized_massflux_updraft_forced.field[
                :
            ],
            "local_mass_entrainment_updraft_upmsemomentumo": locals.mass_entrainment_updraft.field[:],
            "local_mass_detrainment_updraft_upmsemomentumo": locals.mass_detrainment_updraft.field[:],
            "local_mass_entrainment_u_updraft_upmsemomentumo": locals.mass_entrainment_u_updraft.field[:],
            "local_mass_detrainment_u_updraft_upmsemomentumo": locals.mass_detrainment_u_updraft.field[:],
            "mass_detrainment_updraft_forced_upmsemomentumo": state.output.mass_detrainment_updraft_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "mass_entrainment_updraft_forced_upmsemomentumo": state.output.mass_entrainment_updraft_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "u_upmsemomentumo": state.input_output.u.field[:],
            "v_upmsemomentumo": state.input_output.v.field[:],
            "local_u_c_upmsemomentumo": locals.u_c.field[:],
            "local_v_c_upmsemomentumo": locals.v_c.field[:],
            "local_u_cloud_levels_upmsemomentumo": locals.u_cloud_levels.field[:],
            "local_v_cloud_levels_upmsemomentumo": locals.v_cloud_levels.field[:],
            "local_partition_liquid_ice_upmsemomentumo": locals.partition_liquid_ice.field[:],
            "cloud_liquid_after_rain_forced_upmsemomentumo": state.output.cloud_liquid_after_rain_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_vapor_excess_upmsemomentumo": locals.vapor_excess.field[:],
            "local_t_excess_upmsemomentumo": locals.t_excess.field[:],
            "local_add_buoyancy_upmsemomentumo": locals.add_buoyancy.field[:],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_UpdraftMoistStaticEnergyAndMomentumBudget_shallow(
    TranslateFortranData2Py
):
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


class TranslateGF2020_CumulusParameterization_UpdraftMoistStaticEnergyAndMomentumBudget_mid(
    TranslateFortranData2Py
):
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


class TranslateGF2020_CumulusParameterization_UpdraftMoistStaticEnergyAndMomentumBudget_deep(
    TranslateFortranData2Py
):
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
