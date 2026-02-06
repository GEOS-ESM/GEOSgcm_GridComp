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
from pyMoist.convection.GF_2020.cumulus_parameterization.downdraft import (
    downdraft_moist_static_energy_and_buoyancy,
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
            "error_code": {},
            "downdraft_origin_level": {},
            "u": {},
            "local_u_cloud_levels": {},
            "local_u_c_downdraft": {},
            "v": {},
            "local_v_cloud_levels": {},
            "local_v_c_downdraft": {},
            "local_env_moist_static_energy_forced": {},
            "local_env_saturation_moist_static_energy_cloud_levels_forced": {},
            "local_cloud_moist_static_energy": {},
            "local_cloud_moist_static_energy_downdraft_forced": {},
            "local_d_buoyancy_downdraft_forced": {},
            "local_t_wetbulb": {},
            "local_vapor_wetbulb": {},
            "local_geopotential_height_cloud_levels_forced": {},
            "normalized_massflux_downdraft_forced": {},
            "mass_entrainment_downdraft_forced": {},
            "mass_detrainment_downdraft_forced": {},
            "local_mass_entrainment_u_downdraft": {},
            "local_mass_detrainment_u_downdraft": {},
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
        state.output.downdraft_origin_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "downdraft_origin_level"
        ] - 1
        state.input_output.u.data[:] = inputs["u"]
        locals.u_cloud_levels.data[:] = inputs["local_u_cloud_levels"]
        locals.u_c_downdraft.data[:] = inputs["local_u_c_downdraft"]
        state.input_output.v.data[:] = inputs["v"]
        locals.v_cloud_levels.data[:] = inputs["local_v_cloud_levels"]
        locals.v_c_downdraft.data[:] = inputs["local_v_c_downdraft"]
        locals.environment_moist_static_energy_forced.data[:] = inputs["local_env_moist_static_energy_forced"]
        locals.environment_saturation_moist_static_energy_cloud_levels_forced.data[:] = inputs[
            "local_env_saturation_moist_static_energy_cloud_levels_forced"
        ]
        locals.cloud_moist_static_energy.data[:] = inputs["local_cloud_moist_static_energy"]
        locals.cloud_moist_static_energy_downdraft_forced.data[:] = inputs[
            "local_cloud_moist_static_energy_downdraft_forced"
        ]
        locals.d_buoyancy_downdraft_forced.data[:] = inputs["local_d_buoyancy_downdraft_forced"]
        locals.t_wetbulb.data[:] = inputs["local_t_wetbulb"]
        locals.vapor_wetbulb.data[:] = inputs["local_vapor_wetbulb"]
        locals.geopotential_height_cloud_levels_forced.data[:] = inputs[
            "local_geopotential_height_cloud_levels_forced"
        ]
        state.output.normalized_massflux_downdraft_forced.data[
            :, :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["normalized_massflux_downdraft_forced"]
        state.output.mass_entrainment_downdraft_forced.data[
            :, :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["mass_entrainment_downdraft_forced"]
        state.output.mass_detrainment_downdraft_forced.data[
            :, :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["mass_detrainment_downdraft_forced"]
        locals.mass_entrainment_u_downdraft.data[:] = inputs["local_mass_entrainment_u_downdraft"]
        locals.mass_detrainment_u_downdraft.data[:] = inputs["local_mass_detrainment_u_downdraft"]

        # initialize test code
        code = self.stencil_factory.from_dims_halo(
            func=downdraft_moist_static_energy_and_buoyancy,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "USE_WETBULB": cumulus_parameterization_config.USE_WETBULB,
                "PRESSURE_GRADIENT_CONSTANT": cumulus_parameterization_config.PRESSURE_GRADIENT_CONSTANT,
            },
        )

        # call test code
        if plume_dependent_constants.ENABLE_PLUME == 1:
            if cumulus_parameterization_config.FIRST_GUESS_W == 0:
                code(
                    error_code=state.output.error_code,
                    downdraft_origin_level=state.output.downdraft_origin_level,
                    u=state.input_output.u,
                    u_cloud_levels=locals.u_cloud_levels,
                    u_c_downdraft=locals.u_c_downdraft,
                    v=state.input_output.v,
                    v_cloud_levels=locals.v_cloud_levels,
                    v_c_downdraft=locals.v_c_downdraft,
                    environment_moist_static_energy_forced=locals.environment_moist_static_energy_forced,
                    environment_saturation_moist_static_energy_cloud_levels_forced=locals.environment_saturation_moist_static_energy_cloud_levels_forced,
                    cloud_moist_static_energy=locals.cloud_moist_static_energy,
                    cloud_moist_static_energy_downdraft_forced=locals.cloud_moist_static_energy_downdraft_forced,
                    buoyancy_downdraft_forced=locals.d_buoyancy_downdraft_forced,
                    t_wetbulb=locals.t_wetbulb,
                    vapor_wetbulb=locals.vapor_wetbulb,
                    geopotential_height_cloud_levels_forced=locals.geopotential_height_cloud_levels_forced,
                    normalized_massflux_downdraft_forced=state.output.normalized_massflux_downdraft_forced,
                    mass_entrainment_downdraft_forced=state.output.mass_entrainment_downdraft_forced,
                    mass_detrainment_downdraft_forced=state.output.mass_detrainment_downdraft_forced,
                    mass_entrainment_u_downdraft=locals.mass_entrainment_u_downdraft,
                    mass_detrainment_u_downdraft=locals.mass_detrainment_u_downdraft,
                    plume=plume_dependent_constants.PLUME_INDEX,
                )

        # write output
        outputs = {
            "error_code": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "downdraft_origin_level": state.output.downdraft_origin_level.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ] + 1,
            "u": state.input_output.u.field[:],
            "local_u_cloud_levels": locals.u_cloud_levels.field[:],
            "local_u_c_downdraft": locals.u_c_downdraft.field[:],
            "v": state.input_output.v.field[:],
            "local_v_cloud_levels": locals.v_cloud_levels.field[:],
            "local_v_c_downdraft": locals.v_c_downdraft.field[:],
            "local_env_moist_static_energy_forced": locals.environment_moist_static_energy_forced.field[:],
            "local_env_saturation_moist_static_energy_cloud_levels_forced": locals.environment_saturation_moist_static_energy_cloud_levels_forced.field[
                :
            ],
            "local_cloud_moist_static_energy": locals.cloud_moist_static_energy.field[:],
            "local_cloud_moist_static_energy_downdraft_forced": locals.cloud_moist_static_energy_downdraft_forced.field[
                :
            ],
            "local_d_buoyancy_downdraft_forced": locals.d_buoyancy_downdraft_forced.field[:],
            "local_t_wetbulb": locals.t_wetbulb.field[:],
            "local_vapor_wetbulb": locals.vapor_wetbulb.field[:],
            "local_geopotential_height_cloud_levels_forced": locals.geopotential_height_cloud_levels_forced.field[
                :
            ],
            "normalized_massflux_downdraft_forced": state.output.normalized_massflux_downdraft_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "mass_entrainment_downdraft_forced": state.output.mass_entrainment_downdraft_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "mass_detrainment_downdraft_forced": state.output.mass_detrainment_downdraft_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_mass_entrainment_u_downdraft": locals.mass_entrainment_u_downdraft.field[:],
            "local_mass_detrainment_u_downdraft": locals.mass_detrainment_u_downdraft.field[:],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_DowndraftMSEAndBuoyancy_shallow(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_DowndraftMSEAndBuoyancy_mid(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_DowndraftMSEAndBuoyancy_deep(TranslateFortranData2Py):
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
