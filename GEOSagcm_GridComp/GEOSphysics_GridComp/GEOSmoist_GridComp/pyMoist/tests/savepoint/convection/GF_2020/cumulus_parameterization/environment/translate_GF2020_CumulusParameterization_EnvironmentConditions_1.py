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
    environment_conditions,
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
            "local_geopotential_height": {},
            "local_env_saturation_mixing_ratio": {},
            "local_env_moist_static_energy": {},
            "local_env_saturation_moist_static_energy": {},
            "t_old": {},
            "vapor_old": {},
            "p_forced": {},
            "topography_height_no_negative": {},
            "p_surface": {},
            "error_code": {},
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
        locals.geopotential_height.data[:] = inputs["local_geopotential_height"]
        locals.environment_saturation_mixing_ratio.data[:] = inputs["local_env_saturation_mixing_ratio"]
        locals.environment_moist_static_energy.data[:] = inputs["local_env_moist_static_energy"]
        locals.environment_saturation_moist_static_energy.data[:] = inputs[
            "local_env_saturation_moist_static_energy"
        ]
        state.input_output.t_old.data[:] = inputs["t_old"]
        state.input_output.vapor_old.data[:] = inputs["vapor_old"]
        state.input_output.p_forced.data[:] = inputs["p_forced"]
        state.input_output.topography_height_no_negative.data[:] = inputs["topography_height_no_negative"]
        state.input_output.p_surface.data[:] = inputs["p_surface"]
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["error_code"]

        code = self.stencil_factory.from_dims_halo(
            func=environment_conditions,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "SATURATION_CALCULATION_CHOICE": cumulus_parameterization_config.SATURATION_CALCULATION_CHOICE
            },
        )

        if plume_dependent_constants.ENABLE_PLUME == 1:
            code(
                p=state.input_output.p_forced,
                p_surface=state.input_output.p_surface,
                t=state.input_output.t_old,
                vapor=state.input_output.vapor_old,
                topography_height_no_negative=state.input_output.topography_height_no_negative,
                moist_static_energy=locals.environment_moist_static_energy,
                saturation_moist_static_energy=locals.environment_saturation_moist_static_energy,
                saturation_mixing_ratio=locals.environment_saturation_mixing_ratio,
                geopotential_height=locals.geopotential_height,
                error_code=state.output.error_code,
                plume=plume_dependent_constants.PLUME_INDEX,
            )

        outputs = {
            # state fields
            "local_geopotential_height": locals.geopotential_height.field[:],
            "local_env_saturation_mixing_ratio": locals.environment_saturation_mixing_ratio.field[:],
            "local_env_moist_static_energy": locals.environment_moist_static_energy.field[:],
            "local_env_saturation_moist_static_energy": locals.environment_saturation_moist_static_energy.field[
                :
            ],
            "t_old": state.input_output.t_old.field[:],
            "vapor_old": state.input_output.vapor_old.field[:],
            "p_forced": state.input_output.p_forced.field[:],
            "topography_height_no_negative": state.input_output.topography_height_no_negative.field[:],
            "p_surface": state.input_output.p_surface.field[:],
            "error_code": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_EnvironmentConditions_1_shallow(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_EnvironmentConditions_1_mid(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_EnvironmentConditions_1_deep(TranslateFortranData2Py):
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
