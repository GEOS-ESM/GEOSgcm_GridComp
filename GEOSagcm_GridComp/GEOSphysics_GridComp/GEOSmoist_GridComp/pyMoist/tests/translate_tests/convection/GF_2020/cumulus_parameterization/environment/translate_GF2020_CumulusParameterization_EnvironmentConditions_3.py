from f90nml import Namelist
from ndsl import StencilFactory
from ndsl.constants import I_DIM, J_DIM, K_DIM
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py

from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.constants import MAXENS1, MAXENS2, MAXENS3, NUMBER_OF_PLUMES
from pyMoist.convection.GF_2020.cumulus_parameterization.environment import environment_conditions
from pyMoist.convection.GF_2020.cumulus_parameterization.locals import GF2020CumulusParameterizationLocals
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import GF2020PlumeDependentConstants
from pyMoist.convection.GF_2020.cumulus_parameterization.setup.set_constants import set_constants
from pyMoist.convection.GF_2020.cumulus_parameterization.state import GF2020CumulusParameterizationState


class TestCore:
    def __init__(
        self,
        grid: Grid,
        stencil_factory: StencilFactory,
        in_vars: dict,
        out_vars: dict,
    ):
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        in_vars["data_vars"] = {
            "local_geopotential_height_modified": {},
            "local_env_saturation_mixing_ratio_modified": {},
            "local_env_moist_static_energy_modified": {},
            "local_env_saturation_moist_static_energy_modified": {},
            "local_t_modified": {},
            "local_vapor_modified": {},
            "p_forced": {},
            "topography_height_no_negative": {},
            "p_surface": {},
            "error_code": {},
        }

        out_vars.update(in_vars["data_vars"])

    def __call__(self, constants: dict, cu_param_constants: dict, plume: str, **inputs):
        # initialize constants
        config = GF2020Config(**constants)
        cumulus_parameterization_config = GF2020CumulusParameterizationConfig(**cu_param_constants)
        plume_dependent_constants = GF2020PlumeDependentConstants()
        plume_dependent_constants = set_constants(cumulus_parameterization_config, plume_dependent_constants, plume)

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
        locals.geopotential_height_modified.data[:] = inputs["local_geopotential_height_modified"]
        locals.environment_saturation_mixing_ratio_modified.data[:] = inputs["local_env_saturation_mixing_ratio_modified"]
        locals.environment_moist_static_energy_modified.data[:] = inputs["local_env_moist_static_energy_modified"]
        locals.environment_saturation_moist_static_energy_modified.data[:] = inputs["local_env_saturation_moist_static_energy_modified"]
        locals.t_modified.data[:] = inputs["local_t_modified"]
        locals.vapor_modified.data[:] = inputs["local_vapor_modified"]
        state.input_output.p_forced.data[:] = inputs["p_forced"]
        state.input_output.topography_height_no_negative.data[:] = inputs["topography_height_no_negative"]
        state.input_output.p_surface.data[:] = inputs["p_surface"]
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["error_code"]

        code = self.stencil_factory.from_dims_halo(
            func=environment_conditions,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={"SATURATION_CALCULATION_CHOICE": cumulus_parameterization_config.SATURATION_CALCULATION_CHOICE},
        )

        if plume_dependent_constants.ENABLE_PLUME == 1:
            code(
                p=state.input_output.p_forced,
                p_surface=state.input_output.p_surface,
                t=locals.t_modified,
                vapor=locals.vapor_modified,
                topography_height_no_negative=state.input_output.topography_height_no_negative,
                moist_static_energy=locals.environment_moist_static_energy_modified,
                saturation_moist_static_energy=locals.environment_saturation_moist_static_energy_modified,
                saturation_mixing_ratio=locals.environment_saturation_mixing_ratio_modified,
                geopotential_height=locals.geopotential_height_modified,
                error_code=state.output.error_code,
                plume=plume_dependent_constants.PLUME_INDEX,
            )

        outputs = {
            "local_geopotential_height_modified": locals.geopotential_height_modified.field[:],
            "local_env_saturation_mixing_ratio_modified": locals.environment_saturation_mixing_ratio_modified.field[:],
            "local_env_moist_static_energy_modified": locals.environment_moist_static_energy_modified.field[:],
            "local_env_saturation_moist_static_energy_modified": locals.environment_saturation_moist_static_energy_modified.field[:],
            "local_t_modified": locals.t_modified.field[:],
            "local_vapor_modified": locals.vapor_modified.field[:],
            "p_forced": state.input_output.p_forced.field[:],
            "topography_height_no_negative": state.input_output.topography_height_no_negative.field[:],
            "p_surface": state.input_output.p_surface.field[:],
            "error_code": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_EnvironmentConditions_3_shallow(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        _namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)

        self.test_core = TestCore(grid, stencil_factory, self.in_vars, self.out_vars)

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF2020-constants")
        self.cu_param_constants = data_loader.load("GF2020_CumulusParameterization-constants")

    def compute_func(self, **inputs):
        outputs = self.test_core(self.constants, self.cu_param_constants, "shallow", **inputs)

        return outputs


class TranslateGF2020_CumulusParameterization_EnvironmentConditions_3_mid(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        _namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)

        self.test_core = TestCore(grid, stencil_factory, self.in_vars, self.out_vars)

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF2020-constants")
        self.cu_param_constants = data_loader.load("GF2020_CumulusParameterization-constants")

    def compute_func(self, **inputs):
        outputs = self.test_core(self.constants, self.cu_param_constants, "mid", **inputs)

        return outputs


class TranslateGF2020_CumulusParameterization_EnvironmentConditions_3_deep(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        _namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)

        self.test_core = TestCore(grid, stencil_factory, self.in_vars, self.out_vars)

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF2020-constants")
        self.cu_param_constants = data_loader.load("GF2020_CumulusParameterization-constants")

    def compute_func(self, **inputs):
        outputs = self.test_core(self.constants, self.cu_param_constants, "deep", **inputs)

        return outputs
