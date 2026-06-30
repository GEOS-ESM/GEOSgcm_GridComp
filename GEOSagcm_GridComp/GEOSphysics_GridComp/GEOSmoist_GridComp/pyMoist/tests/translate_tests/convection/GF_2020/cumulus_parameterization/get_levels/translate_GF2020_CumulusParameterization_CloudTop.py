from f90nml import Namelist
from ndsl import StencilFactory
from ndsl.constants import I_DIM, J_DIM, K_DIM
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py

from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.constants import MAXENS1, MAXENS2, MAXENS3, NUMBER_OF_PLUMES
from pyMoist.convection.GF_2020.cumulus_parameterization.get_levels import cloud_top_checks, get_cloud_top
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
            "error_code": {},
            "cloud_top_level": {},
            "entrainment_rate": {},
            "local_moist_static_energy_origin_level_forced": {},
            "local_env_moist_static_energy_forced": {},
            "local_env_saturation_moist_static_energy_cloud_levels_forced": {},
            "local_geopotential_height_cloud_levels_forced": {},
            "updraft_lfc_level": {},
            "local_cloud_moist_static_energy_forced_transported": {},
            "p_cloud_levels_forced": {},
            "local_geopotential_height_cloud_levels": {},
            "last_error_code": {},
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
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["error_code"]
        state.output.cloud_top_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["cloud_top_level"]
        state.output.entrainment_rate.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs["entrainment_rate"]
        locals.moist_static_energy_origin_level_forced.data[:] = inputs["local_moist_static_energy_origin_level_forced"]
        locals.environment_moist_static_energy_forced.data[:] = inputs["local_env_moist_static_energy_forced"]
        locals.environment_saturation_moist_static_energy_cloud_levels_forced.data[:] = inputs["local_env_saturation_moist_static_energy_cloud_levels_forced"]
        locals.geopotential_height_cloud_levels_forced.data[:] = inputs["local_geopotential_height_cloud_levels_forced"]
        state.output.updraft_lfc_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["updraft_lfc_level"] - 1
        locals.cloud_moist_static_energy_forced_transported.data[:] = inputs["local_cloud_moist_static_energy_forced_transported"]
        state.output.p_cloud_levels_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs["p_cloud_levels_forced"]
        locals.geopotential_height_cloud_levels.data[:] = inputs["local_geopotential_height_cloud_levels"]
        state.input.last_error_code.data[:] = inputs["last_error_code"]

        code_part_1 = self.stencil_factory.from_dims_halo(
            func=get_cloud_top,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={"OVERSHOOT": cumulus_parameterization_config.OVERSHOOT},
        )

        code_part_2 = self.stencil_factory.from_dims_halo(
            func=cloud_top_checks,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )

        if plume_dependent_constants.ENABLE_PLUME == 1:
            code_part_1(
                entrainment_rate=state.output.entrainment_rate,
                environment_moist_static_energy_forced=locals.environment_moist_static_energy_forced,
                environment_saturation_moist_static_energy_cloud_levels_forced=locals.environment_saturation_moist_static_energy_cloud_levels_forced,
                moist_static_energy_origin_level_forced=locals.moist_static_energy_origin_level_forced,
                updraft_lfc_level=state.output.updraft_lfc_level,
                geopotential_height_cloud_levels_forced=locals.geopotential_height_cloud_levels_forced,
                cloud_moist_static_energy_forced_transported=locals.cloud_moist_static_energy_forced_transported,
                error_code=state.output.error_code,
                cloud_top_level=state.output.cloud_top_level,
                plume=plume_dependent_constants.PLUME_INDEX,
            )

            code_part_2(
                cloud_top_level=state.output.cloud_top_level,
                p_cloud_levels_forced=state.output.p_cloud_levels_forced,
                geopotential_height_cloud_levels=locals.geopotential_height_cloud_levels,
                error_code=state.output.error_code,
                last_error_code=state.input.last_error_code,
                updraft_lfc_level=state.output.updraft_lfc_level,
                MINIMUM_DEPTH=plume_dependent_constants.MINIMUM_DEPTH,
                plume=plume_dependent_constants.PLUME_INDEX,
            )

            state.output.cloud_top_level.field[:] = state.output.cloud_top_level.field[:] + 1

        outputs = {
            "error_code": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "cloud_top_level": state.output.cloud_top_level.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "entrainment_rate": state.output.entrainment_rate.field[:, :, :, plume_dependent_constants.PLUME_INDEX],
            "local_moist_static_energy_origin_level_forced": locals.moist_static_energy_origin_level_forced.field[:],
            "local_env_moist_static_energy_forced": locals.environment_moist_static_energy_forced.field[:],
            "local_env_saturation_moist_static_energy_cloud_levels_forced": locals.environment_saturation_moist_static_energy_cloud_levels_forced.field[:],
            "local_geopotential_height_cloud_levels_forced": locals.geopotential_height_cloud_levels_forced.field[:],
            "updraft_lfc_level": state.output.updraft_lfc_level.field[:, :, plume_dependent_constants.PLUME_INDEX] + 1,
            "local_cloud_moist_static_energy_forced_transported": locals.cloud_moist_static_energy_forced_transported.field[:],
            "p_cloud_levels_forced": state.output.p_cloud_levels_forced.field[:, :, :, plume_dependent_constants.PLUME_INDEX],
            "local_geopotential_height_cloud_levels": locals.geopotential_height_cloud_levels.field[:],
            "last_error_code": state.input.last_error_code.field[:],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_CloudTop_shallow(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_CloudTop_mid(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_CloudTop_deep(TranslateFortranData2Py):
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
