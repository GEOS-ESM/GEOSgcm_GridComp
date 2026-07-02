from f90nml import Namelist
from ndsl import StencilFactory
from ndsl.constants import I_DIM, J_DIM, K_DIM
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py

from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.constants import MAXENS1, MAXENS2, MAXENS3, NUMBER_OF_PLUMES
from pyMoist.convection.GF_2020.cumulus_parameterization.get_levels import get_convective_cloud_base_level, set_start_level
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
            "local_cap_max_increment": {},
            "local_cap_max": {},
            "local_env_moist_static_energy_cloud_levels_forced": {},
            "local_env_saturation_moist_static_energy_cloud_levels_forced": {},
            "local_vapor_cloud_levels_forced": {},
            "local_env_saturation_mixing_ratio_cloud_levels_forced": {},
            "p_forced": {},
            "p_cloud_levels_forced": {},
            "local_geopotential_height_cloud_levels_forced": {},
            "local_env_moist_static_energy_forced": {},
            "local_moist_static_energy_origin_level_forced": {},
            "local_vapor_forced": {},
            "local_env_saturation_mixing_ratio_forced": {},
            "entrainment_rate": {},
            "local_cloud_moist_static_energy_forced_transported": {},
            "updraft_origin_level": {},
            "local_maximum_updraft_origin_level": {},
            "lcl_level": {},
            "updraft_lfc_level": {},
            "cloud_top_level": {},
            "local_negative_buoyancy_depth": {},
            "local_frh_lfc": {},
            "t_perturbation": {},
            "local_start_level": {},
            "local_vapor_excess": {},
            "local_t_excess": {},
            "local_add_buoyancy": {},
            "ocean_fraction": {},
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
        locals.cap_max_increment.data[:] = inputs["local_cap_max_increment"]
        locals.cap_max.data[:] = inputs["local_cap_max"]
        locals.environment_moist_static_energy_cloud_levels_forced.data[:] = inputs["local_env_moist_static_energy_cloud_levels_forced"]
        locals.environment_saturation_moist_static_energy_cloud_levels_forced.data[:] = inputs["local_env_saturation_moist_static_energy_cloud_levels_forced"]
        locals.vapor_cloud_levels_forced.data[:] = inputs["local_vapor_cloud_levels_forced"]
        locals.environment_saturation_mixing_ratio_cloud_levels_forced.data[:] = inputs["local_env_saturation_mixing_ratio_cloud_levels_forced"]
        state.input_output.p_forced.data[:] = inputs["p_forced"]
        state.output.p_cloud_levels_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs["p_cloud_levels_forced"]
        locals.geopotential_height_cloud_levels_forced.data[:] = inputs["local_geopotential_height_cloud_levels_forced"]
        locals.environment_moist_static_energy_forced.data[:] = inputs["local_env_moist_static_energy_forced"]
        locals.moist_static_energy_origin_level_forced.data[:] = inputs["local_moist_static_energy_origin_level_forced"]
        locals.vapor_forced.data[:] = inputs["local_vapor_forced"]
        locals.environment_saturation_mixing_ratio_forced.data[:] = inputs["local_env_saturation_mixing_ratio_forced"]
        state.output.entrainment_rate.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs["entrainment_rate"]
        locals.cloud_moist_static_energy_forced_transported.data[:] = inputs["local_cloud_moist_static_energy_forced_transported"]
        state.output.updraft_origin_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["updraft_origin_level"] - 1
        locals.maximum_updraft_origin_level.data[:] = inputs["local_maximum_updraft_origin_level"] - 1
        state.output.lcl_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["lcl_level"] - 1
        state.output.updraft_lfc_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["updraft_lfc_level"] - 1
        state.output.cloud_top_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["cloud_top_level"] - 1
        locals.negative_buoyancy_depth.data[:] = inputs["local_negative_buoyancy_depth"]
        locals.frh_lfc.data[:] = inputs["local_frh_lfc"]
        state.output.t_perturbation.data[:] = inputs["t_perturbation"]
        locals.start_level.data[:] = inputs["local_start_level"] - 1
        locals.vapor_excess.data[:] = inputs["local_vapor_excess"]
        locals.t_excess.data[:] = inputs["local_t_excess"]
        locals.add_buoyancy.data[:] = inputs["local_add_buoyancy"]
        locals.ocean_fraction.data[:] = inputs["ocean_fraction"]

        code_part_1 = self.stencil_factory.from_dims_halo(
            func=set_start_level,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )
        code_part_2 = self.stencil_factory.from_dims_halo(
            func=get_convective_cloud_base_level,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "OVERSHOOT": cumulus_parameterization_config.OVERSHOOT,
                "ZERO_DIFF": cumulus_parameterization_config.ZERO_DIFF,
                "MOIST_TRIGGER": cumulus_parameterization_config.MOIST_TRIGGER,
                "USE_MEMORY": cumulus_parameterization_config.USE_MEMORY,
                "BOUNDARY_CONDITION_METHOD": cumulus_parameterization_config.BOUNDARY_CONDITION_METHOD,
            },
        )

        if plume_dependent_constants.ENABLE_PLUME == 1:
            code_part_1(
                lcl_level=state.output.lcl_level,
                start_level=locals.start_level,
                plume=plume_dependent_constants.PLUME_INDEX,
            )

            code_part_2(
                error_code=state.output.error_code,
                lcl_level=state.output.lcl_level,
                cloud_moist_static_energy_forced_transported=locals.cloud_moist_static_energy_forced_transported,
                cap_max=locals.cap_max,
                updraft_origin_level=state.output.updraft_origin_level,
                start_level=locals.start_level,
                moist_static_energy_origin_level_forced=locals.moist_static_energy_origin_level_forced,
                updraft_lfc_level=state.output.updraft_lfc_level,
                maximum_updraft_origin_level=locals.maximum_updraft_origin_level,
                negative_buoyancy_depth=locals.negative_buoyancy_depth,
                frh_lfc=locals.frh_lfc,
                geopotential_height_cloud_levels_forced=locals.geopotential_height_cloud_levels_forced,
                entrainment_rate=state.output.entrainment_rate,
                environment_moist_static_energy_forced=locals.environment_moist_static_energy_forced,
                environment_moist_static_energy_cloud_levels_forced=locals.environment_moist_static_energy_cloud_levels_forced,
                environment_saturation_moist_static_energy_cloud_levels_forced=locals.environment_saturation_moist_static_energy_cloud_levels_forced,
                t_excess=locals.t_excess,
                vapor_excess=locals.vapor_excess,
                add_buoyancy=locals.add_buoyancy,
                p_cloud_levels_forced=state.output.p_cloud_levels_forced,
                vapor_forced=locals.vapor_forced,
                environment_saturation_mixing_ratio_forced=locals.environment_saturation_mixing_ratio_forced,
                ocean_fraction=locals.ocean_fraction,
                cap_max_increment=locals.cap_max_increment,
                t_perturbation=state.output.t_perturbation,
                p_forced=state.input_output.p_forced,
                cloud_top_level=state.output.cloud_top_level,
                AVERAGE_LAYER_DEPTH=plume_dependent_constants.AVERAGE_LAYER_DEPTH,
                plume=plume_dependent_constants.PLUME_INDEX,
            )

        outputs = {
            "error_code": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "local_cap_max_increment": locals.cap_max_increment.field[:],
            "local_cap_max": locals.cap_max.field[:],
            "local_env_moist_static_energy_cloud_levels_forced": locals.environment_moist_static_energy_cloud_levels_forced.field[:],
            "local_env_saturation_moist_static_energy_cloud_levels_forced": locals.environment_saturation_moist_static_energy_cloud_levels_forced.field[:],
            "local_vapor_cloud_levels_forced": locals.vapor_cloud_levels_forced.field[:],
            "local_env_saturation_mixing_ratio_cloud_levels_forced": locals.environment_saturation_mixing_ratio_cloud_levels_forced.field[:],
            "p_forced": state.input_output.p_forced.field[:],
            "p_cloud_levels_forced": state.output.p_cloud_levels_forced.field[:, :, :, plume_dependent_constants.PLUME_INDEX],
            "local_geopotential_height_cloud_levels_forced": locals.geopotential_height_cloud_levels_forced.field[:],
            "local_env_moist_static_energy_forced": locals.environment_moist_static_energy_forced.field[:],
            "local_moist_static_energy_origin_level_forced": locals.moist_static_energy_origin_level_forced.field[:],
            "local_vapor_forced": locals.vapor_forced.field[:],
            "local_env_saturation_mixing_ratio_forced": locals.environment_saturation_mixing_ratio_forced.field[:],
            "entrainment_rate": state.output.entrainment_rate.field[:, :, :, plume_dependent_constants.PLUME_INDEX],
            "local_cloud_moist_static_energy_forced_transported": locals.cloud_moist_static_energy_forced_transported.field[:],
            "updraft_origin_level": state.output.updraft_origin_level.field[:, :, plume_dependent_constants.PLUME_INDEX] + 1,
            "local_maximum_updraft_origin_level": locals.maximum_updraft_origin_level.field[:] + 1,
            "lcl_level": state.output.lcl_level.field[:, :, plume_dependent_constants.PLUME_INDEX] + 1,
            "updraft_lfc_level": state.output.updraft_lfc_level.field[:, :, plume_dependent_constants.PLUME_INDEX] + 1,
            "cloud_top_level": state.output.cloud_top_level.field[:, :, plume_dependent_constants.PLUME_INDEX] + 1,
            "local_negative_buoyancy_depth": locals.negative_buoyancy_depth.field[:],
            "local_frh_lfc": locals.frh_lfc.field[:],
            "t_perturbation": state.output.t_perturbation.field[:],
            "local_start_level": locals.start_level.field[:] + 1,
            "local_vapor_excess": locals.vapor_excess.field[:],
            "local_t_excess": locals.t_excess.field[:],
            "local_add_buoyancy": locals.add_buoyancy.field[:],
            "ocean_fraction": locals.ocean_fraction.field[:],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_ConvectiveCloudBaseLevel_shallow(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_ConvectiveCloudBaseLevel_mid(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_ConvectiveCloudBaseLevel_deep(TranslateFortranData2Py):
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
