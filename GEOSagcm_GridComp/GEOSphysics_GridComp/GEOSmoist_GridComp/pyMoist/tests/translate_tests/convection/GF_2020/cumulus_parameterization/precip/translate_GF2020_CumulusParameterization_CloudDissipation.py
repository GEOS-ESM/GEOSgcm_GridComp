from f90nml import Namelist
from ndsl import StencilFactory
from ndsl.constants import I_DIM, J_DIM, K_DIM
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py

from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.constants import MAXENS1, MAXENS2, MAXENS3, NUMBER_OF_PLUMES
from pyMoist.convection.GF_2020.cumulus_parameterization.locals import GF2020CumulusParameterizationLocals
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import GF2020PlumeDependentConstants
from pyMoist.convection.GF_2020.cumulus_parameterization.precip import cloud_dissipation
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
            "updraft_lfc_level": {},
            "cloud_top_level": {},
            "local_hydrostatic_air_density": {},
            "geopotential_height_forced": {},
            "local_t_cloud_levels_forced": {},
            "normalized_massflux_updraft_forced": {},
            "cloud_base_mass_flux_modified": {},
            "local_vapor_cloud_levels_forced": {},
            "cloud_liquid_after_rain_forced": {},
            "local_env_saturation_mixing_ratio_cloud_levels_forced": {},
            "local_env_saturation_moist_static_energy_cloud_levels_forced": {},
            "local_vertical_velocity_3d": {},
            "scale_dependence_factor": {},
            "dtdt": {},
            "dvapordt": {},
            "dcloudicedt": {},
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
        state.output.updraft_lfc_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["updraft_lfc_level"] - 1
        state.output.cloud_top_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["cloud_top_level"] - 1
        locals.hydrostatic_air_density.data[:] = inputs["local_hydrostatic_air_density"]
        state.input_output.geopotential_height_forced.data[:] = inputs["geopotential_height_forced"]
        locals.t_cloud_levels_forced.data[:] = inputs["local_t_cloud_levels_forced"]
        state.output.normalized_massflux_updraft_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs["normalized_massflux_updraft_forced"]
        state.output.cloud_base_mass_flux_modified.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["cloud_base_mass_flux_modified"]
        locals.vapor_cloud_levels_forced.data[:] = inputs["local_vapor_cloud_levels_forced"]
        state.output.cloud_liquid_after_rain_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs["cloud_liquid_after_rain_forced"]
        locals.environment_saturation_mixing_ratio_cloud_levels_forced.data[:] = inputs["local_env_saturation_mixing_ratio_cloud_levels_forced"]
        locals.environment_saturation_moist_static_energy_cloud_levels_forced.data[:] = inputs["local_env_saturation_moist_static_energy_cloud_levels_forced"]
        locals.vertical_velocity_3d.data[:] = inputs["local_vertical_velocity_3d"]
        state.output.scale_dependence_factor.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["scale_dependence_factor"]
        state.output.dtdt.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs["dtdt"]
        state.output.dvapordt.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs["dvapordt"]
        state.output.dcloudicedt.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs["dcloudicedt"]

        code = self.stencil_factory.from_dims_halo(
            func=cloud_dissipation,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "USE_CLOUD_DISSIPATION": cumulus_parameterization_config.USE_CLOUD_DISSIPATION,
                "DT_MOIST": config.DT_MOIST,
            },
        )

        if plume_dependent_constants.ENABLE_PLUME == 1:
            code(
                error_code=state.output.error_code,
                updraft_lfc_level=state.output.updraft_lfc_level,
                cloud_top_level=state.output.cloud_top_level,
                hydrostatic_air_density=locals.hydrostatic_air_density,
                geopotential_height_forced=state.input_output.geopotential_height_forced,
                t_cloud_levels_forced=locals.t_cloud_levels_forced,
                normalized_massflux_updraft_forced=state.output.normalized_massflux_updraft_forced,
                cloud_base_mass_flux_modified=state.output.cloud_base_mass_flux_modified,
                vapor_cloud_levels_forced=locals.vapor_cloud_levels_forced,
                cloud_liquid_after_rain_forced=state.output.cloud_liquid_after_rain_forced,
                environment_saturation_mixing_ratio_cloud_levels_forced=locals.environment_saturation_mixing_ratio_cloud_levels_forced,
                environment_saturation_moist_static_energy_cloud_levels_forced=locals.environment_saturation_moist_static_energy_cloud_levels_forced,
                vertical_velocity_3d=locals.vertical_velocity_3d,
                scale_dependence_factor=state.output.scale_dependence_factor,
                dtdt=state.output.dtdt,
                dvapordt=state.output.dvapordt,
                dcloudicedt=state.output.dcloudicedt,
                plume=plume_dependent_constants.PLUME_INDEX,
            )

        outputs = {
            "error_code": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "updraft_lfc_level": state.output.updraft_lfc_level.field[:, :, plume_dependent_constants.PLUME_INDEX] + 1,
            "cloud_top_level": state.output.cloud_top_level.field[:, :, plume_dependent_constants.PLUME_INDEX] + 1,
            "local_hydrostatic_air_density": locals.hydrostatic_air_density.field[:],
            "geopotential_height_forced": state.input_output.geopotential_height_forced.field[:],
            "local_t_cloud_levels_forced": locals.t_cloud_levels_forced.field[:],
            "normalized_massflux_updraft_forced": state.output.normalized_massflux_updraft_forced.field[:, :, :, plume_dependent_constants.PLUME_INDEX],
            "cloud_base_mass_flux_modified": state.output.cloud_base_mass_flux_modified.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "local_vapor_cloud_levels_forced": locals.vapor_cloud_levels_forced.field[:],
            "cloud_liquid_after_rain_forced": state.output.cloud_liquid_after_rain_forced.field[:, :, :, plume_dependent_constants.PLUME_INDEX],
            "local_env_saturation_mixing_ratio_cloud_levels_forced": locals.environment_saturation_mixing_ratio_cloud_levels_forced.field[:],
            "local_env_saturation_moist_static_energy_cloud_levels_forced": locals.environment_saturation_moist_static_energy_cloud_levels_forced.field[:],
            "local_vertical_velocity_3d": locals.vertical_velocity_3d.field[:],
            "scale_dependence_factor": state.output.scale_dependence_factor.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "dtdt": state.output.dtdt.field[:, :, :, plume_dependent_constants.PLUME_INDEX],
            "dvapordt": state.output.dvapordt.field[:, :, :, plume_dependent_constants.PLUME_INDEX],
            "dcloudicedt": state.output.dcloudicedt.field[:, :, :, plume_dependent_constants.PLUME_INDEX],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_CloudDissipation_shallow(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_CloudDissipation_mid(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_CloudDissipation_deep(TranslateFortranData2Py):
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
