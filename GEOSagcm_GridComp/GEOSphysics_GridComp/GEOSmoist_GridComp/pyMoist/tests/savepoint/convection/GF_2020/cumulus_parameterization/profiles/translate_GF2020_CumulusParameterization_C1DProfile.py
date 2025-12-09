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
from pyMoist.convection.GF_2020.cumulus_parameterization.constants import MAXENS1, MAXENS2, MAXENS3
from pyMoist.convection.GF_2020.cumulus_parameterization.profiles.profiles import C1DProfile
from pyMoist.convection.GF_2020.cumulus_parameterization.setup.set_constants import set_constants


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
            "local_start_level_c1dprofile": {},
            "lcl_level_c1dprofile": {},
            "error_code_c1dprofile": {},
            "local_geopotential_height_cloud_levels_forced_c1dprofile": {},
            "local_cloud_vapor_mixing_ratio_forced_c1dprofile": {},
            "cloud_liquid_after_rain_forced_c1dprofile": {},
            "precipitable_water_updraft_forced_c1dprofile": {},
            "total_normalized_integrated_condensate_forced_c1dprofile": {},
            "local_cloud_moist_static_energy_forced_c1dprofile": {},
            "local_unspecifid_temperature_c1dprofile": {},
            "ocean_fraction_c1dprofile": {},
            "convection_fraction_c1dprofile": {},
            "surface_type_c1dprofile": {},
            "p_forced_c1dprofile": {},
            "local_p_cloud_levels_c1dprofile": {},
            "updraft_lfc_level_c1dprofile": {},
            "cloud_top_level_c1dprofile": {},
            "local_detrainment_function_updraft_c1dprofile": {},
            "local_buoyancy_forced_c1dprofile": {},
            "local_cloud_liquid_before_rain_forced_c1dprofile": {},
            "local_t_cloud_levels_c1dprofile": {},
            "local_vapor_forced_c1dprofile": {},
            "local_gamma_cloud_levels_forced_c1dprofile": {},
            "normalized_massflux_updraft_forced_c1dprofile": {},
            "local_env_saturation_mixing_ratio_cloud_levels_forced_c1dprofile": {},
            "updraft_origin_level_c1dprofile": {},
            "local_vapor_cloud_levels_forced_c1dprofile": {},
            "local_vapor_excess_c1dprofile": {},
            "air_density_c1dprofile": {},
            "local_mass_entrainment_updraft_c1dprofile": {},
            "local_mass_detrainment_updraft_c1dprofile": {},
            "local_psum_c1dprofile": {},
            "local_psumh_c1dprofile": {},
            "local_c1d_c1dprofile": {},
            "local_add_buoyancy_c1dprofile": {},
            "local_vertical_velocity_3d_c1dprofile": {},
            "local_vertical_velocity_2d_c1dprofile": {},
            "convective_scale_velocity_c1dprofile": {},
            "entrainment_rate_c1dprofile": {},
            "geopotential_height_forced_c1dprofile": {},
            "local_t_cloud_levels_forced_c1dprofile": {},
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
        locals.start_level.data[:] = inputs["local_start_level_c1dprofile"]
        state.output.lcl_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "lcl_level_c1dprofile"
        ]
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "error_code_c1dprofile"
        ]
        locals.geopotential_height_cloud_levels_forced.data[:] = inputs[
            "local_geopotential_height_cloud_levels_forced_c1dprofile"
        ]
        locals.cloud_vapor_mixing_ratio_forced.data[:] = inputs[
            "local_cloud_vapor_mixing_ratio_forced_c1dprofile"
        ]
        state.output.cloud_liquid_after_rain_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["cloud_liquid_after_rain_forced_c1dprofile"]
        )
        state.output.precipitable_water_updraft_forced.data[
            :, :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["precipitable_water_updraft_forced_c1dprofile"]
        state.output.total_normalized_integrated_condensate_forced.data[
            :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["total_normalized_integrated_condensate_forced_c1dprofile"]
        locals.cloud_moist_static_energy_forced.data[:] = inputs[
            "local_cloud_moist_static_energy_forced_c1dprofile"
        ]
        locals.unspecifid_temperature.data[:] = inputs["local_unspecifid_temperature_c1dprofile"]
        state.input.ocean_fraction.data[:] = inputs["ocean_fraction_c1dprofile"]
        state.input.convection_fraction.data[:] = inputs["convection_fraction_c1dprofile"]
        state.input.surface_type.data[:] = inputs["surface_type_c1dprofile"]
        state.input_output.p_forced.data[:] = inputs["p_forced_c1dprofile"]
        locals.p_cloud_levels.data[:] = inputs["local_p_cloud_levels_c1dprofile"]
        state.output.updraft_lfc_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "updraft_lfc_level_c1dprofile"
        ]
        state.output.cloud_top_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "cloud_top_level_c1dprofile"
        ]
        locals.detrainment_function_updraft.data[:] = inputs["local_detrainment_function_updraft_c1dprofile"]
        locals.buoyancy_forced.data[:] = inputs["local_buoyancy_forced_c1dprofile"]
        locals.cloud_liquid_before_rain_forced.data[:] = inputs[
            "local_cloud_liquid_before_rain_forced_c1dprofile"
        ]
        locals.t_cloud_levels.data[:] = inputs["local_t_cloud_levels_c1dprofile"]
        locals.vapor_forced.data[:] = inputs["local_vapor_forced_c1dprofile"]
        locals.gamma_cloud_levels_forced.data[:] = inputs["local_gamma_cloud_levels_forced_c1dprofile"]
        state.output.normalized_massflux_updraft_forced.data[
            :, :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["normalized_massflux_updraft_forced_c1dprofile"]
        locals.environment_saturation_mixing_ratio_cloud_levels_forced.data[:] = inputs[
            "local_env_saturation_mixing_ratio_cloud_levels_forced_c1dprofile"
        ]
        state.output.updraft_origin_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "updraft_origin_level_c1dprofile"
        ]
        locals.vapor_cloud_levels_forced.data[:] = inputs["local_vapor_cloud_levels_forced_c1dprofile"]
        locals.vapor_excess.data[:] = inputs["local_vapor_excess_c1dprofile"]
        state.input_output.air_density.data[:] = inputs["air_density_c1dprofile"]
        locals.mass_entrainment_updraft.data[:] = inputs["local_mass_entrainment_updraft_c1dprofile"]
        locals.mass_detrainment_updraft.data[:] = inputs["local_mass_detrainment_updraft_c1dprofile"]
        locals.psum.data[:] = inputs["local_psum_c1dprofile"]
        locals.psumh.data[:] = inputs["local_psumh_c1dprofile"]
        locals.c1d.data[:] = inputs["local_c1d_c1dprofile"]
        locals.add_buoyancy.data[:] = inputs["local_add_buoyancy_c1dprofile"]
        locals.vertical_velocity_3d.data[:] = inputs["local_vertical_velocity_3d_c1dprofile"]
        locals.vertical_velocity_2d.data[:] = inputs["local_vertical_velocity_2d_c1dprofile"]
        state.input_output.convective_scale_velocity.data[:] = inputs["convective_scale_velocity_c1dprofile"]
        state.output.entrainment_rate.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "entrainment_rate_c1dprofile"
        ]
        state.input_output.geopotential_height_forced.data[:] = inputs[
            "geopotential_height_forced_c1dprofile"
        ]
        locals.t_cloud_levels_forced.data[:] = inputs["local_t_cloud_levels_forced_c1dprofile"]

        # initalize test code
        code = C1DProfile(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        # call test code
        if plume_dependent_constants.ENABLE_PLUME == 1:
            code(
                state=state,
                locals=locals,
                plume_dependent_constants=plume_dependent_constants,
            )

        # write output
        outputs = {
            "local_start_level_c1dprofile": locals.start_level.field[:],
            "lcl_level_c1dprofile": state.output.lcl_level.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "error_code_c1dprofile": state.output.error_code.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_geopotential_height_cloud_levels_forced_c1dprofile": locals.geopotential_height_cloud_levels_forced.field[
                :
            ],
            "local_cloud_vapor_mixing_ratio_forced_c1dprofile": locals.cloud_vapor_mixing_ratio_forced.field[
                :
            ],
            "cloud_liquid_after_rain_forced_c1dprofile": state.output.cloud_liquid_after_rain_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "precipitable_water_updraft_forced_c1dprofile": state.output.precipitable_water_updraft_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "total_normalized_integrated_condensate_forced_c1dprofile": state.output.total_normalized_integrated_condensate_forced.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_cloud_moist_static_energy_forced_c1dprofile": locals.cloud_moist_static_energy_forced.field[
                :
            ],
            "local_unspecifid_temperature_c1dprofile": locals.unspecifid_temperature.field[:],
            "ocean_fraction_c1dprofile": state.input.ocean_fraction.field[:],
            "convection_fraction_c1dprofile": state.input.convection_fraction.field[:],
            "surface_type_c1dprofile": state.input.surface_type.field[:],
            "p_forced_c1dprofile": state.input_output.p_forced.field[:],
            "local_p_cloud_levels_c1dprofile": locals.p_cloud_levels.field[:],
            "updraft_lfc_level_c1dprofile": state.output.updraft_lfc_level.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "cloud_top_level_c1dprofile": state.output.cloud_top_level.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_detrainment_function_updraft_c1dprofile": locals.detrainment_function_updraft.field[:],
            "local_buoyancy_forced_c1dprofile": locals.buoyancy_forced.field[:],
            "local_cloud_liquid_before_rain_forced_c1dprofile": locals.cloud_liquid_before_rain_forced.field[
                :
            ],
            "local_t_cloud_levels_c1dprofile": locals.t_cloud_levels.field[:],
            "local_vapor_forced_c1dprofile": locals.vapor_forced.field[:],
            "local_gamma_cloud_levels_forced_c1dprofile": locals.gamma_cloud_levels_forced.field[:],
            "normalized_massflux_updraft_forced_c1dprofile": state.output.normalized_massflux_updraft_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_env_saturation_mixing_ratio_cloud_levels_forced_c1dprofile": locals.environment_saturation_mixing_ratio_cloud_levels_forced.field[
                :
            ],
            "updraft_origin_level_c1dprofile": state.output.updraft_origin_level.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_vapor_cloud_levels_forced_c1dprofile": locals.vapor_cloud_levels_forced.field[:],
            "local_vapor_excess_c1dprofile": locals.vapor_excess.field[:],
            "air_density_c1dprofile": state.input_output.air_density.field[:],
            "local_mass_entrainment_updraft_c1dprofile": locals.mass_entrainment_updraft.field[:],
            "local_mass_detrainment_updraft_c1dprofile": locals.mass_detrainment_updraft.field[:],
            "local_psum_c1dprofile": locals.psum.field[:],
            "local_psumh_c1dprofile": locals.psumh.field[:],
            "local_c1d_c1dprofile": locals.c1d.field[:],
            "local_add_buoyancy_c1dprofile": locals.add_buoyancy.field[:],
            "local_vertical_velocity_3d_c1dprofile": locals.vertical_velocity_3d.field[:],
            "local_vertical_velocity_2d_c1dprofile": locals.vertical_velocity_2d.field[:],
            "convective_scale_velocity_c1dprofile": state.input_output.convective_scale_velocity.field[:],
            "entrainment_rate_c1dprofile": state.output.entrainment_rate.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "geopotential_height_forced_c1dprofile": state.input_output.geopotential_height_forced.field[:],
            "local_t_cloud_levels_forced_c1dprofile": locals.t_cloud_levels_forced.field[:],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_C1DProfile_shallow(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_C1DProfile_mid(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_C1DProfile_deep(TranslateFortranData2Py):
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
