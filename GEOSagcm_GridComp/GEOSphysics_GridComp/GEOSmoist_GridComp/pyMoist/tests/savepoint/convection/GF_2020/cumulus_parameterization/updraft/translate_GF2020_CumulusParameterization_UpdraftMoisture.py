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
from pyMoist.convection.GF_2020.cumulus_parameterization.updraft import updraft_moisture
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
            "local_start_level": {},
            "error_code": {},
            "local_geopotential_height_cloud_levels_forced": {},
            "local_cloud_total_water_after_entrainment_forced": {},
            "cloud_liquid_after_rain_forced": {},
            "condensate_to_fall_forced": {},
            "total_normalized_integrated_condensate_forced": {},
            "local_cloud_moist_static_energy_forced": {},
            "local_updraft_column_temperature_forced": {},
            "ocean_fraction": {},
            "convection_fraction": {},
            "surface_type": {},
            "p_forced": {},
            "cloud_top_level": {},
            "local_d_buoyancy_forced": {},
            "local_cloud_liquid_before_rain_forced": {},
            "local_t_cloud_levels": {},
            "local_vapor_forced": {},
            "local_gamma_cloud_levels_forced": {},
            "normalized_massflux_updraft_forced": {},
            "local_env_saturation_mixing_ratio_cloud_levels_forced": {},
            "updraft_origin_level": {},
            "local_vapor_cloud_levels_forced": {},
            "local_vapor_excess": {},
            "ccn": {},
            "local_mass_entrainment_updraft": {},
            "local_mass_detrainment_updraft": {},
            "local_psum": {},
            "local_psumh": {},
            "local_c1d": {},
            "local_add_buoyancy": {},
            "local_vertical_velocity_3d": {},
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
        locals.start_level.data[:] = inputs["local_start_level"] - 1
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["error_code"]
        locals.geopotential_height_cloud_levels_forced.data[:] = inputs[
            "local_geopotential_height_cloud_levels_forced"
        ]
        locals.cloud_total_water_after_entrainment_forced.data[:] = inputs[
            "local_cloud_total_water_after_entrainment_forced"
        ]
        state.output.cloud_liquid_after_rain_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["cloud_liquid_after_rain_forced"]
        )
        state.output.condensate_to_fall_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "condensate_to_fall_forced"
        ]
        state.output.total_normalized_integrated_condensate_forced.data[
            :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["total_normalized_integrated_condensate_forced"]
        locals.cloud_moist_static_energy_forced.data[:] = inputs["local_cloud_moist_static_energy_forced"]
        locals.updraft_column_temperature_forced.data[:] = inputs["local_updraft_column_temperature_forced"]
        state.input.ocean_fraction.data[:] = inputs["ocean_fraction"]
        state.input.convection_fraction.data[:] = inputs["convection_fraction"]
        state.input.surface_type.data[:] = inputs["surface_type"]
        state.input_output.p_forced.data[:] = inputs["p_forced"]
        state.output.cloud_top_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["cloud_top_level"] - 1
        )
        locals.d_buoyancy_forced.data[:] = inputs["local_d_buoyancy_forced"]
        locals.cloud_liquid_before_rain_forced.data[:] = inputs["local_cloud_liquid_before_rain_forced"]
        locals.t_cloud_levels.data[:] = inputs["local_t_cloud_levels"]
        locals.vapor_forced.data[:] = inputs["local_vapor_forced"]
        locals.gamma_cloud_levels_forced.data[:] = inputs["local_gamma_cloud_levels_forced"]
        state.output.normalized_massflux_updraft_forced.data[
            :, :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["normalized_massflux_updraft_forced"]
        locals.environment_saturation_mixing_ratio_cloud_levels_forced.data[:] = inputs[
            "local_env_saturation_mixing_ratio_cloud_levels_forced"
        ]
        state.output.updraft_origin_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["updraft_origin_level"] - 1
        )
        locals.vapor_cloud_levels_forced.data[:] = inputs["local_vapor_cloud_levels_forced"]
        locals.vapor_excess.data[:] = inputs["local_vapor_excess"]
        state.input_output.ccn.data[:] = inputs["ccn"]
        locals.mass_entrainment_updraft.data[:] = inputs["local_mass_entrainment_updraft"]
        locals.mass_detrainment_updraft.data[:] = inputs["local_mass_detrainment_updraft"]
        locals.psum.data[:] = inputs["local_psum"]
        locals.psumh.data[:] = inputs["local_psumh"]
        locals.c1d.data[:] = inputs["local_c1d"]
        locals.add_buoyancy.data[:] = inputs["local_add_buoyancy"]
        locals.vertical_velocity_3d.data[:] = inputs["local_vertical_velocity_3d"]

        # initialize test code
        code = self.stencil_factory.from_dims_halo(
            func=updraft_moisture,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "BOUNDARY_CONDITION_METHOD": cumulus_parameterization_config.BOUNDARY_CONDITION_METHOD,
                "USE_LINEAR_SUBCLOUD_MOISTURE_FLUXES": cumulus_parameterization_config.USE_LINEAR_SUBCLOUD_MOISTURE_FLUXES,
                "AUTOCONV": config.AUTOCONV,
                "CRITICAL_MIXING_RATIO_OVER_OCEAN": cumulus_parameterization_config.CRITICAL_MIXING_RATIO_OVER_OCEAN,
                "CRITICAL_MIXING_RATIO_OVER_LAND": cumulus_parameterization_config.CRITICAL_MIXING_RATIO_OVER_LAND,
                "FRAC_MODIS": cumulus_parameterization_config.FRAC_MODIS,
                "ZERO_DIFF": cumulus_parameterization_config.ZERO_DIFF,
            },
        )

        # call test code
        if plume_dependent_constants.ENABLE_PLUME == 1:
            if cumulus_parameterization_config.FIRST_GUESS_W == 0:
                code(
                    start_level=locals.start_level,
                    error_code=state.output.error_code,
                    geopotential_height_cloud_levels_forced=locals.geopotential_height_cloud_levels_forced,
                    cloud_total_water_after_entrainment_forced=locals.cloud_total_water_after_entrainment_forced,
                    cloud_liquid_after_rain_forced=state.output.cloud_liquid_after_rain_forced,
                    condensate_to_fall_forced=state.output.condensate_to_fall_forced,
                    total_normalized_integrated_condensate_forced=state.output.total_normalized_integrated_condensate_forced,
                    cloud_moist_static_energy_forced=locals.cloud_moist_static_energy_forced,
                    updraft_column_temperature_forced=locals.updraft_column_temperature_forced,
                    ocean_fraction=state.input.ocean_fraction,
                    convection_fraction=state.input.convection_fraction,
                    surface_type=state.input.surface_type,
                    p_forced=state.input_output.p_forced,
                    cloud_top_level=state.output.cloud_top_level,
                    d_buoyancy_forced=locals.d_buoyancy_forced,
                    cloud_liquid_before_rain_forced=locals.cloud_liquid_before_rain_forced,
                    t_cloud_levels=locals.t_cloud_levels,
                    vapor_forced=locals.vapor_forced,
                    gamma_cloud_levels_forced=locals.gamma_cloud_levels_forced,
                    normalized_massflux_updraft_forced=state.output.normalized_massflux_updraft_forced,
                    environment_saturation_mixing_ratio_cloud_levels_forced=locals.environment_saturation_mixing_ratio_cloud_levels_forced,
                    updraft_origin_level=state.output.updraft_origin_level,
                    vapor_cloud_levels_forced=locals.vapor_cloud_levels_forced,
                    vapor_excess=locals.vapor_excess,
                    ccn=state.input_output.ccn,
                    mass_entrainment_updraft=locals.mass_entrainment_updraft,
                    mass_detrainment_updraft=locals.mass_detrainment_updraft,
                    psum=locals.psum,
                    psumh=locals.psumh,
                    c1d=locals.c1d,
                    add_buoyancy=locals.add_buoyancy,
                    vertical_velocity_3d=locals.vertical_velocity_3d,
                    C0=plume_dependent_constants.C0,
                    AVERAGE_LAYER_DEPTH=plume_dependent_constants.AVERAGE_LAYER_DEPTH,
                    plume=plume_dependent_constants.PLUME_INDEX,
                )

        # write output
        outputs = {
            "local_start_level": locals.start_level.data[:] + 1,
            "error_code": state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX],
            "local_geopotential_height_cloud_levels_forced": locals.geopotential_height_cloud_levels_forced.data[
                :
            ],
            "local_cloud_total_water_after_entrainment_forced": locals.cloud_total_water_after_entrainment_forced.data[
                :
            ],
            "cloud_liquid_after_rain_forced": state.output.cloud_liquid_after_rain_forced.data[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "condensate_to_fall_forced": state.output.condensate_to_fall_forced.data[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "total_normalized_integrated_condensate_forced": state.output.total_normalized_integrated_condensate_forced.data[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_cloud_moist_static_energy_forced": locals.cloud_moist_static_energy_forced.data[:],
            "local_updraft_column_temperature_forced": locals.updraft_column_temperature_forced.data[:],
            "ocean_fraction": state.input.ocean_fraction.data[:],
            "convection_fraction": state.input.convection_fraction.data[:],
            "surface_type": state.input.surface_type.data[:],
            "p_forced": state.input_output.p_forced.data[:],
            "cloud_top_level": state.output.cloud_top_level.data[:, :, plume_dependent_constants.PLUME_INDEX]
            + 1,
            "local_d_buoyancy_forced": locals.d_buoyancy_forced.data[:],
            "local_cloud_liquid_before_rain_forced": locals.cloud_liquid_before_rain_forced.data[:],
            "local_t_cloud_levels": locals.t_cloud_levels.data[:],
            "local_vapor_forced": locals.vapor_forced.data[:],
            "local_gamma_cloud_levels_forced": locals.gamma_cloud_levels_forced.data[:],
            "normalized_massflux_updraft_forced": state.output.normalized_massflux_updraft_forced.data[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_env_saturation_mixing_ratio_cloud_levels_forced": locals.environment_saturation_mixing_ratio_cloud_levels_forced.data[
                :
            ],
            "updraft_origin_level": state.output.updraft_origin_level.data[
                :, :, plume_dependent_constants.PLUME_INDEX
            ]
            + 1,
            "local_vapor_cloud_levels_forced": locals.vapor_cloud_levels_forced.data[:],
            "local_vapor_excess": locals.vapor_excess.data[:],
            "ccn": state.input_output.ccn.data[:],
            "local_mass_entrainment_updraft": locals.mass_entrainment_updraft.data[:],
            "local_mass_detrainment_updraft": locals.mass_detrainment_updraft.data[:],
            "local_psum": locals.psum.data[:],
            "local_psumh": locals.psumh.data[:],
            "local_c1d": locals.c1d.data[:],
            "local_add_buoyancy": locals.add_buoyancy.data[:],
            "local_vertical_velocity_3d": locals.vertical_velocity_3d.data[:],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_UpdraftMoisture_shallow(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_UpdraftMoisture_mid(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_UpdraftMoisture_deep(TranslateFortranData2Py):
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
