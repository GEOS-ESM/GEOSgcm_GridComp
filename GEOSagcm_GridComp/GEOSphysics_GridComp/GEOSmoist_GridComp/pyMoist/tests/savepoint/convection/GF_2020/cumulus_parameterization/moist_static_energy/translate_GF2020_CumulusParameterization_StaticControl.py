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
from pyMoist.convection.GF_2020.cumulus_parameterization.moist_static_energy import StaticControl
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
            "error_code": {},
            "local_start_level": {},
            "lcl_level": {},
            "updraft_lfc_level": {},
            "cloud_top_level": {},
            "local_cloud_moist_static_energy_modified": {},
            "local_moist_static_energy_origin_level_modified": {},
            "local_env_saturation_moist_static_energy_modified": {},
            "local_env_moist_static_energy_cloud_levels_modified": {},
            "local_env_saturation_moist_static_energy_cloud_levels_modified": {},
            "mass_detrainment_updraft_forced": {},
            "mass_entrainment_updraft_forced": {},
            "local_normalized_massflux_updraft_modified": {},
            "local_partition_liquid_ice": {},
            "local_vapor_excess": {},
            "local_t_excess": {},
            "local_add_buoyancy": {},
            "cloud_liquid_after_rain_forced": {},
            "local_d_buoyancy_modified": {},
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
                "tracers": config.NUMBER_OF_TRACERS,
            },
        )

        locals = GF2020CumulusParameterizationLocals.zeros(
            self.quantity_factory,
            data_dimensions={
                "ensemble_1": MAXENS1,
                "ensemble_2": MAXENS2,
                "ensemble_3": MAXENS3,
                "ensemble_members": MAXENS1 * MAXENS2 * MAXENS3,
                "tracers": config.NUMBER_OF_TRACERS,
            },
        )

        # fill relevant parts of dataclasses
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["error_code"]
        locals.start_level.data[:] = inputs["local_start_level"] - 1
        state.output.lcl_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["lcl_level"] - 1
        state.output.updraft_lfc_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["updraft_lfc_level"] - 1
        )
        state.output.cloud_top_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["cloud_top_level"] - 1
        )
        locals.cloud_moist_static_energy_modified.data[:] = inputs["local_cloud_moist_static_energy_modified"]
        locals.moist_static_energy_origin_level_modified.data[:] = inputs[
            "local_moist_static_energy_origin_level_modified"
        ]
        locals.environment_saturation_moist_static_energy_modified.data[:] = inputs[
            "local_env_saturation_moist_static_energy_modified"
        ]
        locals.environment_moist_static_energy_cloud_levels_modified.data[:] = inputs[
            "local_env_moist_static_energy_cloud_levels_modified"
        ]
        locals.environment_saturation_moist_static_energy_cloud_levels_modified.data[:] = inputs[
            "local_env_saturation_moist_static_energy_cloud_levels_modified"
        ]
        state.output.mass_detrainment_updraft_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["mass_detrainment_updraft_forced"]
        )
        state.output.mass_entrainment_updraft_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["mass_entrainment_updraft_forced"]
        )
        locals.normalized_massflux_updraft_modified.data[:] = inputs[
            "local_normalized_massflux_updraft_modified"
        ]
        locals.partition_liquid_ice.data[:] = inputs["local_partition_liquid_ice"]
        locals.vapor_excess.data[:] = inputs["local_vapor_excess"]
        locals.t_excess.data[:] = inputs["local_t_excess"]
        locals.add_buoyancy.data[:] = inputs["local_add_buoyancy"]
        state.output.cloud_liquid_after_rain_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["cloud_liquid_after_rain_forced"]
        )
        locals.d_buoyancy_modified.data[:] = inputs["local_d_buoyancy_modified"]

        # initalize test code
        code = StaticControl(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        # call test code
        if plume_dependent_constants.ENABLE_PLUME == 1:
            code(
                error_code=state.output.error_code,
                start_level=locals.start_level,
                lcl_level=state.output.lcl_level,
                updraft_lfc_level=state.output.updraft_lfc_level,
                cloud_top_level=state.output.cloud_top_level,
                cloud_moist_static_energy_modified=locals.cloud_moist_static_energy_modified,
                moist_static_energy_origin_level_modified=locals.moist_static_energy_origin_level_modified,
                environment_saturation_moist_static_energy_modified=locals.environment_saturation_moist_static_energy_modified,
                environment_moist_static_energy_cloud_levels_modified=locals.environment_moist_static_energy_cloud_levels_modified,
                environment_saturation_moist_static_energy_cloud_levels_modified=locals.environment_saturation_moist_static_energy_cloud_levels_modified,
                mass_detrainment_updraft_forced=state.output.mass_detrainment_updraft_forced,
                mass_entrainment_updraft_forced=state.output.mass_entrainment_updraft_forced,
                normalized_massflux_updraft_modified=locals.normalized_massflux_updraft_modified,
                partition_liquid_ice=locals.partition_liquid_ice,
                vapor_excess=locals.vapor_excess,
                t_excess=locals.t_excess,
                add_buoyancy=locals.add_buoyancy,
                cloud_liquid_after_rain_forced=state.output.cloud_liquid_after_rain_forced,
                d_buoyancy_modified=locals.d_buoyancy_modified,
                plume_dependent_constants=plume_dependent_constants,
            )

        # write output
        outputs = {
            "error_code": state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX],
            "local_start_level": locals.start_level.data[:] + 1,
            "lcl_level": state.output.lcl_level.data[:, :, plume_dependent_constants.PLUME_INDEX] + 1,
            "updraft_lfc_level": state.output.updraft_lfc_level.data[
                :, :, plume_dependent_constants.PLUME_INDEX
            ]
            + 1,
            "cloud_top_level": state.output.cloud_top_level.data[:, :, plume_dependent_constants.PLUME_INDEX]
            + 1,
            "local_cloud_moist_static_energy_modified": locals.cloud_moist_static_energy_modified.data[:],
            "local_moist_static_energy_origin_level_modified": locals.moist_static_energy_origin_level_modified.data[
                :
            ],
            "local_env_saturation_moist_static_energy_modified": locals.environment_saturation_moist_static_energy_modified.data[
                :
            ],
            "local_env_moist_static_energy_cloud_levels_modified": locals.environment_moist_static_energy_cloud_levels_modified.data[
                :
            ],
            "local_env_saturation_moist_static_energy_cloud_levels_modified": locals.environment_saturation_moist_static_energy_cloud_levels_modified.data[
                :
            ],
            "mass_detrainment_updraft_forced": state.output.mass_detrainment_updraft_forced.data[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "mass_entrainment_updraft_forced": state.output.mass_entrainment_updraft_forced.data[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_normalized_massflux_updraft_modified": locals.normalized_massflux_updraft_modified.data[:],
            "local_partition_liquid_ice": locals.partition_liquid_ice.data[:],
            "local_vapor_excess": locals.vapor_excess.data[:],
            "local_t_excess": locals.t_excess.data[:],
            "local_add_buoyancy": locals.add_buoyancy.data[:],
            "cloud_liquid_after_rain_forced": state.output.cloud_liquid_after_rain_forced.data[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_d_buoyancy_modified": locals.d_buoyancy_modified.data[:],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_StaticControl_shallow(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_StaticControl_mid(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_StaticControl_deep(TranslateFortranData2Py):
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
