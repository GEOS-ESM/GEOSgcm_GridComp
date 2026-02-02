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
from pyMoist.convection.GF_2020.cumulus_parameterization.shared_stencils import updraft_vertical_velocity
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
            "local_vertical_velocity_3d": {},
            "local_vertical_velocity_2d": {},
            "convective_scale_velocity": {},
            "entrainment_rate": {},
            "error_code": {},
            "local_detrainment_function_updraft": {},
            "local_geopotential_height_cloud_levels_forced": {},
            "local_t_cloud_levels_forced": {},
            "local_updraft_column_temperature_forced": {},
            "local_cloud_total_water_after_entrainment_forced": {},
            "cloud_liquid_after_rain_forced": {},
            "local_vapor_forced": {},
            "lcl_level": {},
            "cloud_top_level": {},
            "updraft_lfc_level": {},
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
        locals.vertical_velocity_3d.data[:] = inputs["local_vertical_velocity_3d"]
        locals.vertical_velocity_2d.data[:] = inputs["local_vertical_velocity_2d"]
        state.input_output.convective_scale_velocity.data[:] = inputs["convective_scale_velocity"]
        state.output.entrainment_rate.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "entrainment_rate"
        ]
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["error_code"]
        locals.detrainment_function_updraft.data[:] = inputs["local_detrainment_function_updraft"]
        locals.geopotential_height_cloud_levels_forced.data[:] = inputs[
            "local_geopotential_height_cloud_levels_forced"
        ]
        locals.t_cloud_levels_forced.data[:] = inputs["local_t_cloud_levels_forced"]
        locals.updraft_column_temperature_forced.data[:] = inputs["local_updraft_column_temperature_forced"]
        locals.cloud_total_water_after_entrainment_forced.data[:] = inputs[
            "local_cloud_total_water_after_entrainment_forced"
        ]
        state.output.cloud_liquid_after_rain_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["cloud_liquid_after_rain_forced"]
        )
        locals.vapor_forced.data[:] = inputs["local_vapor_forced"]
        state.output.lcl_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["lcl_level"] - 1
        state.output.cloud_top_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["cloud_top_level"] - 1
        )
        state.output.updraft_lfc_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["updraft_lfc_level"] - 1
        )

        # initialize test code
        code = self.stencil_factory.from_dims_halo(
            func=updraft_vertical_velocity,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ZERO_DIFF": cumulus_parameterization_config.ZERO_DIFF},
        )

        # call test code
        if plume_dependent_constants.ENABLE_PLUME == 1:
            if cumulus_parameterization_config.FIRST_GUESS_W == 0:
                code(
                    vertical_velocity_3d=locals.vertical_velocity_3d,
                    vertical_velocity_2d=locals.vertical_velocity_2d,
                    convective_scale_velocity=state.input_output.convective_scale_velocity,
                    entrainment_rate=state.output.entrainment_rate,
                    detrainment_function_updraft=locals.detrainment_function_updraft,
                    geopotential_height_cloud_levels_forced=locals.geopotential_height_cloud_levels_forced,
                    t_cloud_levels_forced=locals.t_cloud_levels_forced,
                    updraft_column_temperature_forced=locals.updraft_column_temperature_forced,
                    cloud_total_water_after_entrainment_forced=locals.cloud_total_water_after_entrainment_forced,
                    cloud_liquid_after_rain_forced=state.output.cloud_liquid_after_rain_forced,
                    vapor_forced=locals.vapor_forced,
                    updraft_lfc_level=state.output.updraft_lfc_level,
                    cloud_top_level=state.output.cloud_top_level,
                    error_code=state.output.error_code,
                    plume=plume_dependent_constants.PLUME_INDEX,
                )

        # write output
        outputs = {
            "local_vertical_velocity_3d": locals.vertical_velocity_3d.field[:],
            "local_vertical_velocity_2d": locals.vertical_velocity_2d.field[:],
            "convective_scale_velocity": state.input_output.convective_scale_velocity.field[:],
            "entrainment_rate": state.output.entrainment_rate.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "error_code": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "local_detrainment_function_updraft": locals.detrainment_function_updraft.field[:],
            "local_geopotential_height_cloud_levels_forced": locals.geopotential_height_cloud_levels_forced.field[
                :
            ],
            "local_t_cloud_levels_forced": locals.t_cloud_levels_forced.field[:],
            "local_updraft_column_temperature_forced": locals.updraft_column_temperature_forced.field[:],
            "local_cloud_total_water_after_entrainment_forced": locals.cloud_total_water_after_entrainment_forced.field[
                :
            ],
            "cloud_liquid_after_rain_forced": state.output.cloud_liquid_after_rain_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_vapor_forced": locals.vapor_forced.field[:],
            "lcl_level": state.output.lcl_level.field[:, :, plume_dependent_constants.PLUME_INDEX] + 1,
            "cloud_top_level": state.output.cloud_top_level.field[:, :, plume_dependent_constants.PLUME_INDEX]
            + 1,
            "updraft_lcf_level": state.output.updraft_lfc_level.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ]
            + 1,
        }

        return outputs


class TranslateGF2020_CumulusParameterization_UpdraftVerticalVelocity_shallow(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_UpdraftVerticalVelocity_mid(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_UpdraftVerticalVelocity_deep(TranslateFortranData2Py):
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
