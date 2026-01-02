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
from pyMoist.convection.GF_2020.cumulus_parameterization.precip import partition_liquid_ice
from pyMoist.convection.GF_2020.cumulus_parameterization.get_levels import (
    find_maximum_updraft_origin_level,
    find_detrainmet_start_level,
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
            "error_code_partition": {},
            "local_t_new_partition": {},
            "topography_height_no_negative_partition": {},
            "local_geopotential_height_cloud_levels_forced_partition": {},
            "p_cloud_levels_forced_partition": {},
            "local_partition_liquid_ice_partition": {},
            "local_melting_layer_partition": {},
            "convection_fraction_partition": {},
            "surface_type_partition": {},
            "local_maximum_updraft_origin_level_partition": {},
            "local_detrainment_start_level_partition": {},
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
            },
        )

        locals = GF2020CumulusParameterizationLocals.zeros(
            self.quantity_factory,
            data_dimensions={
                "ensemble_1": MAXENS1,
                "ensemble_2": MAXENS2,
                "ensemble_3": MAXENS3,
                "ensemble_members": MAXENS1 * MAXENS2 * MAXENS3,
            },
        )

        # fill relevant parts of dataclasses
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "error_code_partition"
        ]
        locals.t_new.data[:] = inputs["local_t_new_partition"]
        state.input_output.topography_height_no_negative.data[:] = inputs[
            "topography_height_no_negative_partition"
        ]
        locals.geopotential_height_cloud_levels_forced.data[:] = inputs[
            "local_geopotential_height_cloud_levels_forced_partition"
        ]
        state.output.p_cloud_levels_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "p_cloud_levels_forced_partition"
        ]
        locals.partition_liquid_ice.data[:] = inputs["local_partition_liquid_ice_partition"]
        locals.melting_layer.data[:] = inputs["local_melting_layer_partition"]
        state.input.convection_fraction.data[:] = inputs["convection_fraction_partition"]
        state.input.surface_type.data[:] = inputs["surface_type_partition"]
        locals.maximum_updraft_origin_level.data[:] = (
            inputs["local_maximum_updraft_origin_level_partition"] - 1
        )
        locals.detrainment_start_level.data[:] = inputs["local_detrainment_start_level_partition"] - 1

        code_part_1 = self.stencil_factory.from_dims_halo(
            func=partition_liquid_ice,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "MELT_GLAC": cumulus_parameterization_config.MELT_GLAC,
                "FRAC_MODIS": cumulus_parameterization_config.FRAC_MODIS,
            },
        )

        code_part_2 = self.stencil_factory.from_dims_halo(
            func=find_maximum_updraft_origin_level,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        code_part_3 = self.stencil_factory.from_dims_halo(
            func=find_detrainmet_start_level,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        if plume_dependent_constants.ENABLE_PLUME == 1:
            code_part_1(
                t=locals.t_new,
                p=state.output.p_cloud_levels_forced,
                geopotential_height=locals.geopotential_height_cloud_levels_forced,
                topography_height_no_negative=state.input_output.topography_height_no_negative,
                surface_type=state.input.surface_type,
                convection_fraction=state.input.convection_fraction,
                error_code=state.output.error_code,
                melting_layer=locals.melting_layer,
                part_liquid_ice=locals.partition_liquid_ice,
                plume=plume_dependent_constants.PLUME_INDEX,
            )

            code_part_2(
                geopotential_height=locals.geopotential_height_cloud_levels_forced,
                topography_height_no_negative=state.input_output.topography_height_no_negative,
                error_code=state.output.error_code,
                maximum_updraft_origin_level=locals.maximum_updraft_origin_level,
                MAX_UPDRAFT_ORIGIN_HEIGHT=plume_dependent_constants.MAX_UPDRAFT_ORIGIN_HEIGHT,
                plume=plume_dependent_constants.PLUME_INDEX,
            )

            code_part_3(
                geopotential_height=locals.geopotential_height_cloud_levels_forced,
                topography_height_no_negative=state.input_output.topography_height_no_negative,
                error_code=state.output.error_code,
                detrainment_start_level=locals.detrainment_start_level,
                DETRAINMENT_CRITICAL_DEPTH=plume_dependent_constants.DETRAINMENT_CRITICAL_DEPTH,
                plume=plume_dependent_constants.PLUME_INDEX,
            )

        outputs = {
            "error_code_partition": state.output.error_code.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_t_new_partition": locals.t_new.field[:],
            "topography_height_no_negative_partition": state.input_output.topography_height_no_negative.field[
                :
            ],
            "local_geopotential_height_cloud_levels_forced_partition": locals.geopotential_height_cloud_levels_forced.field[
                :
            ],
            "p_cloud_levels_forced_partition": state.output.p_cloud_levels_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_partition_liquid_ice_partition": locals.partition_liquid_ice.field[:],
            "local_melting_layer_partition": locals.melting_layer.field[:],
            "convection_fraction_partition": state.input.convection_fraction.field[:],
            "surface_type_partition": state.input.surface_type.field[:],
            "local_maximum_updraft_origin_level_partition": locals.maximum_updraft_origin_level.field[:] + 1,
            "local_detrainment_start_level_partition": locals.detrainment_start_level.field[:] + 1,
        }

        return outputs


class TranslateGF2020_CumulusParameterization_PartitionLiquidIceAndGetLevels_shallow(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_PartitionLiquidIceAndGetLevels_mid(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_PartitionLiquidIceAndGetLevels_deep(TranslateFortranData2Py):
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
