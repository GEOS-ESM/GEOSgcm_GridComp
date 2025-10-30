from ndsl import Namelist, StencilFactory
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
from pyMoist.convection.GF_2020.cumulus_parameterization.precip.precip import PartitionLiquidIce
from pyMoist.convection.GF_2020.cumulus_parameterization.get_levels.get_levels import (
    MaximumUpdraftOriginLevel,
    DowndraftDetrainmentLevel,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.setup.set_constants import set_constants


class TranslateGF2020_CumulusParameterization_PartitionLiquidIce_shallow(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        self.in_vars["data_vars"] = {
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

        self.out_vars = self.in_vars["data_vars"].copy()

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF2020-constants")
        self.cu_param_constants = data_loader.load("GF2020_CumulusParameterization-constants")

    def compute_func(self, **inputs):
        # initalize constants
        config = GF2020Config(SINGLE_COLUMN_MODE=False, **self.constants)
        cumulus_parameterization_config = GF2020CumulusParameterizationConfig(**self.cu_param_constants)
        plume_dependent_constants = GF2020PlumeDependentConstants()
        plume_dependent_constants = set_constants(
            cumulus_parameterization_config, plume_dependent_constants, "shallow"
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
        locals.maximum_updraft_origin_level.data[:] = inputs["local_maximum_updraft_origin_level_partition"]
        locals.detrainment_start_level.data[:] = inputs["local_detrainment_start_level_partition"]

        # initalize test code
        partition_liquid_ice = PartitionLiquidIce(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        maximum_updraft_origin_level = MaximumUpdraftOriginLevel(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        downdraft_detrainment_level = DowndraftDetrainmentLevel(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        # call test code
        if plume_dependent_constants.ENABLE_PLUME == 1:
            partition_liquid_ice(
                state=state,
                locals=locals,
                plume_dependent_constants=plume_dependent_constants,
            )
            maximum_updraft_origin_level(
                state=state,
                locals=locals,
                plume_dependent_constants=plume_dependent_constants,
            )
            downdraft_detrainment_level(
                state=state,
                locals=locals,
                plume_dependent_constants=plume_dependent_constants,
            )

            locals.maximum_updraft_origin_level.field[:] = locals.maximum_updraft_origin_level.field[:] + 1
            locals.detrainment_start_level.field[:] = locals.detrainment_start_level.field[:] + 1

        # write output
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
            "local_maximum_updraft_origin_level_partition": locals.maximum_updraft_origin_level.field[:],
            "local_detrainment_start_level_partition": locals.detrainment_start_level.field[:],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_PartitionLiquidIce_mid(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        self.in_vars["data_vars"] = {
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

        self.out_vars = self.in_vars["data_vars"].copy()

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF2020-constants")
        self.cu_param_constants = data_loader.load("GF2020_CumulusParameterization-constants")

    def compute_func(self, **inputs):
        # initalize constants
        config = GF2020Config(SINGLE_COLUMN_MODE=False, **self.constants)
        cumulus_parameterization_config = GF2020CumulusParameterizationConfig(**self.cu_param_constants)
        plume_dependent_constants = GF2020PlumeDependentConstants()
        plume_dependent_constants = set_constants(
            cumulus_parameterization_config, plume_dependent_constants, "mid"
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
        locals.maximum_updraft_origin_level.data[:] = inputs["local_maximum_updraft_origin_level_partition"]
        locals.detrainment_start_level.data[:] = inputs["local_detrainment_start_level_partition"]

        # initalize test code
        partition_liquid_ice = PartitionLiquidIce(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        maximum_updraft_origin_level = MaximumUpdraftOriginLevel(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        downdraft_detrainment_level = DowndraftDetrainmentLevel(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        # call test code
        if plume_dependent_constants.ENABLE_PLUME == 1:
            partition_liquid_ice(
                state=state,
                locals=locals,
                plume_dependent_constants=plume_dependent_constants,
            )
            maximum_updraft_origin_level(
                state=state,
                locals=locals,
                plume_dependent_constants=plume_dependent_constants,
            )
            downdraft_detrainment_level(
                state=state,
                locals=locals,
                plume_dependent_constants=plume_dependent_constants,
            )

            locals.maximum_updraft_origin_level.field[:] = locals.maximum_updraft_origin_level.field[:] + 1
            locals.detrainment_start_level.field[:] = locals.detrainment_start_level.field[:] + 1

        # write output
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
            "local_maximum_updraft_origin_level_partition": locals.maximum_updraft_origin_level.field[:],
            "local_detrainment_start_level_partition": locals.detrainment_start_level.field[:],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_PartitionLiquidIce_deep(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        self.in_vars["data_vars"] = {
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

        self.out_vars = self.in_vars["data_vars"].copy()

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF2020-constants")
        self.cu_param_constants = data_loader.load("GF2020_CumulusParameterization-constants")

    def compute_func(self, **inputs):
        # initalize constants
        config = GF2020Config(SINGLE_COLUMN_MODE=False, **self.constants)
        cumulus_parameterization_config = GF2020CumulusParameterizationConfig(**self.cu_param_constants)
        plume_dependent_constants = GF2020PlumeDependentConstants()
        plume_dependent_constants = set_constants(
            cumulus_parameterization_config, plume_dependent_constants, "deep"
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
        locals.maximum_updraft_origin_level.data[:] = inputs["local_maximum_updraft_origin_level_partition"]
        locals.detrainment_start_level.data[:] = inputs["local_detrainment_start_level_partition"]

        # initalize test code
        partition_liquid_ice = PartitionLiquidIce(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        maximum_updraft_origin_level = MaximumUpdraftOriginLevel(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        downdraft_detrainment_level = DowndraftDetrainmentLevel(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        # call test code
        if plume_dependent_constants.ENABLE_PLUME == 1:
            partition_liquid_ice(
                state=state,
                locals=locals,
                plume_dependent_constants=plume_dependent_constants,
            )
            maximum_updraft_origin_level(
                state=state,
                locals=locals,
                plume_dependent_constants=plume_dependent_constants,
            )
            downdraft_detrainment_level(
                state=state,
                locals=locals,
                plume_dependent_constants=plume_dependent_constants,
            )

            locals.maximum_updraft_origin_level.field[:] = locals.maximum_updraft_origin_level.field[:] + 1
            locals.detrainment_start_level.field[:] = locals.detrainment_start_level.field[:] + 1

        # write output
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
            "local_maximum_updraft_origin_level_partition": locals.maximum_updraft_origin_level.field[:],
            "local_detrainment_start_level_partition": locals.detrainment_start_level.field[:],
        }

        return outputs
