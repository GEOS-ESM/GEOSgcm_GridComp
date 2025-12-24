from f90nml import Namelist
from ndsl import StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
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
from pyMoist.convection.GF_2020.cumulus_parameterization.profiles import (
    get_melting_profile,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.setup.set_constants import (
    set_constants,
)


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
            "error_code_getmeltpro": {},
            "local_melting_layer_getmeltpro": {},
            "local_partition_liquid_ice_getmeltpro": {},
            "p_cloud_levels_forced_getmeltpro": {},
            "condensate_to_fall_forced_getmeltpro": {},
            "local_melting_getmeltpro": {},
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
                "ensemble_members": MAXENS1 * MAXENS2 * MAXENS3,
            },
        )

        # fill relevant parts of dataclasses
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "error_code_getmeltpro"
        ]

        state.output.p_cloud_levels_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "p_cloud_levels_forced_getmeltpro"
        ]

        state.output.condensate_to_fall_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "condensate_to_fall_forced_getmeltpro"
        ]

        locals.melting_layer.data[:] = inputs["local_melting_layer_getmeltpro"]

        locals.partition_liquid_ice.data[:] = inputs["local_partition_liquid_ice_getmeltpro"]

        locals.melting.data[:] = inputs["local_melting_getmeltpro"]

        # initalize test code
        code = self.stencil_factory.from_dims_halo(
            func=get_melting_profile,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "MELT_GLAC": cumulus_parameterization_config.MELT_GLAC,
            },
        )

        # call test code
        if plume_dependent_constants.ENABLE_PLUME == 1:
            code(
                error_code=state.output.error_code,
                plume=plume_dependent_constants.PLUME_INDEX,
                melting_layer=locals.melting_layer,
                partition_liquid_ice=locals.partition_liquid_ice,
                p_cloud_levels_forced=state.output.p_cloud_levels_forced,
                condensate_to_fall_forced=state.output.condensate_to_fall_forced,
                melting=locals.melting,
            )

        # write output
        outputs = {
            "error_code_getmeltpro": state.output.error_code.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "p_cloud_levels_forced": state.output.p_cloud_levels_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "condensate_to_fall_forced_getmeltpro": state.output.condensate_to_fall_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_melting_layer_getmeltpro": locals.melting_layer.field[:],
            "local_melting_getmeltpro": locals.melting.field[:],
            "local_partition_liquid_ice_getmeltpro": locals.partition_liquid_ice.field[:],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_MeltingProfile_shallow(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_MeltingProfile_mid(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_MeltingProfile_deep(TranslateFortranData2Py):
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
