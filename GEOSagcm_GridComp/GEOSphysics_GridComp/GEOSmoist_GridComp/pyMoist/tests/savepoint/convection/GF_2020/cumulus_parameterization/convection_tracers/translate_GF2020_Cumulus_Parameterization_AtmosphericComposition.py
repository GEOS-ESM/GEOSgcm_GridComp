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
from pyMoist.convection.GF_2020.cumulus_parameterization.convective_tracers import AtmosphericComposition

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
            # "chemistry_tracers": {},
            # "local_chemistry_tracers_cloud_levels": {},
        }

        out_vars.update(in_vars["data_vars"])
        out_vars.update(
            {
                "chemistry_tracers": {},
                "local_chemistry_tracers_cloud_levels": {},
            }
        )

    def __call__(self, constants: dict, cu_param_constants: dict, plume: str, ddim_fields: dict, **inputs):
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
        state.input_output.chemistry_tracers.field[:] = ddim_fields["chemistry_tracers"]
        locals.chemistry_tracers_cloud_levels.field[:] = ddim_fields["local_chemistry_tracers_cloud_levels"]

        code = AtmosphericComposition(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        if plume_dependent_constants.ENABLE_PLUME == 1:
            code(
                error_code=state.output.error_code,
                chemistry_tracers=state.input_output.chemistry_tracers,
                chemistry_tracers_cloud_levels=locals.chemistry_tracers_cloud_levels,
                plume_dependent_constants=plume_dependent_constants,
            )

        outputs = {
            "error_code": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "chemistry_tracers": state.input_output.chemistry_tracers.field[:],
            "local_chemistry_tracers_cloud_levels": locals.chemistry_tracers_cloud_levels.field[:],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_AtmosphericComposition_shallow(TranslateFortranData2Py):
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
        self.ddim_fields = data_loader.load(
            "GF2020_CumulusParameterization_AtmosphericComposition_shallow-In"
        )

    def compute_func(self, **inputs):
        outputs = self.test_core(
            self.constants, self.cu_param_constants, "shallow", ddim_fields=self.ddim_fields, **inputs
        )

        return outputs


class TranslateGF2020_CumulusParameterization_AtmosphericComposition_mid(TranslateFortranData2Py):
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
        self.ddim_fields = data_loader.load("GF2020_CumulusParameterization_AtmosphericComposition_mid-In")

    def compute_func(self, **inputs):
        outputs = self.test_core(
            self.constants, self.cu_param_constants, "mid", ddim_fields=self.ddim_fields, **inputs
        )

        return outputs


class TranslateGF2020_CumulusParameterization_AtmosphericComposition_deep(TranslateFortranData2Py):
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
        self.ddim_fields = data_loader.load("GF2020_CumulusParameterization_AtmosphericComposition_deep-In")

    def compute_func(self, **inputs):
        outputs = self.test_core(
            self.constants, self.cu_param_constants, "deep", ddim_fields=self.ddim_fields, **inputs
        )

        return outputs
