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
)
from pyMoist.convection.GF_2020.cumulus_parameterization.profiles.profiles import (
    InCloudTemperature,
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
            "error_code_incloudtemp": {},
            "local_incloud_air_temp_forced_incloudtemp": {},
            "local_hcdo_incloudtemp": {},
            "local_geopotential_height_cloud_levels_forced_incloudtemp": {},
            "local_incloud_water_vapor_mixing_ratio_forced_incloudtemp": {},
            "local_t_cloud_levels_forced_incloudtemp": {},
        }

        out_vars.update(in_vars["data_vars"])

    def __call__(self, constants: dict, cu_param_constants: dict, plume: str, **inputs):
        # initalize constants
        config = GF2020Config(SINGLE_COLUMN_MODE=False, **constants)
        cumulus_parameterization_config = GF2020CumulusParameterizationConfig(
            **cu_param_constants
        )
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
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["error_code_incloudtemp"]
        )
        locals.incloud_air_temp_forced.data[:] = inputs[
            "local_incloud_air_temp_forced_incloudtemp"
        ]
        locals.hcdo.data[:] = inputs["local_hcdo_incloudtemp"]
        locals.geopotential_height_cloud_levels_forced.data[:] = inputs[
            "local_geopotential_height_cloud_levels_forced_incloudtemp"
        ]
        locals.incloud_water_vapor_mixing_ratio_forced.data[:] = inputs[
            "local_incloud_water_vapor_mixing_ratio_forced_incloudtemp"
        ]
        locals.t_cloud_levels_forced.data[:] = inputs[
            "local_t_cloud_levels_forced_incloudtemp"
        ]

        # initalize test code
        code = InCloudTemperature(
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
            "error_code_incloudtemp": state.output.error_code.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_incloud_air_temp_forced_incloudtemp": locals.incloud_air_temp_forced.field[
                :
            ],
            "local_hcdo_incloudtemp": locals.hcdo.field[:],
            "local_geopotential_height_cloud_levels_forced_incloudtemp": locals.geopotential_height_cloud_levels_forced.field[
                :
            ],
            "local_incloud_water_vapor_mixing_ratio_forced_incloudtemp": locals.incloud_water_vapor_mixing_ratio_forced.field[
                :
            ],
            "local_t_cloud_levels_forced_incloudtemp": locals.t_cloud_levels_forced.field[
                :
            ],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_InCloudTemperature_shallow(
    TranslateFortranData2Py
):
    def __init__(
        self,
        grid: Grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        self.test_core = TestCore(
            grid, namelist, stencil_factory, self.in_vars, self.out_vars
        )

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF2020-constants")
        self.cu_param_constants = data_loader.load(
            "GF2020_CumulusParameterization-constants"
        )

    def compute_func(self, **inputs):
        outputs = self.test_core(
            self.constants, self.cu_param_constants, "shallow", **inputs
        )

        return outputs


class TranslateGF2020_CumulusParameterization_InCloudTemperature_mid(
    TranslateFortranData2Py
):
    def __init__(
        self,
        grid: Grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        self.test_core = TestCore(
            grid, namelist, stencil_factory, self.in_vars, self.out_vars
        )

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF2020-constants")
        self.cu_param_constants = data_loader.load(
            "GF2020_CumulusParameterization-constants"
        )

    def compute_func(self, **inputs):
        outputs = self.test_core(
            self.constants, self.cu_param_constants, "mid", **inputs
        )

        return outputs


class TranslateGF2020_CumulusParameterization_InCloudTemperature_deep(
    TranslateFortranData2Py
):
    def __init__(
        self,
        grid: Grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        self.test_core = TestCore(
            grid, namelist, stencil_factory, self.in_vars, self.out_vars
        )

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF2020-constants")
        self.cu_param_constants = data_loader.load(
            "GF2020_CumulusParameterization-constants"
        )

    def compute_func(self, **inputs):
        outputs = self.test_core(
            self.constants, self.cu_param_constants, "deep", **inputs
        )

        return outputs
