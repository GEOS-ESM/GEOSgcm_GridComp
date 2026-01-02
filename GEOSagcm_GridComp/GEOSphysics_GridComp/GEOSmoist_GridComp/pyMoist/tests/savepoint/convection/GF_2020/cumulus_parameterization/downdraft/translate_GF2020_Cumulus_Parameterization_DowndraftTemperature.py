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
from pyMoist.convection.GF_2020.cumulus_parameterization.downdraft import downdraft_temperature
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
            "error_code_downtemp": {},
            "local_downdraft_column_temperature_forced_downtemp": {},
            "local_cloud_moist_static_energy_downdraft_forced_downtemp": {},
            "local_geopotential_height_cloud_levels_forced_downtemp": {},
            "local_cloud_total_water_after_entrainment_downdraft_forced_downmoisture": {},
            "local_t_cloud_levels_forced_downtemp": {},
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
            "error_code_downtemp"
        ]
        locals.downdraft_column_temperature_forced.data[:] = inputs[
            "local_downdraft_column_temperature_forced_downtemp"
        ]
        locals.cloud_moist_static_energy_downdraft_forced.data[:] = inputs[
            "local_cloud_moist_static_energy_downdraft_forced_downtemp"
        ]
        locals.geopotential_height_cloud_levels_forced.data[:] = inputs[
            "local_geopotential_height_cloud_levels_forced_downtemp"
        ]
        locals.cloud_total_water_after_entrainment_downdraft_forced.data[:] = inputs[
            "local_cloud_total_water_after_entrainment_downdraft_forced_downmoisture"
        ]
        locals.t_cloud_levels_forced.data[:] = inputs["local_t_cloud_levels_forced_downtemp"]

        # initalize test code
        code = self.stencil_factory.from_dims_halo(
            func=downdraft_temperature,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        # call test code
        if plume_dependent_constants.ENABLE_PLUME == 1:
            code(
                error_code=state.output.error_code,
                downdraft_column_temperature_forced=locals.downdraft_column_temperature_forced,
                cloud_moist_static_energy_downdraft_forced=locals.cloud_moist_static_energy_downdraft_forced,
                geopotential_height_cloud_levels_forced=locals.geopotential_height_cloud_levels_forced,
                cloud_total_water_after_entrainment_downdraft_forced=locals.cloud_total_water_after_entrainment_downdraft_forced,
                t_cloud_levels_forced=locals.t_cloud_levels_forced,
                plume=plume_dependent_constants.PLUME_INDEX,
            )

        # write output
        outputs = {
            "error_code_downtemp": state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX],
            "local_downdraft_column_temperature_forced_downtemp": locals.downdraft_column_temperature_forced.data[
                :
            ],
            "local_cloud_moist_static_energy_downdraft_forced_downtemp": locals.cloud_moist_static_energy_downdraft_forced.data[
                :
            ],
            "local_geopotential_height_cloud_levels_forced_downtemp": locals.geopotential_height_cloud_levels_forced.data[
                :
            ],
            "local_cloud_total_water_after_entrainment_downdraft_forced_downmoisture": locals.cloud_total_water_after_entrainment_downdraft_forced.data[
                :
            ],
            "local_t_cloud_levels_forced_downtemp": locals.t_cloud_levels_forced.data[:],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_DowndraftTemperature_shallow(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_DowndraftTemperature_mid(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_DowndraftTemperature_deep(TranslateFortranData2Py):
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
