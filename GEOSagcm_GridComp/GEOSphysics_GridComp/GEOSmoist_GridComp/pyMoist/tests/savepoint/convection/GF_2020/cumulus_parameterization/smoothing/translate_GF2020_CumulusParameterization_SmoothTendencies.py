from f90nml import Namelist

from ndsl import StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.constants import (
    MAXENS1,
    MAXENS2,
    MAXENS3,
    NUMBER_OF_PLUMES,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.locals import GF2020CumulusParameterizationLocals
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.setup.set_constants import set_constants
from pyMoist.convection.GF_2020.cumulus_parameterization.smoothing import smooth_tendencies
from pyMoist.convection.GF_2020.cumulus_parameterization.state import GF2020CumulusParameterizationState


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
            "cloud_top_level": {},
            "p_cloud_levels_forced": {},
            "local_del_moist_static_energy_cloud_ensemble": {},
            "local_del_vapor_cloud_ensemble": {},
            "local_del_cloud_liquid_cloud_ensemble": {},
            "local_del_u_cloud_ensemble": {},
            "local_del_v_cloud_ensemble": {},
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
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["error_code"]
        state.output.cloud_top_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["cloud_top_level"] - 1
        )
        state.output.p_cloud_levels_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "p_cloud_levels_forced"
        ]
        locals.del_moist_static_energy_cloud_ensemble.data[:] = inputs[
            "local_del_moist_static_energy_cloud_ensemble"
        ]
        locals.del_vapor_cloud_ensemble.data[:] = inputs["local_del_vapor_cloud_ensemble"]
        locals.del_cloud_liquid_cloud_ensemble.data[:] = inputs["local_del_cloud_liquid_cloud_ensemble"]
        locals.del_u_cloud_ensemble.data[:] = inputs["local_del_u_cloud_ensemble"]
        locals.del_v_cloud_ensemble.data[:] = inputs["local_del_v_cloud_ensemble"]

        # initialize test code
        code = self.stencil_factory.from_dims_halo(
            func=smooth_tendencies,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"USE_SMOOTH_TENDENCIES": cumulus_parameterization_config.USE_SMOOTH_TENDENCIES},
        )

        # call test code
        if plume_dependent_constants.ENABLE_PLUME == 1:
            code(
                error_code=state.output.error_code,
                cloud_top_level=state.output.cloud_top_level,
                p_cloud_levels_forced=state.output.p_cloud_levels_forced,
                del_moist_static_energy_cloud_ensemble=locals.del_moist_static_energy_cloud_ensemble,
                del_vapor_cloud_ensemble=locals.del_vapor_cloud_ensemble,
                del_cloud_liquid_cloud_ensemble=locals.del_cloud_liquid_cloud_ensemble,
                del_u_cloud_ensemble=locals.del_u_cloud_ensemble,
                del_v_cloud_ensemble=locals.del_v_cloud_ensemble,
                plume=plume_dependent_constants.PLUME_INDEX,
            )

        # write output
        outputs = {
            "error_code": state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX],
            "cloud_top_level": state.output.cloud_top_level.data[:, :, plume_dependent_constants.PLUME_INDEX]
            + 1,
            "p_cloud_levels_forced": state.output.p_cloud_levels_forced.data[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_del_moist_static_energy_cloud_ensemble": locals.del_moist_static_energy_cloud_ensemble.data[
                :
            ],
            "local_del_vapor_cloud_ensemble": locals.del_vapor_cloud_ensemble.data[:],
            "local_del_cloud_liquid_cloud_ensemble": locals.del_cloud_liquid_cloud_ensemble.data[:],
            "local_del_u_cloud_ensemble": locals.del_u_cloud_ensemble.data[:],
            "local_del_v_cloud_ensemble": locals.del_v_cloud_ensemble.data[:],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_SmoothTendencies_shallow(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_SmoothTendencies_mid(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_SmoothTendencies_deep(TranslateFortranData2Py):
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
