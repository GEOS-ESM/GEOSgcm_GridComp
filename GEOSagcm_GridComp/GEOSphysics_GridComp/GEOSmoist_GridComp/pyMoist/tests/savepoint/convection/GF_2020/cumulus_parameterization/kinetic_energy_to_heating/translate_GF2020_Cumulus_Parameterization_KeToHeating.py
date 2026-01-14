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
from pyMoist.convection.GF_2020.cumulus_parameterization.kinetic_energy_to_heating.stencils import (
    ke_to_heating,
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
            "error_code_kth": {},
            "p_cloud_levels_forced": {},
            "u": {},
            "v": {},
            "local_del_u_cloud_ensemble": {},
            "cloud_top_level": {},
            "local_del_v_cloud_ensemble": {},
            "local_del_t_cloud_ensemble": {},
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
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["error_code_kth"]
        state.output.cloud_top_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["cloud_top_level"] - 1
        )
        state.output.p_cloud_levels_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "p_cloud_levels_forced"
        ]
        state.input_output.u.data[:] = inputs["u"]
        state.input_output.v.data[:] = inputs["v"]
        locals.del_u_cloud_ensemble.data[:] = inputs["local_del_u_cloud_ensemble"]
        locals.del_v_cloud_ensemble.data[:] = inputs["local_del_v_cloud_ensemble"]
        locals.del_t_cloud_ensemble.data[:] = inputs["local_del_t_cloud_ensemble"]

        code_part_1 = self.stencil_factory.from_dims_halo(
            func=ke_to_heating,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        if plume_dependent_constants.ENABLE_PLUME == 1:
            code_part_1(
                del_u_cloud_ensemble=locals.del_u_cloud_ensemble,
                del_v_cloud_ensemble=locals.del_v_cloud_ensemble,
                error_code=state.output.error_code,
                plume=plume_dependent_constants.PLUME_INDEX,
                cloud_top_level=state.output.cloud_top_level,
                p_cloud_levels_forced=state.output.p_cloud_levels_forced,
                u=state.input_output.u,
                v=state.input_output.v,
                del_t_cloud_ensemble=locals.del_t_cloud_ensemble,
            )

        outputs = {
            "error_code_kth": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "cloud_top_level": state.output.cloud_top_level.field[:, :, plume_dependent_constants.PLUME_INDEX]
            + 1,
            "p_cloud_levels_forced": state.output.p_cloud_levels_forced.data[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "u": state.input_output.u.data[:],
            "v": state.input_output.v.data[:],
            "local_del_u_cloud_ensemble": locals.del_u_cloud_ensemble.data[:],
            "local_del_v_cloud_ensemble": locals.del_v_cloud_ensemble.data[:],
            "local_del_t_cloud_ensemble": locals.del_t_cloud_ensemble.data[:],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_KeToHeating_shallow(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_KeToHeating_mid(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_KeToHeating_deep(TranslateFortranData2Py):
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
