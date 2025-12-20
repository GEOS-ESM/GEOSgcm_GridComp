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
from pyMoist.convection.GF_2020.cumulus_parameterization.downdraft import (
    DowndraftWindshear,
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
            "error_code_downdraftwindshear": {},
            "u_downdraftwindshear": {},
            "v_downdraftwindshear": {},
            "geopotential_height_forced_downdraftwindshear": {},
            "cloud_top_downdraftwindshear": {},
            "updraft_lfc_level_downdraftwindshear": {},
            "local_epsilon_downdraftwindshear": {},
            "p_forced_downdraftwindshear": {},
            "local_pwavo_downdraftwindshear": {},
            "ccn_downdraftwindshear_in_out": {},
            "local_pwevo_downdraftwindshear": {},
            "local_epsilon_min_downdraftwindshear": {},
            "local_epsilon_max_downdraftwindshear": {},
            "local_psum_downdraftwindshear": {},
            "local_psumh_downdraftwindshear": {},
            "epsilon_downdraftwindshear": {},
            "local_scale_dependence_factor_downdraft_downdraftwindshear": {},
            "local_edtc_downdraftwindshear": {},
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
            "error_code_downdraftwindshear"
        ]
        state.input_output.u.data[:, :, :] = inputs["u_downdraftwindshear"]
        state.input_output.v.data[:, :, :] = inputs["v_downdraftwindshear"]
        state.input_output.ccn.data[:, :] = inputs["ccn_downdraftwindshear_in_out"]
        locals.epsilon_max.data[:] = inputs["local_epsilon_max_downdraftwindshear"]
        locals.epsilon_min.data[:] = inputs["local_epsilon_min_downdraftwindshear"]
        state.output.updraft_lfc_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "updraft_lfc_level_downdraftwindshear"
        ]
        state.output.cloud_top_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "cloud_top_downdraftwindshear"
        ]
        state.input_output.p_forced.data[:] = inputs["p_forced_downdraftwindshear"]
        locals.psum.data[:] = inputs["local_psum_downdraftwindshear"]
        locals.psumh.data[:] = inputs["local_psumh_downdraftwindshear"]
        locals.pwavo.data[:] = inputs["local_pwavo_downdraftwindshear"]
        locals.pwevo.data[:] = inputs["local_pwevo_downdraftwindshear"]
        state.input_output.geopotential_height_forced.data[:] = inputs[
            "geopotential_height_forced_downdraftwindshear"
        ]
        locals.epsilon.data[:] = inputs["local_epsilon_downdraftwindshear"]
        locals.edtc.data[:] = inputs["local_edtc_downdraftwindshear"]
        locals.scale_dependence_factor_downdraft.data[:] = inputs[
            "local_scale_dependence_factor_downdraft_downdraftwindshear"
        ]
        state.output.epsilon.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "epsilon_downdraftwindshear"
        ]

        # initalize test code
        code = DowndraftWindshear(
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
            "error_code_downdraftwindshear": state.output.error_code.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "updraft_lfc_level_downdraftwindshear": state.output.updraft_lfc_level.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "cloud_top_downdraftwindshear": state.output.cloud_top_level.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "u_downdraftwindshear": state.input_output.u.field[:, :, :],
            "v_downdraftwindshear": state.input_output.v.field[:, :, :],
            "geopotential_height_forced_downdraftwindshear": state.input_output.geopotential_height_forced.field[
                :, :, :
            ],
            "local_epsilon_downdraftwindshear": locals.epsilon.data[:],
            "p_forced_downdraftwindshear": state.input_output.p_forced.field[:, :, :],
            "local_pwavo_downdraftwindshear": locals.pwavo.data[:],
            "ccn_downdraftwindshear_in_out": state.input_output.ccn.field[:, :],
            "local_pwevo_downdraftwindshear": locals.pwevo.data[:],
            "local_epsilon_min_downdraftwindshear": locals.epsilon_min.data[:],
            "local_epsilon_max_downdraftwindshear": locals.epsilon_max.data[:],
            "local_psum_downdraftwindshear": locals.psum.data[:],
            "local_psumh_downdraftwindshear": locals.psumh.data[:],
            "epsilon_downdraftwindshear": state.output.epsilon.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_scale_dependence_factor_downdraft_downdraftwindshear": locals.scale_dependence_factor_downdraft.data[
                :
            ],
            "local_edtc_downdraftwindshear": locals.edtc.data[:, :],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_DowndraftWindshear_shallow(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_DowndraftWindshear_mid(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_DowndraftWindshear_deep(TranslateFortranData2Py):
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
