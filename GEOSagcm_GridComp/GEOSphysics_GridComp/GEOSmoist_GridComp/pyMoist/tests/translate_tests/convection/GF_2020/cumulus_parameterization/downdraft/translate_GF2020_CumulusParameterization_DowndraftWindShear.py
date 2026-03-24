from f90nml import Namelist

from ndsl import StencilFactory
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
from pyMoist.convection.GF_2020.cumulus_parameterization.downdraft import DowndraftWindShear
from pyMoist.convection.GF_2020.cumulus_parameterization.locals import GF2020CumulusParameterizationLocals
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.setup.set_constants import set_constants
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
            "updraft_lfc_level": {},
            "cloud_top_level": {},
            "geopotential_height_forced": {},
            "p_forced": {},
            "u": {},
            "v": {},
            "ccn": {},
            "total_normalized_integrated_condensate_forced": {},
            "local_total_normalized_integrated_evaporate_forced": {},
            "local_psum": {},
            "local_psumh": {},
            "local_scale_dependence_factor_downdraft": {},
            "local_epsilon": {},
            "local_epsilon_min": {},
            "local_epsilon_max": {},
            "local_epsilon_computed": {},
            "epsilon_forced": {},
        }

        out_vars.update(in_vars["data_vars"])

    def __call__(self, constants: dict, cu_param_constants: dict, plume: str, **inputs):
        # initialize constants
        config = GF2020Config(**constants)
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
        state.output.updraft_lfc_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["updraft_lfc_level"] - 1
        )
        state.output.cloud_top_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["cloud_top_level"] - 1
        )
        state.input_output.geopotential_height_forced.data[:] = inputs["geopotential_height_forced"]
        state.input_output.p_forced.data[:] = inputs["p_forced"]
        state.input_output.u.data[:] = inputs["u"]
        state.input_output.v.data[:] = inputs["v"]
        state.input_output.ccn.data[:] = inputs["ccn"]
        state.output.total_normalized_integrated_condensate_forced.data[
            :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["total_normalized_integrated_condensate_forced"]
        locals.total_normalized_integrated_evaporate_forced.data[:,] = inputs[
            "local_total_normalized_integrated_evaporate_forced"
        ]
        locals.psum.data[:] = inputs["local_psum"]
        locals.psumh.data[:] = inputs["local_psumh"]
        locals.scale_dependence_factor_downdraft.data[:] = inputs["local_scale_dependence_factor_downdraft"]
        locals.epsilon.data[:] = inputs["local_epsilon"]
        locals.epsilon_min.data[:] = inputs["local_epsilon_min"]
        locals.epsilon_max.data[:] = inputs["local_epsilon_max"]
        size_epsilon_computed = len(inputs["local_epsilon_computed"].shape)
        if size_epsilon_computed == 2:
            import numpy as np

            locals.epsilon_computed.data[:] = inputs["local_epsilon_computed"][:, :, np.newaxis]
        else:
            locals.epsilon_computed.data[:] = inputs["local_epsilon_computed"]
        state.output.epsilon_forced.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "epsilon_forced"
        ]

        # initialize test code
        code = DowndraftWindShear(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        # call test code
        if plume_dependent_constants.ENABLE_PLUME == 1:
            code(
                error_code=state.output.error_code,
                updraft_lfc_level=state.output.updraft_lfc_level,
                cloud_top_level=state.output.cloud_top_level,
                geopotential_height_forced=state.input_output.geopotential_height_forced,
                p_forced=state.input_output.p_forced,
                u=state.input_output.u,
                v=state.input_output.v,
                ccn=state.input_output.ccn,
                psum=locals.psum,
                psumh=locals.psumh,
                total_normalized_integrated_condensate_forced=state.output.total_normalized_integrated_condensate_forced,
                total_normalized_integrated_evaporate_forced=locals.total_normalized_integrated_evaporate_forced,
                scale_dependence_factor_downdraft=locals.scale_dependence_factor_downdraft,
                epsilon=locals.epsilon,
                epsilon_min=locals.epsilon_min,
                epsilon_max=locals.epsilon_max,
                epsilon_computed=locals.epsilon_computed,
                epsilon_forced=state.output.epsilon_forced,
                plume_dependent_constants=plume_dependent_constants,
            )

        # write output
        outputs = {
            "error_code": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "updraft_lfc_level": state.output.updraft_lfc_level.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ]
            + 1,
            "cloud_top_level": state.output.cloud_top_level.field[:, :, plume_dependent_constants.PLUME_INDEX]
            + 1,
            "geopotential_height_forced": state.input_output.geopotential_height_forced.field[:],
            "p_forced": state.input_output.p_forced.field[:],
            "u": state.input_output.u.field[:],
            "v": state.input_output.v.field[:],
            "ccn": state.input_output.ccn.field[:],
            "total_normalized_integrated_condensate_forced": state.output.total_normalized_integrated_condensate_forced.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_total_normalized_integrated_evaporate_forced": locals.total_normalized_integrated_evaporate_forced.field[
                :,
            ],
            "local_psum": locals.psum.field[:],
            "local_psumh": locals.psumh.field[:],
            "local_scale_dependence_factor_downdraft": locals.scale_dependence_factor_downdraft.field[:],
            "local_epsilon": locals.epsilon.field[:],
            "local_epsilon_min": locals.epsilon_min.field[:],
            "local_epsilon_max": locals.epsilon_max.field[:],
            "local_epsilon_computed": locals.epsilon_computed.field[:],
            "epsilon_forced": state.output.epsilon_forced.field[:, :, plume_dependent_constants.PLUME_INDEX],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_DowndraftWindShear_shallow(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_DowndraftWindShear_mid(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_DowndraftWindShear_deep(TranslateFortranData2Py):
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
