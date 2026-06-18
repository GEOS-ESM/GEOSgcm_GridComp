from f90nml import Namelist
from ndsl import StencilFactory
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py

from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.constants import MAXENS1, MAXENS2, MAXENS3, NUMBER_OF_PLUMES
from pyMoist.convection.GF_2020.cumulus_parameterization.locals import GF2020CumulusParameterizationLocals
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import GF2020PlumeDependentConstants
from pyMoist.convection.GF_2020.cumulus_parameterization.setup.set_constants import set_constants
from pyMoist.convection.GF_2020.cumulus_parameterization.state import GF2020CumulusParameterizationState
from pyMoist.convection.GF_2020.cumulus_parameterization.updraft import UpdraftCIN


class TestCore:
    def __init__(
        self,
        grid: Grid,
        stencil_factory: StencilFactory,
        in_vars: dict,
        out_vars: dict,
    ):
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        in_vars["data_vars"] = {
            "error_code": {},
            "updraft_origin_level": {},
            "updraft_lfc_level": {},
            "cloud_top_level": {},
            "local_geopotential_height_cloud_levels": {},
            "local_geopotential_height_cloud_levels_forced": {},
            "local_normalized_massflux_updraft": {},
            "normalized_massflux_updraft_forced": {},
            "local_d_buoyancy": {},
            "local_d_buoyancy_forced": {},
            "local_gamma_cloud_levels": {},
            "local_gamma_cloud_levels_forced": {},
            "local_t_cloud_levels": {},
            "local_t_cloud_levels_forced": {},
            "local_cin_0": {},
            "local_cin_1": {},
        }

        out_vars.update(in_vars["data_vars"])

    def __call__(self, constants: dict, cu_param_constants: dict, plume: str, **inputs):
        # initialize constants
        config = GF2020Config(**constants)
        cumulus_parameterization_config = GF2020CumulusParameterizationConfig(**cu_param_constants)
        plume_dependent_constants = GF2020PlumeDependentConstants()
        plume_dependent_constants = set_constants(cumulus_parameterization_config, plume_dependent_constants, plume)

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
        state.output.updraft_origin_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["updraft_origin_level"] - 1
        state.output.updraft_lfc_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["updraft_lfc_level"] - 1
        state.output.cloud_top_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["cloud_top_level"] - 1
        locals.geopotential_height_cloud_levels.data[:] = inputs["local_geopotential_height_cloud_levels"]
        locals.geopotential_height_cloud_levels_forced.data[:] = inputs["local_geopotential_height_cloud_levels_forced"]
        locals.normalized_massflux_updraft.data[:] = inputs["local_normalized_massflux_updraft"]
        state.output.normalized_massflux_updraft_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs["normalized_massflux_updraft_forced"]
        locals.d_buoyancy.data[:] = inputs["local_d_buoyancy"]
        locals.d_buoyancy_forced.data[:] = inputs["local_d_buoyancy_forced"]
        locals.gamma_cloud_levels.data[:] = inputs["local_gamma_cloud_levels"]
        locals.gamma_cloud_levels_forced.data[:] = inputs["local_gamma_cloud_levels_forced"]
        locals.t_cloud_levels.data[:] = inputs["local_t_cloud_levels"]
        locals.t_cloud_levels_forced.data[:] = inputs["local_t_cloud_levels_forced"]
        locals.cin_0.data[:] = inputs["local_cin_0"]
        locals.cin_1.data[:] = inputs["local_cin_1"]

        # initialize test code
        code = UpdraftCIN(
            stencil_factory=self.stencil_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        # call test code
        if plume_dependent_constants.ENABLE_PLUME == 1:
            code(
                error_code=state.output.error_code,
                updraft_origin_level=state.output.updraft_origin_level,
                updraft_lfc_level=state.output.updraft_lfc_level,
                cloud_top_level=state.output.cloud_top_level,
                geopotential_height_cloud_levels=locals.geopotential_height_cloud_levels,
                geopotential_height_cloud_levels_forced=locals.geopotential_height_cloud_levels_forced,
                normalized_massflux_updraft=locals.normalized_massflux_updraft,
                normalized_massflux_updraft_forced=state.output.normalized_massflux_updraft_forced,
                d_buoyancy=locals.d_buoyancy,
                d_buoyancy_forced=locals.d_buoyancy_forced,
                gamma_cloud_levels=locals.gamma_cloud_levels,
                gamma_cloud_levels_forced=locals.gamma_cloud_levels_forced,
                t_cloud_levels=locals.t_cloud_levels,
                t_cloud_levels_forced=locals.t_cloud_levels_forced,
                cin_0=locals.cin_0,
                cin_1=locals.cin_1,
                plume_dependent_constants=plume_dependent_constants,
            )

        # write output
        outputs = {
            "error_code": state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX],
            "updraft_origin_level": state.output.updraft_origin_level.data[:, :, plume_dependent_constants.PLUME_INDEX] + 1,
            "updraft_lfc_level": state.output.updraft_lfc_level.data[:, :, plume_dependent_constants.PLUME_INDEX] + 1,
            "cloud_top_level": state.output.cloud_top_level.data[:, :, plume_dependent_constants.PLUME_INDEX] + 1,
            "local_geopotential_height_cloud_levels": locals.geopotential_height_cloud_levels.data[:],
            "local_geopotential_height_cloud_levels_forced": locals.geopotential_height_cloud_levels_forced.data[:],
            "local_normalized_massflux_updraft": locals.normalized_massflux_updraft.data[:],
            "normalized_massflux_updraft_forced": state.output.normalized_massflux_updraft_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX],
            "local_d_buoyancy": locals.d_buoyancy.data[:],
            "local_d_buoyancy_forced": locals.d_buoyancy_forced.data[:],
            "local_gamma_cloud_levels": locals.gamma_cloud_levels.data[:],
            "local_gamma_cloud_levels_forced": locals.gamma_cloud_levels_forced.data[:],
            "local_t_cloud_levels": locals.t_cloud_levels.data[:],
            "local_t_cloud_levels_forced": locals.t_cloud_levels_forced.data[:],
            "local_cin_0": locals.cin_0.data[:],
            "local_cin_1": locals.cin_1.data[:],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_UpdraftCIN_shallow(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        _namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)

        self.test_core = TestCore(grid, stencil_factory, self.in_vars, self.out_vars)

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF2020-constants")
        self.cu_param_constants = data_loader.load("GF2020_CumulusParameterization-constants")

    def compute_func(self, **inputs):
        outputs = self.test_core(self.constants, self.cu_param_constants, "shallow", **inputs)

        return outputs


class TranslateGF2020_CumulusParameterization_UpdraftCIN_mid(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        _namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)

        self.test_core = TestCore(grid, stencil_factory, self.in_vars, self.out_vars)

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF2020-constants")
        self.cu_param_constants = data_loader.load("GF2020_CumulusParameterization-constants")

    def compute_func(self, **inputs):
        outputs = self.test_core(self.constants, self.cu_param_constants, "mid", **inputs)

        return outputs


class TranslateGF2020_CumulusParameterization_UpdraftCIN_deep(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        _namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)

        self.test_core = TestCore(grid, stencil_factory, self.in_vars, self.out_vars)

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF2020-constants")
        self.cu_param_constants = data_loader.load("GF2020_CumulusParameterization-constants")

    def compute_func(self, **inputs):
        outputs = self.test_core(self.constants, self.cu_param_constants, "deep", **inputs)

        return outputs
