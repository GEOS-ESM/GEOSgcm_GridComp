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
from pyMoist.convection.GF_2020.cumulus_parameterization.updraft import UpdraftMassFlux
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
            "error_code": {},
            "updraft_origin_level": {},
            "cloud_top_level": {},
            "pbl_level": {},
            "updraft_lfc_level": {},
            "lcl_level": {},
            "p_cloud_levels_forced": {},
            "p_surface": {},
            "ocean_fraction": {},
            "local_normalized_massflux_updraft": {},
            "normalized_massflux_updraft_forced": {},
            "local_normalized_massflux_updraft_modified": {},
            "local_random_number": {},
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
        state.output.updraft_origin_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["updraft_origin_level"] - 1
        )
        state.output.cloud_top_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["cloud_top_level"] - 1
        )
        state.input_output.pbl_level.data[:] = inputs["pbl_level"] - 1
        state.output.updraft_lfc_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["updraft_lfc_level"] - 1
        )
        state.output.lcl_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["lcl_level"] - 1
        state.output.p_cloud_levels_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "p_cloud_levels_forced"
        ]
        state.input_output.p_surface.data[:] = inputs["p_surface"]
        state.input.ocean_fraction.data[:] = inputs["ocean_fraction"]
        locals.normalized_massflux_updraft.data[:] = inputs["local_normalized_massflux_updraft"]
        state.output.normalized_massflux_updraft_forced.data[
            :, :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["normalized_massflux_updraft_forced"]
        locals.normalized_massflux_updraft_modified.data[:] = inputs[
            "local_normalized_massflux_updraft_modified"
        ]
        locals.random_number.data[:] = inputs["local_random_number"]

        # initalize test code
        code = UpdraftMassFlux(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        # call test code
        if plume_dependent_constants.ENABLE_PLUME == 1:
            if cumulus_parameterization_config.FIRST_GUESS_W == 0:
                code(
                    error_code=state.output.error_code,
                    updraft_origin_level=state.output.updraft_origin_level,
                    cloud_top_level=state.output.cloud_top_level,
                    pbl_level=state.input_output.pbl_level,
                    updraft_lfc_level=state.output.updraft_lfc_level,
                    lcl_level=state.output.lcl_level,
                    p_cloud_levels_forced=state.output.p_cloud_levels_forced,
                    p_surface=state.input_output.p_surface,
                    ocean_fraction=state.input.ocean_fraction,
                    normalized_massflux_updraft=locals.normalized_massflux_updraft,
                    normalized_massflux_updraft_forced=state.output.normalized_massflux_updraft_forced,
                    normalized_massflux_updraft_modified=locals.normalized_massflux_updraft_modified,
                    random_number=locals.random_number,
                    plume_dependent_constants=plume_dependent_constants,
                )

        # write output
        outputs = {
            "error_code": state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX],
            "updraft_origin_level": state.output.updraft_origin_level.data[
                :, :, plume_dependent_constants.PLUME_INDEX
            ]
            + 1,
            "cloud_top_level": state.output.cloud_top_level.data[:, :, plume_dependent_constants.PLUME_INDEX]
            + 1,
            "pbl_level": state.input_output.pbl_level.data[:] + 1,
            "updraft_lfc_level": state.output.updraft_lfc_level.data[
                :, :, plume_dependent_constants.PLUME_INDEX
            ]
            + 1,
            "lcl_level": state.output.lcl_level.data[:, :, plume_dependent_constants.PLUME_INDEX] + 1,
            "p_cloud_levels_forced": state.output.p_cloud_levels_forced.data[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "p_surface": state.input_output.p_surface.data[:],
            "ocean_fraction": state.input.ocean_fraction.data[:],
            "local_normalized_massflux_updraft": locals.normalized_massflux_updraft.data[:],
            "normalized_massflux_updraft_forced": state.output.normalized_massflux_updraft_forced.data[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_normalized_massflux_updraft_modified": locals.normalized_massflux_updraft_modified.data[:],
            "local_random_number": locals.random_number.data[:],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_UpdraftMassFlux_shallow(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_UpdraftMassFlux_mid(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_UpdraftMassFlux_deep(TranslateFortranData2Py):
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
