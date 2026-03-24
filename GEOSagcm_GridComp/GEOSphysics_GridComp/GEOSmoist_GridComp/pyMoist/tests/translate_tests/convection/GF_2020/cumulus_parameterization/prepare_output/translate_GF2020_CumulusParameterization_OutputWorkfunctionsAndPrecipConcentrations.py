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
from pyMoist.convection.GF_2020.cumulus_parameterization.locals import GF2020CumulusParameterizationLocals
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.prepare_output import (
    OutputWorkfunctionsAndPrecipConcentrations,
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
            "cloud_top_level": {},
            "convection_fraction": {},
            "surface_type": {},
            "cloud_workfunction_0": {},
            "cloud_workfunction_1": {},
            "local_cloud_workfunction_0": {},
            "local_cloud_workfunction_1": {},
            "air_density": {},
            "local_updraft_column_temperature_forced": {},
            "dcloudicedt": {},
            "dnliquiddt": {},
            "dnicedt": {},
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
        state.output.cloud_top_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["cloud_top_level"] - 1
        )
        state.input.convection_fraction.data[:] = inputs["convection_fraction"]
        state.input.surface_type.data[:] = inputs["surface_type"]
        state.output.cloud_workfunction_0.data[:] = inputs["cloud_workfunction_0"]
        state.output.cloud_workfunction_1.data[:] = inputs["cloud_workfunction_1"]
        locals.cloud_workfunction_0.data[:] = inputs["local_cloud_workfunction_0"]
        locals.cloud_workfunction_1.data[:] = inputs["local_cloud_workfunction_1"]
        state.input_output.air_density.data[:] = inputs["air_density"]
        locals.updraft_column_temperature_forced.data[:] = inputs["local_updraft_column_temperature_forced"]
        state.output.dcloudicedt.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs["dcloudicedt"]
        state.output.dnliquiddt.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs["dnliquiddt"]
        state.output.dnicedt.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs["dnicedt"]

        code = OutputWorkfunctionsAndPrecipConcentrations(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        if plume_dependent_constants.ENABLE_PLUME == 1:
            code(
                error_code=state.output.error_code,
                cloud_top_level=state.output.cloud_top_level,
                convection_fraction=state.input.convection_fraction,
                surface_type=state.input.surface_type,
                cloud_workfunction_0_output=state.output.cloud_workfunction_0,
                cloud_workfunction_1_output=state.output.cloud_workfunction_1,
                cloud_workfunction_0=locals.cloud_workfunction_0,
                cloud_workfunction_1=locals.cloud_workfunction_1,
                air_density=state.input_output.air_density,
                updraft_column_temperature_forced=locals.updraft_column_temperature_forced,
                dcloudicedt=state.output.dcloudicedt,
                dnliquiddt=state.output.dnliquiddt,
                dnicedt=state.output.dnicedt,
                plume_dependent_constants=plume_dependent_constants,
            )

        outputs = {
            "error_code": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "cloud_top_level": state.output.cloud_top_level.field[:, :, plume_dependent_constants.PLUME_INDEX]
            + 1,
            "convection_fraction": state.input.convection_fraction.field[:],
            "surface_type": state.input.surface_type.field[:],
            "cloud_workfunction_0": state.output.cloud_workfunction_0.field[:],
            "cloud_workfunction_1": state.output.cloud_workfunction_1.field[:],
            "local_cloud_workfunction_0": locals.cloud_workfunction_0.field[:],
            "local_cloud_workfunction_1": locals.cloud_workfunction_1.field[:],
            "air_density": state.input_output.air_density.field[:],
            "local_updraft_column_temperature_forced": locals.updraft_column_temperature_forced.field[:],
            "dcloudicedt": state.output.dcloudicedt.field[:, :, :, plume_dependent_constants.PLUME_INDEX],
            "dnliquiddt": state.output.dnliquiddt.field[:, :, :, plume_dependent_constants.PLUME_INDEX],
            "dnicedt": state.output.dnicedt.field[:, :, :, plume_dependent_constants.PLUME_INDEX],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_OutputWorkfunctionsAndPrecipConcentrations_shallow(
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

        self.test_core = TestCore(grid, namelist, stencil_factory, self.in_vars, self.out_vars)

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF2020-constants")
        self.cu_param_constants = data_loader.load("GF2020_CumulusParameterization-constants")

    def compute_func(self, **inputs):
        outputs = self.test_core(self.constants, self.cu_param_constants, "shallow", **inputs)

        return outputs


class TranslateGF2020_CumulusParameterization_OutputWorkfunctionsAndPrecipConcentrations_mid(
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

        self.test_core = TestCore(grid, namelist, stencil_factory, self.in_vars, self.out_vars)

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF2020-constants")
        self.cu_param_constants = data_loader.load("GF2020_CumulusParameterization-constants")

    def compute_func(self, **inputs):
        outputs = self.test_core(self.constants, self.cu_param_constants, "mid", **inputs)

        return outputs


class TranslateGF2020_CumulusParameterization_OutputWorkfunctionsAndPrecipConcentrations_deep(
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

        self.test_core = TestCore(grid, namelist, stencil_factory, self.in_vars, self.out_vars)

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF2020-constants")
        self.cu_param_constants = data_loader.load("GF2020_CumulusParameterization-constants")

    def compute_func(self, **inputs):
        outputs = self.test_core(self.constants, self.cu_param_constants, "deep", **inputs)

        return outputs
