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
from pyMoist.convection.GF_2020.cumulus_parameterization.precip import rain_evaporation_below_cloud_base

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
            "updraft_lfc_level": {},
            "cloud_top_level": {},
            "ocean_fraction": {},
            "p_cloud_levels_forced": {},
            "p_surface": {},
            "local_t_cloud_levels": {},
            "local_vapor_cloud_levels_forced": {},
            "local_env_saturation_mixing_ratio_cloud_levels": {},
            "epsilon_forced": {},
            "cloud_base_mass_flux_modified": {},
            "condensate_to_fall_forced": {},
            "evaporate_in_downdraft_forced": {},
            "precip": {},
            "local_precipitation_flux": {},
            "local_evaporation_flux": {},
            "local_evaporation_below_cloud_base": {},
            "dtdt": {},
            "dvapordt": {},
            "dbuoyancydt": {},
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
        state.output.updraft_lfc_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["updraft_lfc_level"] - 1
        )
        state.output.cloud_top_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["cloud_top_level"] - 1
        )
        state.input.ocean_fraction.data[:] = inputs["ocean_fraction"]
        state.output.p_cloud_levels_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "p_cloud_levels_forced"
        ]
        state.input_output.p_surface.data[:] = inputs["p_surface"]
        locals.t_cloud_levels.data[:] = inputs["local_t_cloud_levels"]
        locals.vapor_cloud_levels_forced.data[:] = inputs["local_vapor_cloud_levels_forced"]
        locals.environment_saturation_mixing_ratio_cloud_levels.data[:] = inputs[
            "local_env_saturation_mixing_ratio_cloud_levels"
        ]
        state.output.epsilon_forced.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "epsilon_forced"
        ]
        state.output.cloud_base_mass_flux_modified.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "cloud_base_mass_flux_modified"
        ]
        state.output.condensate_to_fall_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "condensate_to_fall_forced"
        ]
        state.output.evaporate_in_downdraft_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["evaporate_in_downdraft_forced"]
        )
        state.output.precip.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["precip"]
        locals.precipitation_flux.data[:] = inputs["local_precipitation_flux"]
        locals.evaporation_flux.data[:] = inputs["local_evaporation_flux"]
        locals.evaporation_below_cloud_base.data[:] = inputs["local_evaporation_below_cloud_base"]
        state.output.dtdt.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs["dtdt"]
        state.output.dvapordt.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs["dvapordt"]
        state.output.dbuoyancydt.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs["dbuoyancydt"]

        code = self.stencil_factory.from_dims_halo(
            func=rain_evaporation_below_cloud_base,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        if plume_dependent_constants.ENABLE_PLUME == 1:
            code(
                error_code=state.output.error_code,
                updraft_lfc_level=state.output.updraft_lfc_level,
                cloud_top_level=state.output.cloud_top_level,
                ocean_fraction=state.input.ocean_fraction,
                p_cloud_levels_forced=state.output.p_cloud_levels_forced,
                p_surface=state.input_output.p_surface,
                t_cloud_levels=locals.t_cloud_levels,
                vapor_cloud_levels_forced=locals.vapor_cloud_levels_forced,
                environment_saturation_mixing_ratio_cloud_levels=locals.environment_saturation_mixing_ratio_cloud_levels,
                epsilon_forced=state.output.epsilon_forced,
                cloud_base_mass_flux_modified=state.output.cloud_base_mass_flux_modified,
                condensate_to_fall_forced=state.output.condensate_to_fall_forced,
                evaporate_in_downdraft_forced=state.output.evaporate_in_downdraft_forced,
                precip=state.output.precip,
                precipitation_flux=locals.precipitation_flux,
                evaporation_flux=locals.evaporation_flux,
                evaporation_below_cloud_base=locals.evaporation_below_cloud_base,
                dtdt=state.output.dtdt,
                dvapordt=state.output.dvapordt,
                dbuoyancydt=state.output.dbuoyancydt,
                plume=plume_dependent_constants.PLUME_INDEX,
            )

        outputs = {
            "error_code": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "updraft_lfc_level": state.output.updraft_lfc_level.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ]
            + 1,
            "cloud_top_level": state.output.cloud_top_level.field[:, :, plume_dependent_constants.PLUME_INDEX]
            + 1,
            "ocean_fraction": state.input.ocean_fraction.field[:],
            "p_cloud_levels_forced": state.output.p_cloud_levels_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "p_surface": state.input_output.p_surface.field[:],
            "local_t_cloud_levels": locals.t_cloud_levels.field[:],
            "local_vapor_cloud_levels_forced": locals.vapor_cloud_levels_forced.field[:],
            "local_env_saturation_mixing_ratio_cloud_levels": locals.environment_saturation_mixing_ratio_cloud_levels.field[
                :
            ],
            "epsilon_forced": state.output.epsilon_forced.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "cloud_base_mass_flux_modified": state.output.cloud_base_mass_flux_modified.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "condensate_to_fall_forced": state.output.condensate_to_fall_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "evaporate_in_downdraft_forced": state.output.evaporate_in_downdraft_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "precip": state.output.precip.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "local_precipitation_flux": locals.precipitation_flux.field[:],
            "local_evaporation_flux": locals.evaporation_flux.field[:],
            "local_evaporation_below_cloud_base": locals.evaporation_below_cloud_base.field[:],
            "dtdt": state.output.dtdt.field[:, :, :, plume_dependent_constants.PLUME_INDEX],
            "dvapordt": state.output.dvapordt.field[:, :, :, plume_dependent_constants.PLUME_INDEX],
            "dbuoyancydt": state.output.dbuoyancydt.field[:, :, :, plume_dependent_constants.PLUME_INDEX],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_RainEvaporationBelowCloudBase_shallow(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_RainEvaporationBelowCloudBase_mid(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_RainEvaporationBelowCloudBase_deep(TranslateFortranData2Py):
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
