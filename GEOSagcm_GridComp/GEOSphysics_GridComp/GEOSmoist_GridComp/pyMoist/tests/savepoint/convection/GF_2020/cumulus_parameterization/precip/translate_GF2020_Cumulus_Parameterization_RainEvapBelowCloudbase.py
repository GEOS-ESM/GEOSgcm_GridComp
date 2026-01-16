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
from pyMoist.convection.GF_2020.cumulus_parameterization.precip import rain_evap_below_cloudbase

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
            "error_code_rebc": {},
            "updraft_lfc_level_rebc": {},
            "cloud_top_level": {},
            "p_surface": {},
            "ocean_fraction": {},
            "epsilon_forced_rebc": {},
            "cloud_base_mass_flux_rebc": {},
            "p_cloud_levels_forced": {},
            "local_vapor_cloud_levels_forced": {},
            "local_env_saturation_mixing_ratio_cloud_levels_rebc": {},
            "condensate_to_fall_forced_rebc": {},
            "evaporate_in_downdraft_forced": {},
            "precip_rebc": {},
            "t_rebc": {},
            "dvapordt_rebc": {},
            "dbuoyancydt_rebc": {},
            "prec_flux_rebc": {},
            "evap_flux_rebc": {},
            "local_evap_bcb_rebc": {},
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
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["error_code_rebc"]
        state.output.cloud_top_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["cloud_top_level"] - 1
        )
        state.output.updraft_lfc_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["updraft_lfc_level_rebc"] - 1
        )
        state.input_output.p_surface.data[:] = inputs["p_surface"]
        locals.ocean_fraction.data[:] = inputs["ocean_fraction"]
        state.output.epsilon_forced.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "epsilon_forced_rebc"
        ]
        state.output.cloud_base_mass_flux.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "cloud_base_mass_flux_rebc"
        ]
        state.output.p_cloud_levels_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "p_cloud_levels_forced"
        ]
        locals.vapor_cloud_levels_forced.data[:] = inputs["local_vapor_cloud_levels_forced"]
        locals.environment_saturation_mixing_ratio_cloud_levels.data[:] = inputs[
            "local_env_saturation_mixing_ratio_cloud_levels_rebc"
        ]
        state.output.condensate_to_fall_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "condensate_to_fall_forced_rebc"
        ]
        state.output.evaporate_in_downdraft_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["evaporate_in_downdraft_forced"]
        )
        state.output.precip.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["precip_rebc"]
        state.output.t.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs["t_rebc"]
        state.output.dvapordt.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs["dvapordt_rebc"]
        state.output.dbuoyancydt.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "dbuoyancydt_rebc"
        ]
        locals.prec_flux.data[:] = inputs["prec_flux_rebc"]
        locals.evap_flux.data[:] = inputs["evap_flux_rebc"]
        locals.evap_bcb.data[:] = inputs["local_evap_bcb_rebc"]

        code_part_1 = self.stencil_factory.from_dims_halo(
            func=rain_evap_below_cloudbase,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        if plume_dependent_constants.ENABLE_PLUME == 1:
            code_part_1(
                error_code=state.output.error_code,
                plume=plume_dependent_constants.PLUME_INDEX,
                epsilon_forced=state.output.epsilon_forced,
                updraft_lfc_level=state.output.updraft_lfc_level,
                cloud_top_level=state.output.cloud_top_level,
                p_cloud_levels_forced=state.output.p_cloud_levels_forced,
                p_surface=state.input_output.p_surface,
                evaporate_in_downdraft_forced=state.output.evaporate_in_downdraft_forced,
                condensate_to_fall_forced=state.output.condensate_to_fall_forced,
                local_env_saturation_mixing_ratio_cloud_levels=locals.environment_saturation_mixing_ratio_cloud_levels,
                local_vapor_cloud_levels_forced=locals.vapor_cloud_levels_forced,
                ocean_fraction=locals.ocean_fraction,
                cloud_base_mass_flux=state.output.cloud_base_mass_flux,
                local_evap_bcb=locals.evap_bcb,
                evap_flux=locals.evap_flux,
                dbuoyancydt=state.output.dbuoyancydt,
                dvapordt=state.output.dvapordt,
                t=state.output.t,
                precip=state.output.precip,
                prec_flux=locals.prec_flux,
            )

        outputs = {
            "error_code_rebc": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "cloud_top_level": state.output.cloud_top_level.data[:, :, plume_dependent_constants.PLUME_INDEX]
            + 1,
            "updraft_lfc_level_rebc": state.output.updraft_lfc_level.data[
                :, :, plume_dependent_constants.PLUME_INDEX
            ]
            + 1,
            "p_surface": state.input_output.p_surface.data[:],
            "ocean_fraction": locals.ocean_fraction.data[:],
            "epsilon_forced_rebc": state.output.epsilon_forced.data[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "cloud_base_mass_flux_rebc": state.output.cloud_base_mass_flux.data[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "p_cloud_levels_forced": state.output.p_cloud_levels_forced.data[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_vapor_cloud_levels_forced": locals.vapor_cloud_levels_forced.data[:],
            "local_env_saturation_mixing_ratio_cloud_levels_rebc": locals.environment_saturation_mixing_ratio_cloud_levels.data[
                :
            ],
            "condensate_to_fall_forced_rebc": state.output.condensate_to_fall_forced.data[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "evaporate_in_downdraft_forced": state.output.evaporate_in_downdraft_forced.data[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "precip_rebc": state.output.precip.data[:, :, plume_dependent_constants.PLUME_INDEX],
            "t_rebc": state.output.t.data[:, :, :, plume_dependent_constants.PLUME_INDEX],
            "dvapordt_rebc": state.output.dvapordt.data[:, :, :, plume_dependent_constants.PLUME_INDEX],
            "dbuoyancydt_rebc": state.output.dbuoyancydt.data[:, :, :, plume_dependent_constants.PLUME_INDEX],
            "prec_flux_rebc": locals.prec_flux.data[:],
            "evap_flux_rebc": locals.evap_flux.data[:],
            "local_evap_bcb_rebc": locals.evap_bcb.data[:],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_RainEvapBelowCloudbase_shallow(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_RainEvapBelowCloudbase_mid(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_RainEvapBelowCloudbase_deep(TranslateFortranData2Py):
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
