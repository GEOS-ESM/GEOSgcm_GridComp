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
from pyMoist.convection.GF_2020.cumulus_parameterization.large_scale_forcing import LargeScaleForcing
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
            "local_error_code_2": {},
            "local_error_code_3": {},
            "updraft_origin_level": {},
            "updraft_lfc_level": {},
            "cloud_top_level": {},
            "pbl_level": {},
            "local_ocean_fraction": {},
            "p_cloud_levels_forced": {},
            "local_vapor_forced": {},
            "condensate_to_fall_forced": {},
            "local_effective_condensate_to_fall_forced": {},
            "evaporate_in_downdraft_forced": {},
            "omega": {},
            "convective_scale_velocity": {},
            "normalized_massflux_updraft_forced": {},
            "normalized_massflux_downdraft_forced": {},
            "local_cloud_moist_static_energy": {},
            "local_cloud_moist_static_energy_forced": {},
            "local_env_moist_static_energy_cloud_levels": {},
            "local_env_moist_static_energy_cloud_levels_forced": {},
            "local_dmoist_static_energydt": {},
            "local_cloud_workfunction_0": {},
            "local_cloud_workfunction_0_modified": {},
            "local_cloud_workfunction_1": {},
            "local_cloud_workfunction_1_pbl": {},
            "local_arbitrary_numerical_parameter": {},
            "local_f_dicycle_modified": {},
            "local_cape_removal_time_scale": {},
            "epsilon_forced": {},
            "local_k_x_modified": {},
            "local_mass_flux_ensemble": {},
            "local_precipitation_ensemble": {},
            "local_xff_mid": {},
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
        locals.error_code_2.data[:] = inputs["local_error_code_2"]
        locals.error_code_3.data[:] = inputs["local_error_code_3"]
        state.output.updraft_origin_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["updraft_origin_level"] - 1
        )
        state.output.updraft_lfc_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["updraft_lfc_level"] - 1
        )
        state.output.cloud_top_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["cloud_top_level"] - 1
        )
        state.input_output.pbl_level.data[:] = inputs["pbl_level"] - 1
        locals.ocean_fraction.data[:] = inputs["local_ocean_fraction"]
        state.output.p_cloud_levels_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "p_cloud_levels_forced"
        ]
        locals.vapor_forced.data[:] = inputs["local_vapor_forced"]
        state.output.condensate_to_fall_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "condensate_to_fall_forced"
        ]
        locals.effective_condensate_to_fall_forced.data[:, :, :] = inputs[
            "local_effective_condensate_to_fall_forced"
        ]
        state.output.evaporate_in_downdraft_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["evaporate_in_downdraft_forced"]
        )
        state.input_output.omega.data[:] = inputs["omega"]
        state.input_output.convective_scale_velocity.data[:] = inputs["convective_scale_velocity"]
        state.output.normalized_massflux_updraft_forced.data[
            :, :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["normalized_massflux_updraft_forced"]
        state.output.normalized_massflux_downdraft_forced.data[
            :, :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["normalized_massflux_downdraft_forced"]
        locals.cloud_moist_static_energy.data[:] = inputs["local_cloud_moist_static_energy"]
        locals.cloud_moist_static_energy_forced.data[:] = inputs["local_cloud_moist_static_energy_forced"]
        locals.environment_moist_static_energy_cloud_levels.data[:] = inputs[
            "local_env_moist_static_energy_cloud_levels"
        ]
        locals.environment_moist_static_energy_cloud_levels_forced.data[:] = inputs[
            "local_env_moist_static_energy_cloud_levels_forced"
        ]
        locals.dmoist_static_energydt.data[:] = inputs["local_dmoist_static_energydt"]
        locals.cloud_workfunction_0.data[:] = inputs["local_cloud_workfunction_0"]
        locals.cloud_workfunction_0_modified.data[:] = inputs["local_cloud_workfunction_0_modified"]
        locals.cloud_workfunction_1.data[:] = inputs["local_cloud_workfunction_1"]
        locals.cloud_workfunction_1_pbl.data[:] = inputs["local_cloud_workfunction_1_pbl"]
        locals.arbitrary_numerical_parameter.data[:] = inputs["local_arbitrary_numerical_parameter"]
        locals.f_dicycle_modified.data[:] = inputs["local_f_dicycle_modified"]
        locals.cape_removal_time_scale.data[:] = inputs["local_cape_removal_time_scale"]
        state.output.epsilon_forced.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "epsilon_forced"
        ]
        locals.k_x_modified.data[:] = inputs["local_k_x_modified"]
        locals.mass_flux_ensemble.data[:] = inputs["local_mass_flux_ensemble"][:, :, 0:16]
        locals.precipitation_ensemble.data[:] = inputs["local_precipitation_ensemble"][:, :, 0:16]
        locals.xff_mid.data[:] = inputs["local_xff_mid"][:, :, 0:16]

        code = LargeScaleForcing(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        if plume_dependent_constants.ENABLE_PLUME == 1:
            code(
                error_code=state.output.error_code,
                error_code_2=locals.error_code_2,
                error_code_3=locals.error_code_3,
                updraft_origin_level=state.output.updraft_origin_level,
                updraft_lfc_level=state.output.updraft_lfc_level,
                cloud_top_level=state.output.cloud_top_level,
                pbl_level=state.input_output.pbl_level,
                ocean_fraction=locals.ocean_fraction,
                p_cloud_levels_forced=state.output.p_cloud_levels_forced,
                vapor_forced=locals.vapor_forced,
                condensate_to_fall_forced=state.output.condensate_to_fall_forced,
                effective_condensate_to_fall_forced=locals.effective_condensate_to_fall_forced,
                evaporate_in_downdraft_forced=state.output.evaporate_in_downdraft_forced,
                omega=state.input_output.omega,
                convective_scale_velocity=state.input_output.convective_scale_velocity,
                normalized_massflux_updraft_forced=state.output.normalized_massflux_updraft_forced,
                normalized_massflux_downdraft_forced=state.output.normalized_massflux_downdraft_forced,
                cloud_moist_static_energy=locals.cloud_moist_static_energy,
                cloud_moist_static_energy_forced=locals.cloud_moist_static_energy_forced,
                environment_moist_static_energy_cloud_levels=locals.environment_moist_static_energy_cloud_levels,
                environment_moist_static_energy_cloud_levels_forced=locals.environment_moist_static_energy_cloud_levels_forced,
                dmoist_static_energydt=locals.dmoist_static_energydt,
                cloud_workfunction_0=locals.cloud_workfunction_0,
                cloud_workfunction_0_modified=locals.cloud_workfunction_0_modified,
                cloud_workfunction_1=locals.cloud_workfunction_1,
                cloud_workfunction_1_pbl=locals.cloud_workfunction_1_pbl,
                arbitrary_numerical_parameter=locals.arbitrary_numerical_parameter,
                f_dicycle_modified=locals.f_dicycle_modified,
                cape_removal_time_scale=locals.cape_removal_time_scale,
                epsilon_forced=state.output.epsilon_forced,
                k_x_modified=locals.k_x_modified,
                mass_flux_ensemble=locals.mass_flux_ensemble,
                precipitation_ensemble=locals.precipitation_ensemble,
                xff_mid=locals.xff_mid,
                plume_dependent_constants=plume_dependent_constants,
            )

        outputs = {
            "error_code": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "local_error_code_2": locals.error_code_2.field[:],
            "local_error_code_3": locals.error_code_3.field[:],
            "updraft_origin_level": state.output.updraft_origin_level.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ]
            + 1,
            "updraft_lfc_level": state.output.updraft_lfc_level.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ]
            + 1,
            "cloud_top_level": state.output.cloud_top_level.field[:, :, plume_dependent_constants.PLUME_INDEX]
            + 1,
            "pbl_level": state.input_output.pbl_level.field[:] + 1,
            "local_ocean_fraction": locals.ocean_fraction.field[:],
            "p_cloud_levels_forced": state.output.p_cloud_levels_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_vapor_forced": locals.vapor_forced.field[:],
            "condensate_to_fall_forced": state.output.condensate_to_fall_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_effective_condensate_to_fall_forced": locals.effective_condensate_to_fall_forced.field[
                :, :, :
            ],
            "evaporate_in_downdraft_forced": state.output.evaporate_in_downdraft_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "omega": state.input_output.omega.field[:],
            "convective_scale_velocity": state.input_output.convective_scale_velocity.field[:],
            "normalized_massflux_updraft_forced": state.output.normalized_massflux_updraft_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "normalized_massflux_downdraft_forced": state.output.normalized_massflux_downdraft_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_cloud_moist_static_energy": locals.cloud_moist_static_energy.field[:],
            "local_cloud_moist_static_energy_forced": locals.cloud_moist_static_energy_forced.field[:],
            "local_env_moist_static_energy_cloud_levels": locals.environment_moist_static_energy_cloud_levels.field[
                :
            ],
            "local_env_moist_static_energy_cloud_levels_forced": locals.environment_moist_static_energy_cloud_levels_forced.field[
                :
            ],
            "local_dmoist_static_energydt": locals.dmoist_static_energydt.field[:],
            "local_cloud_workfunction_0": locals.cloud_workfunction_0.field[:],
            "local_cloud_workfunction_0_modified": locals.cloud_workfunction_0_modified.field[:],
            "local_cloud_workfunction_1": locals.cloud_workfunction_1.field[:],
            "local_cloud_workfunction_1_pbl": locals.cloud_workfunction_1_pbl.field[:],
            "local_arbitrary_numerical_parameter": locals.arbitrary_numerical_parameter.field[:],
            "local_f_dicycle_modified": locals.f_dicycle_modified.field[:],
            "local_cape_removal_time_scale": locals.cape_removal_time_scale.field[:],
            "epsilon_forced": state.output.epsilon_forced.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "local_k_x_modified": locals.k_x_modified.field[:],
            "local_mass_flux_ensemble": locals.mass_flux_ensemble.field[:],
            "local_precipitation_ensemble": locals.precipitation_ensemble.field[:],
            "local_xff_mid": locals.xff_mid.field[:],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_LargeScaleForcing_shallow(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_LargeScaleForcing_mid(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_LargeScaleForcing_deep(TranslateFortranData2Py):
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
