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
from pyMoist.convection.GF_2020.cumulus_parameterization.cloud_base_mass_flux.cloud_base_mass_flux import (
    cup_forcing_ens_3d_mid,
    cup_forcing_ens_3d,
    update_condensate_to_fall,
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
            "error_code_cbmf": {},
            "updraft_lfc_level": {},
            "local_cloud_workfunction_0_cbmf": {},
            "local_cloud_workfunction_1_cbmf": {},
            "local_cloud_workfunction_0_modified_cbmf": {},
            "local_cloud_workfunction_1_pbl_cbmf": {},
            "local_arbitrary_numerical_parameter_cbmf": {},
            "local_tau_ecmwf_cbmf": {},
            "local_f_dicycle_modified_cbmf": {},
            "updraft_origin_level": {},
            "pbl_level": {},
            "local_moist_static_energy_cbmf": {},
            "p_cloud_levels_forced": {},
            "local_cloud_moist_static_energy_forced_cbmf": {},
            "local_env_moist_static_energy_cloud_levels_forced_cbmf": {},
            "convective_scale_velocity_cbmf": {},
            "local_ichoice_cbmf": {},
            "local_xff_mid_cbmf": {},
            "local_precipitation_ensemble_cbmf": {},
            "local_xf_ens_cbmf": {},
            "local_moisture_convergence_cbmf": {},
            "local_vapor_forced_cbmf": {},
            "omega_cbmf": {},
            "ocean_fraction": {},
            "error_code2_cbmf": {},
            "error_code3_cbmf": {},
            "effective_condensate_to_fall_forced": {},
            "condensate_to_fall_forced": {},
            "epsilon_forced": {},
            "evaporate_in_downdraft_forced": {},
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
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["error_code_cbmf"]
        state.output.updraft_lfc_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["updraft_lfc_level"] - 1
        )

        locals.cloud_workfunction_0.data[:] = inputs["local_cloud_workfunction_0_cbmf"]
        locals.cloud_workfunction_1.data[:] = inputs["local_cloud_workfunction_1_cbmf"]
        locals.cloud_workfunction_0_modified.data[:] = inputs["local_cloud_workfunction_0_modified_cbmf"]
        locals.cloud_work_function_1_pbl.data[:] = inputs["local_cloud_workfunction_1_pbl_cbmf"]
        locals.arbitrary_numerical_parameter.data[:] = inputs["local_arbitrary_numerical_parameter_cbmf"]
        locals.tau_ecmwf.data[:] = inputs["local_tau_ecmwf_cbmf"]
        locals.f_dicycle_modified.data[:] = inputs["local_f_dicycle_modified_cbmf"]
        state.output.updraft_origin_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["updraft_origin_level"] - 1
        )
        state.input_output.pbl_level.data[:] = inputs["pbl_level"] - 1
        locals.moist_static_energy.data[:] = inputs["local_moist_static_energy_cbmf"]
        state.output.p_cloud_levels_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "p_cloud_levels_forced"
        ]
        locals.cloud_moist_static_energy_forced.data[:] = inputs[
            "local_cloud_moist_static_energy_forced_cbmf"
        ]
        locals.environment_moist_static_energy_cloud_levels_forced.data[:] = inputs[
            "local_env_moist_static_energy_cloud_levels_forced_cbmf"
        ]
        state.input_output.convective_scale_velocity.data[:] = inputs["convective_scale_velocity_cbmf"]
        locals.ichoice.data[:] = inputs["local_ichoice_cbmf"]

        locals.precipitation_ensemble.data[:] = inputs["local_precipitation_ensemble_cbmf"][:, :, 0:16]
        locals.xf_ens.data[:] = inputs["local_xf_ens_cbmf"][:, :, 0:16]
        locals.xff_mid.data[:] = inputs["local_xff_mid_cbmf"][:, :, 0:16]
        locals.moisture_convergence.data[:] = inputs["local_moisture_convergence_cbmf"]
        locals.vapor_forced.data[:] = inputs["local_vapor_forced_cbmf"]
        state.input_output.omega.data[:] = inputs["omega_cbmf"]
        locals.ocean_fraction.data[:] = inputs["ocean_fraction"]
        state.output.error_code2.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "error_code2_cbmf"
        ]
        state.output.error_code3.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "error_code3_cbmf"
        ]
        state.output.effective_condensate_to_fall_forced.data[
            :, :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["effective_condensate_to_fall_forced"]
        state.output.condensate_to_fall_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "condensate_to_fall_forced"
        ]
        state.output.epsilon_forced.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "epsilon_forced"
        ]
        state.output.evaporate_in_downdraft_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["evaporate_in_downdraft_forced"]
        )

        code = self.stencil_factory.from_dims_halo(
            func=cup_forcing_ens_3d_mid,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "DIURNAL_CYCLE": cumulus_parameterization_config.DIURNAL_CYCLE,
                "DTIME": cumulus_parameterization_config.DTIME,
                # "ENSEMBLE_MEMBERS": cumulus_parameterization_config.ENSEMBLE_MEMBERS,
            },
        )
        code2 = self.stencil_factory.from_dims_halo(
            func=cup_forcing_ens_3d,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "DIURNAL_CYCLE": cumulus_parameterization_config.DIURNAL_CYCLE,
                "DTIME": cumulus_parameterization_config.DTIME,
                "ZERO_DIFF": cumulus_parameterization_config.ZERO_DIFF,
                # "ENSEMBLE_MEMBERS": cumulus_parameterization_config.ENSEMBLE_MEMBERS,
            },
        )
        code3 = self.stencil_factory.from_dims_halo(
            func=update_condensate_to_fall,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        if plume_dependent_constants.ENABLE_PLUME == 1:

            if plume_dependent_constants.PLUME_INDEX == 2:
                code2(
                    xf_ens=locals.xf_ens,
                    error_code=state.output.error_code,
                    error_code2=state.output.error_code2,
                    error_code3=state.output.error_code3,
                    plume=plume_dependent_constants.PLUME_INDEX,
                    cloud_workfunction_1=locals.cloud_workfunction_1,
                    cloud_workfunction_0=locals.cloud_workfunction_0,
                    cloud_workfunction_1_pbl=locals.cloud_work_function_1_pbl,
                    xff_ens3=locals.xff_ens3,
                    updraft_lfc_level=state.output.updraft_lfc_level,
                    omega=state.input_output.omega,
                    moisture_convergence=locals.moisture_convergence,
                    vapor_forced=locals.vapor_forced,
                    tau_ecmwf=locals.tau_ecmwf,
                    ichoice=locals.ichoice,
                    cloud_workfunction_0_modified=locals.cloud_workfunction_0_modified,
                    arbitrary_numerical_parameter=locals.arbitrary_numerical_parameter,
                    ocean_fraction=locals.ocean_fraction,
                    precipitation_ensemble=locals.precipitation_ensemble,
                    f_dicycle_modified=locals.f_dicycle_modified,
                )

            if plume_dependent_constants.PLUME_INDEX == 1:
                code(
                    error_code=state.output.error_code,
                    plume=plume_dependent_constants.PLUME_INDEX,
                    k_x=locals.k_x,
                    cloud_workfunction_1=locals.cloud_workfunction_1,
                    cloud_workfunction_0=locals.cloud_workfunction_0,
                    cloud_workfunction_0_modified=locals.cloud_workfunction_0_modified,
                    cloud_workfunction_1_pbl=locals.cloud_work_function_1_pbl,
                    arbitrary_numerical_parameter=locals.arbitrary_numerical_parameter,
                    tau_ecmwf=locals.tau_ecmwf,
                    f_dicycle_modified=locals.f_dicycle_modified,
                    xff_mid=locals.xff_mid,
                    updraft_origin_level=state.output.updraft_origin_level,
                    pbl_level=state.input_output.pbl_level,
                    updraft_lfc_level=state.output.updraft_lfc_level,
                    moist_static_energy=locals.moist_static_energy,
                    p_cloud_levels_forced=state.output.p_cloud_levels_forced,
                    cloud_moist_static_energy_forced=locals.cloud_moist_static_energy_forced,
                    env_moist_static_energy_cloud_levels_forced=locals.environment_moist_static_energy_cloud_levels_forced,
                    convective_scale_velocity=state.input_output.convective_scale_velocity,
                    ichoice=locals.ichoice,
                )
            if plume_dependent_constants.PLUME_INDEX == 0:
                raise NotImplementedError("Plume == 0 (shallow) not expected!!")

            code3(
                error_code=state.output.error_code,
                plume=plume_dependent_constants.PLUME_INDEX,
                effective_condensate_to_fall_forced=state.output.effective_condensate_to_fall_forced,
                condensate_to_fall_forced=state.output.condensate_to_fall_forced,
                epsilon_forced=state.output.epsilon_forced,
                evaporate_in_downdraft_forced=state.output.evaporate_in_downdraft_forced,
            )

        outputs = {
            "error_code_cbmf": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "updraft_lfc_level": state.output.updraft_lfc_level.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ]
            + 1,
            "local_cloud_workfunction_0_cbmf": locals.cloud_workfunction_0.data[:],
            "local_cloud_workfunction_1_cbmf": locals.cloud_workfunction_1.data[:],
            "local_cloud_workfunction_0_modified_cbmf": locals.cloud_workfunction_0_modified.data[:],
            "local_cloud_workfunction_1_pbl_cbmf": locals.cloud_work_function_1_pbl.data[:],
            "local_arbitrary_numerical_parameter_cbmf": locals.arbitrary_numerical_parameter.data[:],
            "local_tau_ecmwf_cbmf": locals.tau_ecmwf.data[:],
            "local_f_dicycle_modified_cbmf": locals.f_dicycle_modified.data[:],
            "updraft_origin_level": state.output.updraft_origin_level.data[
                :, :, plume_dependent_constants.PLUME_INDEX
            ]
            + 1,
            "pbl_level": state.input_output.pbl_level.data[:] + 1,
            "local_moist_static_energy_cbmf": locals.moist_static_energy.data[:],
            "p_cloud_levels_forced": state.output.p_cloud_levels_forced.data[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_cloud_moist_static_energy_forced_cbmf": locals.cloud_moist_static_energy_forced.data[:],
            "local_env_moist_static_energy_cloud_levels_forced_cbmf": locals.environment_moist_static_energy_cloud_levels_forced.data[
                :
            ],
            "convective_scale_velocity_cbmf": state.input_output.convective_scale_velocity.data[:],
            "local_ichoice_cbmf": locals.ichoice.data[:],
            "local_xff_mid_cbmf": locals.xff_mid.data[:],
            "local_precipitation_ensemble_cbmf": locals.precipitation_ensemble.data[:],
            "local_xf_ens_cbmf": locals.xf_ens.data[:],
            "local_moisture_convergence_cbmf": locals.moisture_convergence.data[:],
            "local_vapor_forced_cbmf": locals.vapor_forced.data[:],
            "omega_cbmf": state.input_output.omega.data[:],
            "ocean_fraction": locals.ocean_fraction.data[:],
            "error_code2_cbmf": state.output.error_code2.data[:, :, plume_dependent_constants.PLUME_INDEX],
            "error_code3_cbmf": state.output.error_code3.data[:, :, plume_dependent_constants.PLUME_INDEX],
            "effective_condensate_to_fall_forced": state.output.effective_condensate_to_fall_forced.data[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "condensate_to_fall_forced": state.output.condensate_to_fall_forced.data[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "epsilon_forced": state.output.epsilon_forced.data[:, :, plume_dependent_constants.PLUME_INDEX],
            "evaporate_in_downdraft_forced": state.output.evaporate_in_downdraft_forced.data[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_CloudBaseMassFlux_shallow(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_CloudBaseMassFlux_mid(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_CloudBaseMassFlux_deep(TranslateFortranData2Py):
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
