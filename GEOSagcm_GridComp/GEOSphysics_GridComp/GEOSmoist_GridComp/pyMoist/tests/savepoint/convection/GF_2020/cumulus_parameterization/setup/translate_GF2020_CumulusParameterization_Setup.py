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
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.convection.GF_2020.cumulus_parameterization.setup.setup import Setup
from pyMoist.convection.GF_2020.cumulus_parameterization.constants import (
    MAXENS1,
    MAXENS2,
    MAXENS3,
)
import numpy as np


class TranslateGF2020_CumulusParameterization_Setup_shallow(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, namelist, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        self.in_vars["data_vars"] = {
            "t_excess_cu_param_setup": {},
            "vapor_excess_cu_param_setup": {},
            "ocean_fraction_cu_param_setup": {},
            "t_old_cu_param_setup": {},
            "vapor_old_cu_param_setup": {},
            "grid_scale_forcing_t_cu_param_setup": {},
            "grid_scale_forcing_vapor_cu_param_setup": {},
            "subgrid_scale_forcing_t_cu_param_setup": {},
            "subgrid_scale_forcing_vapor_cu_param_setup": {},
            "geopotential_height_forced_cu_param_setup": {},
            "epsilon_cu_param_setup": {},
            "precip_cu_param_setup": {},
            "scale_dependence_factor_cu_param_setup": {},
            "lightning_density_cu_param_setup": {},
            "seed_convection_cu_param_setup": {},
            "error_code_cu_param_setup": {},
            "grid_length_cu_param_setup": {},
            "lateral_entrainment_rate_cu_param_setup": {},
            "entrainment_rate_cu_param_setup": {},
        }

        self.out_vars = self.in_vars["data_vars"].copy()
        del (
            self.out_vars["t_excess_cu_param_setup"],
            self.out_vars["vapor_excess_cu_param_setup"],
            self.out_vars["ocean_fraction_cu_param_setup"],
            self.out_vars["t_old_cu_param_setup"],
            self.out_vars["vapor_old_cu_param_setup"],
            self.out_vars["grid_scale_forcing_t_cu_param_setup"],
            self.out_vars["grid_scale_forcing_vapor_cu_param_setup"],
            self.out_vars["subgrid_scale_forcing_t_cu_param_setup"],
            self.out_vars["subgrid_scale_forcing_vapor_cu_param_setup"],
            self.out_vars["geopotential_height_forced_cu_param_setup"],
            self.out_vars["seed_convection_cu_param_setup"],
        )
        self.out_vars.update(
            {
                "local_t_excess_cu_param_setup": {},
                "local_vapor_excess_cu_param_setup": {},
                "local_t_new_cu_param_setup": {},
                "local_vapor_forced_cu_param_setup": {},
                "local_t_new_pbl_cu_param_setup": {},
                "local_vapor_forced_pbl_cu_param_setup": {},
                "local_moist_static_energy_cu_param_setup": {},
                "local_maximum_updraft_origin_level_cu_param_setup": {},
                "local_kstabm_cu_param_setup": {},
                "local_ocean_fraction_cu_param_setup": {},
                "local_cap_max_cu_param_setup": {},
                "local_error_code_2_cu_param_setup": {},
                "local_error_code_3_cu_param_setup": {},
                "local_cap_max_increment_cu_param_setup": {},
                "local_geopotential_height_cu_param_setup": {},
                "local_geopotential_height_modified_cu_param_setup": {},
                "local_cloud_work_function_0_cu_param_setup": {},
                "local_cloud_work_function_1_cu_param_setup": {},
                "local_cloud_work_function_2_cu_param_setup": {},
                "local_cloud_work_function_3_cu_param_setup": {},
                "local_cloud_work_function_0_pbl_cu_param_setup": {},
                "local_cloud_work_function_1_pbl_cu_param_setup": {},
                "local_cloud_work_function_1_fa_cu_param_setup": {},
                "local_cin1_cu_param_setup": {},
                "local_k_x_modified_cu_param_setup": {},
                "local_epsilon_cu_param_setup": {},
                "local_pbl_time_scale_cu_param_setup": {},
                "local_t_wetbulb_cu_param_setup": {},
                "local_vapor_wetbulb_cu_param_setup": {},
                "local_tau_ecmwf_cu_param_setup": {},
                "local_f_dicycle_modified_cu_param_setup": {},
                "local_add_buoy_cu_param_setup": {},
                "local_hcdo_cu_param_setup": {},
                "local_cupclw_cu_param_setup": {},
                "local_qrcdo_cu_param_setup": {},
                "local_cloud_moist_static_energy_forced_t_cu_param_setup": {},
                "local_c1d_cu_param_setup": {},
                "local_evap_bcb_cu_param_setup": {},
                "local_random_number_cu_param_setup": {},
                "local_updraft_detrainment_function_cu_param_setup": {},
                "local_epsilon_min_cu_param_setup": {},
                "local_epsilon_max_cu_param_setup": {},
                "local_arbitrary_numerical_parameter_cu_param_setup": {},
            }
        )

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF2020-constants")
        self.cu_param_constants = data_loader.load("GF2020_CumulusParameterization-constants")

    def compute_func(self, **inputs):
        # initalize constants
        config = GF2020Config(SINGLE_COLUMN_MODE=False, **self.constants)
        cumulus_parameterization_config = GF2020CumulusParameterizationConfig(**self.cu_param_constants)
        plume_dependent_constants = GF2020PlumeDependentConstants()

        # initalize saturation tables
        saturation_tables = SaturationVaporPressureTable(self.stencil_factory.backend)

        # initalize dataclasses
        state = GF2020CumulusParameterizationState.zeros(
            self.quantity_factory,
            data_dimensions={
                "plumes": 3,
            },
        )

        locals = GF2020CumulusParameterizationLocals.zeros(
            self.quantity_factory,
            data_dimensions={
                "ensemble_members": MAXENS1 * MAXENS2 * MAXENS3,
            },
        )

        # fill relevant parts of dataclasses
        state.input.t_excess.data[:] = inputs["t_excess_cu_param_setup"]
        state.input.vapor_excess.data[:] = inputs["vapor_excess_cu_param_setup"]
        state.input.ocean_fraction.data[:] = inputs["ocean_fraction_cu_param_setup"]
        state.input_output.t_old.data[:] = inputs["t_old_cu_param_setup"]
        state.input_output.vapor_old.data[:] = inputs["vapor_old_cu_param_setup"]
        state.input.grid_scale_forcing_t.data[:] = inputs["grid_scale_forcing_t_cu_param_setup"]
        state.input.grid_scale_forcing_vapor.data[:] = inputs["grid_scale_forcing_vapor_cu_param_setup"]
        state.input.subgrid_scale_forcing_t.data[:] = inputs["subgrid_scale_forcing_t_cu_param_setup"]
        state.input.subgrid_scale_forcing_vapor.data[:] = inputs["subgrid_scale_forcing_vapor_cu_param_setup"]
        state.input_output.geopotential_height_forced.data[:] = inputs[
            "geopotential_height_forced_cu_param_setup"
        ]
        state.output.epsilon.data[:, :, 0] = inputs["epsilon_cu_param_setup"]  # plume dependent
        state.output.precip.data[:, :, 0] = inputs["precip_cu_param_setup"]  # plume dependent
        state.output.scale_dependence_factor.data[:, :, 0] = inputs[
            "scale_dependence_factor_cu_param_setup"
        ]  # plume dependent
        state.output.lightning_density.data[:] = inputs["lightning_density_cu_param_setup"]
        state.input.seed_convection.data[:] = inputs["seed_convection_cu_param_setup"]
        state.output.error_code.data[:, :, 0] = inputs["error_code_cu_param_setup"]  # plume dependent
        state.input_output.grid_length.data[:] = inputs["grid_length_cu_param_setup"]
        state.input.lateral_entrainment_rate.data[:] = inputs["lateral_entrainment_rate_cu_param_setup"]
        state.output.entrainment_rate.data[:, :, :, 0] = inputs[
            "entrainment_rate_cu_param_setup"
        ]  # plume dependent

        setup = Setup(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        setup(
            state=state,
            locals=locals,
            saturation_tables=saturation_tables,
            plume_dependent_constants=plume_dependent_constants,
            plume="shallow",
        )

        # fortran initalized to nan and never touches the top level of some variables
        # set those rows to nan here so the tests pass - even though they are zero in python
        # done like this to avoid introduce nans into the python
        locals.t_new.field[:, :, -1] = np.nan
        locals.vapor_forced.field[:, :, -1] = np.nan
        locals.t_new_pbl.field[:, :, -1] = np.nan
        locals.vapor_forced_pbl.field[:, :, -1] = np.nan
        locals.moist_static_energy.field[:, :, -1] = np.nan
        locals.detrainment_function_updraft.field[:, :, -1] = np.nan

        outputs = {
            # state fields
            "epsilon_cu_param_setup": state.output.epsilon.field[:, :, 0],
            "precip_cu_param_setup": state.output.precip.field[:, :, 0],
            "scale_dependence_factor_cu_param_setup": state.output.scale_dependence_factor.field[:, :, 0],
            "lightning_density_cu_param_setup": state.output.lightning_density.field[:],
            "error_code_cu_param_setup": state.output.error_code.field[:, :, 0],
            "grid_length_cu_param_setup": state.input_output.grid_length.field[:],
            "lateral_entrainment_rate_cu_param_setup": state.input.lateral_entrainment_rate.field[:],
            "entrainment_rate_cu_param_setup": state.output.entrainment_rate.field[:, :, :, 0],
            # local fields
            "local_t_excess_cu_param_setup": locals.t_excess.field[:],
            "local_vapor_excess_cu_param_setup": locals.vapor_excess.field[:],
            "local_t_new_cu_param_setup": locals.t_new.field[:],
            "local_vapor_forced_cu_param_setup": locals.vapor_forced.field[:],
            "local_t_new_pbl_cu_param_setup": locals.t_new_pbl.field[:],
            "local_vapor_forced_pbl_cu_param_setup": locals.vapor_forced_pbl.field[:],
            "local_moist_static_energy_cu_param_setup": locals.moist_static_energy.field[:],
            "local_maximum_updraft_origin_level_cu_param_setup": locals.maximum_updraft_origin_level.field[:],
            "local_kstabm_cu_param_setup": locals.kstabm.field[:] + 1,  # +1 b/c python counts from 0
            "local_ocean_fraction_cu_param_setup": locals.ocean_fraction.field[:],
            "local_cap_max_cu_param_setup": locals.cap_max.field[:],
            "local_error_code_2_cu_param_setup": locals.error_code_2.field[:],
            "local_error_code_3_cu_param_setup": locals.error_code_3.field[:],
            "local_cap_max_increment_cu_param_setup": locals.cap_max_increment.field[:],
            "local_geopotential_height_cu_param_setup": locals.geopotential_height.field[:],
            "local_geopotential_height_modified_cu_param_setup": locals.geopotential_height_modified.field[:],
            "local_cloud_work_function_0_cu_param_setup": locals.cloud_work_function_0.field[:],
            "local_cloud_work_function_1_cu_param_setup": locals.cloud_work_function_1.field[:],
            "local_cloud_work_function_2_cu_param_setup": locals.cloud_work_function_2.field[:],
            "local_cloud_work_function_3_cu_param_setup": locals.cloud_work_function_3.field[:],
            "local_cloud_work_function_0_pbl_cu_param_setup": locals.cloud_work_function_0_pbl.field[:],
            "local_cloud_work_function_1_pbl_cu_param_setup": locals.cloud_work_function_1_pbl.field[:],
            "local_cloud_work_function_1_fa_cu_param_setup": locals.cloud_work_function_1_fa.field[:],
            "local_cin1_cu_param_setup": locals.cin1.field[:],
            "local_k_x_modified_cu_param_setup": locals.k_x_modified.field[:],
            "local_epsilon_cu_param_setup": locals.epsilon.field[:],
            "local_pbl_time_scale_cu_param_setup": locals.pbl_time_scale.field[:],
            "local_t_wetbulb_cu_param_setup": locals.t_wetbulb.field[:],
            "local_vapor_wetbulb_cu_param_setup": locals.vapor_wetbulb.field[:],
            "local_tau_ecmwf_cu_param_setup": locals.tau_ecmwf.field[:],
            "local_f_dicycle_modified_cu_param_setup": locals.f_dicycle_modified.field[:],
            "local_add_buoy_cu_param_setup": locals.add_buoyancy.field[:],
            "local_hcdo_cu_param_setup": locals.hcdo.field[:],
            "local_cupclw_cu_param_setup": locals.cupclw.field[:],
            "local_qrcdo_cu_param_setup": locals.qrcdo.field[:],
            "local_cloud_moist_static_energy_forced_t_cu_param_setup": locals.cloud_moist_static_energy_forced_transported.field[
                :
            ],
            "local_c1d_cu_param_setup": locals.c1d.field[:],
            "local_evap_bcb_cu_param_setup": locals.evap_bcb.field[:],
            "local_mass_flux_ensemble_cu_param_setup": locals.mass_flux_ensemble.field[:],
            "local_precipitation_ensemble_cu_param_setup": locals.precipitation_ensemble.field[:],
            "local_downdraft_scale_dependence_factor_cu_param_setup": locals.scale_dependence_factor_downdraft.field[
                :
            ],
            "local_random_number_cu_param_setup": locals.random_number.field[:],
            "local_updraft_detrainment_function_cu_param_setup": locals.detrainment_function_updraft.field[:],
            "local_epsilon_min_cu_param_setup": locals.epsilon_min.field[:],
            "local_epsilon_max_cu_param_setup": locals.epsilon_max.field[:],
            "local_arbitrary_numerical_parameter_cu_param_setup": locals.arbitrary_numerical_parameter.field[
                :
            ],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_Setup_mid(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, namelist, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        self.in_vars["data_vars"] = {
            "t_excess_cu_param_setup": {},
            "vapor_excess_cu_param_setup": {},
            "ocean_fraction_cu_param_setup": {},
            "t_old_cu_param_setup": {},
            "vapor_old_cu_param_setup": {},
            "grid_scale_forcing_t_cu_param_setup": {},
            "grid_scale_forcing_vapor_cu_param_setup": {},
            "subgrid_scale_forcing_t_cu_param_setup": {},
            "subgrid_scale_forcing_vapor_cu_param_setup": {},
            "geopotential_height_forced_cu_param_setup": {},
            "epsilon_cu_param_setup": {},
            "precip_cu_param_setup": {},
            "scale_dependence_factor_cu_param_setup": {},
            "lightning_density_cu_param_setup": {},
            "seed_convection_cu_param_setup": {},
            "error_code_cu_param_setup": {},
            "grid_length_cu_param_setup": {},
            "lateral_entrainment_rate_cu_param_setup": {},
            "entrainment_rate_cu_param_setup": {},
        }

        self.out_vars = self.in_vars["data_vars"].copy()
        del (
            self.out_vars["t_excess_cu_param_setup"],
            self.out_vars["vapor_excess_cu_param_setup"],
            self.out_vars["ocean_fraction_cu_param_setup"],
            self.out_vars["t_old_cu_param_setup"],
            self.out_vars["vapor_old_cu_param_setup"],
            self.out_vars["grid_scale_forcing_t_cu_param_setup"],
            self.out_vars["grid_scale_forcing_vapor_cu_param_setup"],
            self.out_vars["subgrid_scale_forcing_t_cu_param_setup"],
            self.out_vars["subgrid_scale_forcing_vapor_cu_param_setup"],
            self.out_vars["geopotential_height_forced_cu_param_setup"],
            self.out_vars["seed_convection_cu_param_setup"],
        )
        self.out_vars.update(
            {
                "local_t_excess_cu_param_setup": {},
                "local_vapor_excess_cu_param_setup": {},
                "local_t_new_cu_param_setup": {},
                "local_vapor_forced_cu_param_setup": {},
                "local_t_new_pbl_cu_param_setup": {},
                "local_vapor_forced_pbl_cu_param_setup": {},
                "local_moist_static_energy_cu_param_setup": {},
                "local_maximum_updraft_origin_level_cu_param_setup": {},
                "local_kstabm_cu_param_setup": {},
                "local_ocean_fraction_cu_param_setup": {},
                "local_cap_max_cu_param_setup": {},
                "local_error_code_2_cu_param_setup": {},
                "local_error_code_3_cu_param_setup": {},
                "local_cap_max_increment_cu_param_setup": {},
                "local_geopotential_height_cu_param_setup": {},
                "local_geopotential_height_modified_cu_param_setup": {},
                "local_cloud_work_function_0_cu_param_setup": {},
                "local_cloud_work_function_1_cu_param_setup": {},
                "local_cloud_work_function_2_cu_param_setup": {},
                "local_cloud_work_function_3_cu_param_setup": {},
                "local_cloud_work_function_0_pbl_cu_param_setup": {},
                "local_cloud_work_function_1_pbl_cu_param_setup": {},
                "local_cloud_work_function_1_fa_cu_param_setup": {},
                "local_cin1_cu_param_setup": {},
                "local_k_x_modified_cu_param_setup": {},
                "local_epsilon_cu_param_setup": {},
                "local_pbl_time_scale_cu_param_setup": {},
                "local_t_wetbulb_cu_param_setup": {},
                "local_vapor_wetbulb_cu_param_setup": {},
                "local_tau_ecmwf_cu_param_setup": {},
                "local_f_dicycle_modified_cu_param_setup": {},
                "local_add_buoy_cu_param_setup": {},
                "local_hcdo_cu_param_setup": {},
                "local_cupclw_cu_param_setup": {},
                "local_qrcdo_cu_param_setup": {},
                "local_cloud_moist_static_energy_forced_t_cu_param_setup": {},
                "local_c1d_cu_param_setup": {},
                "local_evap_bcb_cu_param_setup": {},
                "local_random_number_cu_param_setup": {},
                "local_updraft_detrainment_function_cu_param_setup": {},
                "local_epsilon_min_cu_param_setup": {},
                "local_epsilon_max_cu_param_setup": {},
                "local_arbitrary_numerical_parameter_cu_param_setup": {},
            }
        )

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF2020-constants")
        self.cu_param_constants = data_loader.load("GF2020_CumulusParameterization-constants")

    def compute_func(self, **inputs):
        # initalize constants
        config = GF2020Config(SINGLE_COLUMN_MODE=False, **self.constants)
        cumulus_parameterization_config = GF2020CumulusParameterizationConfig(**self.cu_param_constants)
        plume_dependent_constants = GF2020PlumeDependentConstants()

        # initalize saturation tables
        saturation_tables = SaturationVaporPressureTable(self.stencil_factory.backend)

        # initalize dataclasses
        state = GF2020CumulusParameterizationState.zeros(
            self.quantity_factory,
            data_dimensions={
                "plumes": 3,
            },
        )

        locals = GF2020CumulusParameterizationLocals.zeros(
            self.quantity_factory,
            data_dimensions={
                "ensemble_members": MAXENS1 * MAXENS2 * MAXENS3,
            },
        )

        # fill relevant parts of dataclasses
        state.input.t_excess.data[:] = inputs["t_excess_cu_param_setup"]
        state.input.vapor_excess.data[:] = inputs["vapor_excess_cu_param_setup"]
        state.input.ocean_fraction.data[:] = inputs["ocean_fraction_cu_param_setup"]
        state.input_output.t_old.data[:] = inputs["t_old_cu_param_setup"]
        state.input_output.vapor_old.data[:] = inputs["vapor_old_cu_param_setup"]
        state.input.grid_scale_forcing_t.data[:] = inputs["grid_scale_forcing_t_cu_param_setup"]
        state.input.grid_scale_forcing_vapor.data[:] = inputs["grid_scale_forcing_vapor_cu_param_setup"]
        state.input.subgrid_scale_forcing_t.data[:] = inputs["subgrid_scale_forcing_t_cu_param_setup"]
        state.input.subgrid_scale_forcing_vapor.data[:] = inputs["subgrid_scale_forcing_vapor_cu_param_setup"]
        state.input_output.geopotential_height_forced.data[:] = inputs[
            "geopotential_height_forced_cu_param_setup"
        ]
        state.output.epsilon.data[:, :, 1] = inputs["epsilon_cu_param_setup"]  # plume dependent
        state.output.precip.data[:, :, 1] = inputs["precip_cu_param_setup"]  # plume dependent
        state.output.scale_dependence_factor.data[:, :, 1] = inputs[
            "scale_dependence_factor_cu_param_setup"
        ]  # plume dependent
        state.output.lightning_density.data[:] = inputs["lightning_density_cu_param_setup"]
        state.input.seed_convection.data[:] = inputs["seed_convection_cu_param_setup"]
        state.output.error_code.data[:, :, 1] = inputs["error_code_cu_param_setup"]  # plume dependent
        state.input_output.grid_length.data[:] = inputs["grid_length_cu_param_setup"]
        state.input.lateral_entrainment_rate.data[:] = inputs["lateral_entrainment_rate_cu_param_setup"]
        state.output.entrainment_rate.data[:, :, :, 1] = inputs[
            "entrainment_rate_cu_param_setup"
        ]  # plume dependent

        setup = Setup(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        setup(
            state=state,
            locals=locals,
            saturation_tables=saturation_tables,
            plume_dependent_constants=plume_dependent_constants,
            plume="mid",
        )

        # fortran initalized to nan and never touches the top level of some variables
        # set those rows to nan here so the tests pass - even though they are zero in python
        # done like this to avoid introduce nans into the python
        locals.t_new.field[:, :, -1] = np.nan
        locals.vapor_forced.field[:, :, -1] = np.nan
        locals.t_new_pbl.field[:, :, -1] = np.nan
        locals.vapor_forced_pbl.field[:, :, -1] = np.nan
        locals.moist_static_energy.field[:, :, -1] = np.nan
        locals.detrainment_function_updraft.field[:, :, -1] = np.nan

        outputs = {
            # state fields
            "epsilon_cu_param_setup": state.output.epsilon.field[:, :, 1],
            "precip_cu_param_setup": state.output.precip.field[:, :, 1],
            "scale_dependence_factor_cu_param_setup": state.output.scale_dependence_factor.field[:, :, 1],
            "lightning_density_cu_param_setup": state.output.lightning_density.field[:],
            "error_code_cu_param_setup": state.output.error_code.field[:, :, 1],
            "grid_length_cu_param_setup": state.input_output.grid_length.field[:],
            "lateral_entrainment_rate_cu_param_setup": state.input.lateral_entrainment_rate.field[:],
            "entrainment_rate_cu_param_setup": state.output.entrainment_rate.field[:, :, :, 1],
            # local fields
            "local_t_excess_cu_param_setup": locals.t_excess.field[:],
            "local_vapor_excess_cu_param_setup": locals.vapor_excess.field[:],
            "local_t_new_cu_param_setup": locals.t_new.field[:],
            "local_vapor_forced_cu_param_setup": locals.vapor_forced.field[:],
            "local_t_new_pbl_cu_param_setup": locals.t_new_pbl.field[:],
            "local_vapor_forced_pbl_cu_param_setup": locals.vapor_forced_pbl.field[:],
            "local_moist_static_energy_cu_param_setup": locals.moist_static_energy.field[:],
            "local_maximum_updraft_origin_level_cu_param_setup": locals.maximum_updraft_origin_level.field[:],
            "local_kstabm_cu_param_setup": locals.kstabm.field[:] + 1,  # +1 b/c python counts from 0
            "local_ocean_fraction_cu_param_setup": locals.ocean_fraction.field[:],
            "local_cap_max_cu_param_setup": locals.cap_max.field[:],
            "local_error_code_2_cu_param_setup": locals.error_code_2.field[:],
            "local_error_code_3_cu_param_setup": locals.error_code_3.field[:],
            "local_cap_max_increment_cu_param_setup": locals.cap_max_increment.field[:],
            "local_geopotential_height_cu_param_setup": locals.geopotential_height.field[:],
            "local_geopotential_height_modified_cu_param_setup": locals.geopotential_height_modified.field[:],
            "local_cloud_work_function_0_cu_param_setup": locals.cloud_work_function_0.field[:],
            "local_cloud_work_function_1_cu_param_setup": locals.cloud_work_function_1.field[:],
            "local_cloud_work_function_2_cu_param_setup": locals.cloud_work_function_2.field[:],
            "local_cloud_work_function_3_cu_param_setup": locals.cloud_work_function_3.field[:],
            "local_cloud_work_function_0_pbl_cu_param_setup": locals.cloud_work_function_0_pbl.field[:],
            "local_cloud_work_function_1_pbl_cu_param_setup": locals.cloud_work_function_1_pbl.field[:],
            "local_cloud_work_function_1_fa_cu_param_setup": locals.cloud_work_function_1_fa.field[:],
            "local_cin1_cu_param_setup": locals.cin1.field[:],
            "local_k_x_modified_cu_param_setup": locals.k_x_modified.field[:],
            "local_epsilon_cu_param_setup": locals.epsilon.field[:],
            "local_pbl_time_scale_cu_param_setup": locals.pbl_time_scale.field[:],
            "local_t_wetbulb_cu_param_setup": locals.t_wetbulb.field[:],
            "local_vapor_wetbulb_cu_param_setup": locals.vapor_wetbulb.field[:],
            "local_tau_ecmwf_cu_param_setup": locals.tau_ecmwf.field[:],
            "local_f_dicycle_modified_cu_param_setup": locals.f_dicycle_modified.field[:],
            "local_add_buoy_cu_param_setup": locals.add_buoyancy.field[:],
            "local_hcdo_cu_param_setup": locals.hcdo.field[:],
            "local_cupclw_cu_param_setup": locals.cupclw.field[:],
            "local_qrcdo_cu_param_setup": locals.qrcdo.field[:],
            "local_cloud_moist_static_energy_forced_t_cu_param_setup": locals.cloud_moist_static_energy_forced_transported.field[
                :
            ],
            "local_c1d_cu_param_setup": locals.c1d.field[:],
            "local_evap_bcb_cu_param_setup": locals.evap_bcb.field[:],
            "local_mass_flux_ensemble_cu_param_setup": locals.mass_flux_ensemble.field[:],
            "local_precipitation_ensemble_cu_param_setup": locals.precipitation_ensemble.field[:],
            "local_downdraft_scale_dependence_factor_cu_param_setup": locals.scale_dependence_factor_downdraft.field[
                :
            ],
            "local_random_number_cu_param_setup": locals.random_number.field[:],
            "local_updraft_detrainment_function_cu_param_setup": locals.detrainment_function_updraft.field[:],
            "local_epsilon_min_cu_param_setup": locals.epsilon_min.field[:],
            "local_epsilon_max_cu_param_setup": locals.epsilon_max.field[:],
            "local_arbitrary_numerical_parameter_cu_param_setup": locals.arbitrary_numerical_parameter.field[
                :
            ],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_Setup_deep(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, namelist, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        self.in_vars["data_vars"] = {
            "t_excess_cu_param_setup": {},
            "vapor_excess_cu_param_setup": {},
            "ocean_fraction_cu_param_setup": {},
            "t_old_cu_param_setup": {},
            "vapor_old_cu_param_setup": {},
            "grid_scale_forcing_t_cu_param_setup": {},
            "grid_scale_forcing_vapor_cu_param_setup": {},
            "subgrid_scale_forcing_t_cu_param_setup": {},
            "subgrid_scale_forcing_vapor_cu_param_setup": {},
            "geopotential_height_forced_cu_param_setup": {},
            "epsilon_cu_param_setup": {},
            "precip_cu_param_setup": {},
            "scale_dependence_factor_cu_param_setup": {},
            "lightning_density_cu_param_setup": {},
            "seed_convection_cu_param_setup": {},
            "error_code_cu_param_setup": {},
            "grid_length_cu_param_setup": {},
            "lateral_entrainment_rate_cu_param_setup": {},
            "entrainment_rate_cu_param_setup": {},
        }

        self.out_vars = self.in_vars["data_vars"].copy()
        del (
            self.out_vars["t_excess_cu_param_setup"],
            self.out_vars["vapor_excess_cu_param_setup"],
            self.out_vars["ocean_fraction_cu_param_setup"],
            self.out_vars["t_old_cu_param_setup"],
            self.out_vars["vapor_old_cu_param_setup"],
            self.out_vars["grid_scale_forcing_t_cu_param_setup"],
            self.out_vars["grid_scale_forcing_vapor_cu_param_setup"],
            self.out_vars["subgrid_scale_forcing_t_cu_param_setup"],
            self.out_vars["subgrid_scale_forcing_vapor_cu_param_setup"],
            self.out_vars["geopotential_height_forced_cu_param_setup"],
            self.out_vars["seed_convection_cu_param_setup"],
        )
        self.out_vars.update(
            {
                "local_t_excess_cu_param_setup": {},
                "local_vapor_excess_cu_param_setup": {},
                "local_t_new_cu_param_setup": {},
                "local_vapor_forced_cu_param_setup": {},
                "local_t_new_pbl_cu_param_setup": {},
                "local_vapor_forced_pbl_cu_param_setup": {},
                "local_moist_static_energy_cu_param_setup": {},
                "local_maximum_updraft_origin_level_cu_param_setup": {},
                "local_kstabm_cu_param_setup": {},
                "local_ocean_fraction_cu_param_setup": {},
                "local_cap_max_cu_param_setup": {},
                "local_error_code_2_cu_param_setup": {},
                "local_error_code_3_cu_param_setup": {},
                "local_cap_max_increment_cu_param_setup": {},
                "local_geopotential_height_cu_param_setup": {},
                "local_geopotential_height_modified_cu_param_setup": {},
                "local_cloud_work_function_0_cu_param_setup": {},
                "local_cloud_work_function_1_cu_param_setup": {},
                "local_cloud_work_function_2_cu_param_setup": {},
                "local_cloud_work_function_3_cu_param_setup": {},
                "local_cloud_work_function_0_pbl_cu_param_setup": {},
                "local_cloud_work_function_1_pbl_cu_param_setup": {},
                "local_cloud_work_function_1_fa_cu_param_setup": {},
                "local_cin1_cu_param_setup": {},
                "local_k_x_modified_cu_param_setup": {},
                "local_epsilon_cu_param_setup": {},
                "local_pbl_time_scale_cu_param_setup": {},
                "local_t_wetbulb_cu_param_setup": {},
                "local_vapor_wetbulb_cu_param_setup": {},
                "local_tau_ecmwf_cu_param_setup": {},
                "local_f_dicycle_modified_cu_param_setup": {},
                "local_add_buoy_cu_param_setup": {},
                "local_hcdo_cu_param_setup": {},
                "local_cupclw_cu_param_setup": {},
                "local_qrcdo_cu_param_setup": {},
                "local_cloud_moist_static_energy_forced_t_cu_param_setup": {},
                "local_c1d_cu_param_setup": {},
                "local_evap_bcb_cu_param_setup": {},
                "local_random_number_cu_param_setup": {},
                "local_updraft_detrainment_function_cu_param_setup": {},
                "local_epsilon_min_cu_param_setup": {},
                "local_epsilon_max_cu_param_setup": {},
                "local_arbitrary_numerical_parameter_cu_param_setup": {},
            }
        )

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF2020-constants")
        self.cu_param_constants = data_loader.load("GF2020_CumulusParameterization-constants")

    def compute_func(self, **inputs):
        # initalize constants
        config = GF2020Config(SINGLE_COLUMN_MODE=False, **self.constants)
        cumulus_parameterization_config = GF2020CumulusParameterizationConfig(**self.cu_param_constants)
        plume_dependent_constants = GF2020PlumeDependentConstants()

        # initalize saturation tables
        saturation_tables = SaturationVaporPressureTable(self.stencil_factory.backend)

        # initalize dataclasses
        state = GF2020CumulusParameterizationState.zeros(
            self.quantity_factory,
            data_dimensions={
                "plumes": 3,
            },
        )

        locals = GF2020CumulusParameterizationLocals.zeros(
            self.quantity_factory,
            data_dimensions={
                "ensemble_members": MAXENS1 * MAXENS2 * MAXENS3,
            },
        )

        # fill relevant parts of dataclasses
        state.input.t_excess.data[:] = inputs["t_excess_cu_param_setup"]
        state.input.vapor_excess.data[:] = inputs["vapor_excess_cu_param_setup"]
        state.input.ocean_fraction.data[:] = inputs["ocean_fraction_cu_param_setup"]
        state.input_output.t_old.data[:] = inputs["t_old_cu_param_setup"]
        state.input_output.vapor_old.data[:] = inputs["vapor_old_cu_param_setup"]
        state.input.grid_scale_forcing_t.data[:] = inputs["grid_scale_forcing_t_cu_param_setup"]
        state.input.grid_scale_forcing_vapor.data[:] = inputs["grid_scale_forcing_vapor_cu_param_setup"]
        state.input.subgrid_scale_forcing_t.data[:] = inputs["subgrid_scale_forcing_t_cu_param_setup"]
        state.input.subgrid_scale_forcing_vapor.data[:] = inputs["subgrid_scale_forcing_vapor_cu_param_setup"]
        state.input_output.geopotential_height_forced.data[:] = inputs[
            "geopotential_height_forced_cu_param_setup"
        ]
        state.output.epsilon.data[:, :, 2] = inputs["epsilon_cu_param_setup"]  # plume dependent
        state.output.precip.data[:, :, 2] = inputs["precip_cu_param_setup"]  # plume dependent
        state.output.scale_dependence_factor.data[:, :, 2] = inputs[
            "scale_dependence_factor_cu_param_setup"
        ]  # plume dependent
        state.output.lightning_density.data[:] = inputs["lightning_density_cu_param_setup"]
        state.input.seed_convection.data[:] = inputs["seed_convection_cu_param_setup"]
        state.output.error_code.data[:, :, 2] = inputs["error_code_cu_param_setup"]  # plume dependent
        state.input_output.grid_length.data[:] = inputs["grid_length_cu_param_setup"]
        state.input.lateral_entrainment_rate.data[:] = inputs["lateral_entrainment_rate_cu_param_setup"]
        state.output.entrainment_rate.data[:, :, :, 2] = inputs[
            "entrainment_rate_cu_param_setup"
        ]  # plume dependent

        setup = Setup(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        setup(
            state=state,
            locals=locals,
            saturation_tables=saturation_tables,
            plume_dependent_constants=plume_dependent_constants,
            plume="deep",
        )

        # fortran initalized to nan and never touches the top level of some variables
        # set those rows to nan here so the tests pass - even though they are zero in python
        # done like this to avoid introduce nans into the python
        locals.t_new.field[:, :, -1] = np.nan
        locals.vapor_forced.field[:, :, -1] = np.nan
        locals.t_new_pbl.field[:, :, -1] = np.nan
        locals.vapor_forced_pbl.field[:, :, -1] = np.nan
        locals.moist_static_energy.field[:, :, -1] = np.nan
        locals.detrainment_function_updraft.field[:, :, -1] = np.nan

        outputs = {
            # state fields
            "epsilon_cu_param_setup": state.output.epsilon.field[:, :, 2],
            "precip_cu_param_setup": state.output.precip.field[:, :, 2],
            "scale_dependence_factor_cu_param_setup": state.output.scale_dependence_factor.field[:, :, 2],
            "lightning_density_cu_param_setup": state.output.lightning_density.field[:],
            "error_code_cu_param_setup": state.output.error_code.field[:, :, 2],
            "grid_length_cu_param_setup": state.input_output.grid_length.field[:],
            "lateral_entrainment_rate_cu_param_setup": state.input.lateral_entrainment_rate.field[:],
            "entrainment_rate_cu_param_setup": state.output.entrainment_rate.field[:, :, :, 2],
            # local fields
            "local_t_excess_cu_param_setup": locals.t_excess.field[:],
            "local_vapor_excess_cu_param_setup": locals.vapor_excess.field[:],
            "local_t_new_cu_param_setup": locals.t_new.field[:],
            "local_vapor_forced_cu_param_setup": locals.vapor_forced.field[:],
            "local_t_new_pbl_cu_param_setup": locals.t_new_pbl.field[:],
            "local_vapor_forced_pbl_cu_param_setup": locals.vapor_forced_pbl.field[:],
            "local_moist_static_energy_cu_param_setup": locals.moist_static_energy.field[:],
            "local_maximum_updraft_origin_level_cu_param_setup": locals.maximum_updraft_origin_level.field[:],
            "local_kstabm_cu_param_setup": locals.kstabm.field[:] + 1,  # +1 b/c python counts from 0
            "local_ocean_fraction_cu_param_setup": locals.ocean_fraction.field[:],
            "local_cap_max_cu_param_setup": locals.cap_max.field[:],
            "local_error_code_2_cu_param_setup": locals.error_code_2.field[:],
            "local_error_code_3_cu_param_setup": locals.error_code_3.field[:],
            "local_cap_max_increment_cu_param_setup": locals.cap_max_increment.field[:],
            "local_geopotential_height_cu_param_setup": locals.geopotential_height.field[:],
            "local_geopotential_height_modified_cu_param_setup": locals.geopotential_height_modified.field[:],
            "local_cloud_work_function_0_cu_param_setup": locals.cloud_work_function_0.field[:],
            "local_cloud_work_function_1_cu_param_setup": locals.cloud_work_function_1.field[:],
            "local_cloud_work_function_2_cu_param_setup": locals.cloud_work_function_2.field[:],
            "local_cloud_work_function_3_cu_param_setup": locals.cloud_work_function_3.field[:],
            "local_cloud_work_function_0_pbl_cu_param_setup": locals.cloud_work_function_0_pbl.field[:],
            "local_cloud_work_function_1_pbl_cu_param_setup": locals.cloud_work_function_1_pbl.field[:],
            "local_cloud_work_function_1_fa_cu_param_setup": locals.cloud_work_function_1_fa.field[:],
            "local_cin1_cu_param_setup": locals.cin1.field[:],
            "local_k_x_modified_cu_param_setup": locals.k_x_modified.field[:],
            "local_epsilon_cu_param_setup": locals.epsilon.field[:],
            "local_pbl_time_scale_cu_param_setup": locals.pbl_time_scale.field[:],
            "local_t_wetbulb_cu_param_setup": locals.t_wetbulb.field[:],
            "local_vapor_wetbulb_cu_param_setup": locals.vapor_wetbulb.field[:],
            "local_tau_ecmwf_cu_param_setup": locals.tau_ecmwf.field[:],
            "local_f_dicycle_modified_cu_param_setup": locals.f_dicycle_modified.field[:],
            "local_add_buoy_cu_param_setup": locals.add_buoyancy.field[:],
            "local_hcdo_cu_param_setup": locals.hcdo.field[:],
            "local_cupclw_cu_param_setup": locals.cupclw.field[:],
            "local_qrcdo_cu_param_setup": locals.qrcdo.field[:],
            "local_cloud_moist_static_energy_forced_t_cu_param_setup": locals.cloud_moist_static_energy_forced_transported.field[
                :
            ],
            "local_c1d_cu_param_setup": locals.c1d.field[:],
            "local_evap_bcb_cu_param_setup": locals.evap_bcb.field[:],
            "local_mass_flux_ensemble_cu_param_setup": locals.mass_flux_ensemble.field[:],
            "local_precipitation_ensemble_cu_param_setup": locals.precipitation_ensemble.field[:],
            "local_downdraft_scale_dependence_factor_cu_param_setup": locals.scale_dependence_factor_downdraft.field[
                :
            ],
            "local_random_number_cu_param_setup": locals.random_number.field[:],
            "local_updraft_detrainment_function_cu_param_setup": locals.detrainment_function_updraft.field[:],
            "local_epsilon_min_cu_param_setup": locals.epsilon_min.field[:],
            "local_epsilon_max_cu_param_setup": locals.epsilon_max.field[:],
            "local_arbitrary_numerical_parameter_cu_param_setup": locals.arbitrary_numerical_parameter.field[
                :
            ],
        }

        return outputs
