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
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.convection.GF_2020.cumulus_parameterization.setup.setup import Setup
from pyMoist.convection.GF_2020.cumulus_parameterization.constants import MAXENS1, MAXENS2, MAXENS3


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
            "t_excess": {} | {"serialname": "t_excess_cu_param_setup"},
            "vapor_excess": {} | {"serialname": "vapor_excess_cu_param_setup"},
            "ocean_fraction": {} | {"serialname": "ocean_fraction_cu_param_setup"},
            "t_old": {} | {"serialname": "t_old_cu_param_setup"},
            "vapor_old": {} | {"serialname": "vapor_old_cu_param_setup"},
            "grid_scale_forcing_t": {} | {"serialname": "grid_scale_forcing_t_cu_param_setup"},
            "grid_scale_forcing_vapor": {} | {"serialname": "grid_scale_forcing_vapor_cu_param_setup"},
            "subgrid_scale_forcing_t": {} | {"serialname": "subgrid_scale_forcing_t_cu_param_setup"},
            "subgrid_scale_forcing_vapor": {} | {"serialname": "subgrid_scale_forcing_vapor_cu_param_setup"},
            "geopotential_height": {} | {"serialname": "geopotential_height_cu_param_setup"},
            "epsilon": {} | {"serialname": "epsilon_cu_param_setup"},
            "precip": {} | {"serialname": "precip_cu_param_setup"},
            "scale_dependence_factor": {} | {"serialname": "scale_dependence_factor_cu_param_setup"},
            "lightning_density": {} | {"serialname": "lightning_density_cu_param_setup"},
            "seed_convection": {} | {"serialname": "seed_convection_cu_param_setup"},
            "error_code": {} | {"serialname": "error_code_cu_param_setup"},
            "grid_length": {} | {"serialname": "grid_length_cu_param_setup"},
            "lateral_entrainment_rate": {} | {"serialname": "lateral_entrainment_rate_cu_param_setup"},
            "entrainment_rate": {} | {"serialname": "entrainment_rate_cu_param_setup"},
        }

        self.out_vars = self.in_vars["data_vars"].copy()
        del (
            self.out_vars["t_excess"],
            self.out_vars["vapor_excess"],
            self.out_vars["ocean_fraction"],
            self.out_vars["t_old"],
            self.out_vars["vapor_old"],
            self.out_vars["grid_scale_forcing_t"],
            self.out_vars["grid_scale_forcing_vapor"],
            self.out_vars["subgrid_scale_forcing_t"],
            self.out_vars["subgrid_scale_forcing_vapor"],
            self.out_vars["geopotential_height"],
            self.out_vars["seed_convection"],
        )
        self.out_vars.update(
            {
                "local_t_excess": {"serialname": "local_t_excess_cu_param_setup"},
                "local_vapor_excess": {"serialname": "local_vapor_excess_cu_param_setup"},
                "local_t_new": {"serialname": "local_t_new_cu_param_setup"},
                "local_vapor_new": {"serialname": "local_vapor_new_cu_param_setup"},
                "local_t_new_pbl": {"serialname": "local_t_new_pbl_cu_param_setup"},
                "local_vapor_new_pbl": {"serialname": "local_vapor_new_pbl_cu_param_setup"},
                "local_moist_static_energy": {"serialname": "local_moist_static_energy_cu_param_setup"},
                "local_maximum_updraft_origin_level": {
                    "serialname": "local_maximum_updraft_origin_level_cu_param_setup"
                },
                "local_kstabm": {"serialname": "local_kstabm_cu_param_setup"},
                "local_ocean_fraction": {"serialname": "local_ocean_fraction_cu_param_setup"},
                "local_cap_max": {"serialname": "local_cap_max_cu_param_setup"},
                "local_error_code_2": {"serialname": "local_error_code_2_cu_param_setup"},
                "local_error_code_3": {"serialname": "local_error_code_3_cu_param_setup"},
                "local_cap_max_increment": {"serialname": "local_cap_max_increment_cu_param_setup"},
                "local_geopotential_height": {"serialname": "local_geopotential_height_cu_param_setup"},
                "local_geopotential_height_modified": {
                    "serialname": "local_geopotential_height_modified_cu_param_setup"
                },
                "local_cloud_work_function_0": {"serialname": "local_cloud_work_function_0_cu_param_setup"},
                "local_cloud_work_function_1": {"serialname": "local_cloud_work_function_1_cu_param_setup"},
                "local_cloud_work_function_2": {"serialname": "local_cloud_work_function_2_cu_param_setup"},
                "local_cloud_work_function_3": {"serialname": "local_cloud_work_function_3_cu_param_setup"},
                "local_cloud_work_function_0_pbl": {
                    "serialname": "local_cloud_work_function_0_pbl_cu_param_setup"
                },
                "local_cloud_work_function_1_pbl": {
                    "serialname": "local_cloud_work_function_1_pbl_cu_param_setup"
                },
                "local_cloud_work_function_1_fa": {
                    "serialname": "local_cloud_work_function_1_fa_cu_param_setup"
                },
                "local_cin1": {"serialname": "local_cin1_cu_param_setup"},
                "local_k_x_modified": {"serialname": "local_k_x_modified_cu_param_setup"},
                "local_epsilon": {"serialname": "local_epsilon_cu_param_setup"},
                "local_pbl_time_scale": {"serialname": "local_pbl_time_scale_cu_param_setup"},
                "local_t_wetbulb": {"serialname": "local_t_wetbulb_cu_param_setup"},
                "local_vapor_wetbulb": {"serialname": "local_vapor_wetbulb_cu_param_setup"},
                "local_tau_ecmwf": {"serialname": "local_tau_ecmwf_cu_param_setup"},
                "local_f_dicycle_modified": {"serialname": "local_f_dicycle_modified_cu_param_setup"},
                "local_add_buoyancy": {"serialname": "local_add_buoyancy_cu_param_setup"},
                "local_hcdo": {"serialname": "local_hcdo_cu_param_setup"},
                "local_cupclw": {"serialname": "local_cupclw_cu_param_setup"},
                "local_qrcdo": {"serialname": "local_qrcdo_cu_param_setup"},
                "local_hcot": {"serialname": "local_hcot_cu_param_setup"},
                "local_c1d": {"serialname": "local_c1d_cu_param_setup"},
                "local_evap_bcb": {"serialname": "local_evap_bcb_cu_param_setup"},
                "local_scale_dependence_factor": {
                    "serialname": "local_scale_dependence_factor_cu_param_setup"
                },
                "local_random_number": {"serialname": "local_random_number_cu_param_setup"},
                "local_updraft_detrainment_function": {
                    "serialname": "local_updraft_detrainment_function_cu_param_setup"
                },
                "local_epsilon_min": {"serialname": "local_epsilon_min_cu_param_setup"},
                "local_epsilon_max": {"serialname": "local_epsilon_max_cu_param_setup"},
                "local_arbitrary_numerical_parameter": {
                    "serialname": "local_arbitrary_numerical_parameter_cu_param_setup"
                },
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
        state.input.t_excess.data[:] = inputs["t_excess"]
        state.input.vapor_excess.data[:] = inputs["vapor_excess"]
        state.input.ocean_fraction.data[:] = inputs["ocean_fraction"]
        state.input_output.t_old.data[:] = inputs["t_old"]
        state.input_output.vapor_old.data[:] = inputs["vapor_old"]
        state.input.grid_scale_forcing_t.data[:] = inputs["grid_scale_forcing_t"]
        state.input.grid_scale_forcing_vapor.data[:] = inputs["grid_scale_forcing_vapor"]
        state.input.subgrid_scale_forcing_t.data[:] = inputs["subgrid_scale_forcing_t"]
        state.input.subgrid_scale_forcing_vapor.data[:] = inputs["subgrid_scale_forcing_vapor"]
        state.input_output.geopotential_height_forced.data[:] = inputs["geopotential_height"]
        state.output.epsilon.data[:, :, 0] = inputs["epsilon"]  # plume dependent
        state.output.precip.data[:, :, 0] = inputs["precip"]  # plume dependent
        state.output.scale_dependence_factor.data[:, :, 0] = inputs[
            "scale_dependence_factor"
        ]  # plume dependent
        state.output.lightning_density.data[:] = inputs["lightning_density"]
        state.input.seed_convection.data[:] = inputs["seed_convection"]
        state.output.error_code.data[:, :, 0] = inputs["error_code"]  # plume dependent
        state.input_output.grid_length.data[:] = inputs["grid_length"]
        state.input.lateral_entrainment_rate.data[:] = inputs["lateral_entrainment_rate"]
        state.output.entrainment_rate.data[:, :, :, 0] = inputs["entrainment_rate"]  # plume dependent

        setup = Setup(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        if cumulus_parameterization_config.ENABLE_SHALLOW == 1:
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
        import numpy as np

        locals.t_new.field[:, :, -1] = np.nan
        locals.vapor_new.field[:, :, -1] = np.nan
        locals.t_new_pbl.field[:, :, -1] = np.nan
        locals.vapor_new_pbl.field[:, :, -1] = np.nan
        locals.moist_static_energy.field[:, :, -1] = np.nan
        locals.updraft_detrainment_function.field[:, :, -1] = np.nan

        outputs = {
            # state fields
            "epsilon": state.output.epsilon.field[:, :, 0],
            "precip": state.output.precip.field[:, :, 0],
            "scale_dependence_factor": state.output.scale_dependence_factor.field[:, :, 0],
            "lightning_density": state.output.lightning_density.field[:],
            "error_code": state.output.error_code.field[:, :, 0],
            "grid_length": state.input_output.grid_length.field[:],
            "lateral_entrainment_rate": state.input.lateral_entrainment_rate.field[:],
            "entrainment_rate": state.output.entrainment_rate.field[:, :, :, 0],
            # local fields
            "local_t_excess": locals.t_excess.field[:],
            "local_vapor_excess": locals.vapor_excess.field[:],
            "local_t_new": locals.t_new.field[:],
            "local_vapor_new": locals.vapor_new.field[:],
            "local_t_new_pbl": locals.t_new_pbl.field[:],
            "local_vapor_new_pbl": locals.vapor_new_pbl.field[:],
            "local_moist_static_energy": locals.moist_static_energy.field[:],
            "local_maximum_updraft_origin_level": locals.maximum_updraft_origin_level.field[:],
            "local_kstabm": locals.kstabm.field[:] + 1,  # +1 b/c python counts from 0
            "local_ocean_fraction": locals.ocean_fraction.field[:],
            "local_cap_max": locals.cap_max.field[:],
            "local_error_code_2": locals.error_code_2.field[:],
            "local_error_code_3": locals.error_code_3.field[:],
            "local_error_code_string": locals.error_code_string.field[:],
            "local_cap_max_increment": locals.cap_max_increment.field[:],
            "local_geopotential_height": locals.geopotential_height.field[:],
            "local_geopotential_height_modified": locals.geopotential_height_modified.field[:],
            "local_cloud_work_function_0": locals.cloud_work_function_0.field[:],
            "local_cloud_work_function_1": locals.cloud_work_function_1.field[:],
            "local_cloud_work_function_2": locals.cloud_work_function_2.field[:],
            "local_cloud_work_function_3": locals.cloud_work_function_3.field[:],
            "local_cloud_work_function_0_pbl": locals.cloud_work_function_0_pbl.field[:],
            "local_cloud_work_function_1_pbl": locals.cloud_work_function_1_pbl.field[:],
            "local_cloud_work_function_1_fa": locals.cloud_work_function_1_fa.field[:],
            "local_cin1": locals.cin1.field[:],
            "local_k_x_modified": locals.k_x_modified.field[:],
            "local_epsilon": locals.epsilon.field[:],
            "local_pbl_time_scale": locals.pbl_time_scale.field[:],
            "local_t_wetbulb": locals.t_wetbulb.field[:],
            "local_vapor_wetbulb": locals.vapor_wetbulb.field[:],
            "local_tau_ecmwf": locals.tau_ecmwf.field[:],
            "local_f_dicycle_modified": locals.f_dicycle_modified.field[:],
            "local_add_buoyancy": locals.add_buoyancy.field[:],
            "local_hcdo": locals.hcdo.field[:],
            "local_cupclw": locals.cupclw.field[:],
            "local_qrcdo": locals.qrcdo.field[:],
            "local_hcot": locals.hcot.field[:],
            "local_c1d": locals.c1d.field[:],
            "local_evap_bcb": locals.evap_bcb.field[:],
            "local_mass_flux_ensemble": locals.mass_flux_ensemble.field[:],
            "local_precipitation_ensemble": locals.precipitation_ensemble.field[:],
            "local_scale_dependence_factor": locals.scale_dependence_factor.field[:],
            "local_random_number": locals.random_number.field[:],
            "local_updraft_detrainment_function": locals.updraft_detrainment_function.field[:],
            "local_epsilon_min": locals.epsilon_min.field[:],
            "local_epsilon_max": locals.epsilon_max.field[:],
            "local_arbitrary_numerical_parameter": locals.arbitrary_numerical_parameter.field[:],
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
            "t_excess": {} | {"serialname": "t_excess_cu_param_setup"},
            "vapor_excess": {} | {"serialname": "vapor_excess_cu_param_setup"},
            "ocean_fraction": {} | {"serialname": "ocean_fraction_cu_param_setup"},
            "t_old": {} | {"serialname": "t_old_cu_param_setup"},
            "vapor_old": {} | {"serialname": "vapor_old_cu_param_setup"},
            "grid_scale_forcing_t": {} | {"serialname": "grid_scale_forcing_t_cu_param_setup"},
            "grid_scale_forcing_vapor": {} | {"serialname": "grid_scale_forcing_vapor_cu_param_setup"},
            "subgrid_scale_forcing_t": {} | {"serialname": "subgrid_scale_forcing_t_cu_param_setup"},
            "subgrid_scale_forcing_vapor": {} | {"serialname": "subgrid_scale_forcing_vapor_cu_param_setup"},
            "geopotential_height": {} | {"serialname": "geopotential_height_cu_param_setup"},
            "epsilon": {} | {"serialname": "epsilon_cu_param_setup"},
            "precip": {} | {"serialname": "precip_cu_param_setup"},
            "scale_dependence_factor": {} | {"serialname": "scale_dependence_factor_cu_param_setup"},
            "lightning_density": {} | {"serialname": "lightning_density_cu_param_setup"},
            "seed_convection": {} | {"serialname": "seed_convection_cu_param_setup"},
            "error_code": {} | {"serialname": "error_code_cu_param_setup"},
            "grid_length": {} | {"serialname": "grid_length_cu_param_setup"},
            "lateral_entrainment_rate": {} | {"serialname": "lateral_entrainment_rate_cu_param_setup"},
            "entrainment_rate": {} | {"serialname": "entrainment_rate_cu_param_setup"},
        }

        self.out_vars = self.in_vars["data_vars"].copy()
        del (
            self.out_vars["t_excess"],
            self.out_vars["vapor_excess"],
            self.out_vars["ocean_fraction"],
            self.out_vars["t_old"],
            self.out_vars["vapor_old"],
            self.out_vars["grid_scale_forcing_t"],
            self.out_vars["grid_scale_forcing_vapor"],
            self.out_vars["subgrid_scale_forcing_t"],
            self.out_vars["subgrid_scale_forcing_vapor"],
            self.out_vars["geopotential_height"],
            self.out_vars["seed_convection"],
        )
        self.out_vars.update(
            {
                "local_t_excess": {"serialname": "local_t_excess_cu_param_setup"},
                "local_vapor_excess": {"serialname": "local_vapor_excess_cu_param_setup"},
                "local_t_new": {"serialname": "local_t_new_cu_param_setup"},
                "local_vapor_new": {"serialname": "local_vapor_new_cu_param_setup"},
                "local_t_new_pbl": {"serialname": "local_t_new_pbl_cu_param_setup"},
                "local_vapor_new_pbl": {"serialname": "local_vapor_new_pbl_cu_param_setup"},
                "local_moist_static_energy": {"serialname": "local_moist_static_energy_cu_param_setup"},
                "local_maximum_updraft_origin_level": {
                    "serialname": "local_maximum_updraft_origin_level_cu_param_setup"
                },
                "local_kstabm": {"serialname": "local_kstabm_cu_param_setup"},
                "local_ocean_fraction": {"serialname": "local_ocean_fraction_cu_param_setup"},
                "local_cap_max": {"serialname": "local_cap_max_cu_param_setup"},
                "local_error_code_2": {"serialname": "local_error_code_2_cu_param_setup"},
                "local_error_code_3": {"serialname": "local_error_code_3_cu_param_setup"},
                "local_cap_max_increment": {"serialname": "local_cap_max_increment_cu_param_setup"},
                "local_geopotential_height": {"serialname": "local_geopotential_height_cu_param_setup"},
                "local_geopotential_height_modified": {
                    "serialname": "local_geopotential_height_modified_cu_param_setup"
                },
                "local_cloud_work_function_0": {"serialname": "local_cloud_work_function_0_cu_param_setup"},
                "local_cloud_work_function_1": {"serialname": "local_cloud_work_function_1_cu_param_setup"},
                "local_cloud_work_function_2": {"serialname": "local_cloud_work_function_2_cu_param_setup"},
                "local_cloud_work_function_3": {"serialname": "local_cloud_work_function_3_cu_param_setup"},
                "local_cloud_work_function_0_pbl": {
                    "serialname": "local_cloud_work_function_0_pbl_cu_param_setup"
                },
                "local_cloud_work_function_1_pbl": {
                    "serialname": "local_cloud_work_function_1_pbl_cu_param_setup"
                },
                "local_cloud_work_function_1_fa": {
                    "serialname": "local_cloud_work_function_1_fa_cu_param_setup"
                },
                "local_cin1": {"serialname": "local_cin1_cu_param_setup"},
                "local_k_x_modified": {"serialname": "local_k_x_modified_cu_param_setup"},
                "local_epsilon": {"serialname": "local_epsilon_cu_param_setup"},
                "local_pbl_time_scale": {"serialname": "local_pbl_time_scale_cu_param_setup"},
                "local_t_wetbulb": {"serialname": "local_t_wetbulb_cu_param_setup"},
                "local_vapor_wetbulb": {"serialname": "local_vapor_wetbulb_cu_param_setup"},
                "local_tau_ecmwf": {"serialname": "local_tau_ecmwf_cu_param_setup"},
                "local_f_dicycle_modified": {"serialname": "local_f_dicycle_modified_cu_param_setup"},
                "local_add_buoyancy": {"serialname": "local_add_buoy_cu_param_setup"},
                "local_hcdo": {"serialname": "local_hcdo_cu_param_setup"},
                "local_cupclw": {"serialname": "local_cupclw_cu_param_setup"},
                "local_qrcdo": {"serialname": "local_qrcdo_cu_param_setup"},
                "local_hcot": {"serialname": "local_hcot_cu_param_setup"},
                "local_c1d": {"serialname": "local_c1d_cu_param_setup"},
                "local_evap_bcb": {"serialname": "local_evap_bcb_cu_param_setup"},
                "local_scale_dependence_factor": {
                    "serialname": "local_scale_dependence_factor_cu_param_setup"
                },
                "local_random_number": {"serialname": "local_random_number_cu_param_setup"},
                "local_updraft_detrainment_function": {
                    "serialname": "local_updraft_detrainment_function_cu_param_setup"
                },
                "local_epsilon_min": {"serialname": "local_epsilon_min_cu_param_setup"},
                "local_epsilon_max": {"serialname": "local_epsilon_max_cu_param_setup"},
                "local_arbitrary_numerical_parameter": {
                    "serialname": "local_arbitrary_numerical_parameter_cu_param_setup"
                },
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
        state.input.t_excess.data[:] = inputs["t_excess"]
        state.input.vapor_excess.data[:] = inputs["vapor_excess"]
        state.input.ocean_fraction.data[:] = inputs["ocean_fraction"]
        state.input_output.t_old.data[:] = inputs["t_old"]
        state.input_output.vapor_old.data[:] = inputs["vapor_old"]
        state.input.grid_scale_forcing_t.data[:] = inputs["grid_scale_forcing_t"]
        state.input.grid_scale_forcing_vapor.data[:] = inputs["grid_scale_forcing_vapor"]
        state.input.subgrid_scale_forcing_t.data[:] = inputs["subgrid_scale_forcing_t"]
        state.input.subgrid_scale_forcing_vapor.data[:] = inputs["subgrid_scale_forcing_vapor"]
        state.input_output.geopotential_height_forced.data[:] = inputs["geopotential_height"]
        state.output.epsilon.data[:, :, 1] = inputs["epsilon"]  # plume dependent
        state.output.precip.data[:, :, 1] = inputs["precip"]  # plume dependent
        state.output.scale_dependence_factor.data[:, :, 1] = inputs[
            "scale_dependence_factor"
        ]  # plume dependent
        state.output.lightning_density.data[:] = inputs["lightning_density"]
        state.input.seed_convection.data[:] = inputs["seed_convection"]
        state.output.error_code.data[:, :, 1] = inputs["error_code"]  # plume dependent
        state.input_output.grid_length.data[:] = inputs["grid_length"]
        state.input.lateral_entrainment_rate.data[:] = inputs["lateral_entrainment_rate"]
        state.output.entrainment_rate.data[:, :, :, 1] = inputs["entrainment_rate"]  # plume dependent

        setup = Setup(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        if cumulus_parameterization_config.ENABLE_MID == 1:
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
        import numpy as np

        locals.t_new.field[:, :, -1] = np.nan
        locals.vapor_new.field[:, :, -1] = np.nan
        locals.t_new_pbl.field[:, :, -1] = np.nan
        locals.vapor_new_pbl.field[:, :, -1] = np.nan
        locals.moist_static_energy.field[:, :, -1] = np.nan
        locals.updraft_detrainment_function.field[:, :, -1] = np.nan

        outputs = {
            # state fields
            "epsilon": state.output.epsilon.field[:, :, 1],
            "precip": state.output.precip.field[:, :, 1],
            "scale_dependence_factor": state.output.scale_dependence_factor.field[:, :, 1],
            "lightning_density": state.output.lightning_density.field[:],
            "error_code": state.output.error_code.field[:, :, 1],
            "grid_length": state.input_output.grid_length.field[:],
            "lateral_entrainment_rate": state.input.lateral_entrainment_rate.field[:],
            "entrainment_rate": state.output.entrainment_rate.field[:, :, :, 1],
            # local fields
            "local_t_excess": locals.t_excess.field[:],
            "local_vapor_excess": locals.vapor_excess.field[:],
            "local_t_new": locals.t_new.field[:],
            "local_vapor_new": locals.vapor_new.field[:],
            "local_t_new_pbl": locals.t_new_pbl.field[:],
            "local_vapor_new_pbl": locals.vapor_new_pbl.field[:],
            "local_moist_static_energy": locals.moist_static_energy.field[:],
            "local_maximum_updraft_origin_level": locals.maximum_updraft_origin_level.field[:],
            "local_kstabm": locals.kstabm.field[:] + 1,  # +1 b/c python counts from 0
            "local_ocean_fraction": locals.ocean_fraction.field[:],
            "local_cap_max": locals.cap_max.field[:],
            "local_error_code_2": locals.error_code_2.field[:],
            "local_error_code_3": locals.error_code_3.field[:],
            "local_error_code_string": locals.error_code_string.field[:],
            "local_cap_max_increment": locals.cap_max_increment.field[:],
            "local_geopotential_height": locals.geopotential_height.field[:],
            "local_geopotential_height_modified": locals.geopotential_height_modified.field[:],
            "local_cloud_work_function_0": locals.cloud_work_function_0.field[:],
            "local_cloud_work_function_1": locals.cloud_work_function_1.field[:],
            "local_cloud_work_function_2": locals.cloud_work_function_2.field[:],
            "local_cloud_work_function_3": locals.cloud_work_function_3.field[:],
            "local_cloud_work_function_0_pbl": locals.cloud_work_function_0_pbl.field[:],
            "local_cloud_work_function_1_pbl": locals.cloud_work_function_1_pbl.field[:],
            "local_cloud_work_function_1_fa": locals.cloud_work_function_1_fa.field[:],
            "local_cin1": locals.cin1.field[:],
            "local_k_x_modified": locals.k_x_modified.field[:],
            "local_epsilon": locals.epsilon.field[:],
            "local_pbl_time_scale": locals.pbl_time_scale.field[:],
            "local_t_wetbulb": locals.t_wetbulb.field[:],
            "local_vapor_wetbulb": locals.vapor_wetbulb.field[:],
            "local_tau_ecmwf": locals.tau_ecmwf.field[:],
            "local_f_dicycle_modified": locals.f_dicycle_modified.field[:],
            "local_add_buoyancy": locals.add_buoyancy_modified.field[:],
            "local_hcdo": locals.hcdo.field[:],
            "local_cupclw": locals.cupclw.field[:],
            "local_qrcdo": locals.qrcdo.field[:],
            "local_hcot": locals.hcot.field[:],
            "local_c1d": locals.c1d.field[:],
            "local_evap_bcb": locals.evap_bcb.field[:],
            "local_mass_flux_ensemble": locals.mass_flux_ensemble.field[:],
            "local_precipitation_ensemble": locals.precipitation_ensemble.field[:],
            "local_scale_dependence_factor": locals.scale_dependence_factor.field[:],
            "local_random_number": locals.random_number.field[:],
            "local_updraft_detrainment_function": locals.updraft_detrainment_function.field[:],
            "local_epsilon_min": locals.epsilon_min.field[:],
            "local_epsilon_max": locals.epsilon_max.field[:],
            "local_arbitrary_numerical_parameter": locals.arbitrary_numerical_parameter.field[:],
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
            "t_excess": {} | {"serialname": "t_excess_cu_param_setup"},
            "vapor_excess": {} | {"serialname": "vapor_excess_cu_param_setup"},
            "ocean_fraction": {} | {"serialname": "ocean_fraction_cu_param_setup"},
            "t_old": {} | {"serialname": "t_old_cu_param_setup"},
            "vapor_old": {} | {"serialname": "vapor_old_cu_param_setup"},
            "grid_scale_forcing_t": {} | {"serialname": "grid_scale_forcing_t_cu_param_setup"},
            "grid_scale_forcing_vapor": {} | {"serialname": "grid_scale_forcing_vapor_cu_param_setup"},
            "subgrid_scale_forcing_t": {} | {"serialname": "subgrid_scale_forcing_t_cu_param_setup"},
            "subgrid_scale_forcing_vapor": {} | {"serialname": "subgrid_scale_forcing_vapor_cu_param_setup"},
            "geopotential_height": {} | {"serialname": "geopotential_height_cu_param_setup"},
            "epsilon": {} | {"serialname": "epsilon_cu_param_setup"},
            "precip": {} | {"serialname": "precip_cu_param_setup"},
            "scale_dependence_factor": {} | {"serialname": "scale_dependence_factor_cu_param_setup"},
            "lightning_density": {} | {"serialname": "lightning_density_cu_param_setup"},
            "seed_convection": {} | {"serialname": "seed_convection_cu_param_setup"},
            "error_code": {} | {"serialname": "error_code_cu_param_setup"},
            "grid_length": {} | {"serialname": "grid_length_cu_param_setup"},
            "lateral_entrainment_rate": {} | {"serialname": "lateral_entrainment_rate_cu_param_setup"},
            "entrainment_rate": {} | {"serialname": "entrainment_rate_cu_param_setup"},
        }

        self.out_vars = self.in_vars["data_vars"].copy()
        del (
            self.out_vars["t_excess"],
            self.out_vars["vapor_excess"],
            self.out_vars["ocean_fraction"],
            self.out_vars["t_old"],
            self.out_vars["vapor_old"],
            self.out_vars["grid_scale_forcing_t"],
            self.out_vars["grid_scale_forcing_vapor"],
            self.out_vars["subgrid_scale_forcing_t"],
            self.out_vars["subgrid_scale_forcing_vapor"],
            self.out_vars["geopotential_height"],
            self.out_vars["seed_convection"],
        )
        self.out_vars.update(
            {
                "local_t_excess": {"serialname": "local_t_excess_cu_param_setup"},
                "local_vapor_excess": {"serialname": "local_vapor_excess_cu_param_setup"},
                "local_t_new": {"serialname": "local_t_new_cu_param_setup"},
                "local_vapor_new": {"serialname": "local_vapor_new_cu_param_setup"},
                "local_t_new_pbl": {"serialname": "local_t_new_pbl_cu_param_setup"},
                "local_vapor_new_pbl": {"serialname": "local_vapor_new_pbl_cu_param_setup"},
                "local_moist_static_energy": {"serialname": "local_moist_static_energy_cu_param_setup"},
                "local_maximum_updraft_origin_level": {
                    "serialname": "local_maximum_updraft_origin_level_cu_param_setup"
                },
                "local_kstabm": {"serialname": "local_kstabm_cu_param_setup"},
                "local_ocean_fraction": {"serialname": "local_ocean_fraction_cu_param_setup"},
                "local_cap_max": {"serialname": "local_cap_max_cu_param_setup"},
                "local_error_code_2": {"serialname": "local_error_code_2_cu_param_setup"},
                "local_error_code_3": {"serialname": "local_error_code_3_cu_param_setup"},
                "local_cap_max_increment": {"serialname": "local_cap_max_increment_cu_param_setup"},
                "local_geopotential_height": {"serialname": "local_geopotential_height_cu_param_setup"},
                "local_geopotential_height_modified": {
                    "serialname": "local_geopotential_height_modified_cu_param_setup"
                },
                "local_cloud_work_function_0": {"serialname": "local_cloud_work_function_0_cu_param_setup"},
                "local_cloud_work_function_1": {"serialname": "local_cloud_work_function_1_cu_param_setup"},
                "local_cloud_work_function_2": {"serialname": "local_cloud_work_function_2_cu_param_setup"},
                "local_cloud_work_function_3": {"serialname": "local_cloud_work_function_3_cu_param_setup"},
                "local_cloud_work_function_0_pbl": {
                    "serialname": "local_cloud_work_function_0_pbl_cu_param_setup"
                },
                "local_cloud_work_function_1_pbl": {
                    "serialname": "local_cloud_work_function_1_pbl_cu_param_setup"
                },
                "local_cloud_work_function_1_fa": {
                    "serialname": "local_cloud_work_function_1_fa_cu_param_setup"
                },
                "local_cin1": {"serialname": "local_cin1_cu_param_setup"},
                "local_k_x_modified": {"serialname": "local_k_x_modified_cu_param_setup"},
                "local_epsilon": {"serialname": "local_epsilon_cu_param_setup"},
                "local_pbl_time_scale": {"serialname": "local_pbl_time_scale_cu_param_setup"},
                "local_t_wetbulb": {"serialname": "local_t_wetbulb_cu_param_setup"},
                "local_vapor_wetbulb": {"serialname": "local_vapor_wetbulb_cu_param_setup"},
                "local_tau_ecmwf": {"serialname": "local_tau_ecmwf_cu_param_setup"},
                "local_f_dicycle_modified": {"serialname": "local_f_dicycle_modified_cu_param_setup"},
                "local_add_buoyancy": {"serialname": "local_add_buoy_cu_param_setup"},
                "local_hcdo": {"serialname": "local_hcdo_cu_param_setup"},
                "local_cupclw": {"serialname": "local_cupclw_cu_param_setup"},
                "local_qrcdo": {"serialname": "local_qrcdo_cu_param_setup"},
                "local_hcot": {"serialname": "local_hcot_cu_param_setup"},
                "local_c1d": {"serialname": "local_c1d_cu_param_setup"},
                "local_evap_bcb": {"serialname": "local_evap_bcb_cu_param_setup"},
                "local_scale_dependence_factor": {
                    "serialname": "local_scale_dependence_factor_cu_param_setup"
                },
                "local_random_number": {"serialname": "local_random_number_cu_param_setup"},
                "local_updraft_detrainment_function": {
                    "serialname": "local_updraft_detrainment_function_cu_param_setup"
                },
                "local_epsilon_min": {"serialname": "local_epsilon_min_cu_param_setup"},
                "local_epsilon_max": {"serialname": "local_epsilon_max_cu_param_setup"},
                "local_arbitrary_numerical_parameter": {
                    "serialname": "local_arbitrary_numerical_parameter_cu_param_setup"
                },
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
        state.input.t_excess.data[:] = inputs["t_excess"]
        state.input.vapor_excess.data[:] = inputs["vapor_excess"]
        state.input.ocean_fraction.data[:] = inputs["ocean_fraction"]
        state.input_output.t_old.data[:] = inputs["t_old"]
        state.input_output.vapor_old.data[:] = inputs["vapor_old"]
        state.input.grid_scale_forcing_t.data[:] = inputs["grid_scale_forcing_t"]
        state.input.grid_scale_forcing_vapor.data[:] = inputs["grid_scale_forcing_vapor"]
        state.input.subgrid_scale_forcing_t.data[:] = inputs["subgrid_scale_forcing_t"]
        state.input.subgrid_scale_forcing_vapor.data[:] = inputs["subgrid_scale_forcing_vapor"]
        state.input_output.geopotential_height_forced.data[:] = inputs["geopotential_height"]
        state.output.epsilon.data[:, :, 2] = inputs["epsilon"]  # plume dependent
        state.output.precip.data[:, :, 2] = inputs["precip"]  # plume dependent
        state.output.scale_dependence_factor.data[:, :, 2] = inputs[
            "scale_dependence_factor"
        ]  # plume dependent
        state.output.lightning_density.data[:] = inputs["lightning_density"]
        state.input.seed_convection.data[:] = inputs["seed_convection"]
        state.output.error_code.data[:, :, 2] = inputs["error_code"]  # plume dependent
        state.input_output.grid_length.data[:] = inputs["grid_length"]
        state.input.lateral_entrainment_rate.data[:] = inputs["lateral_entrainment_rate"]
        state.output.entrainment_rate.data[:, :, :, 2] = inputs["entrainment_rate"]  # plume dependent

        setup = Setup(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        if cumulus_parameterization_config.ENABLE_DEEP == 1:
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
        import numpy as np

        locals.t_new.field[:, :, -1] = np.nan
        locals.vapor_new.field[:, :, -1] = np.nan
        locals.t_new_pbl.field[:, :, -1] = np.nan
        locals.vapor_new_pbl.field[:, :, -1] = np.nan
        locals.moist_static_energy.field[:, :, -1] = np.nan
        locals.updraft_detrainment_function.field[:, :, -1] = np.nan

        outputs = {
            # state fields
            "epsilon": state.output.epsilon.field[:, :, 2],
            "precip": state.output.precip.field[:, :, 2],
            "scale_dependence_factor": state.output.scale_dependence_factor.field[:, :, 2],
            "lightning_density": state.output.lightning_density.field[:],
            "error_code": state.output.error_code.field[:, :, 2],
            "grid_length": state.input_output.grid_length.field[:],
            "lateral_entrainment_rate": state.input.lateral_entrainment_rate.field[:],
            "entrainment_rate": state.output.entrainment_rate.field[:, :, :, 2],
            # local fields
            "local_t_excess": locals.t_excess.field[:],
            "local_vapor_excess": locals.vapor_excess.field[:],
            "local_t_new": locals.t_new.field[:],
            "local_vapor_new": locals.vapor_new.field[:],
            "local_t_new_pbl": locals.t_new_pbl.field[:],
            "local_vapor_new_pbl": locals.vapor_new_pbl.field[:],
            "local_moist_static_energy": locals.moist_static_energy.field[:],
            "local_maximum_updraft_origin_level": locals.maximum_updraft_origin_level.field[:],
            "local_kstabm": locals.kstabm.field[:] + 1,  # +1 b/c python counts from 0
            "local_ocean_fraction": locals.ocean_fraction.field[:],
            "local_cap_max": locals.cap_max.field[:],
            "local_error_code_2": locals.error_code_2.field[:],
            "local_error_code_3": locals.error_code_3.field[:],
            "local_error_code_string": locals.error_code_string.field[:],
            "local_cap_max_increment": locals.cap_max_increment.field[:],
            "local_geopotential_height": locals.geopotential_height.field[:],
            "local_geopotential_height_modified": locals.geopotential_height_modified.field[:],
            "local_cloud_work_function_0": locals.cloud_work_function_0.field[:],
            "local_cloud_work_function_1": locals.cloud_work_function_1.field[:],
            "local_cloud_work_function_2": locals.cloud_work_function_2.field[:],
            "local_cloud_work_function_3": locals.cloud_work_function_3.field[:],
            "local_cloud_work_function_0_pbl": locals.cloud_work_function_0_pbl.field[:],
            "local_cloud_work_function_1_pbl": locals.cloud_work_function_1_pbl.field[:],
            "local_cloud_work_function_1_fa": locals.cloud_work_function_1_fa.field[:],
            "local_cin1": locals.cin1.field[:],
            "local_k_x_modified": locals.k_x_modified.field[:],
            "local_epsilon": locals.epsilon.field[:],
            "local_pbl_time_scale": locals.pbl_time_scale.field[:],
            "local_t_wetbulb": locals.t_wetbulb.field[:],
            "local_vapor_wetbulb": locals.vapor_wetbulb.field[:],
            "local_tau_ecmwf": locals.tau_ecmwf.field[:],
            "local_f_dicycle_modified": locals.f_dicycle_modified.field[:],
            "local_add_buoyancy": locals.add_buoyancy.field[:],
            "local_hcdo": locals.hcdo.field[:],
            "local_cupclw": locals.cupclw.field[:],
            "local_qrcdo": locals.qrcdo.field[:],
            "local_hcot": locals.hcot.field[:],
            "local_c1d": locals.c1d.field[:],
            "local_evap_bcb": locals.evap_bcb.field[:],
            "local_mass_flux_ensemble": locals.mass_flux_ensemble.field[:],
            "local_precipitation_ensemble": locals.precipitation_ensemble.field[:],
            "local_scale_dependence_factor": locals.scale_dependence_factor.field[:],
            "local_random_number": locals.random_number.field[:],
            "local_updraft_detrainment_function": locals.updraft_detrainment_function.field[:],
            "local_epsilon_min": locals.epsilon_min.field[:],
            "local_epsilon_max": locals.epsilon_max.field[:],
            "local_arbitrary_numerical_parameter": locals.arbitrary_numerical_parameter.field[:],
        }

        return outputs
