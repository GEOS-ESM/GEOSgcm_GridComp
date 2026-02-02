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
from pyMoist.convection.GF_2020.cumulus_parameterization.setup.setup import Setup
import numpy as np


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
            "t_excess": {},
            "vapor_excess": {},
            "ocean_fraction": {},
            "t_old": {},
            "vapor_old": {},
            "grid_scale_forcing_t": {},
            "grid_scale_forcing_vapor": {},
            "subgrid_scale_forcing_t": {},
            "subgrid_scale_forcing_vapor": {},
            "geopotential_height_forced": {},
            "epsilon_forced": {},
            "precip": {},
            "scale_dependence_factor": {},
            "lightning_density": {},
            "seed_convection": {},
            "error_code": {},
            "grid_length": {},
            "lateral_entrainment_rate": {},
            "entrainment_rate": {},
        }

        out_vars.update(in_vars["data_vars"])
        del (
            out_vars["t_excess"],
            out_vars["vapor_excess"],
            out_vars["ocean_fraction"],
            out_vars["t_old"],
            out_vars["vapor_old"],
            out_vars["grid_scale_forcing_t"],
            out_vars["grid_scale_forcing_vapor"],
            out_vars["subgrid_scale_forcing_t"],
            out_vars["subgrid_scale_forcing_vapor"],
            out_vars["seed_convection"],
        )
        out_vars.update(
            {
                "local_t_excess": {},
                "local_vapor_excess": {},
                "local_t_new": {},
                "local_vapor_forced": {},
                "local_t_new_pbl": {},
                "local_vapor_forced_pbl": {},
                "local_dmoist_static_energydt": {},
                "local_maximum_updraft_origin_level": {},
                "local_kstabm": {},
                "local_ocean_fraction": {},
                "local_error_code_2": {},
                "local_error_code_3": {},
                "local_cap_max": {},
                "local_cap_max_increment": {},
                "local_geopotential_height": {},
                "local_geopotential_height_modified": {},
                "local_cloud_workfunction_0": {},
                "local_cloud_workfunction_0_pbl": {},
                "local_cloud_workfunction_1": {},
                "local_cloud_workfunction_1_pbl": {},
                "local_cloud_workfunction_1_fa": {},
                "local_cloud_workfunction_2": {},
                "local_cloud_workfunction_3": {},
                "local_cin_1": {},
                "local_k_x_modified": {},
                "local_epsilon": {},
                "local_epsilon_min": {},
                "local_epsilon_max": {},
                "local_pbl_time_scale": {},
                "local_t_wetbulb": {},
                "local_vapor_wetbulb": {},
                "local_cape_removal_time_scale": {},
                "local_f_dicycle_modified": {},
                "local_add_buoyancy": {},
                "local_cloud_moist_static_energy_downdraft_forced": {},
                "local_downdraft_saturation_vapor_forced": {},
                "local_cloud_moist_static_energy_forced_transported": {},
                "local_c1d": {},
                "local_evaporation_below_cloud_base": {},
                "local_mass_flux_ensemble": {},
                "local_precipitation_ensemble": {},
                "local_scale_dependence_factor_downdraft": {},
                "local_random_number": {},
                "local_detrainment_function_updraft": {},
                "local_arbitrary_numerical_parameter": {},
            }
        )

    def __call__(self, constants: dict, cu_param_constants: dict, plume: str, **inputs):
        # initialize constants
        config = GF2020Config(SINGLE_COLUMN_MODE=False, **constants)
        cumulus_parameterization_config = GF2020CumulusParameterizationConfig(**cu_param_constants)
        plume_dependent_constants = GF2020PlumeDependentConstants()

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
        state.input.t_excess.data[:] = inputs["t_excess"]
        state.input.vapor_excess.data[:] = inputs["vapor_excess"]
        state.input.ocean_fraction.data[:] = inputs["ocean_fraction"]
        state.input_output.t_old.data[:] = inputs["t_old"]
        state.input_output.vapor_old.data[:] = inputs["vapor_old"]
        state.input.grid_scale_forcing_t.data[:] = inputs["grid_scale_forcing_t"]
        state.input.grid_scale_forcing_vapor.data[:] = inputs["grid_scale_forcing_vapor"]
        state.input.subgrid_scale_forcing_t.data[:] = inputs["subgrid_scale_forcing_t"]
        state.input.subgrid_scale_forcing_vapor.data[:] = inputs["subgrid_scale_forcing_vapor"]
        state.input_output.geopotential_height_forced.data[:] = inputs["geopotential_height_forced"]
        state.output.epsilon_forced.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "epsilon_forced"
        ]
        state.output.precip.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["precip"]
        state.output.scale_dependence_factor.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "scale_dependence_factor"
        ]
        state.output.lightning_density.data[:] = inputs["lightning_density"]
        state.input.seed_convection.data[:] = inputs["seed_convection"]
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["error_code"]
        state.input_output.grid_length.data[:] = inputs["grid_length"]
        state.input.lateral_entrainment_rate.data[:] = inputs["lateral_entrainment_rate"]
        state.output.entrainment_rate.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "entrainment_rate"
        ]

        # prefill locals with nans to replicate fortran initalization
        locals.t_excess.field[:] = np.nan
        locals.vapor_excess.field[:] = np.nan
        locals.t_new.field[:] = np.nan
        locals.vapor_forced.field[:] = np.nan
        locals.t_new_pbl.field[:] = np.nan
        locals.vapor_forced_pbl.field[:] = np.nan
        locals.dmoist_static_energydt.field[:] = np.nan
        # locals.maximum_updraft_origin_level.field[:] = np.nan
        # locals.kstabm.field[:] = np.nan
        locals.ocean_fraction.field[:] = np.nan
        # locals.error_code_2.field[:] = np.nan
        # locals.error_code_3.field[:] = np.nan
        locals.cap_max.field[:] = np.nan
        locals.cap_max_increment.field[:] = np.nan
        locals.geopotential_height.field[:] = np.nan
        locals.geopotential_height_modified.field[:] = np.nan
        locals.cloud_workfunction_0.field[:] = np.nan
        locals.cloud_workfunction_0_pbl.field[:] = np.nan
        locals.cloud_workfunction_1.field[:] = np.nan
        locals.cloud_workfunction_1_pbl.field[:] = np.nan
        locals.cloud_workfunction_1_fa.field[:] = np.nan
        locals.cloud_workfunction_2.field[:] = np.nan
        locals.cloud_workfunction_3.field[:] = np.nan
        locals.cin_1.field[:] = np.nan
        locals.k_x_modified.field[:] = np.nan
        locals.epsilon.field[:] = np.nan
        locals.epsilon_min.field[:] = np.nan
        locals.epsilon_max.field[:] = np.nan
        locals.pbl_time_scale.field[:] = np.nan
        locals.t_wetbulb.field[:] = np.nan
        locals.vapor_wetbulb.field[:] = np.nan
        locals.cape_removal_time_scale.field[:] = np.nan
        locals.f_dicycle_modified.field[:] = np.nan
        locals.add_buoyancy.field[:] = np.nan
        locals.cloud_moist_static_energy_downdraft_forced.field[:] = np.nan
        locals.downdraft_saturation_vapor_forced.field[:] = np.nan
        locals.cloud_moist_static_energy_forced_transported.field[:] = np.nan
        locals.c1d.field[:] = np.nan
        locals.evaporation_below_cloud_base.field[:] = np.nan
        locals.mass_flux_ensemble.field[:] = np.nan
        locals.precipitation_ensemble.field[:] = np.nan
        locals.scale_dependence_factor_downdraft.field[:] = np.nan
        locals.random_number.field[:] = np.nan
        locals.detrainment_function_updraft.field[:] = np.nan
        locals.arbitrary_numerical_parameter.field[:] = np.nan

        # initialize test code
        code = Setup(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        # call test code
        code(
            error_code=state.output.error_code,
            error_code_2=locals.error_code_2,
            error_code_3=locals.error_code_3,
            maximum_updraft_origin_level=locals.maximum_updraft_origin_level,
            kstabm=locals.kstabm,
            t_excess=state.input.t_excess,
            t_excess_local=locals.t_excess,
            vapor_excess=state.input.vapor_excess,
            vapor_excess_local=locals.vapor_excess,
            ocean_fraction=state.input.ocean_fraction,
            ocean_fraction_local=locals.ocean_fraction,
            t_old=state.input_output.t_old,
            t_new=locals.t_new,
            t_new_pbl=locals.t_new_pbl,
            vapor_old=state.input_output.vapor_old,
            vapor_forced=locals.vapor_forced,
            vapor_forced_pbl=locals.vapor_forced_pbl,
            downdraft_saturation_vapor_forced=locals.downdraft_saturation_vapor_forced,
            grid_scale_forcing_t=state.input.grid_scale_forcing_t,
            grid_scale_forcing_vapor=state.input.grid_scale_forcing_vapor,
            subgrid_scale_forcing_t=state.input.subgrid_scale_forcing_t,
            subgrid_scale_forcing_vapor=state.input.subgrid_scale_forcing_vapor,
            dmoist_static_energydt=locals.dmoist_static_energydt,
            cloud_moist_static_energy_downdraft_forced=locals.cloud_moist_static_energy_downdraft_forced,
            cloud_moist_static_energy_forced_transported=locals.cloud_moist_static_energy_forced_transported,
            cap_max=locals.cap_max,
            cap_max_increment=locals.cap_max_increment,
            geopotential_height=state.input_output.geopotential_height_forced,
            geopotential_height_local=locals.geopotential_height,
            geopotential_height_modified_local=locals.geopotential_height_modified,
            cloud_workfunction_0=locals.cloud_workfunction_0,
            cloud_workfunction_1=locals.cloud_workfunction_1,
            cloud_workfunction_2=locals.cloud_workfunction_2,
            cloud_workfunction_3=locals.cloud_workfunction_3,
            cloud_workfunction_0_pbl=locals.cloud_workfunction_0_pbl,
            cloud_workfunction_1_pbl=locals.cloud_workfunction_1_pbl,
            cloud_workfunction_1_fa=locals.cloud_workfunction_1_fa,
            cin_1=locals.cin_1,
            k_x_modified=locals.k_x_modified,
            epsilon_forced=state.output.epsilon_forced,
            epsilon_local=locals.epsilon,
            epsilon_min=locals.epsilon_min,
            epsilon_max=locals.epsilon_max,
            pbl_time_scale=locals.pbl_time_scale,
            t_wetbulb=locals.t_wetbulb,
            vapor_wetbulb=locals.vapor_wetbulb,
            cape_removal_time_scale=locals.cape_removal_time_scale,
            f_dicycle_modified=locals.f_dicycle_modified,
            add_buoyancy=locals.add_buoyancy,
            scale_dependence_factor=state.output.scale_dependence_factor,
            scale_dependence_factor_downdraft=locals.scale_dependence_factor_downdraft,
            c1d=locals.c1d,
            evaporation_below_cloud_base=locals.evaporation_below_cloud_base,
            mass_flux_ensemble=locals.mass_flux_ensemble,
            precipitation_ensemble=locals.precipitation_ensemble,
            precip=state.output.precip,
            lightning_density=state.output.lightning_density,
            seed_convection=state.input.seed_convection,
            grid_length=state.input_output.grid_length,
            random_number=locals.random_number,
            lateral_entrainment_rate=state.input.lateral_entrainment_rate,
            entrainment_rate=state.output.entrainment_rate,
            detrainment_function_updraft=locals.detrainment_function_updraft,
            arbitrary_numerical_parameter=locals.arbitrary_numerical_parameter,
            plume_dependent_constants=plume_dependent_constants,
            plume=plume,
        )

        # write output
        outputs = {
            "geopotential_height_forced": state.input_output.geopotential_height_forced.field[:],
            "epsilon_forced": state.output.epsilon_forced.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "precip": state.output.precip.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "scale_dependence_factor": state.output.scale_dependence_factor.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "lightning_density": state.output.lightning_density.field[:],
            "error_code": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "grid_length": state.input_output.grid_length.field[:],
            "lateral_entrainment_rate": state.input.lateral_entrainment_rate.field[:],
            "entrainment_rate": state.output.entrainment_rate.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_t_excess": locals.t_excess.field[:],
            "local_vapor_excess": locals.vapor_excess.field[:],
            "local_t_new": locals.t_new.field[:],
            "local_vapor_forced": locals.vapor_forced.field[:],
            "local_t_new_pbl": locals.t_new_pbl.field[:],
            "local_vapor_forced_pbl": locals.vapor_forced_pbl.field[:],
            "local_dmoist_static_energydt": locals.dmoist_static_energydt.field[:],
            "local_maximum_updraft_origin_level": locals.maximum_updraft_origin_level.field[:] + 1,
            "local_kstabm": locals.kstabm.field[:] + 1,
            "local_ocean_fraction": locals.ocean_fraction.field[:],
            "local_error_code_2": locals.error_code_2.field[:],
            "local_error_code_3": locals.error_code_3.field[:],
            "local_cap_max": locals.cap_max.field[:],
            "local_cap_max_increment": locals.cap_max_increment.field[:],
            "local_geopotential_height": locals.geopotential_height.field[:],
            "local_geopotential_height_modified": locals.geopotential_height_modified.field[:],
            "local_cloud_workfunction_0": locals.cloud_workfunction_0.field[:],
            "local_cloud_workfunction_0_pbl": locals.cloud_workfunction_0_pbl.field[:],
            "local_cloud_workfunction_1": locals.cloud_workfunction_1.field[:],
            "local_cloud_workfunction_1_pbl": locals.cloud_workfunction_1_pbl.field[:],
            "local_cloud_workfunction_1_fa": locals.cloud_workfunction_1_fa.field[:],
            "local_cloud_workfunction_2": locals.cloud_workfunction_2.field[:],
            "local_cloud_workfunction_3": locals.cloud_workfunction_3.field[:],
            "local_cin_1": locals.cin_1.field[:],
            "local_k_x_modified": locals.k_x_modified.field[:],
            "local_epsilon": locals.epsilon.field[:],
            "local_epsilon_min": locals.epsilon_min.field[:],
            "local_epsilon_max": locals.epsilon_max.field[:],
            "local_pbl_time_scale": locals.pbl_time_scale.field[:],
            "local_t_wetbulb": locals.t_wetbulb.field[:],
            "local_vapor_wetbulb": locals.vapor_wetbulb.field[:],
            "local_cape_removal_time_scale": locals.cape_removal_time_scale.field[:],
            "local_f_dicycle_modified": locals.f_dicycle_modified.field[:],
            "local_add_buoyancy": locals.add_buoyancy.field[:],
            "local_cloud_moist_static_energy_downdraft_forced": locals.cloud_moist_static_energy_downdraft_forced.field[
                :
            ],
            "local_downdraft_saturation_vapor_forced": locals.downdraft_saturation_vapor_forced.field[:],
            "local_cloud_moist_static_energy_forced_transported": locals.cloud_moist_static_energy_forced_transported.field[
                :
            ],
            "local_c1d": locals.c1d.field[:],
            "local_evaporation_below_cloud_base": locals.evaporation_below_cloud_base.field[:],
            "local_mass_flux_ensemble": locals.mass_flux_ensemble.field[:],
            "local_precipitation_ensemble": locals.precipitation_ensemble.field[:],
            "local_scale_dependence_factor_downdraft": locals.scale_dependence_factor_downdraft.field[:],
            "local_random_number": locals.random_number.field[:],
            "local_detrainment_function_updraft": locals.detrainment_function_updraft.field[:],
            "local_arbitrary_numerical_parameter": locals.arbitrary_numerical_parameter.field[:],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_Setup_shallow(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_Setup_mid(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_Setup_deep(TranslateFortranData2Py):
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
