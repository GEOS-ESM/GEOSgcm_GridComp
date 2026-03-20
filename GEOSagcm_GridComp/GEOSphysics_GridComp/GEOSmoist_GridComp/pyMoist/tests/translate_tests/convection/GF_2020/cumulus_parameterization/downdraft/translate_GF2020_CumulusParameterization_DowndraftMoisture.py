from f90nml import Namelist

from ndsl import StencilFactory
from ndsl.constants import I_DIM, J_DIM, K_DIM
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
from pyMoist.convection.GF_2020.cumulus_parameterization.downdraft import downdraft_moisture
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
            "downdraft_origin_level": {},
            "local_t_cloud_levels_forced": {},
            "local_t_wetbulb": {},
            "local_vapor_forced": {},
            "local_vapor_cloud_levels_forced": {},
            "local_env_saturation_mixing_ratio_cloud_levels_forced": {},
            "local_cloud_total_water_after_entrainment_forced": {},
            "local_cloud_total_water_after_entrainment_downdraft_forced": {},
            "local_downdraft_saturation_vapor_forced": {},
            "local_vapor_wetbulb": {},
            "normalized_massflux_downdraft_forced": {},
            "local_env_moist_static_energy_forced": {},
            "local_env_saturation_moist_static_energy_cloud_levels_forced": {},
            "local_cloud_moist_static_energy_downdraft_forced": {},
            "evaporate_in_downdraft_forced": {},
            "local_geopotential_height_cloud_levels_forced": {},
            "mass_entrainment_downdraft_forced": {},
            "mass_detrainment_downdraft_forced": {},
            "local_gamma_cloud_levels_forced": {},
            "total_normalized_integrated_condensate_forced": {},
            "local_total_normalized_integrated_evaporate_forced": {},
            "local_buoyancy": {},
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
        state.output.downdraft_origin_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["downdraft_origin_level"] - 1
        )
        locals.t_cloud_levels_forced.data[:] = inputs["local_t_cloud_levels_forced"]
        locals.t_wetbulb.data[:] = inputs["local_t_wetbulb"]
        locals.vapor_forced.data[:] = inputs["local_vapor_forced"]
        locals.vapor_cloud_levels_forced.data[:] = inputs["local_vapor_cloud_levels_forced"]
        locals.environment_saturation_mixing_ratio_cloud_levels_forced.data[:] = inputs[
            "local_env_saturation_mixing_ratio_cloud_levels_forced"
        ]
        locals.cloud_total_water_after_entrainment_forced.data[:] = inputs[
            "local_cloud_total_water_after_entrainment_forced"
        ]
        locals.cloud_total_water_after_entrainment_downdraft_forced.data[:] = inputs[
            "local_cloud_total_water_after_entrainment_downdraft_forced"
        ]
        locals.downdraft_saturation_vapor_forced.data[:] = inputs["local_downdraft_saturation_vapor_forced"]
        locals.vapor_wetbulb.data[:] = inputs["local_vapor_wetbulb"]
        state.output.normalized_massflux_downdraft_forced.data[
            :, :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["normalized_massflux_downdraft_forced"]
        locals.environment_moist_static_energy_forced.data[:] = inputs["local_env_moist_static_energy_forced"]
        locals.environment_saturation_moist_static_energy_cloud_levels_forced.data[:] = inputs[
            "local_env_saturation_moist_static_energy_cloud_levels_forced"
        ]
        locals.cloud_moist_static_energy_downdraft_forced.data[:] = inputs[
            "local_cloud_moist_static_energy_downdraft_forced"
        ]
        state.output.evaporate_in_downdraft_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["evaporate_in_downdraft_forced"]
        )
        locals.geopotential_height_cloud_levels_forced.data[:] = inputs[
            "local_geopotential_height_cloud_levels_forced"
        ]
        state.output.mass_entrainment_downdraft_forced.data[
            :, :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["mass_entrainment_downdraft_forced"]
        state.output.mass_detrainment_downdraft_forced.data[
            :, :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["mass_detrainment_downdraft_forced"]
        locals.gamma_cloud_levels_forced.data[:] = inputs["local_gamma_cloud_levels_forced"]
        state.output.total_normalized_integrated_condensate_forced.data[
            :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["total_normalized_integrated_condensate_forced"]
        locals.total_normalized_integrated_evaporate_forced.data[:,] = inputs[
            "local_total_normalized_integrated_evaporate_forced"
        ]
        locals.buoyancy.data[:] = inputs["local_buoyancy"]

        # initialize test code
        code = self.stencil_factory.from_dims_halo(
            func=downdraft_moisture,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "USE_WETBULB": cumulus_parameterization_config.USE_WETBULB,
                "ZERO_DIFF": cumulus_parameterization_config.ZERO_DIFF,
                "EVAP_FIX": cumulus_parameterization_config.EVAP_FIX,
            },
        )

        # call test code
        if plume_dependent_constants.ENABLE_PLUME == 1:
            if cumulus_parameterization_config.FIRST_GUESS_W == 0:
                code(
                    error_code=state.output.error_code,
                    downdraft_origin_level=state.output.downdraft_origin_level,
                    t_cloud_levels_forced=locals.t_cloud_levels_forced,
                    t_wetbulb=locals.t_wetbulb,
                    vapor_forced=locals.vapor_forced,
                    vapor_cloud_levels_forced=locals.vapor_cloud_levels_forced,
                    environment_saturation_mixing_ratio_cloud_levels_forced=locals.environment_saturation_mixing_ratio_cloud_levels_forced,
                    cloud_total_water_after_entrainment_forced=locals.cloud_total_water_after_entrainment_forced,
                    cloud_total_water_after_entrainment_downdraft_forced=locals.cloud_total_water_after_entrainment_downdraft_forced,
                    downdraft_saturation_vapor_forced=locals.downdraft_saturation_vapor_forced,
                    vapor_wetbulb=locals.vapor_wetbulb,
                    normalized_massflux_downdraft_forced=state.output.normalized_massflux_downdraft_forced,
                    environment_moist_static_energy_forced=locals.environment_moist_static_energy_forced,
                    environment_saturation_moist_static_energy_cloud_levels_forced=locals.environment_saturation_moist_static_energy_cloud_levels_forced,
                    cloud_moist_static_energy_downdraft_forced=locals.cloud_moist_static_energy_downdraft_forced,
                    evaporate_in_downdraft_forced=state.output.evaporate_in_downdraft_forced,
                    geopotential_height_cloud_levels_forced=locals.geopotential_height_cloud_levels_forced,
                    mass_entrainment_downdraft_forced=state.output.mass_entrainment_downdraft_forced,
                    mass_detrainment_downdraft_forced=state.output.mass_detrainment_downdraft_forced,
                    gamma_cloud_levels_forced=locals.gamma_cloud_levels_forced,
                    total_normalized_integrated_condensate_forced=state.output.total_normalized_integrated_condensate_forced,
                    total_normalized_integrated_evaporate_forced=locals.total_normalized_integrated_evaporate_forced,
                    buoyancy=locals.buoyancy,
                    plume=plume_dependent_constants.PLUME_INDEX,
                )

        # write output
        outputs = {
            "error_code": state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX],
            "downdraft_origin_level": state.output.downdraft_origin_level.data[
                :, :, plume_dependent_constants.PLUME_INDEX
            ]
            + 1,
            "local_t_cloud_levels_forced": locals.t_cloud_levels_forced.data[:],
            "local_t_wetbulb": locals.t_wetbulb.data[:],
            "local_vapor_forced": locals.vapor_forced.data[:],
            "local_vapor_cloud_levels_forced": locals.vapor_cloud_levels_forced.data[:],
            "local_env_saturation_mixing_ratio_cloud_levels_forced": locals.environment_saturation_mixing_ratio_cloud_levels_forced.data[
                :
            ],
            "local_cloud_total_water_after_entrainment_forced": locals.cloud_total_water_after_entrainment_forced.data[
                :
            ],
            "local_cloud_total_water_after_entrainment_downdraft_forced": locals.cloud_total_water_after_entrainment_downdraft_forced.data[
                :
            ],
            "local_downdraft_saturation_vapor_forced": locals.downdraft_saturation_vapor_forced.data[:],
            "local_vapor_wetbulb": locals.vapor_wetbulb.data[:],
            "normalized_massflux_downdraft_forced": state.output.normalized_massflux_downdraft_forced.data[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_env_moist_static_energy_forced": locals.environment_moist_static_energy_forced.data[:],
            "local_env_saturation_moist_static_energy_cloud_levels_forced": locals.environment_saturation_moist_static_energy_cloud_levels_forced.data[
                :
            ],
            "local_cloud_moist_static_energy_downdraft_forced": locals.cloud_moist_static_energy_downdraft_forced.data[
                :
            ],
            "evaporate_in_downdraft_forced": state.output.evaporate_in_downdraft_forced.data[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_geopotential_height_cloud_levels_forced": locals.geopotential_height_cloud_levels_forced.data[
                :
            ],
            "mass_entrainment_downdraft_forced": state.output.mass_entrainment_downdraft_forced.data[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "mass_detrainment_downdraft_forced": state.output.mass_detrainment_downdraft_forced.data[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_gamma_cloud_levels_forced": locals.gamma_cloud_levels_forced.data[:],
            "total_normalized_integrated_condensate_forced": state.output.total_normalized_integrated_condensate_forced.data[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_total_normalized_integrated_evaporate_forced": locals.total_normalized_integrated_evaporate_forced.data[
                :
            ],
            "local_buoyancy": locals.buoyancy.data[:],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_DowndraftMoisture_shallow(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_DowndraftMoisture_mid(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_DowndraftMoisture_deep(TranslateFortranData2Py):
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
