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
from pyMoist.convection.GF_2020.cumulus_parameterization.constants import MAXENS1, MAXENS2, MAXENS3
from pyMoist.convection.GF_2020.cumulus_parameterization.get_levels.get_levels import (
    ConvectiveCloudBaseLevel,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.setup.set_constants import set_constants


class TranslateGF2020_CumulusParameterization_ConvectiveCloudBaseLevel_deep(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        self.in_vars["data_vars"] = {
            "error_code_ccb": {},
            "local_cap_max_increment_ccb": {},
            "local_cap_max_ccb": {},
            "local_env_moist_static_energy_cloud_levels_forced_ccb": {},
            "local_env_saturation_moist_static_energy_cloud_levels_forced_ccb": {},
            "local_vapor_cloud_levels_forced_ccb": {},
            "local_env_saturation_mixing_ratio_cloud_levels_forced_ccb": {},
            "p_forced_ccb": {},
            "p_cloud_levels_forced_ccb": {},
            "local_geopotential_height_cloud_levels_forced_ccb": {},
            "local_env_moist_static_energy_forced_ccb": {},
            "local_moist_static_energy_origin_level_forced_ccb": {},
            "local_vapor_forced_ccb": {},
            "local_env_saturation_mixing_ratio_forced_ccb": {},
            "entrainment_rate_ccb": {},
            "local_cloud_moist_static_energy_forced_transported_ccb": {},
            "updraft_origin_level_ccb": {},
            "local_maximum_updraft_origin_level_ccb": {},
            "lcl_level_ccb": {},
            "updraft_lfc_level_ccb": {},
            "cloud_top_ccb": {},
            "local_negative_buoyancy_depth_ccb": {},
            "local_frh_lfc_ccb": {},
            "t_perturbation_ccb": {},
            "local_start_level_ccb": {},
            "local_vapor_excess_ccb": {},
            "local_t_excess_ccb": {},
            "local_add_buoyancy_ccb": {},
            "ocean_fraction_ccb": {},
        }

        self.out_vars = self.in_vars["data_vars"].copy()

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF2020-constants")
        self.cu_param_constants = data_loader.load("GF2020_CumulusParameterization-constants")

    def compute_func(self, **inputs):
        # initalize constants
        config = GF2020Config(SINGLE_COLUMN_MODE=False, **self.constants)
        cumulus_parameterization_config = GF2020CumulusParameterizationConfig(**self.cu_param_constants)
        plume_dependent_constants = GF2020PlumeDependentConstants()
        plume_dependent_constants = set_constants(
            cumulus_parameterization_config, plume_dependent_constants, "deep"
        )

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
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["error_code_ccb"]
        locals.cap_max_increment.data[:] = inputs["local_cap_max_increment_ccb"]
        locals.cap_max.data[:] = inputs["local_cap_max_ccb"]
        locals.environment_moist_static_energy_cloud_levels_forced.data[:] = inputs[
            "local_env_moist_static_energy_cloud_levels_forced_ccb"
        ]
        locals.environment_saturation_moist_static_energy_cloud_levels_forced.data[:] = inputs[
            "local_env_saturation_moist_static_energy_cloud_levels_forced_ccb"
        ]
        locals.vapor_cloud_levels_forced.data[:] = inputs["local_vapor_cloud_levels_forced_ccb"]
        locals.environment_saturation_mixing_ratio_cloud_levels_forced.data[:] = inputs[
            "local_env_saturation_mixing_ratio_cloud_levels_forced_ccb"
        ]
        state.input_output.p_forced.data[:] = inputs["p_forced_ccb"]
        state.output.p_cloud_levels_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "p_cloud_levels_forced_ccb"
        ]
        locals.geopotential_height_cloud_levels_forced.data[:] = inputs[
            "local_geopotential_height_cloud_levels_forced_ccb"
        ]
        locals.environment_moist_static_energy_forced.data[:] = inputs[
            "local_env_moist_static_energy_forced_ccb"
        ]
        locals.moist_static_energy_origin_level_forced.data[:] = inputs[
            "local_moist_static_energy_origin_level_forced_ccb"
        ]
        locals.vapor_forced.data[:] = inputs["local_vapor_forced_ccb"]
        locals.environment_saturation_mixing_ratio_forced.data[:] = inputs[
            "local_env_saturation_mixing_ratio_forced_ccb"
        ]
        state.output.entrainment_rate.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "entrainment_rate_ccb"
        ]
        locals.cloud_moist_static_energy_forced_transported.data[:] = inputs[
            "local_cloud_moist_static_energy_forced_transported_ccb"
        ]
        state.output.updraft_origin_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "updraft_origin_level_ccb"
        ]
        locals.maximum_updraft_origin_level.data[:] = inputs["local_maximum_updraft_origin_level_ccb"]
        state.output.lcl_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["lcl_level_ccb"]
        state.output.updraft_lfc_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "updraft_lfc_level_ccb"
        ]
        state.output.cloud_top.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["cloud_top_ccb"]
        locals.negative_buoyancy_depth.data[:] = inputs["local_negative_buoyancy_depth_ccb"]
        locals.frh_lfc.data[:] = inputs["local_frh_lfc_ccb"]
        state.output.t_perturbation.data[:] = inputs["t_perturbation_ccb"]
        locals.start_level.data[:] = inputs["local_start_level_ccb"]
        locals.vapor_excess.data[:] = inputs["local_vapor_excess_ccb"]
        locals.t_excess.data[:] = inputs["local_t_excess_ccb"]
        locals.add_buoyancy.data[:] = inputs["local_add_buoyancy_ccb"]
        locals.ocean_fraction.data[:] = inputs["ocean_fraction_ccb"]

        # initalize test code
        code = ConvectiveCloudBaseLevel(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        # call test code
        if plume_dependent_constants.ENABLE_PLUME == 1:
            code(
                state=state,
                locals=locals,
                plume_dependent_constants=plume_dependent_constants,
            )

            state.output.updraft_origin_level.field[:, :, plume_dependent_constants.PLUME_INDEX] += 1
            state.output.lcl_level.field[:, :, plume_dependent_constants.PLUME_INDEX] += 1

        # write output
        outputs = {
            "error_code_ccb": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "local_cap_max_increment_ccb": locals.cap_max_increment.field[:],
            "local_cap_max_ccb": locals.cap_max.field[:],
            "local_env_moist_static_energy_cloud_levels_forced_ccb": locals.environment_moist_static_energy_cloud_levels_forced.field[
                :
            ],
            "local_env_saturation_moist_static_energy_cloud_levels_forced_ccb": locals.environment_saturation_moist_static_energy_cloud_levels_forced.field[
                :
            ],
            "local_vapor_cloud_levels_forced_ccb": locals.vapor_cloud_levels_forced.field[:],
            "local_env_saturation_mixing_ratio_cloud_levels_forced_ccb": locals.environment_saturation_mixing_ratio_cloud_levels_forced.field[
                :
            ],
            "p_forced_ccb": state.input_output.p_forced.field[:],
            "p_cloud_levels_forced_ccb": state.output.p_cloud_levels_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_geopotential_height_cloud_levels_forced_ccb": locals.geopotential_height_cloud_levels_forced.field[
                :
            ],
            "local_env_moist_static_energy_forced_ccb": locals.environment_moist_static_energy_forced.field[
                :
            ],
            "local_moist_static_energy_origin_level_forced_ccb": locals.moist_static_energy_origin_level_forced.field[
                :
            ],
            "local_vapor_forced_ccb": locals.vapor_forced.field[:],
            "local_env_saturation_mixing_ratio_forced_ccb": locals.environment_saturation_mixing_ratio_forced.field[
                :
            ],
            "entrainment_rate_ccb": state.output.entrainment_rate.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_env_saturation_mixing_ratio_forced_ccb": locals.cloud_moist_static_energy_forced_transported.field[
                :
            ],
            "updraft_origin_level_ccb": state.output.updraft_origin_level.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_maximum_updraft_origin_level_ccb": locals.maximum_updraft_origin_level.field[:],
            "lcl_level_ccb": state.output.lcl_level.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "updraft_lfc_level_ccb": state.output.updraft_lfc_level.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "cloud_top_ccb": state.output.cloud_top.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "local_negative_buoyancy_depth_ccb": locals.negative_buoyancy_depth.field[:],
            "local_frh_lfc_ccb": locals.frh_lfc.field[:],
            "t_perturbation_ccb": state.output.t_perturbation.field[:],
            "local_start_level_ccb": locals.start_level.field[:],
            "local_vapor_excess_ccb": locals.vapor_excess.field[:],
            "local_t_excess_ccb": locals.t_excess.field[:],
            "local_add_buoyancy_ccb": locals.add_buoyancy.field[:],
            "ocean_fraction_ccb": locals.ocean_fraction.field[:],
        }

        return outputs
