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
from pyMoist.convection.GF_2020.cumulus_parameterization.convective_tracers import AtmosphericComposition

from pyMoist.convection.GF_2020.cumulus_parameterization.setup.set_constants import set_constants
from pyMoist.convection_tracers import ConvectionTracers
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
            "error_code": {},
            "cloud_top_level": {},
            "updraft_origin_level": {},
            "downdraft_origin_level": {},
            "ocean_fraction": {},
            "p_forced": {},
            "p_cloud_levels_forced": {},
            "local_geopotential_height_cloud_levels": {},
            "local_environment_massflux": {},
            "normalized_massflux_updraft_forced": {},
            "normalized_massflux_downdraft_forced": {},
            "mass_entrainment_updraft_forced": {},
            "mass_detrainment_updraft_forced": {},
            "mass_entrainment_downdraft_forced": {},
            "mass_detrainment_downdraft_forced": {},
            "local_vertical_velocity_3d": {},
            "total_normalized_integrated_condensate_forced": {},
            "local_total_normalized_integrated_evaporate_forced": {},
            "evaporate_in_downdraft_forced": {},
            "epsilon_forced": {},
            # # THESE BREAK BECAUSE TRANSLATE TEST CANNOT HANDLE 4D FIELDS
            # "chemistry_tracers": {},
            # "chemistry_tracers_output": {},
            # "local_chemistry_tracers_cloud_levels": {},
            # "local_chemistry_tracers_sc_updraft": {},
            # "local_chemistry_tracers_sc_downdraft": {},
            # "local_chemistry_tracers_pw_updraft": {},
            # "local_chemistry_tracers_pw_downdraft": {},
            # "local_chemistry_tracers_total_pw_updraft": {},
            # "local_chemistry_tracers_total_pw_downdraft": {},
        }

        out_vars.update(in_vars["data_vars"])
        out_vars.update(
            {
                "chemistry_tracers": {},
                "chemistry_tracers_output": {},
                "local_chemistry_tracers_cloud_levels": {},
                "local_chemistry_tracers_sc_updraft": {},
                "local_chemistry_tracers_sc_downdraft": {},
                "local_chemistry_tracers_pw_updraft": {},
                "local_chemistry_tracers_pw_downdraft": {},
                "local_chemistry_tracers_total_pw_updraft": {},
                "local_chemistry_tracers_total_pw_downdraft": {},
            }
        )

    def __call__(
        self,
        constants: dict,
        cu_param_constants: dict,
        convection_tracers_input: dict,
        plume: str,
        ddim_fields: dict,
        **inputs,
    ):
        # initalize constants
        config = GF2020Config(SINGLE_COLUMN_MODE=False, **constants)
        cumulus_parameterization_config = GF2020CumulusParameterizationConfig(**cu_param_constants)
        plume_dependent_constants = GF2020PlumeDependentConstants()
        plume_dependent_constants = set_constants(
            cumulus_parameterization_config, plume_dependent_constants, plume
        )

        # initalize convection tracers
        convection_tracers = ConvectionTracers.ones(
            self.quantity_factory,
            data_dimensions={
                "tracers": config.NUMBER_OF_TRACERS,
                "size_three_dimension": 3,
                "size_four_dimension": 4,
            },
        )

        convection_tracers.tracers.field[:] = np.moveaxis(convection_tracers_input["tracers"], 0, 3)
        convection_tracers.vect_hcts.field[:] = convection_tracers_input["vect_hcts"]
        convection_tracers.kc_scal.field[:] = convection_tracers_input["kc_scal"]
        convection_tracers.fscav.field[:] = convection_tracers_input["fscav"]
        convection_tracers.convfaci2g.field[:] = convection_tracers_input["convfaci2g"]
        convection_tracers.retfactor.field[:] = convection_tracers_input["retfactor"]
        convection_tracers.liq_and_gas.field[:] = convection_tracers_input["liq_and_gas"]
        convection_tracers.online_cldliq.field[:] = convection_tracers_input["online_cldliq"]
        convection_tracers.online_vud.field[:] = convection_tracers_input["online_vud"]
        convection_tracers.ftemp_threshold.field[:] = convection_tracers_input["ftemp_threshold"]
        convection_tracers.use_gcc_washout.field[:] = convection_tracers_input["use_gcc_washout"]
        convection_tracers.use_gocart.field[:] = convection_tracers_input["use_gocart"]
        convection_tracers.is_wetdep.field[:] = convection_tracers_input["is_wetdep"]

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
        state.output.cloud_top_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["cloud_top_level"] - 1
        )
        state.output.updraft_origin_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["updraft_origin_level"] - 1
        )
        state.output.downdraft_origin_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["downdraft_origin_level"] - 1
        )
        state.input.ocean_fraction.data[:] = inputs["ocean_fraction"]
        state.input_output.p_forced.data[:] = inputs["p_forced"]
        state.output.p_cloud_levels_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "p_cloud_levels_forced"
        ]
        locals.geopotential_height_cloud_levels.data[:] = inputs["local_geopotential_height_cloud_levels"]
        locals.environment_massflux.data[:] = inputs["local_environment_massflux"]
        state.output.normalized_massflux_updraft_forced.data[
            :, :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["normalized_massflux_updraft_forced"]
        state.output.normalized_massflux_downdraft_forced.data[
            :, :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["normalized_massflux_downdraft_forced"]
        state.output.mass_entrainment_updraft_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["mass_entrainment_updraft_forced"]
        )
        state.output.mass_detrainment_updraft_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["mass_detrainment_updraft_forced"]
        )
        state.output.mass_entrainment_downdraft_forced.data[
            :, :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["mass_entrainment_downdraft_forced"]
        state.output.mass_detrainment_downdraft_forced.data[
            :, :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["mass_detrainment_downdraft_forced"]
        locals.vertical_velocity_3d.data[:] = inputs["local_vertical_velocity_3d"]
        state.output.total_normalized_integrated_condensate_forced.data[
            :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["total_normalized_integrated_condensate_forced"]
        locals.total_normalized_integrated_evaporate_forced.data[:] = inputs[
            "local_total_normalized_integrated_evaporate_forced"
        ]
        state.output.evaporate_in_downdraft_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["evaporate_in_downdraft_forced"]
        )
        state.output.epsilon_forced.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "epsilon_forced"
        ]
        state.input_output.chemistry_tracers.field[:] = ddim_fields["chemistry_tracers"]
        state.input_output.chemistry_tracers_output.field[
            :, :, :, plume_dependent_constants.PLUME_INDEX, :
        ] = ddim_fields["chemistry_tracers_output"]
        locals.chemistry_tracers_cloud_levels.field[:] = ddim_fields["local_chemistry_tracers_cloud_levels"]
        locals.chemistry_tracers_sc_updraft.field[:] = ddim_fields["local_chemistry_tracers_sc_updraft"]
        locals.chemistry_tracers_sc_downdraft.field[:] = ddim_fields["local_chemistry_tracers_sc_downdraft"]
        locals.chemistry_tracers_pw_updraft.field[:] = ddim_fields["local_chemistry_tracers_pw_updraft"]
        locals.chemistry_tracers_pw_downdraft.field[:] = ddim_fields["local_chemistry_tracers_pw_downdraft"]
        locals.chemistry_tracers_total_pw_updraft.field[:] = ddim_fields[
            "local_chemistry_tracers_total_pw_updraft"
        ]
        locals.chemistry_tracers_total_pw_downdraft.field[:] = ddim_fields[
            "local_chemistry_tracers_total_pw_downdraft"
        ]

        code = AtmosphericComposition(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        if plume_dependent_constants.ENABLE_PLUME == 1:
            code(
                error_code=state.output.error_code,
                cloud_top_level=state.output.cloud_top_level,
                updraft_origin_level=state.output.updraft_origin_level,
                downdraft_origin_level=state.output.downdraft_origin_level,
                ocean_fraction=state.input.ocean_fraction,
                p_forced=state.input_output.p_forced,
                p_cloud_levels_forced=state.output.p_cloud_levels_forced,
                geopotential_height_cloud_levels=locals.geopotential_height_cloud_levels,
                environment_massflux=locals.environment_massflux,
                normalized_massflux_updraft_forced=state.output.normalized_massflux_updraft_forced,
                normalized_massflux_downdraft_forced=state.output.normalized_massflux_downdraft_forced,
                mass_entrainment_updraft_forced=state.output.mass_entrainment_updraft_forced,
                mass_detrainment_updraft_forced=state.output.mass_detrainment_updraft_forced,
                mass_entrainment_downdraft_forced=state.output.mass_entrainment_downdraft_forced,
                mass_detrainment_downdraft_forced=state.output.mass_detrainment_downdraft_forced,
                vertical_velocity_3d=locals.vertical_velocity_3d,
                total_normalized_integrated_condensate_forced=state.output.total_normalized_integrated_condensate_forced,
                total_normalized_integrated_evaporate_forced=locals.total_normalized_integrated_evaporate_forced,
                evaporate_in_downdraft_forced=state.output.evaporate_in_downdraft_forced,
                epsilon_forced=state.output.epsilon_forced,
                chemistry_tracers=state.input_output.chemistry_tracers,
                chemistry_tracers_output=state.input_output.chemistry_tracers_output,
                chemistry_tracers_cloud_levels=locals.chemistry_tracers_cloud_levels,
                chemistry_tracers_sc_updraft=locals.chemistry_tracers_sc_updraft,
                chemistry_tracers_sc_downdraft=locals.chemistry_tracers_sc_downdraft,
                chemistry_tracers_pw_updraft=locals.chemistry_tracers_pw_updraft,
                chemistry_tracers_pw_downdraft=locals.chemistry_tracers_pw_downdraft,
                chemistry_tracers_total_pw_updraft=locals.chemistry_tracers_total_pw_updraft,
                chemistry_tracers_total_pw_downdraft=locals.chemistry_tracers_total_pw_downdraft,
                convection_tracers=convection_tracers,
                plume_dependent_constants=plume_dependent_constants,
            )

        outputs = {
            "error_code": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "cloud_top_level": state.output.cloud_top_level.field[:, :, plume_dependent_constants.PLUME_INDEX]
            + 1,
            "updraft_origin_level": state.output.updraft_origin_level.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ]
            + 1,
            "downdraft_origin_level": state.output.downdraft_origin_level.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ]
            + 1,
            "ocean_fraction": state.input.ocean_fraction.field[:],
            "p_forced": state.input_output.p_forced.field[:],
            "p_cloud_levels_forced": state.output.p_cloud_levels_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_geopotential_height_cloud_levels": locals.geopotential_height_cloud_levels.field[:],
            "local_environment_massflux": locals.environment_massflux.field[:],
            "normalized_massflux_updraft_forced": state.output.normalized_massflux_updraft_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "normalized_massflux_downdraft_forced": state.output.normalized_massflux_downdraft_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "mass_entrainment_updraft_forced": state.output.mass_entrainment_updraft_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "mass_detrainment_updraft_forced": state.output.mass_detrainment_updraft_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "mass_entrainment_downdraft_forced": state.output.mass_entrainment_downdraft_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "mass_detrainment_downdraft_forced": state.output.mass_detrainment_downdraft_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_vertical_velocity_3d": locals.vertical_velocity_3d.field[:],
            "total_normalized_integrated_condensate_forced": state.output.total_normalized_integrated_condensate_forced.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_total_normalized_integrated_evaporate_forced": locals.total_normalized_integrated_evaporate_forced.field[
                :
            ],
            "evaporate_in_downdraft_forced": state.output.evaporate_in_downdraft_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "epsilon_forced": state.output.epsilon_forced.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "chemistry_tracers": state.input_output.chemistry_tracers.field[:],
            "chemistry_tracers_output": state.input_output.chemistry_tracers_output.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX, :
            ],
            "local_chemistry_tracers_cloud_levels": locals.chemistry_tracers_cloud_levels.field[:],
            "local_chemistry_tracers_sc_updraft": locals.chemistry_tracers_sc_updraft.field[:],
            "local_chemistry_tracers_sc_downdraft": locals.chemistry_tracers_sc_downdraft.field[:],
            "local_chemistry_tracers_pw_updraft": locals.chemistry_tracers_pw_updraft.field[:],
            "local_chemistry_tracers_pw_downdraft": locals.chemistry_tracers_pw_downdraft.field[:],
            "local_chemistry_tracers_total_pw_updraft": locals.chemistry_tracers_total_pw_updraft.field[:],
            "local_chemistry_tracers_total_pw_downdraft": locals.chemistry_tracers_total_pw_downdraft.field[
                :
            ],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_AtmosphericComposition_shallow(TranslateFortranData2Py):
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
        self.convection_tracers = data_loader.load("GF2020_ConvectionTracers")
        self.ddim_fields = data_loader.load(
            "GF2020_CumulusParameterization_AtmosphericComposition_shallow-In"
        )

    def compute_func(self, **inputs):
        outputs = self.test_core(
            self.constants,
            self.cu_param_constants,
            self.convection_tracers,
            "shallow",
            ddim_fields=self.ddim_fields,
            **inputs,
        )

        return outputs


class TranslateGF2020_CumulusParameterization_AtmosphericComposition_mid(TranslateFortranData2Py):
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
        self.convection_tracers = data_loader.load("GF2020_ConvectionTracers")
        self.ddim_fields = data_loader.load("GF2020_CumulusParameterization_AtmosphericComposition_mid-In")

    def compute_func(self, **inputs):
        outputs = self.test_core(
            self.constants,
            self.cu_param_constants,
            self.convection_tracers,
            "mid",
            ddim_fields=self.ddim_fields,
            **inputs,
        )

        return outputs


class TranslateGF2020_CumulusParameterization_AtmosphericComposition_deep(TranslateFortranData2Py):
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
        self.convection_tracers = data_loader.load("GF2020_ConvectionTracers")
        self.ddim_fields = data_loader.load("GF2020_CumulusParameterization_AtmosphericComposition_deep-In")

    def compute_func(self, **inputs):
        outputs = self.test_core(
            self.constants,
            self.cu_param_constants,
            self.convection_tracers,
            "deep",
            ddim_fields=self.ddim_fields,
            **inputs,
        )

        return outputs
