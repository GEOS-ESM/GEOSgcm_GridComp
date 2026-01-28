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
            "cloud_top_level": {},
            "ocean_fraction": {},
            "p_forced_atmoscomp": {},
            "updraft_origin_level_atmoscomp": {},
            "normalized_massflux_updraft_forced_atmoscomp": {},
            "mass_entrainment_updraft_forced_atmoscomp": {},
            "mass_detrainment_updraft_forced_atmoscomp": {},
            "local_geopotential_height_cloud_levels_atmoscomp": {},
            "local_vertical_velocity_3d_atmoscomp": {},
            "p_cloud_levels_forced_atmoscomp": {},
            "downdraft_origin_level": {},
            "evaporate_in_downdraft_forced_atmoscomp": {},
            "normalized_massflux_downdraft_forced_atmoscomp": {},
            "total_normalized_integrated_condensate_forced_atmoscomp": {},
            "total_normalized_integrated_evaporate_forced_atmoscomp": {},
            "mass_entrainment_downdraft_forced_atmoscomp": {},
            "mass_detrainment_downdraft_forced_atmoscomp": {},
            "epsilon_forced_atmoscomp": {},
            "local_environment_massflux_atmoscomp": {},
        }

        out_vars.update(in_vars["data_vars"])
        out_vars.update(
            {
                "chemistry_tracers": {},
                "local_chemistry_tracers_cloud_levels": {},
                "local_sc_up_chem": {},
                "local_pw_up_chem": {},
                "local_tot_pw_up_chem": {},
                "CNV_Tracers_fscav": {},
                "CNV_Tracers_Vect_hcts": {},
                "local_pw_dn": {},
                "local_sc_dn": {},
                "local_tot_pw_dn_chem": {},
                "local_out_chem": {},
            }
        )

    def __call__(self, constants: dict, cu_param_constants: dict, plume: str, ddim_fields: dict, **inputs):
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
        state.input_output.chemistry_tracers.field[:] = ddim_fields["chemistry_tracers"]
        locals.chemistry_tracers_cloud_levels.field[:] = ddim_fields["local_chemistry_tracers_cloud_levels"]
        locals.sc_up_chem.field[:] = ddim_fields["local_sc_up_chem"]
        locals.pw_up_chem.field[:] = ddim_fields["local_pw_up_chem"]
        locals.tot_pw_up_chem.field[:] = ddim_fields["local_tot_pw_up_chem"]
        locals.CNV_Tracers_fscav.field[:] = ddim_fields["CNV_Tracers_fscav"]
        locals.CNV_Tracers_Vect_hcts.field[:] = ddim_fields["CNV_Tracers_Vect_hcts"]
        locals.pw_dn.field[:] = ddim_fields["local_pw_dn"]
        locals.sc_dn.field[:] = ddim_fields["local_sc_dn"]
        locals.tot_pw_dn_chem.field[:] = ddim_fields["local_tot_pw_dn_chem"]
        state.output.cloud_top_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["cloud_top_level"] - 1
        )
        state.input.ocean_fraction.data[:] = inputs["ocean_fraction"]
        state.input_output.p_forced.data[:] = inputs["p_forced_atmoscomp"]
        state.output.updraft_origin_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["updraft_origin_level_atmoscomp"] - 1
        )
        state.output.normalized_massflux_updraft_forced.data[
            :, :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["normalized_massflux_updraft_forced_atmoscomp"]
        state.output.mass_entrainment_updraft_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["mass_entrainment_updraft_forced_atmoscomp"]
        )
        state.output.mass_detrainment_updraft_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["mass_detrainment_updraft_forced_atmoscomp"]
        )
        locals.geopotential_height_cloud_levels.data[:] = inputs[
            "local_geopotential_height_cloud_levels_atmoscomp"
        ]
        locals.vertical_velocity_3d.data[:] = inputs["local_vertical_velocity_3d_atmoscomp"]
        state.output.p_cloud_levels_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "p_cloud_levels_forced_atmoscomp"
        ]
        state.output.downdraft_origin_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["downdraft_origin_level"] - 1
        )
        state.output.evaporate_in_downdraft_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["evaporate_in_downdraft_forced_atmoscomp"]
        )
        state.output.normalized_massflux_downdraft_forced.data[
            :, :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["normalized_massflux_downdraft_forced_atmoscomp"]
        state.output.total_normalized_integrated_condensate_forced.data[
            :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["total_normalized_integrated_condensate_forced_atmoscomp"]
        state.output.total_normalized_integrated_evaporate_forced.data[
            :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["total_normalized_integrated_evaporate_forced_atmoscomp"]
        state.output.mass_entrainment_downdraft_forced.data[
            :, :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["mass_entrainment_downdraft_forced_atmoscomp"]
        state.output.mass_detrainment_downdraft_forced.data[
            :, :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["mass_detrainment_downdraft_forced_atmoscomp"]
        state.output.epsilon_forced.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "epsilon_forced_atmoscomp"
        ]
        locals.environment_massflux.data[:] = inputs["local_environment_massflux_atmoscomp"]
        locals.out_chem.field[:] = ddim_fields["local_out_chem"]

        code = AtmosphericComposition(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        if plume_dependent_constants.ENABLE_PLUME == 1:
            code(
                error_code=state.output.error_code,
                chemistry_tracers=state.input_output.chemistry_tracers,
                chemistry_tracers_cloud_levels=locals.chemistry_tracers_cloud_levels,
                plume_dependent_constants=plume_dependent_constants,
                sc_up_chem=locals.sc_up_chem,
                pw_up_chem=locals.pw_up_chem,
                tot_pw_up_chem=locals.tot_pw_up_chem,
                ocean_fraction=state.input.ocean_fraction,
                AVERAGE_LAYER_DEPTH=plume_dependent_constants.AVERAGE_LAYER_DEPTH,
                updraft_origin_level=state.output.updraft_origin_level,
                po=state.input_output.p_forced,
                sc_b=locals.sc_b,
                cloud_top_level=state.output.cloud_top_level,
                normalized_massflux_updraft_forced=state.output.normalized_massflux_updraft_forced,
                mass_detrainment_updraft_forced=state.output.mass_detrainment_updraft_forced,
                mass_entrainment_updraft_forced=state.output.mass_entrainment_updraft_forced,
                geopotential_height_cloud_levels=locals.geopotential_height_cloud_levels,
                vertical_velocity_3d=locals.vertical_velocity_3d,
                CNV_Tracers_fscav=locals.CNV_Tracers_fscav,
                CNV_Tracers_Vect_hcts=locals.CNV_Tracers_Vect_hcts,
                p_cloud_levels_forced=state.output.p_cloud_levels_forced,
                sc_dn=locals.sc_dn,
                pw_dn=locals.pw_dn,
                tot_pw_dn_chem=locals.tot_pw_dn_chem,
                total_normalized_integrated_evaporate_forced=state.output.total_normalized_integrated_evaporate_forced,
                total_normalized_integrated_condensate_forced=state.output.total_normalized_integrated_condensate_forced,
                downdraft_origin_level=state.output.downdraft_origin_level,
                evaporate_in_downdraft_forced=state.output.evaporate_in_downdraft_forced,
                normalized_massflux_downdraft_forced=state.output.normalized_massflux_downdraft_forced,
                mass_detrainment_downdraft_forced=state.output.mass_detrainment_downdraft_forced,
                mass_entrainment_downdraft_forced=state.output.mass_entrainment_downdraft_forced,
                environment_massflux=locals.environment_massflux,
                ddtr=locals.ddtr,
                epsilon_forced=state.output.epsilon_forced,
                out_chem=locals.out_chem,
                trash_=locals.trash_,
                trash2_=locals.trash2_,
                evap_=locals.evap_,
                wetdep_=locals.wetdep_,
                residu_=locals.residu_,
            )

        outputs = {
            "error_code": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "chemistry_tracers": state.input_output.chemistry_tracers.field[:],
            "local_chemistry_tracers_cloud_levels": locals.chemistry_tracers_cloud_levels.field[:],
            "local_sc_up_chem": locals.sc_up_chem.field[:],
            "local_pw_up_chem": locals.pw_up_chem.field[:],
            "local_tot_pw_up_chem": locals.tot_pw_up_chem.field[:],
            "CNV_Tracers_fscav": locals.CNV_Tracers_fscav.field[:],
            "CNV_Tracers_Vect_hcts": locals.CNV_Tracers_Vect_hcts.field[:],
            "local_pw_dn": locals.pw_dn.field[:],
            "local_sc_dn": locals.sc_dn.field[:],
            "local_tot_pw_dn_chem": locals.tot_pw_dn_chem.field[:],
            "cloud_top_level": state.output.cloud_top_level.data[:, :, plume_dependent_constants.PLUME_INDEX]
            + 1,
            "ocean_fraction": state.input.ocean_fraction.data[:],
            "p_forced_atmoscomp": state.input_output.p_forced.data[:],
            "updraft_origin_level_atmoscomp": state.output.updraft_origin_level.data[
                :, :, plume_dependent_constants.PLUME_INDEX
            ]
            + 1,
            "normalized_massflux_updraft_forced_atmoscomp": state.output.normalized_massflux_updraft_forced.data[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "mass_entrainment_updraft_forced_atmoscomp": state.output.mass_entrainment_updraft_forced.data[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "mass_detrainment_updraft_forced_atmoscomp": state.output.mass_detrainment_updraft_forced.data[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_geopotential_height_cloud_levels_atmoscomp": locals.geopotential_height_cloud_levels.data[
                :
            ],
            "local_vertical_velocity_3d_atmoscomp": locals.vertical_velocity_3d.data[:],
            "p_cloud_levels_forced_atmoscomp": state.output.p_cloud_levels_forced.data[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "locals_sc_dn": locals.sc_dn.field[:],
            "locals_pw_dn": locals.pw_dn.field[:],
            "locals_tot_pw_dn_chem": locals.tot_pw_dn_chem.field[:],
            "total_normalized_integrated_evaporate_forced": state.output.total_normalized_integrated_evaporate_forced.data[
                :
            ],
            "total_normalized_integrated_condensate_forced": state.output.total_normalized_integrated_condensate_forced.data[
                :
            ],
            "downdraft_origin_level": state.output.downdraft_origin_level.data[
                :, :, plume_dependent_constants.PLUME_INDEX
            ]
            + 1,
            "evaporate_in_downdraft_forced": state.output.evaporate_in_downdraft_forced.data[:],
            "normalized_massflux_downdraft_forced": state.output.normalized_massflux_downdraft_forced.data[:],
            "mass_entrainment_downdraft_forced_atmoscomp": state.output.mass_entrainment_downdraft_forced.data[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "mass_detrainment_downdraft_forced_atmoscomp": state.output.mass_detrainment_downdraft_forced.data[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "epsilon_forced_atmoscomp": state.output.epsilon_forced.data[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_environment_massflux_atmoscomp": locals.environment_massflux.data[:],
            "local_out_chem": locals.out_chem.field[:],
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
        self.ddim_fields = data_loader.load(
            "GF2020_CumulusParameterization_AtmosphericComposition_shallow-In"
        )

    def compute_func(self, **inputs):
        outputs = self.test_core(
            self.constants, self.cu_param_constants, "shallow", ddim_fields=self.ddim_fields, **inputs
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
        self.ddim_fields = data_loader.load("GF2020_CumulusParameterization_AtmosphericComposition_mid-In")

    def compute_func(self, **inputs):
        outputs = self.test_core(
            self.constants, self.cu_param_constants, "mid", ddim_fields=self.ddim_fields, **inputs
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
        self.ddim_fields = data_loader.load("GF2020_CumulusParameterization_AtmosphericComposition_deep-In")

    def compute_func(self, **inputs):
        outputs = self.test_core(
            self.constants, self.cu_param_constants, "deep", ddim_fields=self.ddim_fields, **inputs
        )

        return outputs
