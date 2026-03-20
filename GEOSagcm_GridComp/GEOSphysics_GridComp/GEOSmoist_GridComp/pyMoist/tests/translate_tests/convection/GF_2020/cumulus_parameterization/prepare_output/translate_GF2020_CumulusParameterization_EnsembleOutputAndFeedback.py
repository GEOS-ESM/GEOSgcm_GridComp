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
from pyMoist.convection.GF_2020.cumulus_parameterization.locals import GF2020CumulusParameterizationLocals
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.prepare_output import ensemble_output_and_feedback
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
            "cloud_top_level": {},
            "updraft_lfc_level": {},
            "p_cloud_levels_forced": {},
            "normalized_massflux_updraft_forced": {},
            "precip": {},
            "local_effective_condensate_to_fall_forced": {},
            "cloud_base_mass_flux_modified": {},
            "scale_dependence_factor": {},
            "local_ocean_fraction": {},
            "local_f_dicycle_modified": {},
            "local_del_u_cloud_ensemble": {},
            "local_del_v_cloud_ensemble": {},
            "local_del_t_cloud_ensemble": {},
            "local_del_vapor_cloud_ensemble": {},
            "local_del_cloud_liquid_cloud_ensemble": {},
            "local_del_buoyancy_cloud_ensemble": {},
            "local_del_convective_ice_cloud_ensemble": {},
            "local_del_large_scale_ice_cloud_ensemble": {},
            "local_del_convective_liquid_cloud_ensemble": {},
            "local_del_large_scale_liquid_cloud_ensemble": {},
            "local_del_convective_cloud_fraction_cloud_ensemble": {},
            "local_del_large_scale_cloud_fraction_cloud_ensemble": {},
            "dtdt": {},
            "dvapordt": {},
            "dcloudicedt": {},
            "dudt": {},
            "dvdt": {},
            "dbuoyancydt": {},
            "dconvectiveicedt": {},
            "dlargescaleicedt": {},
            "dconvectiveliquiddt": {},
            "dlargescaleliquiddt": {},
            "dconvectivecloudfractiondt": {},
            "dlargescalecloudfractiondt": {},
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
        state.output.cloud_top_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["cloud_top_level"] - 1
        )
        state.output.updraft_lfc_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["updraft_lfc_level"] - 1
        )
        state.output.p_cloud_levels_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "p_cloud_levels_forced"
        ]
        state.output.normalized_massflux_updraft_forced.data[
            :, :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["normalized_massflux_updraft_forced"]
        state.output.precip.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["precip"]
        locals.effective_condensate_to_fall_forced.data[:] = inputs[
            "local_effective_condensate_to_fall_forced"
        ]
        state.output.cloud_base_mass_flux_modified.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "cloud_base_mass_flux_modified"
        ]
        state.output.scale_dependence_factor.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "scale_dependence_factor"
        ]
        locals.ocean_fraction.data[:] = inputs["local_ocean_fraction"]
        locals.f_dicycle_modified.data[:] = inputs["local_f_dicycle_modified"]
        locals.del_u_cloud_ensemble.data[:] = inputs["local_del_u_cloud_ensemble"]
        locals.del_v_cloud_ensemble.data[:] = inputs["local_del_v_cloud_ensemble"]
        locals.del_t_cloud_ensemble.data[:] = inputs["local_del_t_cloud_ensemble"]
        locals.del_vapor_cloud_ensemble.data[:] = inputs["local_del_vapor_cloud_ensemble"]
        locals.del_cloud_liquid_cloud_ensemble.data[:] = inputs["local_del_cloud_liquid_cloud_ensemble"]
        locals.del_buoyancy_cloud_ensemble.data[:] = inputs["local_del_buoyancy_cloud_ensemble"]
        locals.del_convective_ice_cloud_ensemble.data[:] = inputs["local_del_convective_ice_cloud_ensemble"]
        locals.del_large_scale_ice_cloud_ensemble.data[:] = inputs["local_del_large_scale_ice_cloud_ensemble"]
        locals.del_convective_liquid_cloud_ensemble.data[:] = inputs[
            "local_del_convective_liquid_cloud_ensemble"
        ]
        locals.del_large_scale_liquid_cloud_ensemble.data[:] = inputs[
            "local_del_large_scale_liquid_cloud_ensemble"
        ]
        locals.del_convective_cloud_fraction_cloud_ensemble.data[:] = inputs[
            "local_del_convective_cloud_fraction_cloud_ensemble"
        ]
        locals.del_large_scale_cloud_fraction_cloud_ensemble.data[:] = inputs[
            "local_del_large_scale_cloud_fraction_cloud_ensemble"
        ]
        state.output.dtdt.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs["dtdt"]
        state.output.dvapordt.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs["dvapordt"]
        state.output.dcloudicedt.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs["dcloudicedt"]
        state.output.dudt.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs["dudt"]
        state.output.dvdt.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs["dvdt"]
        state.output.dbuoyancydt.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs["dbuoyancydt"]
        state.output.dconvectiveicedt.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "dconvectiveicedt"
        ]
        state.output.dlargescaleicedt.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "dlargescaleicedt"
        ]
        state.output.dconvectiveliquiddt.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "dconvectiveliquiddt"
        ]
        state.output.dlargescaleliquiddt.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "dlargescaleliquiddt"
        ]
        state.output.dconvectivecloudfractiondt.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "dconvectivecloudfractiondt"
        ]
        state.output.dlargescalecloudfractiondt.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "dlargescalecloudfractiondt"
        ]
        locals.mass_flux_ensemble.data[:] = inputs["local_mass_flux_ensemble"][:, :, 0:16]
        locals.precipitation_ensemble.data[:] = inputs["local_precipitation_ensemble"][:, :, 0:16]
        locals.xff_mid.data[:] = inputs["local_xff_mid"][:, :, 0:16]

        code = self.stencil_factory.from_dims_halo(
            func=ensemble_output_and_feedback,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "DTIME": cumulus_parameterization_config.DTIME,
                "MAX_TEMP_VAPOR_TENDENCY": cumulus_parameterization_config.MAX_TEMP_VAPOR_TENDENCY,
                "APPLY_SUBSIDENCE_MICROPHYSICS": config.APPLY_SUBSIDENCE_MICROPHYSICS,
                "USE_SMOOTH_TENDENCIES": cumulus_parameterization_config.USE_SMOOTH_TENDENCIES,
            },
        )

        if plume_dependent_constants.ENABLE_PLUME == 1:
            code(
                error_code=state.output.error_code,
                error_code_2=locals.error_code_2,
                error_code_3=locals.error_code_3,
                cloud_top_level=state.output.cloud_top_level,
                updraft_lfc_level=state.output.updraft_lfc_level,
                p_cloud_levels_forced=state.output.p_cloud_levels_forced,
                normalized_massflux_updraft_forced=state.output.normalized_massflux_updraft_forced,
                precip=state.output.precip,
                effective_condensate_to_fall_forced=locals.effective_condensate_to_fall_forced,
                cloud_base_mass_flux_modified=state.output.cloud_base_mass_flux_modified,
                scale_dependence_factor=state.output.scale_dependence_factor,
                ocean_fraction=locals.ocean_fraction,
                f_dicycle_modified=locals.f_dicycle_modified,
                del_u_cloud_ensemble=locals.del_u_cloud_ensemble,
                del_v_cloud_ensemble=locals.del_v_cloud_ensemble,
                del_t_cloud_ensemble=locals.del_t_cloud_ensemble,
                del_vapor_cloud_ensemble=locals.del_vapor_cloud_ensemble,
                del_cloud_liquid_cloud_ensemble=locals.del_cloud_liquid_cloud_ensemble,
                del_buoyancy_cloud_ensemble=locals.del_buoyancy_cloud_ensemble,
                del_convective_ice_cloud_ensemble=locals.del_convective_ice_cloud_ensemble,
                del_large_scale_ice_cloud_ensemble=locals.del_large_scale_ice_cloud_ensemble,
                del_convective_liquid_cloud_ensemble=locals.del_convective_liquid_cloud_ensemble,
                del_large_scale_liquid_cloud_ensemble=locals.del_large_scale_liquid_cloud_ensemble,
                del_convective_cloud_fraction_cloud_ensemble=locals.del_convective_cloud_fraction_cloud_ensemble,
                del_large_scale_cloud_fraction_cloud_ensemble=locals.del_large_scale_cloud_fraction_cloud_ensemble,
                dtdt=state.output.dtdt,
                dvapordt=state.output.dvapordt,
                dcloudicedt=state.output.dcloudicedt,
                dudt=state.output.dudt,
                dvdt=state.output.dvdt,
                dbuoyancydt=state.output.dbuoyancydt,
                dconvectiveicedt=state.output.dconvectiveicedt,
                dlargescaleicedt=state.output.dlargescaleicedt,
                dconvectiveliquiddt=state.output.dconvectiveliquiddt,
                dlargescaleliquiddt=state.output.dlargescaleliquiddt,
                dconvectivecloudfractiondt=state.output.dconvectivecloudfractiondt,
                dlargescalecloudfractiondt=state.output.dlargescalecloudfractiondt,
                mass_flux_ensemble=locals.mass_flux_ensemble,
                precipitation_ensemble=locals.precipitation_ensemble,
                xff_mid=locals.xff_mid,
                CLOSURE_CHOICE=plume_dependent_constants.CLOSURE_CHOICE,
                CLOUD_BASE_MASS_FLUX_FACTOR=plume_dependent_constants.CLOUD_BASE_MASS_FLUX_FACTOR,
                plume=plume_dependent_constants.PLUME_INDEX,
            )

        outputs = {
            "error_code": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "local_error_code_2": locals.error_code_2.field[:],
            "local_error_code_3": locals.error_code_3.field[:],
            "cloud_top_level": state.output.cloud_top_level.field[:, :, plume_dependent_constants.PLUME_INDEX]
            + 1,
            "updraft_lfc_level": state.output.updraft_lfc_level.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ]
            + 1,
            "p_cloud_levels_forced": state.output.p_cloud_levels_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "normalized_massflux_updraft_forced": state.output.normalized_massflux_updraft_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "precip": state.output.precip.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "local_effective_condensate_to_fall_forced": locals.effective_condensate_to_fall_forced.field[:],
            "cloud_base_mass_flux_modified": state.output.cloud_base_mass_flux_modified.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "scale_dependence_factor": state.output.scale_dependence_factor.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_ocean_fraction": locals.ocean_fraction.field[:],
            "local_f_dicycle_modified": locals.f_dicycle_modified.field[:],
            "local_del_u_cloud_ensemble": locals.del_u_cloud_ensemble.field[:],
            "local_del_v_cloud_ensemble": locals.del_v_cloud_ensemble.field[:],
            "local_del_t_cloud_ensemble": locals.del_t_cloud_ensemble.field[:],
            "local_del_vapor_cloud_ensemble": locals.del_vapor_cloud_ensemble.field[:],
            "local_del_cloud_liquid_cloud_ensemble": locals.del_cloud_liquid_cloud_ensemble.field[:],
            "local_del_buoyancy_cloud_ensemble": locals.del_buoyancy_cloud_ensemble.field[:],
            "local_del_convective_ice_cloud_ensemble": locals.del_convective_ice_cloud_ensemble.field[:],
            "local_del_large_scale_ice_cloud_ensemble": locals.del_large_scale_ice_cloud_ensemble.field[:],
            "local_del_convective_liquid_cloud_ensemble": locals.del_convective_liquid_cloud_ensemble.field[
                :
            ],
            "local_del_large_scale_liquid_cloud_ensemble": locals.del_large_scale_liquid_cloud_ensemble.field[
                :
            ],
            "local_del_convective_cloud_fraction_cloud_ensemble": locals.del_convective_cloud_fraction_cloud_ensemble.field[
                :
            ],
            "local_del_large_scale_cloud_fraction_cloud_ensemble": locals.del_large_scale_cloud_fraction_cloud_ensemble.field[
                :
            ],
            "dtdt": state.output.dtdt.field[:, :, :, plume_dependent_constants.PLUME_INDEX],
            "dvapordt": state.output.dvapordt.field[:, :, :, plume_dependent_constants.PLUME_INDEX],
            "dcloudicedt": state.output.dcloudicedt.field[:, :, :, plume_dependent_constants.PLUME_INDEX],
            "dudt": state.output.dudt.field[:, :, :, plume_dependent_constants.PLUME_INDEX],
            "dvdt": state.output.dvdt.field[:, :, :, plume_dependent_constants.PLUME_INDEX],
            "dbuoyancydt": state.output.dbuoyancydt.field[:, :, :, plume_dependent_constants.PLUME_INDEX],
            "dconvectiveicedt": state.output.dconvectiveicedt.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "dlargescaleicedt": state.output.dlargescaleicedt.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "dconvectiveliquiddt": state.output.dconvectiveliquiddt.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "dlargescaleliquiddt": state.output.dlargescaleliquiddt.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "dconvectivecloudfractiondt": state.output.dconvectivecloudfractiondt.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "dlargescalecloudfractiondt": state.output.dlargescalecloudfractiondt.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_mass_flux_ensemble": locals.mass_flux_ensemble.field[:],
            "local_precipitation_ensemble": locals.precipitation_ensemble.field[:],
            "local_xff_mid": locals.xff_mid.field[:],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_EnsembleOutputAndFeedback_shallow(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_EnsembleOutputAndFeedback_mid(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_EnsembleOutputAndFeedback_deep(TranslateFortranData2Py):
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
