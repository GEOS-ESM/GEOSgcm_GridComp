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
from pyMoist.convection.GF_2020.cumulus_parameterization.prepare_output import prepare_output

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
            "error_code_po": {},
            "cloud_base_mass_flux_modified_po": {},
            "total_normalized_integrated_condensate_forced_po": {},
            "total_normalized_integrated_evaporate_forced_po": {},
            "normalized_massflux_updraft_forced": {},
            "normalized_massflux_downdraft_forced": {},
            "condensate_to_fall_forced_po": {},
            "evaporate_in_downdraft_forced_po": {},
            "mass_entrainment_updraft_forced": {},
            "mass_detrainment_updraft_forced": {},
            "mass_entrainment_downdraft_forced_po": {},
            "mass_detrainment_downdraft_forced_po": {},
            "local_environment_massflux_po": {},
            "local_vapor_tendency_from_environmental_subsidence": {},
            "local_moist_static_energy_tendency_from_environmental_subsidence": {},
            "local_t_tendency_from_environmental_subsidence": {},
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
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["error_code_po"]
        state.output.cloud_base_mass_flux_modified.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "cloud_base_mass_flux_modified_po"
        ]
        state.output.total_normalized_integrated_condensate_forced.data[
            :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["total_normalized_integrated_condensate_forced_po"]
        state.output.total_normalized_integrated_evaporate_forced.data[
            :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["total_normalized_integrated_evaporate_forced_po"]
        state.output.normalized_massflux_updraft_forced.data[
            :, :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["normalized_massflux_updraft_forced"]
        state.output.normalized_massflux_downdraft_forced.data[
            :, :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["normalized_massflux_downdraft_forced"]
        state.output.condensate_to_fall_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "condensate_to_fall_forced_po"
        ]
        state.output.evaporate_in_downdraft_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["evaporate_in_downdraft_forced_po"]
        )
        state.output.mass_entrainment_updraft_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["mass_entrainment_updraft_forced"]
        )
        state.output.mass_detrainment_updraft_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["mass_detrainment_updraft_forced"]
        )
        state.output.mass_detrainment_downdraft_forced.data[
            :, :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["mass_detrainment_downdraft_forced_po"]
        state.output.mass_entrainment_downdraft_forced.data[
            :, :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["mass_entrainment_downdraft_forced_po"]
        locals.environment_massflux.data[:] = inputs["local_environment_massflux_po"]
        locals.vapor_tendency_from_environmental_subsidence.data[:] = inputs[
            "local_vapor_tendency_from_environmental_subsidence"
        ]
        locals.moist_static_energy_tendency_from_environmental_subsidence.data[:] = inputs[
            "local_moist_static_energy_tendency_from_environmental_subsidence"
        ]
        locals.t_tendency_from_environmental_subsidence.data[:] = inputs[
            "local_t_tendency_from_environmental_subsidence"
        ]

        code_part_1 = self.stencil_factory.from_dims_halo(
            func=prepare_output,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        if plume_dependent_constants.ENABLE_PLUME == 1:
            code_part_1(
                error_code=state.output.error_code,
                plume=plume_dependent_constants.PLUME_INDEX,
                cloud_base_mass_flux_modified=state.output.cloud_base_mass_flux_modified,
                total_normalized_integrated_condensate_forced=state.output.total_normalized_integrated_condensate_forced,
                total_normalized_integrated_evaporate_forced=state.output.total_normalized_integrated_evaporate_forced,
                normalized_massflux_updraft_forced=state.output.normalized_massflux_updraft_forced,
                normalized_massflux_downdraft_forced=state.output.normalized_massflux_downdraft_forced,
                condensate_to_fall_forced=state.output.condensate_to_fall_forced,
                evaporate_in_downdraft_forced=state.output.evaporate_in_downdraft_forced,
                mass_entrainment_updraft_forced=state.output.mass_entrainment_updraft_forced,
                mass_detrainment_updraft_forced=state.output.mass_detrainment_updraft_forced,
                mass_entrainment_downdraft_forced=state.output.mass_entrainment_downdraft_forced,
                mass_detrainment_downdraft_forced=state.output.mass_detrainment_downdraft_forced,
                environment_massflux=locals.environment_massflux,
                vapor_tendency_from_environmental_subsidence=locals.vapor_tendency_from_environmental_subsidence,
                moist_static_energy_tendency_from_environmental_subsidence=locals.moist_static_energy_tendency_from_environmental_subsidence,
                t_tendency_from_environmental_subsidence=locals.t_tendency_from_environmental_subsidence,
            )

        outputs = {
            "error_code_po": state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX],
            "cloud_base_mass_flux_modified_po": state.output.cloud_base_mass_flux_modified.data[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "total_normalized_integrated_condensate_forced_po": state.output.total_normalized_integrated_condensate_forced.data[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "total_normalized_integrated_evaporate_forced_po": state.output.total_normalized_integrated_evaporate_forced.data[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "normalized_massflux_updraft_forced": state.output.normalized_massflux_updraft_forced.data[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "normalized_massflux_downdraft_forced": state.output.normalized_massflux_downdraft_forced.data[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "condensate_to_fall_forced_po": state.output.condensate_to_fall_forced.data[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "evaporate_in_downdraft_forced_po": state.output.evaporate_in_downdraft_forced.data[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "mass_entrainment_updraft_forced": state.output.mass_entrainment_updraft_forced.data[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "mass_detrainment_updraft_forced": state.output.mass_detrainment_updraft_forced.data[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "mass_entrainment_downdraft_forced_po": state.output.mass_entrainment_downdraft_forced.data[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "mass_detrainment_downdraft_forced_po": state.output.mass_detrainment_downdraft_forced.data[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_environment_massflux_po": locals.environment_massflux.data[:],
            "local_vapor_tendency_from_environmental_subsidence": locals.vapor_tendency_from_environmental_subsidence.data[
                :
            ],
            "local_moist_static_energy_tendency_from_environmental_subsidence": locals.moist_static_energy_tendency_from_environmental_subsidence.data[
                :
            ],
            "local_t_tendency_from_environmental_subsidence": locals.t_tendency_from_environmental_subsidence.data[
                :
            ],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_PrepareOutput_shallow(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_PrepareOutput_mid(TranslateFortranData2Py):
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


class TranslateGF2020_CumulusParameterization_PrepareOutput_deep(TranslateFortranData2Py):
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
