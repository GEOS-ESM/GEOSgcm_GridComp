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
from pyMoist.convection.GF_2020.cumulus_parameterization.moist_static_energy.moist_static_energy import (
    ParcelMoistStaticEnergy,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.setup.set_constants import set_constants


class TranslateGF2020_CumulusParameterization_ParcelMoistStaticEnergy_deep(TranslateFortranData2Py):
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
            "error_code_parcelmse": {},
            "ocean_fraction_parcelmse": {},
            "local_vapor_excess_parcelmse": {},
            "local_t_excess_parcelmse": {},
            "local_add_buoy_parcelmse": {},
            "p_forced_parcelmse": {},
            "local_env_moist_static_energy_cloud_levels_parcelmse": {},
            "local_env_moist_static_energy_cloud_levels_forced_parcelmse": {},
            "local_moist_static_energy_origin_level_parcelmse": {},
            "local_moist_static_energy_origin_level_forced_parcelmse": {},
            "local_updraft_origin_level_parcelmse": {},
            "t_perturbation_parcelmse": {},
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
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "error_code_parcelmse"
        ]
        state.input.ocean_fraction.data[:] = inputs["ocean_fraction_parcelmse"]
        locals.vapor_excess.data[:] = inputs["local_vapor_excess_parcelmse"]
        locals.t_excess.data[:] = inputs["local_t_excess_parcelmse"]
        locals.add_buoyancy.data[:] = inputs["local_add_buoy_parcelmse"]
        state.input_output.p_forced.data[:] = inputs["p_forced_parcelmse"]
        locals.environment_moist_static_energy_cloud_levels.data[:] = inputs[
            "local_env_moist_static_energy_cloud_levels_parcelmse"
        ]
        locals.environment_moist_static_energy_cloud_levels_forced.data[:] = inputs[
            "local_env_moist_static_energy_cloud_levels_forced_parcelmse"
        ]
        locals.moist_static_energy_origin_level.data[:] = inputs[
            "local_moist_static_energy_origin_level_parcelmse"
        ]
        locals.moist_static_energy_origin_level_forced.data[:] = inputs[
            "local_moist_static_energy_origin_level_forced_parcelmse"
        ]
        locals.updraft_origin_level.data[:] = inputs["local_updraft_origin_level_parcelmse"]
        state.output.t_perturbation.data[:] = inputs["t_perturbation_parcelmse"]

        # initalize test code
        code = ParcelMoistStaticEnergy(
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

        # write output
        outputs = {
            "error_code_parcelmse": state.output.error_code.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "ocean_fraction_parcelmse": state.input.ocean_fraction.field[:],
            "local_vapor_excess_parcelmse": locals.vapor_excess.field[:],
            "local_t_excess_parcelmse": locals.t_excess.field[:],
            "local_add_buoy_parcelmse": locals.add_buoyancy.field[:],
            "p_forced_parcelmse": state.input_output.p_forced.field[:],
            "local_env_moist_static_energy_cloud_levels_parcelmse": locals.environment_moist_static_energy_cloud_levels.field[
                :
            ],
            "local_env_moist_static_energy_cloud_levels_forced_parcelmse": locals.environment_moist_static_energy_cloud_levels_forced.field[
                :
            ],
            "local_moist_static_energy_origin_level_parcelmse": locals.moist_static_energy_origin_level.field[
                :
            ],
            "local_moist_static_energy_origin_level_forced_parcelmse": locals.moist_static_energy_origin_level_forced.field[
                :
            ],
            "local_updraft_origin_level_parcelmse": locals.updraft_origin_level.field[:],
            "t_perturbation_parcelmse": state.output.t_perturbation.field[:],
        }

        return outputs
