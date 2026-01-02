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
from pyMoist.convection.GF_2020.cumulus_parameterization.entrainment import (
    compute_lateral_massflux,
    compute_uc_vc,
)
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
            "error_code_massentdet": {},
            "cloud_top_level_massentdet": {},
            "local_geopotential_height_cloud_levels_forced_massentdet": {},
            "normalized_massflux_updraft_forced_massentdet": {},
            "local_detrainment_function_updraft_massentdet": {},
            "entrainment_rate_massentdet": {},
            "p_cloud_levels_forced_massentdet": {},
            "mass_entrainment_updraft_forced_massentdet": {},
            "mass_detrainment_updraft_forced_massentdet": {},
            "local_mass_entrainment_updraft_massentdet": {},
            "local_mass_detrainment_updraft_massentdet": {},
            "updraft_lfc_level_massentdet": {},
            "updraft_origin_level_massentdet": {},
            "pbl_level_massentdet": {},
            "local_mass_entrainment_u_updraft_massentdet": {},
            "local_mass_detrainment_u_updraft_massentdet": {},
            "local_start_level_massentdet": {},
            "ocean_fraction_massentdet": {},
            "p_forced_massentdet": {},
            "local_u_cloud_levels_massentdet": {},
            "local_v_cloud_levels_massentdet": {},
            "local_u_c_massentdet": {},
            "local_v_c_massentdet": {},
            "local_moist_static_energy_origin_level_massentdet": {},
            "local_moist_static_energy_origin_level_forced_massentdet": {},
            "local_cloud_moist_static_energy_massentdet": {},
            "local_cloud_moist_static_energy_forced_massentdet": {},
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
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "error_code_massentdet"
        ]
        state.output.cloud_top_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "cloud_top_level_massentdet"
        ]
        locals.geopotential_height_cloud_levels_forced.data[:] = inputs[
            "local_geopotential_height_cloud_levels_forced_massentdet"
        ]
        state.output.normalized_massflux_updraft_forced.data[
            :, :, :, plume_dependent_constants.PLUME_INDEX
        ] = inputs["normalized_massflux_updraft_forced_massentdet"]
        locals.detrainment_function_updraft.data[:] = inputs["local_detrainment_function_updraft_massentdet"]
        state.output.entrainment_rate.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "entrainment_rate_massentdet"
        ]
        state.output.p_cloud_levels_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "p_cloud_levels_forced_massentdet"
        ]
        state.output.mass_entrainment_updraft_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["mass_entrainment_updraft_forced_massentdet"]
        )
        state.output.mass_detrainment_updraft_forced.data[:, :, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["mass_detrainment_updraft_forced_massentdet"]
        )
        locals.mass_entrainment_updraft.data[:] = inputs["local_mass_entrainment_updraft_massentdet"]
        locals.mass_detrainment_updraft.data[:] = inputs["local_mass_detrainment_updraft_massentdet"]
        state.output.updraft_lfc_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs[
            "updraft_lfc_level_massentdet"
        ]
        state.output.updraft_origin_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = (
            inputs["updraft_origin_level_massentdet"] - 1
        )
        state.input_output.pbl_level.data[:] = inputs["pbl_level_massentdet"]
        locals.mass_entrainment_u_updraft.data[:] = inputs["local_mass_entrainment_u_updraft_massentdet"]
        locals.mass_detrainment_u_updraft.data[:] = inputs["local_mass_detrainment_u_updraft_massentdet"]
        locals.start_level.data[:] = inputs["local_start_level_massentdet"] - 1
        state.input.ocean_fraction.data[:] = inputs["ocean_fraction_massentdet"]
        state.input_output.p_forced.data[:] = inputs["p_forced_massentdet"]
        locals.u_cloud_levels.data[:] = inputs["local_u_cloud_levels_massentdet"]
        locals.v_cloud_levels.data[:] = inputs["local_v_cloud_levels_massentdet"]
        locals.u_c.data[:] = inputs["local_u_c_massentdet"]
        locals.v_c.data[:] = inputs["local_v_c_massentdet"]
        locals.moist_static_energy_origin_level.data[:] = inputs[
            "local_moist_static_energy_origin_level_massentdet"
        ]
        locals.moist_static_energy_origin_level_forced.data[:] = inputs[
            "local_moist_static_energy_origin_level_forced_massentdet"
        ]
        locals.cloud_moist_static_energy.data[:] = inputs["local_cloud_moist_static_energy_massentdet"]
        locals.cloud_moist_static_energy_forced.data[:] = inputs[
            "local_cloud_moist_static_energy_forced_massentdet"
        ]

        code_part_1 = self.stencil_factory.from_dims_halo(
            func=compute_lateral_massflux,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        code_part_2 = self.stencil_factory.from_dims_halo(
            func=compute_uc_vc,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "BOUNDARY_CONDITION_METHOD": cumulus_parameterization_config.BOUNDARY_CONDITION_METHOD
            },
        )

        if plume_dependent_constants.ENABLE_PLUME == 1:
            code_part_1(
                error_code=state.output.error_code,
                cloud_top_level=state.output.cloud_top_level,
                geopotential_height=locals.geopotential_height_cloud_levels_forced,
                normalized_massflux_updraft=state.output.normalized_massflux_updraft_forced,
                detrainment_function_updraft=locals.detrainment_function_updraft,
                entrainment_rate=state.output.entrainment_rate,
                p_cloud_levels_forced=state.output.p_cloud_levels_forced,
                mass_entrainment_updraft_forced=state.output.mass_entrainment_updraft_forced,
                mass_detrainment_updraft_forced=state.output.mass_detrainment_updraft_forced,
                mass_entrainment_updraft=locals.mass_entrainment_updraft,
                mass_detrainment_updraft=locals.mass_detrainment_updraft,
                updraft_lfc_level=state.output.updraft_lfc_level,
                updraft_origin_level=state.output.updraft_origin_level,
                pbl_level=state.input_output.pbl_level,
                mass_entrainment_u_updraft=locals.mass_entrainment_u_updraft,
                mass_detrainment_u_updraft=locals.mass_detrainment_u_updraft,
                LAMBDA_DEEP=plume_dependent_constants.LAMBDA_DEEP,
                plume=plume_dependent_constants.PLUME_INDEX,
            )

            code_part_2(
                u_c=locals.u_c,
                v_c=locals.v_c,
                cloud_moist_static_energy=locals.cloud_moist_static_energy,
                cloud_moist_static_energy_forced=locals.cloud_moist_static_energy_forced,
                error_code=state.output.error_code,
                start_level=locals.start_level,
                moist_static_energy_origin_level=locals.moist_static_energy_origin_level,
                moist_static_energy_origin_level_forced=locals.moist_static_energy_origin_level_forced,
                u_cloud_levels=locals.u_cloud_levels,
                v_cloud_levels=locals.v_cloud_levels,
                p=state.input_output.p_forced,
                updraft_origin_level=state.output.updraft_origin_level,
                ocean_fraction=state.input.ocean_fraction,
                AVERAGE_LAYER_DEPTH=plume_dependent_constants.AVERAGE_LAYER_DEPTH,
                plume=plume_dependent_constants.PLUME_INDEX,
            )

        outputs = {
            "error_code_massentdet": state.output.error_code.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "cloud_top_level_massentdet": state.output.cloud_top_level.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_geopotential_height_cloud_levels_forced_massentdet": locals.geopotential_height_cloud_levels_forced.field[
                :
            ],
            "normalized_massflux_updraft_forced_massentdet": state.output.normalized_massflux_updraft_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_detrainment_function_updraft_massentdet": locals.detrainment_function_updraft.field[:],
            "entrainment_rate_massentdet": state.output.entrainment_rate.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "p_cloud_levels_forced_massentdet": state.output.p_cloud_levels_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "mass_entrainment_updraft_forced_massentdet": state.output.mass_entrainment_updraft_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "mass_detrainment_updraft_forced_massentdet": state.output.mass_detrainment_updraft_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "local_mass_entrainment_updraft_massentdet": locals.mass_entrainment_updraft.field[:],
            "local_mass_detrainment_updraft_massentdet": locals.mass_detrainment_updraft.field[:],
            "updraft_lfc_level_massentdet": state.output.updraft_lfc_level.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ],
            "updraft_origin_level_massentdet": state.output.updraft_origin_level.field[
                :, :, plume_dependent_constants.PLUME_INDEX
            ]
            + 1,
            "pbl_level_massentdet": state.input_output.pbl_level.field[:],
            "local_mass_entrainment_u_updraft_massentdet": locals.mass_entrainment_u_updraft.field[:],
            "local_mass_detrainment_u_updraft_massentdet": locals.mass_detrainment_u_updraft.field[:],
            "local_start_level_massentdet": locals.start_level.field[:] + 1,
            "ocean_fraction_massentdet": state.input.ocean_fraction.field[:],
            "p_forced_massentdet": state.input_output.p_forced.field[:],
            "local_u_cloud_levels_massentdet": locals.u_cloud_levels.field[:],
            "local_v_cloud_levels_massentdet": locals.v_cloud_levels.field[:],
            "local_u_c_massentdet": locals.u_c.field[:],
            "local_v_c_massentdet": locals.v_c.field[:],
            "local_moist_static_energy_origin_level_massentdet": locals.moist_static_energy_origin_level.field[
                :
            ],
            "local_moist_static_energy_origin_level_forced_massentdet": locals.moist_static_energy_origin_level_forced.field[
                :
            ],
            "local_cloud_moist_static_energy_massentdet": locals.cloud_moist_static_energy.field[:],
            "local_cloud_moist_static_energy_forced_massentdet": locals.cloud_moist_static_energy_forced.field[
                :
            ],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_CalculateMassEntrainmentDetrainment_shallow(
    TranslateFortranData2Py
):
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


class TranslateGF2020_CumulusParameterization_CalculateMassEntrainmentDetrainment_mid(
    TranslateFortranData2Py
):
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


class TranslateGF2020_CumulusParameterization_CalculateMassEntrainmentDetrainment_deep(
    TranslateFortranData2Py
):
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
