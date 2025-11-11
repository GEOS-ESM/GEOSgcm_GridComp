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
from pyMoist.convection.GF_2020.cumulus_parameterization.entrainment.entrainment import (
    DowndraftEntrainmentProfiles,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.setup.set_constants import set_constants


class TranslateGF2020_CumulusParameterization_DowndraftEntrainmentProfiles_shallow(TranslateFortranData2Py):
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
            "lateral_entrainment_rate_downentrainmentprofile": {},
            "local_entrainment_rate_downdraft_downentrainmntprofile": {},
            "local_detrainment_function_downdraft_downentrainmntprofile": {},
            "local_scale_dependence_factor_downdraft_downentrainmentprofile": {},
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
            cumulus_parameterization_config, plume_dependent_constants, "shallow"
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
        state.input.lateral_entrainment_rate.data[:] = inputs[
            "lateral_entrainment_rate_downentrainmentprofile"
        ]
        locals.entrainment_rate_downdraft.data[:] = inputs[
            "local_entrainment_rate_downdraft_downentrainmntprofile"
        ]
        locals.detrainment_function_downdraft.data[:] = inputs[
            "local_detrainment_function_downdraft_downentrainmntprofile"
        ]
        locals.scale_dependence_factor_downdraft.data[:] = inputs[
            "local_scale_dependence_factor_downdraft_downentrainmentprofile"
        ]

        # initalize test code
        code = DowndraftEntrainmentProfiles(
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
            "lateral_entrainment_rate_downentrainmentprofile": state.input.lateral_entrainment_rate.field[:],
            "local_entrainment_rate_downdraft_downentrainmntprofile": locals.entrainment_rate_downdraft.field[
                :
            ],
            "local_detrainment_function_downdraft_downentrainmntprofile": locals.detrainment_function_downdraft.field[
                :
            ],
            "local_scale_dependence_factor_downdraft_downentrainmentprofile": locals.scale_dependence_factor_downdraft.field[
                :
            ],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_DowndraftEntrainmentProfiles_mid(TranslateFortranData2Py):
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
            "lateral_entrainment_rate_downentrainmentprofile": {},
            "local_entrainment_rate_downdraft_downentrainmntprofile": {},
            "local_detrainment_function_downdraft_downentrainmntprofile": {},
            "local_scale_dependence_factor_downdraft_downentrainmentprofile": {},
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
            cumulus_parameterization_config, plume_dependent_constants, "mid"
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
        state.input.lateral_entrainment_rate.data[:] = inputs[
            "lateral_entrainment_rate_downentrainmentprofile"
        ]
        locals.entrainment_rate_downdraft.data[:] = inputs[
            "local_entrainment_rate_downdraft_downentrainmntprofile"
        ]
        locals.detrainment_function_downdraft.data[:] = inputs[
            "local_detrainment_function_downdraft_downentrainmntprofile"
        ]
        locals.scale_dependence_factor_downdraft.data[:] = inputs[
            "local_scale_dependence_factor_downdraft_downentrainmentprofile"
        ]

        # initalize test code
        code = DowndraftEntrainmentProfiles(
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
            "lateral_entrainment_rate_downentrainmentprofile": state.input.lateral_entrainment_rate.field[:],
            "local_entrainment_rate_downdraft_downentrainmntprofile": locals.entrainment_rate_downdraft.field[
                :
            ],
            "local_detrainment_function_downdraft_downentrainmntprofile": locals.detrainment_function_downdraft.field[
                :
            ],
            "local_scale_dependence_factor_downdraft_downentrainmentprofile": locals.scale_dependence_factor_downdraft.field[
                :
            ],
        }

        return outputs


class TranslateGF2020_CumulusParameterization_DowndraftEntrainmentProfiles_deep(TranslateFortranData2Py):
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
            "lateral_entrainment_rate_downentrainmentprofile": {},
            "local_entrainment_rate_downdraft_downentrainmntprofile": {},
            "local_detrainment_function_downdraft_downentrainmntprofile": {},
            "local_scale_dependence_factor_downdraft_downentrainmentprofile": {},
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
        state.input.lateral_entrainment_rate.data[:] = inputs[
            "lateral_entrainment_rate_downentrainmentprofile"
        ]
        locals.entrainment_rate_downdraft.data[:] = inputs[
            "local_entrainment_rate_downdraft_downentrainmntprofile"
        ]
        locals.detrainment_function_downdraft.data[:] = inputs[
            "local_detrainment_function_downdraft_downentrainmntprofile"
        ]
        locals.scale_dependence_factor_downdraft.data[:] = inputs[
            "local_scale_dependence_factor_downdraft_downentrainmentprofile"
        ]

        # initalize test code
        code = DowndraftEntrainmentProfiles(
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
            "lateral_entrainment_rate_downentrainmentprofile": state.input.lateral_entrainment_rate.field[:],
            "local_entrainment_rate_downdraft_downentrainmntprofile": locals.entrainment_rate_downdraft.field[
                :
            ],
            "local_detrainment_function_downdraft_downentrainmntprofile": locals.detrainment_function_downdraft.field[
                :
            ],
            "local_scale_dependence_factor_downdraft_downentrainmentprofile": locals.scale_dependence_factor_downdraft.field[
                :
            ],
        }

        return outputs
