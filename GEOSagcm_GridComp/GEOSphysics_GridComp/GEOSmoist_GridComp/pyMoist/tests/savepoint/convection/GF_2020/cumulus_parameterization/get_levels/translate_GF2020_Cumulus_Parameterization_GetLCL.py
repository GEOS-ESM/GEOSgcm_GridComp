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
    GetLCL,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.setup.set_constants import set_constants


class TranslateGF2020_CumulusParameterization_GetLCL_deep(TranslateFortranData2Py):
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
            "p_forced_env_clev": {},
            "local_p_cloud_levels": {},
            "local_t_excess": {},
            "local_t_cloud_levels_env_clev": {},
            "t_perturbation": {},
            "local_vapor_excess": {},
            "local_vapor_cloud_levels_forced": {},
            "omega": {},
            "air_density": {},
            "local_geopotential_height_cloud_levels": {},
            "topography_height_no_negative": {},
            "ocean_fraction": {},
            "local_updraft_origin_level": {},
            "grid_length": {},
            "lcl_level": {},
            "error_code": {},
            "vapor_source": {},
            "t_source": {},
            "p_source": {},
            "z_source": {},
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
        state.input_output.p_forced.data[:] = inputs["p_forced_env_clev"]
        locals.p_cloud_levels.data[:] = inputs["local_p_cloud_levels"]
        locals.t_excess.data[:] = inputs["local_t_excess"]
        locals.t_cloud_levels.data[:] = inputs["local_t_cloud_levels_env_clev"]
        state.output.t_perturbation.data[:] = inputs["t_perturbation"]
        locals.vapor_excess.data[:] = inputs["local_vapor_excess"]
        locals.vapor_cloud_levels.data[:] = inputs["local_vapor_cloud_levels_forced"]
        state.input_output.omega.data[:] = inputs["omega"]
        state.input_output.air_density.data[:] = inputs["air_density"]
        locals.geopotential_height_cloud_levels.data[:] = inputs["local_geopotential_height_cloud_levels"]
        state.input_output.topography_height_no_negative.data[:] = inputs["topography_height_no_negative"]
        state.input.ocean_fraction.data[:] = inputs["ocean_fraction"]
        locals.updraft_origin_level.data[:] = inputs["local_updraft_origin_level"] - 1
        state.input_output.grid_length.data[:] = inputs["grid_length"]
        state.output.lcl_level.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["lcl_level"]
        state.output.error_code.data[:, :, plume_dependent_constants.PLUME_INDEX] = inputs["error_code"]

        # initalize test code
        code = GetLCL(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        # call test code
        from ndsl.constants import X_DIM, Y_DIM

        vapor_source = self.quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        t_source = self.quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        p_source = self.quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        z_source = self.quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        if plume_dependent_constants.ENABLE_PLUME == 1:
            code(
                state=state,
                locals=locals,
                plume_dependent_constants=plume_dependent_constants,
                vapor_source=vapor_source,
                t_source=t_source,
                p_source=p_source,
                z_source=z_source,
            )

            locals.updraft_origin_level.field[:] += 1
            state.output.lcl_level.field[:, :, plume_dependent_constants.PLUME_INDEX] += 1

        # write output
        outputs = {
            "p_forced_env_clev": state.input_output.p_forced.field[:],
            "local_p_cloud_levels": locals.p_cloud_levels.field[:],
            "local_t_excess": locals.t_excess.field[:],
            "local_t_cloud_levels_env_clev": locals.t_cloud_levels.field[:],
            "t_perturbation": state.output.t_perturbation.field[:],
            "local_vapor_excess": locals.vapor_excess.field[:],
            "local_vapor_cloud_levels_forced": locals.vapor_cloud_levels.field[:],
            "omega": state.input_output.omega.field[:],
            "air_density": state.input_output.air_density.field[:],
            "local_geopotential_height_cloud_levels": locals.geopotential_height_cloud_levels.field[:],
            "topography_height_no_negative": state.input_output.topography_height_no_negative.field[:],
            "ocean_fraction": state.input.ocean_fraction.field[:],
            "local_updraft_origin_level": locals.updraft_origin_level.field[:],
            "grid_length": state.input_output.grid_length.field[:],
            "lcl_level": state.output.lcl_level.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "error_code": state.output.error_code.field[:, :, plume_dependent_constants.PLUME_INDEX],
            "vapor_source": vapor_source.field[:],
            "t_source": t_source.field[:],
            "p_source": p_source.field[:],
            "z_source": z_source.field[:],
        }

        return outputs
