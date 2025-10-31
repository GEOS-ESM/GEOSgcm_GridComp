from f90nml import Namelist
from ndsl import StencilFactory
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.convection.GF_2020.state import GF2020State
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.cumulus_parameterization import (
    CumulusParameterization,
)
from pyMoist.convection.GF_2020.locals import GF2020Locals
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable


class TranslateGF2020_CumulusParameterization(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, namelist, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        self.in_vars["data_vars"] = {
            "t_excess_cu_param_input": {},
            "vapor_excess_cu_param_input": {},
            "ocean_fraction": {},
            "t_cu_param_input": {},
            "vapor_timestep_start_cu_param_input": {},
            "grid_scale_forcing_t_cu_param_input_data_buffered": {},
            "grid_scale_forcing_vapor_cu_param_input_data_buffered": {},
            "subgrid_scale_forcing_t_cu_param_input_data_buffered": {},
            "subgrid_scale_forcing_vapor_cu_param_input_data_buffered": {},
        }

        self.out_vars = self.in_vars["data_vars"].copy()
        del (
            self.out_vars["t_excess_cu_param_input"],
            self.out_vars["vapor_excess_cu_param_input"],
            self.out_vars["ocean_fraction"],
            self.out_vars["t_cu_param_input"],
            self.out_vars["vapor_timestep_start_cu_param_input"],
            self.out_vars["grid_scale_forcing_t_cu_param_input_data_buffered"],
            self.out_vars["grid_scale_forcing_vapor_cu_param_input_data_buffered"],
            self.out_vars["subgrid_scale_forcing_t_cu_param_input_data_buffered"],
            self.out_vars["subgrid_scale_forcing_vapor_cu_param_input_data_buffered"],
        )
        self.out_vars.update(
            {
                "t_excess_cu_param_internal": {},
                "vapor_excess_cu_param_internal": {},
                "t_cu_param_internal": {},
                "vapor_cu_param_internal": {},
                "moist_static_energy_cu_param_internal": {},
                "t_pbl_cu_param_internal": {},
                "vapor_pbl_cu_param_internal": {},
            }
        )

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF2020-constants")
        self.cu_param_constants = data_loader.load("GF2020_CumulusParameterization-constants")

    def compute(self, inputs):
        config = GF2020Config(SINGLE_COLUMN_MODE=False, **self.constants)
        cumulus_parameterization_config = GF2020CumulusParameterizationConfig(**self.cu_param_constants)

        state = GF2020State.zeros(self.quantity_factory)

        locals = GF2020Locals.zeros(
            self.quantity_factory,
            data_dimensions={
                "plumes": 3,
            },
        )

        import numpy as np

        locals.cumulus_parameterization_input.t_excess.field[:] = inputs.pop("t_excess_cu_param_input")
        locals.cumulus_parameterization_input.vapor_excess.field[:] = inputs.pop(
            "vapor_excess_cu_param_input"
        )
        locals.cumulus_parameterization_input.ocean_fraction.field[:] = inputs.pop("ocean_fraction")
        locals.cumulus_parameterization_input.t.field[:] = inputs.pop("t_cu_param_input")
        locals.cumulus_parameterization_input.vapor_timestep_start.field[:] = inputs.pop(
            "vapor_timestep_start_cu_param_input"
        )
        locals.cumulus_parameterization_input.grid_scale_forcing_t.field[:] = inputs.pop(
            "grid_scale_forcing_t_cu_param_input_data_buffered"
        )
        locals.cumulus_parameterization_input.grid_scale_forcing_vapor.field[:] = inputs.pop(
            "grid_scale_forcing_vapor_cu_param_input_data_buffered"
        )
        locals.cumulus_parameterization_input.subgrid_scale_forcing_t.field[:] = inputs.pop(
            "subgrid_scale_forcing_t_cu_param_input_data_buffered"
        )
        locals.cumulus_parameterization_input.subgrid_scale_forcing_vapor.field[:] = inputs.pop(
            "subgrid_scale_forcing_vapor_cu_param_input_data_buffered"
        )

        saturation_tables = SaturationVaporPressureTable(self.stencil_factory.backend)

        cumulus_parameterization = CumulusParameterization(
            self.stencil_factory, self.quantity_factory, config, cumulus_parameterization_config
        )

        cumulus_parameterization(state, locals, saturation_tables)

        # top rows are not computed in Fortran, retains initalized value (nan)
        # Python initalizes with zero, fill top row with nan for test passage
        locals.cumulus_parameterization_internal.t.field[:, :, -1] = np.nan
        locals.cumulus_parameterization_internal.vapor.field[:, :, -1] = np.nan
        locals.cumulus_parameterization_internal.t_pbl.field[:, :, -1] = np.nan
        locals.cumulus_parameterization_internal.vapor_pbl.field[:, :, -1] = np.nan
        locals.cumulus_parameterization_internal.moist_static_energy.field[:, :, -1] = np.nan

        inputs.update(
            {
                "t_excess_cu_param_internal": locals.cumulus_parameterization_internal.t_excess.field[:],
                "vapor_excess_cu_param_internal": locals.cumulus_parameterization_internal.vapor_excess.field[
                    :
                ],
                "t_cu_param_internal": locals.cumulus_parameterization_internal.t.field[:],
                "vapor_cu_param_internal": locals.cumulus_parameterization_internal.vapor.field[:],
                "t_pbl_cu_param_internal": locals.cumulus_parameterization_internal.t_pbl.field[:],
                "vapor_pbl_cu_param_internal": locals.cumulus_parameterization_internal.vapor_pbl.field[:],
                "moist_static_energy_cu_param_internal": locals.cumulus_parameterization_internal.moist_static_energy.field[
                    :
                ],
            }
        )

        return inputs
