from ndsl import StencilFactory, QuantityFactory, ndsl_log
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.state import GF2020State
from pyMoist.convection.GF_2020.locals import GF2020Locals
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.convection.GF_2020.cumulus_parameterization.setup.stencils import (
    set_plume_dependent_fields,
)
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)


class Setup:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GF2020Config,
        cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
    ):
        self.config = config
        self.cu_param_config = cumulus_parameterization_config

        self._set_plume_dependent_fields = stencil_factory.from_dims_halo(
            func=set_plume_dependent_fields,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"DT_MOIST": config.DT_MOIST},
        )

    def __call__(
        self,
        state: GF2020State,
        locals: GF2020Locals,
        saturation_tables: SaturationVaporPressureTable,
        plume: str,
        plume_depenedent_constants: GF2020PlumeDependentConstants,
    ):
        if plume == "shallow":
            plume_depenedent_constants.DOWNDRAFT_MAX_HEIGHT_LAND = (
                self.cu_param_config.DOWNDRAFT_MAX_HEIGHT_LAND_SHALLOW
            )
            plume_depenedent_constants.DOWNDRAFT_MAX_HEIGHT_OCEAN = (
                self.cu_param_config.DOWNDRAFT_MAX_HEIGHT_OCEAN_SHALLOW
            )
            plume_depenedent_constants.UPDRAFT_MAX_HEIGHT_LAND = (
                self.cu_param_config.UPDRAFT_MAX_HEIGHT_LAND_SHALLOW
            )
            plume_depenedent_constants.UPDRAFT_MAX_HEIGHT_OCEAN = (
                self.cu_param_config.UPDRAFT_MAX_HEIGHT_OCEAN_SHALLOW
            )
            plume_depenedent_constants.MINIMUM_EVAP_FRACTION_LAND = (
                self.cu_param_config.MINIMUM_EVAP_FRACTION_LAND_SHALLOW
            )
            plume_depenedent_constants.MINIMUM_EVAP_FRACTION_OCEAN = (
                self.cu_param_config.MINIMUM_EVAP_FRACTION_OCEAN_SHALLOW
            )
            plume_depenedent_constants.MAXIMUM_EVAP_FRACTION_LAND = (
                self.cu_param_config.MAXIMUM_EVAP_FRACTION_LAND_SHALLOW
            )
            plume_depenedent_constants.MAXIMUM_EVAP_FRACTION_OCEAN = (
                self.cu_param_config.MAXIMUM_EVAP_FRACTION_OCEAN_SHALLOW
            )
            plume_depenedent_constants.CLOUD_BASE_MASS_FLUX_FACTOR = (
                self.cu_param_config.CLOUD_BASE_MASS_FLUX_FACTOR_SHALLOW
            )
            plume_depenedent_constants.USE_EXCESS = self.cu_param_config.USE_EXCESS_SHALLOW
            plume_depenedent_constants.ENABLE_PLUME = self.cu_param_config.ENABLE_SHALLOW
            if self.cu_param_config.ZERO_DIFF == 1 or self.cu_param_config.MOIST_TRIGGER == 0:
                plume_depenedent_constants.CAP_MAX_INC = 25.0
            else:
                plume_depenedent_constants.CAP_MAX_INC = 10.0
            if self.cu_param_config.PRESSURE_GRADIENT_CONSTANT != 0.0:
                plume_depenedent_constants.LAMBDA_DEEP = 0.0
                plume_depenedent_constants.LAMBDA_DOWN = 0.0
            else:
                plume_depenedent_constants.LAMBDA_DEEP = self.cu_param_config.LAMBDA_DEEP
                plume_depenedent_constants.LAMBDA_DOWN = self.cu_param_config.LAMBDA_SHALLOW_DOWN
        if plume == "mid":
            plume_depenedent_constants.DOWNDRAFT_MAX_HEIGHT_LAND = (
                self.cu_param_config.DOWNDRAFT_MAX_HEIGHT_LAND_MID
            )
            plume_depenedent_constants.DOWNDRAFT_MAX_HEIGHT_OCEAN = (
                self.cu_param_config.DOWNDRAFT_MAX_HEIGHT_OCEAN_MID
            )
            plume_depenedent_constants.UPDRAFT_MAX_HEIGHT_LAND = (
                self.cu_param_config.UPDRAFT_MAX_HEIGHT_LAND_MID
            )
            plume_depenedent_constants.UPDRAFT_MAX_HEIGHT_OCEAN = (
                self.cu_param_config.UPDRAFT_MAX_HEIGHT_OCEAN_MID
            )
            plume_depenedent_constants.MINIMUM_EVAP_FRACTION_LAND = (
                self.cu_param_config.MINIMUM_EVAP_FRACTION_LAND_MID
            )
            plume_depenedent_constants.MINIMUM_EVAP_FRACTION_OCEAN = (
                self.cu_param_config.MINIMUM_EVAP_FRACTION_OCEAN_MID
            )
            plume_depenedent_constants.MAXIMUM_EVAP_FRACTION_LAND = (
                self.cu_param_config.MAXIMUM_EVAP_FRACTION_LAND_MID
            )
            plume_depenedent_constants.MAXIMUM_EVAP_FRACTION_OCEAN = (
                self.cu_param_config.MAXIMUM_EVAP_FRACTION_OCEAN_MID
            )
            plume_depenedent_constants.CLOUD_BASE_MASS_FLUX_FACTOR = (
                self.cu_param_config.CLOUD_BASE_MASS_FLUX_FACTOR_MID
            )
            plume_depenedent_constants.USE_EXCESS = self.cu_param_config.USE_EXCESS_MID
            plume_depenedent_constants.ENABLE_PLUME = self.cu_param_config.ENABLE_MID
            if self.cu_param_config.ZERO_DIFF == 1 or self.cu_param_config.MOIST_TRIGGER == 0:
                plume_depenedent_constants.CAP_MAX_INC = 10.0
            else:
                plume_depenedent_constants.CAP_MAX_INC = 90.0
            if self.cu_param_config.PRESSURE_GRADIENT_CONSTANT != 0.0:
                plume_depenedent_constants.LAMBDA_DEEP = 0.0
                plume_depenedent_constants.LAMBDA_DOWN = 0.0
            else:
                plume_depenedent_constants.LAMBDA_DEEP = self.cu_param_config.LAMBDA_SHALLOW_DOWN
                plume_depenedent_constants.LAMBDA_DOWN = self.cu_param_config.LAMBDA_SHALLOW_DOWN
        if plume == "deep":
            plume_depenedent_constants.DOWNDRAFT_MAX_HEIGHT_LAND = (
                self.cu_param_config.DOWNDRAFT_MAX_HEIGHT_LAND_DEEP
            )
            plume_depenedent_constants.DOWNDRAFT_MAX_HEIGHT_OCEAN = (
                self.cu_param_config.DOWNDRAFT_MAX_HEIGHT_OCEAN_DEEP
            )
            plume_depenedent_constants.UPDRAFT_MAX_HEIGHT_LAND = (
                self.cu_param_config.UPDRAFT_MAX_HEIGHT_LAND_DEEP
            )
            plume_depenedent_constants.UPDRAFT_MAX_HEIGHT_OCEAN = (
                self.cu_param_config.UPDRAFT_MAX_HEIGHT_OCEAN_DEEP
            )
            plume_depenedent_constants.MINIMUM_EVAP_FRACTION_LAND = (
                self.cu_param_config.MINIMUM_EVAP_FRACTION_LAND_DEEP
            )
            plume_depenedent_constants.MINIMUM_EVAP_FRACTION_OCEAN = (
                self.cu_param_config.MINIMUM_EVAP_FRACTION_OCEAN_DEEP
            )
            plume_depenedent_constants.MAXIMUM_EVAP_FRACTION_LAND = (
                self.cu_param_config.MAXIMUM_EVAP_FRACTION_LAND_DEEP
            )
            plume_depenedent_constants.MAXIMUM_EVAP_FRACTION_OCEAN = (
                self.cu_param_config.MAXIMUM_EVAP_FRACTION_OCEAN_DEEP
            )
            plume_depenedent_constants.CLOUD_BASE_MASS_FLUX_FACTOR = (
                self.cu_param_config.CLOUD_BASE_MASS_FLUX_FACTOR_DEEP
            )
            plume_depenedent_constants.USE_EXCESS = self.cu_param_config.USE_EXCESS_DEEP
            plume_depenedent_constants.ENABLE_PLUME = self.cu_param_config.ENABLE_DEEP
            if self.cu_param_config.ZERO_DIFF == 1 or self.cu_param_config.MOIST_TRIGGER == 0:
                plume_depenedent_constants.CAP_MAX_INC = 20.0
            else:
                plume_depenedent_constants.CAP_MAX_INC = 90.0
            if self.cu_param_config.PRESSURE_GRADIENT_CONSTANT != 0.0:
                plume_depenedent_constants.LAMBDA_DEEP = 0.0
                plume_depenedent_constants.LAMBDA_DOWN = 0.0
            else:
                plume_depenedent_constants.LAMBDA_DEEP = self.cu_param_config.LAMBDA_SHALLOW_DOWN
                plume_depenedent_constants.LAMBDA_DOWN = self.cu_param_config.LAMBDA_SHALLOW_DOWN

        self._set_plume_dependent_fields(
            t_excess_input=locals.cumulus_parameterization_input.t_excess,
            t_excess_internal=locals.cumulus_parameterization_internal.t_excess,
            vapor_excess_input=locals.cumulus_parameterization_input.vapor_excess,
            vapor_excess_internal=locals.cumulus_parameterization_internal.vapor_excess,
            ocean_fraction=locals.cumulus_parameterization_input.ocean_fraction,
            use_excess=plume_depenedent_constants.USE_EXCESS,
            t_input=locals.cumulus_parameterization_input.t,
            vapor_input=locals.cumulus_parameterization_input.vapor_timestep_start,
            grid_scale_forcing_t=locals.cumulus_parameterization_input.grid_scale_forcing_t,
            grid_scale_forcing_vapor=locals.cumulus_parameterization_input.grid_scale_forcing_vapor,
            subgrid_scale_forcing_t=locals.cumulus_parameterization_input.subgrid_scale_forcing_t,
            subgrid_scale_forcing_vapor=locals.cumulus_parameterization_input.subgrid_scale_forcing_vapor,
            t_internal=locals.cumulus_parameterization_internal.t,
            vapor_internal=locals.cumulus_parameterization_internal.vapor,
            t_pbl_internal=locals.cumulus_parameterization_internal.t_pbl,
            vapor_pbl_internal=locals.cumulus_parameterization_internal.vapor_pbl,
            moist_static_energy=locals.cumulus_parameterization_internal.moist_static_energy,
        )
