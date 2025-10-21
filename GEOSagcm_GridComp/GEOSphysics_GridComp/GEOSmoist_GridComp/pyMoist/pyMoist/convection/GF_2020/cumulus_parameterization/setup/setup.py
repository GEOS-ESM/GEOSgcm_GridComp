from ndsl import StencilFactory, QuantityFactory, ndsl_log
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.state import GF2020CumulusParameterizationState
from pyMoist.convection.GF_2020.cumulus_parameterization.locals import GF2020CumulusParameterizationLocals
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.convection.GF_2020.cumulus_parameterization.setup.stencils import (
    set_plume_dependent_fields,
    prefil_internal_fields,
    compute_scale_dependence_factor,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.constants import MAXENS1, MAXENS2, MAXENS3
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

        self._prefil_internal_fields = stencil_factory.from_dims_halo(
            func=prefil_internal_fields,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "CAP_MAXS": cumulus_parameterization_config.CAP_MAXS,
                "ENSEMBLE_MEMBERS": MAXENS1 * MAXENS2 * MAXENS3,
            },
        )

        self._compute_scale_dependence_factor = stencil_factory.from_dims_halo(
            func=compute_scale_dependence_factor,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"USE_SCALE_DEP": cumulus_parameterization_config.USE_SCALE_DEP},
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        saturation_tables: SaturationVaporPressureTable,
        plume: str,
        plume_depenedent_constants: GF2020PlumeDependentConstants,
    ):
        if plume == "shallow":
            # set a number of plume dependent constants
            plume_depenedent_constants.PLUME_INDEX = 1
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

            # maximum depth (mb) of capping inversion (larger cap = no convection)
            if self.cu_param_config.ZERO_DIFF == 1 or self.cu_param_config.MOIST_TRIGGER == 0:
                plume_depenedent_constants.CAP_MAX_INC = 25.0
            else:
                plume_depenedent_constants.CAP_MAX_INC = 10.0

            # lambda_U parameter for momentum transport
            if self.cu_param_config.PRESSURE_GRADIENT_CONSTANT != 0.0:
                plume_depenedent_constants.LAMBDA_DEEP = 0.0
                plume_depenedent_constants.LAMBDA_DOWN = 0.0
            else:
                plume_depenedent_constants.LAMBDA_DEEP = self.cu_param_config.LAMBDA_DEEP
                plume_depenedent_constants.LAMBDA_DOWN = self.cu_param_config.LAMBDA_SHALLOW_DOWN
        if plume == "mid":
            # set a number of plume dependent constants
            plume_depenedent_constants.PLUME_INDEX = 2
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

            # maximum depth (mb) of capping inversion (larger cap = no convection)
            if self.cu_param_config.ZERO_DIFF == 1 or self.cu_param_config.MOIST_TRIGGER == 0:
                plume_depenedent_constants.CAP_MAX_INC = 10.0
            else:
                plume_depenedent_constants.CAP_MAX_INC = 90.0

            # lambda_U parameter for momentum transport
            if self.cu_param_config.PRESSURE_GRADIENT_CONSTANT != 0.0:
                plume_depenedent_constants.LAMBDA_DEEP = 0.0
                plume_depenedent_constants.LAMBDA_DOWN = 0.0
            else:
                plume_depenedent_constants.LAMBDA_DEEP = self.cu_param_config.LAMBDA_SHALLOW_DOWN
                plume_depenedent_constants.LAMBDA_DOWN = self.cu_param_config.LAMBDA_SHALLOW_DOWN
        if plume == "deep":
            # set a number of plume dependent constants
            plume_depenedent_constants.PLUME_INDEX = 3
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

            # maximum depth (mb) of capping inversion (larger cap = no convection)
            if self.cu_param_config.ZERO_DIFF == 1 or self.cu_param_config.MOIST_TRIGGER == 0:
                plume_depenedent_constants.CAP_MAX_INC = 20.0
            else:
                plume_depenedent_constants.CAP_MAX_INC = 90.0

            # lambda_U parameter for momentum transport
            if self.cu_param_config.PRESSURE_GRADIENT_CONSTANT != 0.0:
                plume_depenedent_constants.LAMBDA_DEEP = 0.0
                plume_depenedent_constants.LAMBDA_DOWN = 0.0
            else:
                plume_depenedent_constants.LAMBDA_DEEP = self.cu_param_config.LAMBDA_SHALLOW_DOWN
                plume_depenedent_constants.LAMBDA_DOWN = self.cu_param_config.LAMBDA_SHALLOW_DOWN

        # compute/prefil the last few fields needed for the rest of the scheme
        self._set_plume_dependent_fields(
            t_excess=state.input.t_excess,
            t_excess_local=locals.t_excess,
            vapor_excess=state.input.vapor_excess,
            vapor_excess_local=locals.vapor_excess,
            ocean_fraction=state.input.ocean_fraction,
            use_excess=plume_depenedent_constants.USE_EXCESS,
            t_old=state.input_output.t_old,
            vapor_old=state.input_output.vapor_old,
            grid_scale_forcing_t=state.input.grid_scale_forcing_t,
            grid_scale_forcing_vapor=state.input.grid_scale_forcing_vapor,
            subgrid_scale_forcing_t=state.input.subgrid_scale_forcing_t,
            subgrid_scale_forcing_vapor=state.input.subgrid_scale_forcing_vapor,
            t_new=locals.t_new,
            vapor_new=locals.vapor_new,
            t_new_pbl=locals.t_new_pbl,
            vapor_new_pbl=locals.vapor_new_pbl,
            moist_static_energy=locals.moist_static_energy,
        )

        self._prefil_internal_fields(
            kbmax=locals.kbmax,
            kstamb=locals.kstamb,
            ocean_fraction=locals.ocean_fraction,
            ocean_fraction_local=locals.ocean_fraction,
            cap_max=locals.cap_max,
            ierr2=locals.ierr2,
            ierr3=locals.ierr3,
            ierrc=locals.ierrc,
            CAP_MAX_INC=plume_depenedent_constants.CAP_MAX_INC,
            max_increment=locals.max_increment,
            geopotential_height=locals.geopotential_height,
            geopotential_height_local=locals.geopotential_height,
            geopotential_height_modified_local=locals.geopotential_height_modified,
            cloud_work_function_0=locals.cloud_work_function_0,
            cloud_work_function_1=locals.cloud_work_function_1,
            cloud_work_function_2=locals.cloud_work_function_2,
            cloud_work_function_3=locals.cloud_work_function_3,
            cloud_work_function_0_pbl=locals.cloud_work_function_0_pbl,
            cloud_work_function_1_pbl=locals.cloud_work_function_1_pbl,
            cloud_work_function_1_fa=locals.cloud_work_function_1_fa,
            cin1=locals.cin1,
            k_x_modified=locals.k_x_modified,
            epsilon_local=locals.epsilon,
            pbl_time_scale=locals.pbl_time_scale,
            plume=plume_depenedent_constants.PLUME_INDEX,
            t_wetbulb=locals.t_wetbulb,
            vapor_wetbulb=locals.vapor_wetbulb,
            tau_ecmwf=locals.tau_ecmwf,
            f_dicycle_modified=locals.f_dicycle_modified,
            add_buoy_modified=locals.add_buoy_modified,
            hcdo=locals.hcdo,
            cupclw=locals.cupclw,
            qrcdo=locals.qrcdo,
            hcot=locals.hcot,
            c1d=locals.c1d,
            evap_bcb=locals.evap_bcb,
            mass_flux_ensemble=locals.mass_flux_ensemble,
            precipitation_ensemble=locals.precipitation_ensemble,
            epsilon=state.output.epsilon,
            precip=state.output.precip,
            scale_dependence_factor=state.output.scale_dependence_factor,
            lightning_density=state.output.lightning_density,
        )

        # scale dependence factor (sig), version new
        self._compute_scale_dependence_factor(
            plume=plume_depenedent_constants.PLUME_INDEX,
            scale_dependence_factor=locals.scale_dependence_factor,
            seed_convection=state.input.seed_convection,
            ierr=state.output.ierr,
            ierrc=locals.ierrc,
            grid_length=state.input_output.grid_length,
        )

        # create a real random number in the interval [-use_random_num, +use_random_num]
