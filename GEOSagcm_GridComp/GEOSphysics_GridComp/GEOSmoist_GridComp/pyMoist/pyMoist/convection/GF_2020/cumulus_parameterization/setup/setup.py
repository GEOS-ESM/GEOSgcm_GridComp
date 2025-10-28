from ndsl import StencilFactory, QuantityFactory, ndsl_log
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.state import GF2020CumulusParameterizationState
from pyMoist.convection.GF_2020.cumulus_parameterization.locals import GF2020CumulusParameterizationLocals
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.convection.GF_2020.cumulus_parameterization.setup.stencils import (
    set_plume_dependent_fields,
    prefil_internal_fields,
    compute_scale_dependence_factor,
    get_random_number,
    initial_entrainment_detrainment,
    epsilon_min_max,
    calculate_arbitrary_numerical_parameter,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.setup.set_constants import set_constants
from pyMoist.convection.GF_2020.cumulus_parameterization.constants import MAXENS1, MAXENS2, MAXENS3
from ndsl.constants import X_DIM, Y_DIM, Z_DIM

from ndsl.dsl.typing import Int, Float


class Setup:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GF2020Config,
        cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
    ):
        # make configuration visible at runtime
        self.config = config
        self.cu_param_config = cumulus_parameterization_config

        # construct stencils and functions
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

        self._get_random_number = stencil_factory.from_dims_halo(
            func=get_random_number,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"USE_RANDOM_NUMBER": cumulus_parameterization_config.USE_RANDOM_NUMBER},
        )

        self._initial_entrainment_detrainment = stencil_factory.from_dims_halo(
            func=initial_entrainment_detrainment,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._epsilon_min_max = stencil_factory.from_dims_halo(
            func=epsilon_min_max,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._calculate_arbitrary_numerical_parameter = stencil_factory.from_dims_halo(
            func=calculate_arbitrary_numerical_parameter,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        saturation_tables: SaturationVaporPressureTable,
        plume_dependent_constants: GF2020PlumeDependentConstants,
        plume: str,
    ):
        plume_dependent_constants = set_constants(self.cu_param_config, plume_dependent_constants, plume)

        if plume_dependent_constants.ENABLE_PLUME == 1:
            # compute/prefil the last few fields needed for the rest of the scheme
            self._set_plume_dependent_fields(
                t_excess=state.input.t_excess,
                t_excess_local=locals.t_excess,
                vapor_excess=state.input.vapor_excess,
                vapor_excess_local=locals.vapor_excess,
                ocean_fraction=state.input.ocean_fraction,
                use_excess=plume_dependent_constants.USE_EXCESS,
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
                plume=plume_dependent_constants.PLUME_INDEX,
                maximum_updraft_origin_level=locals.maximum_updraft_origin_level,
                kstabm=locals.kstabm,
                ocean_fraction=state.input.ocean_fraction,
                ocean_fraction_local=locals.ocean_fraction,
                cap_max=locals.cap_max,
                error_code_2=locals.error_code_2,
                error_code_3=locals.error_code_3,
                error_code_string=locals.error_code_string,
                CAP_MAX_INC=plume_dependent_constants.CAP_MAX_INC,
                cap_max_increment=locals.cap_max_increment,
                geopotential_height=state.input_output.geopotential_height,
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
                t_wetbulb=locals.t_wetbulb,
                vapor_wetbulb=locals.vapor_wetbulb,
                tau_ecmwf=locals.tau_ecmwf,
                f_dicycle_modified=locals.f_dicycle_modified,
                add_buoy_modified=locals.add_buoy_modified,
                scale_dependence_factor_downdraft=locals.scale_dependence_factor_downdraft,
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
                plume=plume_dependent_constants.PLUME_INDEX,
                scale_dependence_factor=state.output.scale_dependence_factor,
                seed_convection=state.input.seed_convection,
                error_code=state.output.error_code,
                error_code_string=locals.error_code_string,
                grid_length=state.input_output.grid_length,
            )

            # create a real random number in the interval [-use_random_num, +use_random_num]
            self._get_random_number(
                plume=plume_dependent_constants.PLUME_INDEX,
                random_number=locals.random_number,
            )

            # define entrainment/detrainment profiles for updrafts
            self._initial_entrainment_detrainment(
                plume=plume_dependent_constants.PLUME_INDEX,
                lateral_entrainment_rate=state.input.lateral_entrainment_rate,
                current_plume_rate=plume_dependent_constants.ENTRAINMENT_RATE,
                entrainment_rate=state.output.entrainment_rate,
                updraft_detrainment_function=locals.updraft_detrainment_function,
            )

            # max/min allowed value for epsilon (ratio downdraft base mass flux/updraft base mass flux
            # note : to make the evaporation stronger => increase "epsilon_min"
            self._epsilon_min_max(
                ocean_fraction=state.input.ocean_fraction,
                epsilon_min=locals.epsilon_min,
                epsilon_max=locals.epsilon_max,
                MINIMUM_EVAP_FRACTION_OCEAN=plume_dependent_constants.MINIMUM_EVAP_FRACTION_OCEAN,
                MAXIMUM_EVAP_FRACTION_OCEAN=plume_dependent_constants.MAXIMUM_EVAP_FRACTION_OCEAN,
                MINIMUM_EVAP_FRACTION_LAND=plume_dependent_constants.MINIMUM_EVAP_FRACTION_LAND,
                MAXIMUM_EVAP_FRACTION_LAND=plume_dependent_constants.MAXIMUM_EVAP_FRACTION_LAND,
            )

            # calculate arbitrary numerical parameter
            self._calculate_arbitrary_numerical_parameter(
                arbitrary_numerical_parameter=locals.arbitrary_numerical_parameter,
            )
