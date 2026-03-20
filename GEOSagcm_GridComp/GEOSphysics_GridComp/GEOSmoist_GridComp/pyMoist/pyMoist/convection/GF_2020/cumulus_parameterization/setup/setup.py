import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
from ndsl import NDSLRuntime, Quantity, QuantityFactory, StencilFactory
from ndsl.constants import I_DIM, J_DIM, K_DIM
from ndsl.dsl.gt4py import FORWARD, PARALLEL, computation, interval
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ, Int, IntFieldIJ
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.constants import MAXENS1, MAXENS2, MAXENS3
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import (
    FloatField_Plume,
    FloatFieldIJ_Ensemble,
    FloatFieldIJ_Plume,
    IntFieldIJ_Plume,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.setup.set_constants import set_constants
from pyMoist.shared_generic_math import sigma


def set_plume_dependent_fields(
    t_excess: FloatFieldIJ,
    t_excess_local: FloatFieldIJ,
    vapor_excess: FloatFieldIJ,
    vapor_excess_local: FloatFieldIJ,
    ocean_fraction: FloatFieldIJ,
    use_excess: Int,
    t_old: FloatField,
    vapor_old: FloatField,
    grid_scale_forcing_t: FloatField,
    grid_scale_forcing_vapor: FloatField,
    subgrid_scale_forcing_t: FloatField,
    subgrid_scale_forcing_vapor: FloatField,
    t_new: FloatField,
    vapor_forced: FloatField,
    t_new_pbl: FloatField,
    vapor_forced_pbl: FloatField,
    dmoist_static_energydt: FloatField,
):
    """Fill or modify existing values of temperature/vapor excess before cumulus parameterization

    Args:
        t_excess (FloatFieldIJ)
        t_excess_local (FloatFieldIJ)
        vapor_excess (FloatFieldIJ)
        vapor_excess_local (FloatFieldIJ)
        ocean_fraction (FloatFieldIJ)
        use_excess (Int): trigger to control behavior
            0: fill with zero
            1: no change (retain existing values)
            2: enforce min/max everywhere
            other: enforce min/max only over ocean
        t_old (FloatField)
        vapor_old (FloatField)
        grid_scale_forcing_t (FloatField)
        grid_scale_forcing_vapor (FloatField)
        subgrid_scale_forcing_t (FloatField)
        subgrid_scale_forcing_vapor (FloatField)
        t_new (FloatField)
        vapor_forced (FloatField)
        t_new_pbl (FloatField)
        vapor_forced_pbl (FloatField)
        dmoist_static_energydt (FloatField)
    """
    from __externals__ import DT_MOIST

    with computation(FORWARD), interval(0, 1):
        if use_excess == 0:
            t_excess_local = 0.0
            vapor_excess_local = 0.0
        elif use_excess == 2:
            t_excess_local = min(0.5, max(0.2, t_excess))  # Kelvin
            vapor_excess_local = min(5.0e-4, max(1.0e-4, vapor_excess))  # kg kg^-1
        else:
            if ocean_fraction > 0.98:  # ocean
                t_excess_local = min(0.5, max(0.2, t_excess))  # Kelvin
                vapor_excess_local = min(5.0e-4, max(1.0e-4, vapor_excess))  # kg kg^-1

    with computation(PARALLEL), interval(0, -1):
        t_new = t_old + (subgrid_scale_forcing_t + grid_scale_forcing_t) * DT_MOIST
        vapor_forced = vapor_old + (subgrid_scale_forcing_vapor + grid_scale_forcing_vapor) * DT_MOIST
        vapor_forced = max(cumulus_parameterization_constants.smaller_qv, vapor_forced)

        # temp/water vapor modified only by bl processes
        t_new_pbl = t_old + (subgrid_scale_forcing_t) * DT_MOIST
        vapor_forced_pbl = vapor_old + (subgrid_scale_forcing_vapor) * DT_MOIST

        # moist static energy
        dmoist_static_energydt = cumulus_parameterization_constants.CP * (
            subgrid_scale_forcing_t + grid_scale_forcing_t
        ) + cumulus_parameterization_constants.XLV * (subgrid_scale_forcing_vapor + grid_scale_forcing_vapor)


def prefil_internal_fields(
    plume: Int,
    maximum_updraft_origin_level: IntFieldIJ,
    kstabm: IntFieldIJ_Plume,
    ocean_fraction: FloatFieldIJ,
    ocean_fraction_local: FloatFieldIJ,
    cap_max: FloatFieldIJ,
    error_code_2: IntFieldIJ,
    error_code_3: IntFieldIJ,
    CAP_MAX_INC: Float,
    cap_max_increment: FloatFieldIJ,
    geopotential_height: FloatField,
    geopotential_height_local: FloatField,
    geopotential_height_modified_local: FloatField,
    cloud_workfunction_0: FloatFieldIJ,
    cloud_workfunction_1: FloatFieldIJ,
    cloud_workfunction_2: FloatFieldIJ,
    cloud_workfunction_3: FloatFieldIJ,
    cloud_workfunction_0_pbl: FloatFieldIJ,
    cloud_workfunction_1_pbl: FloatFieldIJ,
    cloud_workfunction_1_fa: FloatFieldIJ,
    cin_1: FloatFieldIJ,
    scale_dependence_factor: FloatFieldIJ_Plume,
    scale_dependence_factor_downdraft: FloatFieldIJ,
    epsilon_forced: FloatFieldIJ_Plume,
    epsilon_local: FloatFieldIJ,
    cloud_moist_static_energy_downdraft_forced: FloatField,
    cloud_moist_static_energy_forced_transported: FloatField,
    k_x_modified: FloatFieldIJ,
    pbl_time_scale: FloatFieldIJ,
    t_wetbulb: FloatFieldIJ,
    vapor_wetbulb: FloatFieldIJ,
    cape_removal_time_scale: FloatFieldIJ,
    f_dicycle_modified: FloatFieldIJ,
    add_buoyancy: FloatFieldIJ,
    downdraft_saturation_vapor_forced: FloatField,
    c1d: FloatField,
    evaporation_below_cloud_base: FloatField,
    mass_flux_ensemble: FloatFieldIJ_Ensemble,
    precipitation_ensemble: FloatFieldIJ_Ensemble,
    precip: FloatFieldIJ_Plume,
    lightning_density: FloatFieldIJ,
):
    """Fill internal fields (from the cumulus parameterization state) to ensure no data is left over from
    the previous call. Most fields are filled with zero, but not all.

    Args:
        plume (Int): _description_
        maximum_updraft_origin_level (IntFieldIJ): _description_
        kstabm (IntFieldIJ_Plume): _description_
        ocean_fraction (FloatFieldIJ): _description_
        ocean_fraction_local (FloatFieldIJ): _description_
        cap_max (FloatFieldIJ): _description_
        error_code_2 (IntFieldIJ): _description_
        error_code_3 (IntFieldIJ): _description_
        CAP_MAX_INC (Float): _description_
        cap_max_increment (FloatFieldIJ): _description_
        geopotential_height (FloatField): _description_
        geopotential_height_local (FloatField): _description_
        geopotential_height_modified_local (FloatField): _description_
        cloud_workfunction_0 (FloatFieldIJ): _description_
        cloud_workfunction_1 (FloatFieldIJ): _description_
        cloud_workfunction_2 (FloatFieldIJ): _description_
        cloud_workfunction_3 (FloatFieldIJ): _description_
        cloud_workfunction_0_pbl (FloatFieldIJ): _description_
        cloud_workfunction_1_pbl (FloatFieldIJ): _description_
        cloud_workfunction_1_fa (FloatFieldIJ): _description_
        cin_1 (FloatFieldIJ): _description_
        scale_dependence_factor (FloatFieldIJ_Plume): _description_
        scale_dependence_factor_downdraft (FloatFieldIJ): _description_
        epsilon_forced (FloatFieldIJ_Plume): _description_
        epsilon_local (FloatFieldIJ): _description_
        cloud_moist_static_energy_downdraft_forced (FloatField): _description_
        cloud_moist_static_energy_forced_transported (FloatField): _description_
        k_x_modified (FloatFieldIJ): _description_
        pbl_time_scale (FloatFieldIJ): _description_
        t_wetbulb (FloatFieldIJ): _description_
        vapor_wetbulb (FloatFieldIJ): _description_
        cape_removal_time_scale (FloatFieldIJ): _description_
        f_dicycle_modified (FloatFieldIJ): _description_
        add_buoyancy (FloatFieldIJ): _description_
        downdraft_saturation_vapor_forced (FloatField): _description_
        c1d (FloatField): _description_
        evaporation_below_cloud_base (FloatField): _description_
        mass_flux_ensemble (FloatFieldIJ_Ensemble): _description_
        precipitation_ensemble (FloatFieldIJ_Ensemble): _description_
        precip (FloatFieldIJ_Plume): _description_
        lightning_density (FloatFieldIJ): _description_
    """
    from __externals__ import CAP_MAXS, ENSEMBLE_MEMBERS, k_end

    # reset to zero manually. cannot use locals.fill(0) because not all fields are reset to zero b/t plumes
    # internal fields
    with computation(FORWARD), interval(0, 1):
        maximum_updraft_origin_level = 0
        kstabm[0, 0][plume] = k_end - 2
        ocean_fraction_local = ocean_fraction
        cap_max = CAP_MAXS
        error_code_2 = 0
        error_code_3 = 0
        cap_max_increment = CAP_MAX_INC
        cloud_workfunction_0 = 0.0
        cloud_workfunction_1 = 0.0
        cloud_workfunction_2 = 0.0
        cloud_workfunction_3 = 0.0
        cloud_workfunction_0_pbl = 0.0
        cloud_workfunction_1_pbl = 0.0
        cloud_workfunction_1_fa = 0.0
        cin_1 = 0.0
        k_x_modified = 0.0
        epsilon_local = 0.0
        pbl_time_scale = 0.0
        t_wetbulb = 0.0
        vapor_wetbulb = 0.0
        cape_removal_time_scale = 0.0
        f_dicycle_modified = 0.0
        add_buoyancy = 0.0
        scale_dependence_factor_downdraft = 0.0

    # ensure locals are zero, no lingering data from previous call
    with computation(PARALLEL), interval(...):
        geopotential_height_local = geopotential_height
        geopotential_height_modified_local = geopotential_height
        cloud_moist_static_energy_downdraft_forced = 0.0
        downdraft_saturation_vapor_forced = 0.0
        cloud_moist_static_energy_forced_transported = 0.0
        c1d = 0.0
        evaporation_below_cloud_base = 0.0

    # ensure locals are zero, no lingering data from previous call
    # reset all ensemble members
    with computation(FORWARD), interval(0, 1):
        member = 0
        while member < ENSEMBLE_MEMBERS:
            mass_flux_ensemble[0, 0][member] = 0.0
            precipitation_ensemble[0, 0][member] = 0.0
            member = member + 1

    # reset state fields to zero
    with computation(FORWARD), interval(0, 1):
        epsilon_forced[0, 0][plume] = 0.0
        precip[0, 0][plume] = 0.0
        scale_dependence_factor[0, 0][plume] = 0.0
        lightning_density = 0.0


def compute_scale_dependence_factor(
    plume: Int,
    scale_dependence_factor: FloatFieldIJ_Plume,
    seed_convection: FloatFieldIJ,
    error_code: IntFieldIJ_Plume,
    grid_length: FloatFieldIJ,
):
    """Compute scale dependence factor for use later.

    Args:
        plume (Int)
        scale_dependence_factor (FloatFieldIJ_Plume)
        seed_convection (FloatFieldIJ)
        error_code (IntFieldIJ_Plume)
        grid_length (FloatFieldIJ)
    """
    from __externals__ import USE_SCALE_DEP

    # prepare mask to stop loop in next computation
    with computation(FORWARD), interval(0, 1):
        error_at_point = False

    with computation(FORWARD), interval(0, 1):
        if USE_SCALE_DEP == 0:
            scale_dependence_factor[0, 0][plume] = 1.0
        elif USE_SCALE_DEP == 1:
            if plume == cumulus_parameterization_constants.SHALLOW:
                scale_dependence_factor[0, 0][plume] = 1.0
            else:
                if seed_convection < 0.0:
                    error_code[0, 0][plume] = 1
                    error_at_point = True
                if error_at_point == False:
                    scale_dependence_factor[0, 0][plume] = sigma(grid_length)
                if seed_convection != 1.0:
                    scale_dependence_factor[0, 0][plume] = scale_dependence_factor[0, 0][plume] ** (
                        seed_convection * max(1.0, scale_dependence_factor[0, 0][plume])
                    )
                scale_dependence_factor[0, 0][plume] = max(
                    0.1, min(scale_dependence_factor[0, 0][plume], 1.0)
                )
                if scale_dependence_factor[0, 0][plume] <= 0.1:
                    error_code[0, 0][plume] = 1
                    error_at_point = True


def get_random_number(
    plume: Int,
    random_number: FloatFieldIJ,
):
    """Generate a random number for each column.

    THIS IS UNFINISHED, RIGHT NOW IT USES SERIALIZED FORTRAN DATA.
    THIS WILL NOT WORK WITH A PROPER INTEGRATION, NEED A BETTER SOLUTION.

    Args:
        plume (Int)
        random_number (FloatFieldIJ)
    """
    from __externals__ import USE_RANDOM_NUMBER

    with computation(FORWARD), interval(0, 1):
        if plume == cumulus_parameterization_constants.DEEP and USE_RANDOM_NUMBER > 1.0e-6:
            # need to figure out how to get system clock data
            random_number = random_number  # keep input data from fortran for now
        else:
            random_number = 0.0


def initial_entrainment_detrainment(
    plume: Int,
    lateral_entrainment_rate: FloatField,
    current_plume_rate: Float,
    entrainment_rate: FloatField_Plume,
    detrainment_function_updraft: FloatField,
):
    """Get initial entrainment/detrainment estimates based on data from the overarching model.

    Args:
        plume (Int)
        lateral_entrainment_rate (FloatField)
        current_plume_rate (Float)
        entrainment_rate (FloatField_Plume)
        detrainment_function_updraft (FloatField)
    """
    with computation(PARALLEL), interval(0, -1):
        entrainment_rate[0, 0, 0][plume] = lateral_entrainment_rate * current_plume_rate
        detrainment_function_updraft = lateral_entrainment_rate * current_plume_rate


def epsilon_min_max(
    ocean_fraction: FloatFieldIJ,
    epsilon_min: FloatFieldIJ,
    epsilon_max: FloatFieldIJ,
    MINIMUM_EVAP_FRACTION_OCEAN: Float,
    MAXIMUM_EVAP_FRACTION_OCEAN: Float,
    MINIMUM_EVAP_FRACTION_LAND: Float,
    MAXIMUM_EVAP_FRACTION_LAND: Float,
):
    """Set min/max epsilon for the current plume.

    Args:
        ocean_fraction (FloatFieldIJ)
        epsilon_min (FloatFieldIJ)
        epsilon_max (FloatFieldIJ)
        MINIMUM_EVAP_FRACTION_OCEAN (Float)
        MAXIMUM_EVAP_FRACTION_OCEAN (Float)
        MINIMUM_EVAP_FRACTION_LAND (Float)
        MAXIMUM_EVAP_FRACTION_LAND (Float)
    """
    with computation(FORWARD), interval(0, 1):
        if ocean_fraction > 0.99:  # water
            epsilon_min = MINIMUM_EVAP_FRACTION_OCEAN
            epsilon_max = MAXIMUM_EVAP_FRACTION_OCEAN
        else:  # land
            epsilon_min = MINIMUM_EVAP_FRACTION_LAND
            epsilon_max = MAXIMUM_EVAP_FRACTION_LAND


def calculate_arbitrary_numerical_parameter(
    arbitrary_numerical_parameter: FloatFieldIJ,
):
    """Set a scaling factor, used primarially (but not exclusively) for mass flux calculations.

    Args:
        arbitrary_numerical_parameter (FloatFieldIJ)
    """
    with computation(FORWARD), interval(0, 1):
        arbitrary_numerical_parameter = 0.1
        # approx xmb * timescale
        # other options (from fortran):
        # 0.1 * dtime*xmb_nm1(i)
        # 100.*(p_cup(i,kbcon(i))-p_cup(i,kbcon(i)+1))/(g*dtime)
        # 0.1*mbdt(i)


class Setup(NDSLRuntime):
    """Setup the GF2020 cumulus parameterization core. Prefill locals with appropriate values based on
    configuration, set plume dependent constants for the correct plume, and

    """

    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GF2020Config,
        cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
    ):
        # init NDSLRuntime
        super().__init__(stencil_factory)

        # make configuration visible at runtime
        self.config = config
        self.cu_param_config = cumulus_parameterization_config

        # construct stencils and functions
        self._set_plume_dependent_fields = stencil_factory.from_dims_halo(
            func=set_plume_dependent_fields,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={"DT_MOIST": config.DT_MOIST},
        )

        self._prefil_internal_fields = stencil_factory.from_dims_halo(
            func=prefil_internal_fields,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "CAP_MAXS": cumulus_parameterization_config.CAP_MAXS,
                "ENSEMBLE_MEMBERS": MAXENS1 * MAXENS2 * MAXENS3,
            },
        )

        self._compute_scale_dependence_factor = stencil_factory.from_dims_halo(
            func=compute_scale_dependence_factor,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={"USE_SCALE_DEP": cumulus_parameterization_config.USE_SCALE_DEP},
        )

        self._get_random_number = stencil_factory.from_dims_halo(
            func=get_random_number,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={"USE_RANDOM_NUMBER": cumulus_parameterization_config.USE_RANDOM_NUMBER},
        )

        self._initial_entrainment_detrainment = stencil_factory.from_dims_halo(
            func=initial_entrainment_detrainment,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )

        self._epsilon_min_max = stencil_factory.from_dims_halo(
            func=epsilon_min_max,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )

        self._calculate_arbitrary_numerical_parameter = stencil_factory.from_dims_halo(
            func=calculate_arbitrary_numerical_parameter,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )

    def __call__(
        self,
        error_code: Quantity,
        error_code_2: Quantity,
        error_code_3: Quantity,
        maximum_updraft_origin_level: Quantity,
        kstabm: Quantity,
        t_excess: Quantity,
        t_excess_local: Quantity,
        vapor_excess: Quantity,
        vapor_excess_local: Quantity,
        ocean_fraction: Quantity,
        ocean_fraction_local: Quantity,
        t_old: Quantity,
        t_new: Quantity,
        t_new_pbl: Quantity,
        vapor_old: Quantity,
        vapor_forced: Quantity,
        vapor_forced_pbl: Quantity,
        downdraft_saturation_vapor_forced: Quantity,
        grid_scale_forcing_t: Quantity,
        grid_scale_forcing_vapor: Quantity,
        subgrid_scale_forcing_t: Quantity,
        subgrid_scale_forcing_vapor: Quantity,
        dmoist_static_energydt: Quantity,
        cloud_moist_static_energy_downdraft_forced: Quantity,
        cloud_moist_static_energy_forced_transported: Quantity,
        cap_max: Quantity,
        cap_max_increment: Quantity,
        geopotential_height: Quantity,
        geopotential_height_local: Quantity,
        geopotential_height_modified_local: Quantity,
        cloud_workfunction_0: Quantity,
        cloud_workfunction_1: Quantity,
        cloud_workfunction_2: Quantity,
        cloud_workfunction_3: Quantity,
        cloud_workfunction_0_pbl: Quantity,
        cloud_workfunction_1_pbl: Quantity,
        cloud_workfunction_1_fa: Quantity,
        cin_1: Quantity,
        k_x_modified: Quantity,
        epsilon_forced: Quantity,
        epsilon_local: Quantity,
        epsilon_min: Quantity,
        epsilon_max: Quantity,
        pbl_time_scale: Quantity,
        t_wetbulb: Quantity,
        vapor_wetbulb: Quantity,
        cape_removal_time_scale: Quantity,
        f_dicycle_modified: Quantity,
        add_buoyancy: Quantity,
        scale_dependence_factor: Quantity,
        scale_dependence_factor_downdraft: Quantity,
        c1d: Quantity,
        evaporation_below_cloud_base: Quantity,
        mass_flux_ensemble: Quantity,
        precipitation_ensemble: Quantity,
        precip: Quantity,
        lightning_density: Quantity,
        seed_convection: Quantity,
        grid_length: Quantity,
        random_number: Quantity,
        lateral_entrainment_rate: Quantity,
        entrainment_rate: Quantity,
        detrainment_function_updraft: Quantity,
        arbitrary_numerical_parameter: Quantity,
        plume_dependent_constants: GF2020PlumeDependentConstants,
        plume: str,
    ):
        plume_dependent_constants = set_constants(self.cu_param_config, plume_dependent_constants, plume)

        if plume_dependent_constants.ENABLE_PLUME == 1:
            # compute/prefil the last few fields needed for the rest of the scheme
            self._set_plume_dependent_fields(
                t_excess=t_excess,
                t_excess_local=t_excess_local,
                vapor_excess=vapor_excess,
                vapor_excess_local=vapor_excess_local,
                ocean_fraction=ocean_fraction,
                use_excess=plume_dependent_constants.USE_EXCESS,
                t_old=t_old,
                vapor_old=vapor_old,
                grid_scale_forcing_t=grid_scale_forcing_t,
                grid_scale_forcing_vapor=grid_scale_forcing_vapor,
                subgrid_scale_forcing_t=subgrid_scale_forcing_t,
                subgrid_scale_forcing_vapor=subgrid_scale_forcing_vapor,
                t_new=t_new,
                vapor_forced=vapor_forced,
                t_new_pbl=t_new_pbl,
                vapor_forced_pbl=vapor_forced_pbl,
                dmoist_static_energydt=dmoist_static_energydt,
            )

            self._prefil_internal_fields(
                plume=plume_dependent_constants.PLUME_INDEX,
                maximum_updraft_origin_level=maximum_updraft_origin_level,
                kstabm=kstabm,
                ocean_fraction=ocean_fraction,
                ocean_fraction_local=ocean_fraction_local,
                cap_max=cap_max,
                error_code_2=error_code_2,
                error_code_3=error_code_3,
                CAP_MAX_INC=plume_dependent_constants.CAP_MAX_INC,
                cap_max_increment=cap_max_increment,
                geopotential_height=geopotential_height,
                geopotential_height_local=geopotential_height_local,
                geopotential_height_modified_local=geopotential_height_modified_local,
                cloud_workfunction_0=cloud_workfunction_0,
                cloud_workfunction_1=cloud_workfunction_1,
                cloud_workfunction_2=cloud_workfunction_2,
                cloud_workfunction_3=cloud_workfunction_3,
                cloud_workfunction_0_pbl=cloud_workfunction_0_pbl,
                cloud_workfunction_1_pbl=cloud_workfunction_1_pbl,
                cloud_workfunction_1_fa=cloud_workfunction_1_fa,
                cin_1=cin_1,
                k_x_modified=k_x_modified,
                epsilon_forced=epsilon_forced,
                epsilon_local=epsilon_local,
                pbl_time_scale=pbl_time_scale,
                t_wetbulb=t_wetbulb,
                vapor_wetbulb=vapor_wetbulb,
                cape_removal_time_scale=cape_removal_time_scale,
                f_dicycle_modified=f_dicycle_modified,
                add_buoyancy=add_buoyancy,
                scale_dependence_factor=scale_dependence_factor,
                scale_dependence_factor_downdraft=scale_dependence_factor_downdraft,
                cloud_moist_static_energy_downdraft_forced=cloud_moist_static_energy_downdraft_forced,
                downdraft_saturation_vapor_forced=downdraft_saturation_vapor_forced,
                cloud_moist_static_energy_forced_transported=cloud_moist_static_energy_forced_transported,
                c1d=c1d,
                evaporation_below_cloud_base=evaporation_below_cloud_base,
                mass_flux_ensemble=mass_flux_ensemble,
                precipitation_ensemble=precipitation_ensemble,
                precip=precip,
                lightning_density=lightning_density,
            )

            # scale dependence factor (sig), version new
            self._compute_scale_dependence_factor(
                plume=plume_dependent_constants.PLUME_INDEX,
                scale_dependence_factor=scale_dependence_factor,
                seed_convection=seed_convection,
                error_code=error_code,
                grid_length=grid_length,
            )

            # create a real random number in the interval [-use_random_num, +use_random_num]
            self._get_random_number(
                plume=plume_dependent_constants.PLUME_INDEX,
                random_number=random_number,
            )

            # define entrainment/detrainment profiles for updrafts
            self._initial_entrainment_detrainment(
                plume=plume_dependent_constants.PLUME_INDEX,
                lateral_entrainment_rate=lateral_entrainment_rate,
                current_plume_rate=plume_dependent_constants.ENTRAINMENT_RATE,
                entrainment_rate=entrainment_rate,
                detrainment_function_updraft=detrainment_function_updraft,
            )

            # max/min allowed value for epsilon (ratio downdraft base mass flux/updraft base mass flux
            # note : to make the evaporation stronger => increase "epsilon_min"
            self._epsilon_min_max(
                ocean_fraction=ocean_fraction,
                epsilon_min=epsilon_min,
                epsilon_max=epsilon_max,
                MINIMUM_EVAP_FRACTION_OCEAN=plume_dependent_constants.MINIMUM_EVAP_FRACTION_OCEAN,
                MAXIMUM_EVAP_FRACTION_OCEAN=plume_dependent_constants.MAXIMUM_EVAP_FRACTION_OCEAN,
                MINIMUM_EVAP_FRACTION_LAND=plume_dependent_constants.MINIMUM_EVAP_FRACTION_LAND,
                MAXIMUM_EVAP_FRACTION_LAND=plume_dependent_constants.MAXIMUM_EVAP_FRACTION_LAND,
            )

            # calculate arbitrary numerical parameter
            self._calculate_arbitrary_numerical_parameter(
                arbitrary_numerical_parameter=arbitrary_numerical_parameter,
            )
