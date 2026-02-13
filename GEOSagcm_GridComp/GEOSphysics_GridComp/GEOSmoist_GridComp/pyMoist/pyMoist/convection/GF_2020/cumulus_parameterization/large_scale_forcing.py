from ndsl import StencilFactory, QuantityFactory, Quantity, Local, NDSLRuntime
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)
from ndsl.dsl.gt4py import computation, interval, FORWARD, PARALLEL, K
from ndsl.dsl.typing import FloatField, FloatFieldIJ, IntFieldIJ, Int
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import (
    FloatField_Plume,
    IntFieldIJ_Plume,
    FloatFieldIJ_Ensemble,
    FloatFieldIJ_Plume,
)
import pyMoist.constants as constants
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants


def ddim_copy(
    field_in: IntFieldIJ_Plume,
    field_out: IntFieldIJ,
    plume: Int,
):
    with computation(FORWARD), interval(0, 1):
        field_out = field_in[0, 0][plume]


def ensemble_forcing(
    error_code: IntFieldIJ_Plume,
    error_code_2: IntFieldIJ,
    error_code_3: IntFieldIJ,
    updraft_lfc_level: IntFieldIJ_Plume,
    vapor_forced: FloatField,
    ocean_fraction: FloatFieldIJ,
    omega: FloatField,
    cape_removal_time_scale: FloatFieldIJ,
    arbitrary_numerical_parameter: FloatFieldIJ,
    f_dicycle_modified: FloatFieldIJ,
    cloud_workfunction_0: FloatFieldIJ,
    cloud_workfunction_0_modified: FloatFieldIJ,
    cloud_workfunction_1: FloatFieldIJ,
    cloud_workfunction_1_pbl: FloatFieldIJ,
    mass_flux_ensemble: FloatFieldIJ_Ensemble,
    internal_mass_flux_ensemble: FloatFieldIJ_Ensemble,
    precipitation_ensemble: FloatFieldIJ_Ensemble,
    CLOSURE_CHOICE: Int,
    plume: Int,
):
    from __externals__ import ENSEMBLE_MEMBERS, DTIME, ZERO_DIFF, DIURNAL_CYCLE

    with computation(FORWARD), interval(0, 1):
        ensemble_adjustment: FloatFieldIJ = 1.0
        member = 0
        while member < ENSEMBLE_MEMBERS:
            mass_flux_ensemble[0, 0][member] = 0.0
            member += 1

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            workfunction_diff_1: FloatFieldIJ = (cloud_workfunction_1 - cloud_workfunction_0) / DTIME

            internal_mass_flux_ensemble[0, 0][0] = max(
                0.0, (cloud_workfunction_1 - cloud_workfunction_0) / DTIME
            )

            internal_mass_flux_ensemble[0, 0][1] = internal_mass_flux_ensemble[0, 0][0]
            internal_mass_flux_ensemble[0, 0][2] = internal_mass_flux_ensemble[0, 0][0]
            internal_mass_flux_ensemble[0, 0][15] = internal_mass_flux_ensemble[0, 0][0]

            # more like Brown (1979), or Frank-Cohen (199?)
            # omega is in Pa/s
            omega_modified: FloatFieldIJ = 0.0

            level_counter: IntFieldIJ = 0

            internal_mass_flux_ensemble[0, 0][3] = 0.0

    with computation(FORWARD), interval(...):
        if (
            error_code[0, 0][plume] == 0
            and K >= max(0, updraft_lfc_level[0, 0][plume] - 1)
            and K <= updraft_lfc_level[0, 0][plume] + 1
        ):
            beta = 1.0
            omega_modified = omega_modified - omega / constants.MAPL_GRAV / beta
            level_counter = level_counter + 1

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            if level_counter > 0:
                internal_mass_flux_ensemble[0, 0][3] = omega_modified / level_counter  # kg[air]/m^3 * m/s
            internal_mass_flux_ensemble[0, 0][3] = max(0.0, internal_mass_flux_ensemble[0, 0][3])
            internal_mass_flux_ensemble[0, 0][4] = internal_mass_flux_ensemble[0, 0][3]
            internal_mass_flux_ensemble[0, 0][5] = internal_mass_flux_ensemble[0, 0][3]
            internal_mass_flux_ensemble[0, 0][13] = internal_mass_flux_ensemble[0, 0][3]

            # more like Krishnamurti et al.;
            # assuming that omega.at(K=cloud_top_level)*vapor_forced.at(K=cloud_top_level)<<
            # omega.at(K=updraft_lfc_level)*vapor_forced.at(K=updraft_lfc_level)
            moisture_convergence: FloatFieldIJ = (
                -omega.at(K=updraft_lfc_level[0, 0][plume])
                * vapor_forced.at(K=updraft_lfc_level[0, 0][plume])
                / constants.MAPL_GRAV
            )

            moisture_convergence = max(0.0, moisture_convergence)
            internal_mass_flux_ensemble[0, 0][6] = moisture_convergence
            internal_mass_flux_ensemble[0, 0][7] = internal_mass_flux_ensemble[0, 0][6]
            internal_mass_flux_ensemble[0, 0][8] = internal_mass_flux_ensemble[0, 0][6]
            internal_mass_flux_ensemble[0, 0][14] = internal_mass_flux_ensemble[0, 0][6]

            # more like  Betchold et al (2014). Note that AA1 already includes the forcings tendencies
            internal_mass_flux_ensemble[0, 0][9] = cloud_workfunction_1 / cape_removal_time_scale

            internal_mass_flux_ensemble[0, 0][10] = internal_mass_flux_ensemble[0, 0][9]
            internal_mass_flux_ensemble[0, 0][11] = internal_mass_flux_ensemble[0, 0][9]
            internal_mass_flux_ensemble[0, 0][12] = internal_mass_flux_ensemble[0, 0][9]

            if CLOSURE_CHOICE == 0 and workfunction_diff_1 < 0.0:
                internal_mass_flux_ensemble[0, 0][0] = 0.0
                internal_mass_flux_ensemble[0, 0][1] = 0.0
                internal_mass_flux_ensemble[0, 0][2] = 0.0
                internal_mass_flux_ensemble[0, 0][15] = 0.0
                internal_mass_flux_ensemble[0, 0][9] = 0.0
                internal_mass_flux_ensemble[0, 0][10] = 0.0
                internal_mass_flux_ensemble[0, 0][11] = 0.0
                internal_mass_flux_ensemble[0, 0][12] = 0.0

            workfunction_diff_2: FloatFieldIJ = (
                cloud_workfunction_0_modified - (cloud_workfunction_1)
            ) / arbitrary_numerical_parameter

            if workfunction_diff_2 <= 0.0 and workfunction_diff_2 > -0.1 * arbitrary_numerical_parameter:
                workfunction_diff_2 = -0.1 * arbitrary_numerical_parameter

            if workfunction_diff_2 > 0.0 and workfunction_diff_2 < 1.0e-2:
                workfunction_diff_2 = 1.0e-2

            # over water, enforce small cap for some of the closures
            if ocean_fraction < 0.1 and (error_code_2 > 0.0 or error_code_3 > 0):
                member = 0
                while member < 16:
                    internal_mass_flux_ensemble[0, 0][member] = (
                        ensemble_adjustment * internal_mass_flux_ensemble[0, 0][member]
                    )
                    member += 1

            # special treatment for stability closures
            if workfunction_diff_2 < 0.0:
                if internal_mass_flux_ensemble[0, 0][0] > 0.0:
                    mass_flux_ensemble[0, 0][0] = max(
                        0.0, -internal_mass_flux_ensemble[0, 0][0] / workfunction_diff_2
                    )

                if internal_mass_flux_ensemble[0, 0][1] > 0.0:
                    mass_flux_ensemble[0, 0][1] = max(
                        0.0, -internal_mass_flux_ensemble[0, 0][1] / workfunction_diff_2
                    )

                if internal_mass_flux_ensemble[0, 0][2] > 0.0:
                    mass_flux_ensemble[0, 0][2] = max(
                        0.0, -internal_mass_flux_ensemble[0, 0][2] / workfunction_diff_2
                    )

                if internal_mass_flux_ensemble[0, 0][15] > 0.0:
                    mass_flux_ensemble[0, 0][15] = max(
                        0.0, -internal_mass_flux_ensemble[0, 0][15] / workfunction_diff_2
                    )
            else:
                internal_mass_flux_ensemble[0, 0][0] = 0.0
                internal_mass_flux_ensemble[0, 0][1] = 0.0
                internal_mass_flux_ensemble[0, 0][2] = 0.0
                internal_mass_flux_ensemble[0, 0][15] = 0.0

            mass_flux_ensemble[0, 0][3] = max(0.0, internal_mass_flux_ensemble[0, 0][3])
            mass_flux_ensemble[0, 0][4] = max(0.0, internal_mass_flux_ensemble[0, 0][4])
            mass_flux_ensemble[0, 0][5] = max(0.0, internal_mass_flux_ensemble[0, 0][5])
            mass_flux_ensemble[0, 0][13] = max(0.0, internal_mass_flux_ensemble[0, 0][13])

            adjustment = max(1.0e-3, precipitation_ensemble[0, 0][6])
            mass_flux_ensemble[0, 0][6] = max(0.0, internal_mass_flux_ensemble[0, 0][6] / adjustment)

            adjustment = max(1.0e-3, precipitation_ensemble[0, 0][7])
            mass_flux_ensemble[0, 0][7] = max(0.0, internal_mass_flux_ensemble[0, 0][7] / adjustment)

            adjustment = max(1.0e-3, precipitation_ensemble[0, 0][8])
            mass_flux_ensemble[0, 0][8] = max(0.0, internal_mass_flux_ensemble[0, 0][8] / adjustment)

            adjustment = max(1.0e-3, precipitation_ensemble[0, 0][14])
            mass_flux_ensemble[0, 0][14] = max(0.0, internal_mass_flux_ensemble[0, 0][14] / adjustment)

            if workfunction_diff_2 < 0.0:
                mass_flux_ensemble[0, 0][9] = max(
                    0.0, -internal_mass_flux_ensemble[0, 0][9] / workfunction_diff_2
                )
                mass_flux_ensemble[0, 0][10] = max(
                    0.0, -internal_mass_flux_ensemble[0, 0][10] / workfunction_diff_2
                )
                mass_flux_ensemble[0, 0][11] = max(
                    0.0, -internal_mass_flux_ensemble[0, 0][11] / workfunction_diff_2
                )
                mass_flux_ensemble[0, 0][12] = max(
                    0.0, -internal_mass_flux_ensemble[0, 0][12] / workfunction_diff_2
                )
            else:
                mass_flux_ensemble[0, 0][9] = 0.0
                mass_flux_ensemble[0, 0][10] = 0.0
                mass_flux_ensemble[0, 0][11] = 0.0
                mass_flux_ensemble[0, 0][12] = 0.0

            if CLOSURE_CHOICE >= 1:
                member = 0
                while member < 16:
                    mass_flux_ensemble[0, 0][member] = mass_flux_ensemble[0, 0][CLOSURE_CHOICE]
                    member += 1

            # over the land, only applies closure 9
            if ZERO_DIFF == 0 and CLOSURE_CHOICE == 0:
                member = 0
                while member < 16:
                    mass_flux_ensemble[0, 0][member] = (1.0 - ocean_fraction) * mass_flux_ensemble[0, 0][
                        9
                    ] + ocean_fraction * mass_flux_ensemble[0, 0][member]
                    member += 1

    with computation(FORWARD), interval(0, 1):
        if DIURNAL_CYCLE == 1 or DIURNAL_CYCLE == 6:
            f_dicycle_modified = 0.0
            if error_code[0, 0][plume] == 0:
                workfunction_diff_3 = (
                    cloud_workfunction_1 - cloud_workfunction_1_pbl
                ) / cape_removal_time_scale
                if workfunction_diff_2 < 0:
                    f_dicycle_modified = max(0.0, -workfunction_diff_3 / workfunction_diff_2)
                f_dicycle_modified = mass_flux_ensemble[0, 0][9] - f_dicycle_modified


def ensemble_forcing_mid_plume(
    error_code: IntFieldIJ_Plume,
    updraft_origin_level: IntFieldIJ_Plume,
    updraft_lfc_level: IntFieldIJ_Plume,
    pbl_level: IntFieldIJ,
    p_cloud_levels_forced: FloatField_Plume,
    cloud_moist_static_energy_forced: FloatField,
    environment_moist_static_energy_cloud_levels_forced: FloatField,
    dmoist_static_energydt: FloatField,
    convective_scale_velocity: FloatFieldIJ,
    cape_removal_time_scale: FloatFieldIJ,
    arbitrary_numerical_parameter: FloatFieldIJ,
    f_dicycle_modified: FloatFieldIJ,
    cloud_workfunction_0: FloatFieldIJ,
    cloud_workfunction_0_modified: FloatFieldIJ,
    cloud_workfunction_1: FloatFieldIJ,
    cloud_workfunction_1_pbl: FloatFieldIJ,
    xff_mid: FloatFieldIJ_Ensemble,
    CLOSURE_CHOICE: Int,
    plume: Int,
):
    from __externals__ import DIURNAL_CYCLE, DTIME

    with computation(FORWARD), interval(0, 1):
        # initialization
        member = 0
        while member < 16:  # 16 should come from config file
            xff_mid[0, 0][member] = 0.0
            member += 1
        f_dicycle_modified = 0.0

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            workfunction_diff_1: FloatFieldIJ = (
                cloud_workfunction_0_modified - cloud_workfunction_1
            ) / arbitrary_numerical_parameter

            if workfunction_diff_1 <= 0.0 and workfunction_diff_1 > -0.1 * arbitrary_numerical_parameter:
                workfunction_diff_1 = -0.1 * arbitrary_numerical_parameter

            if workfunction_diff_1 > 0.0 and workfunction_diff_1 < 0.01:
                workfunction_diff_1 = 0.01

            # diurnal cycle mass flux
            if DIURNAL_CYCLE == 1 or DIURNAL_CYCLE == 6 or DIURNAL_CYCLE == 0:
                internal_workfunction_modified = cloud_workfunction_1_pbl / cape_removal_time_scale
                f_dicycle_modified = max(0.0, -internal_workfunction_modified / workfunction_diff_1)

            # closures 3 and 4 for mid
            if workfunction_diff_1 < 0.0:
                xff_mid[0, 0][2] = max(
                    0.0, -(cloud_workfunction_1 / cape_removal_time_scale) / workfunction_diff_1
                )
                xff_mid[0, 0][3] = f_dicycle_modified

    with computation(FORWARD), interval(0, 1):
        # boundary layer quasi-equilibrium (Raymond 1995)
        if error_code[0, 0][plume] == 0 and updraft_origin_level[0, 0][plume] < pbl_level + 1:
            blqe: FloatFieldIJ = 0.0

    with computation(FORWARD), interval(0, -1):
        if (
            error_code[0, 0][plume] == 0
            and updraft_origin_level[0, 0][plume] < pbl_level + 1
            and K < updraft_lfc_level[0, 0][plume]
        ):
            blqe = (
                blqe
                + 100.0
                * dmoist_static_energydt
                * (p_cloud_levels_forced[0, 0, 0][plume] - p_cloud_levels_forced[0, 0, 1][plume])
                / constants.MAPL_GRAV
            )

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            if updraft_origin_level[0, 0][plume] < pbl_level + 1:
                trash = max(
                    (
                        cloud_moist_static_energy_forced.at(K=updraft_lfc_level)
                        - environment_moist_static_energy_cloud_levels_forced.at(K=updraft_lfc_level)
                    ),
                    1.0e1,
                )
                xff_mid[0, 0][1] = max(0.0, blqe / trash)

            # W* closure (Grant,2001)
            xff_mid[0, 0][0] = 0.03 * convective_scale_velocity

    with computation(FORWARD), interval(0, 1):
        # set f_dicycle_modified=0 in case the CLOSURE_CHOICE is 4 and DIURNAL_CYCLE closure will be applied
        if DIURNAL_CYCLE > 0 and CLOSURE_CHOICE == 4:
            f_dicycle_modified = 0.0

    with computation(FORWARD), interval(...):
        if CLOSURE_CHOICE == 5:
            if DIURNAL_CYCLE > 0:
                if error_code[0, 0][plume] == 0:
                    xff_mid[0, 0][4] = 0.05 * xff_mid[0, 0][3]
                    f_dicycle_modified = 0.0

            elif DIURNAL_CYCLE == 0:
                if error_code[0, 0][plume] == 0:
                    xff_mid[0, 0][4] = 0.05 * xff_mid[0, 0][3]
                    f_dicycle_modified = 0.0


def effective_precipitation(
    error_code: IntFieldIJ_Plume,
    epsilon_forced: FloatFieldIJ_Plume,
    condensate_to_fall_forced: FloatField_Plume,
    evaporate_in_downdraft_forced: FloatField_Plume,
    effective_condensate_to_fall_forced: FloatField,
    plume: Int,
):
    with computation(FORWARD), interval(0, -1):
        if error_code[0, 0][plume] == 0:
            effective_condensate_to_fall_forced = (
                condensate_to_fall_forced[0, 0, 0][plume]
                + epsilon_forced[0, 0][plume] * evaporate_in_downdraft_forced[0, 0, 0][plume]
            )


class LargeScaleForcing(NDSLRuntime):
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GF2020Config,
        cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
    ):
        # init NDSLRuntime
        super().__init__(stencil_factory)
        
        self.config = config
        self.cumulus_parameterization_config = cumulus_parameterization_config

        # initialize locals
        quantity_factory.add_data_dimensions(
            {
                "ensemble_members": cumulus_parameterization_constants.MAXENS1
                * cumulus_parameterization_constants.MAXENS2
                * cumulus_parameterization_constants.MAXENS3
            }
        )
        self._internal_mass_flux_ensemble: Local = quantity_factory.zeros(
            [X_DIM, Y_DIM, "ensemble_members"], "n/a"
        )

        # construct stencils
        self._copy = stencil_factory.from_dims_halo(
            func=ddim_copy,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self.ensemble_forcing = stencil_factory.from_dims_halo(
            func=ensemble_forcing,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "DIURNAL_CYCLE": cumulus_parameterization_config.DIURNAL_CYCLE,
                "DTIME": cumulus_parameterization_config.DTIME,
                "ZERO_DIFF": cumulus_parameterization_config.ZERO_DIFF,
                "ENSEMBLE_MEMBERS": cumulus_parameterization_constants.MAXENS1
                * cumulus_parameterization_constants.MAXENS2
                * cumulus_parameterization_constants.MAXENS3,
            },
        )

        self._ensemble_forcing_mid_plume = stencil_factory.from_dims_halo(
            func=ensemble_forcing_mid_plume,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "DIURNAL_CYCLE": cumulus_parameterization_config.DIURNAL_CYCLE,
                "DTIME": cumulus_parameterization_config.DTIME,
                # "ENSEMBLE_MEMBERS": cumulus_parameterization_config.ENSEMBLE_MEMBERS,
            },
        )

        self._effective_precipitation = stencil_factory.from_dims_halo(
            func=effective_precipitation,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        error_code: Quantity,
        error_code_2: Quantity,
        error_code_3: Quantity,
        updraft_origin_level: Quantity,
        updraft_lfc_level: Quantity,
        cloud_top_level: Quantity,
        pbl_level: Quantity,
        ocean_fraction: Quantity,
        p_cloud_levels_forced: Quantity,
        vapor_forced: Quantity,
        condensate_to_fall_forced: Quantity,
        effective_condensate_to_fall_forced: Quantity,
        evaporate_in_downdraft_forced: Quantity,
        omega: Quantity,
        convective_scale_velocity: Quantity,
        normalized_massflux_updraft_forced: Quantity,
        normalized_massflux_downdraft_forced: Quantity,
        cloud_moist_static_energy: Quantity,
        cloud_moist_static_energy_forced: Quantity,
        environment_moist_static_energy_cloud_levels: Quantity,
        environment_moist_static_energy_cloud_levels_forced: Quantity,
        dmoist_static_energydt: Quantity,
        cloud_workfunction_0: Quantity,
        cloud_workfunction_0_modified: Quantity,
        cloud_workfunction_1: Quantity,
        cloud_workfunction_1_pbl: Quantity,
        arbitrary_numerical_parameter: Quantity,
        f_dicycle_modified: Quantity,
        cape_removal_time_scale: Quantity,
        epsilon_forced: Quantity,
        k_x_modified: Quantity,
        mass_flux_ensemble: Quantity,
        precipitation_ensemble: Quantity,
        xff_mid: Quantity,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        # copy error codes
        self._copy(field_in=error_code, field_out=error_code_2, plume=plume_dependent_constants.PLUME_INDEX)
        self._copy(field_in=error_code, field_out=error_code_3, plume=plume_dependent_constants.PLUME_INDEX)

        if plume_dependent_constants.PLUME_INDEX == cumulus_parameterization_constants.DEEP:

            self.ensemble_forcing(
                error_code=error_code,
                error_code_2=error_code_2,
                error_code_3=error_code_3,
                updraft_lfc_level=updraft_lfc_level,
                vapor_forced=vapor_forced,
                ocean_fraction=ocean_fraction,
                omega=omega,
                cape_removal_time_scale=cape_removal_time_scale,
                arbitrary_numerical_parameter=arbitrary_numerical_parameter,
                f_dicycle_modified=f_dicycle_modified,
                cloud_workfunction_0=cloud_workfunction_0,
                cloud_workfunction_0_modified=cloud_workfunction_0_modified,
                cloud_workfunction_1=cloud_workfunction_1,
                cloud_workfunction_1_pbl=cloud_workfunction_1_pbl,
                mass_flux_ensemble=mass_flux_ensemble,
                internal_mass_flux_ensemble=self._internal_mass_flux_ensemble,
                precipitation_ensemble=precipitation_ensemble,
                CLOSURE_CHOICE=plume_dependent_constants.CLOSURE_CHOICE,
                plume=plume_dependent_constants.PLUME_INDEX,
            )

        if plume_dependent_constants.PLUME_INDEX == cumulus_parameterization_constants.MID:

            self._ensemble_forcing_mid_plume(
                error_code=error_code,
                updraft_origin_level=updraft_origin_level,
                updraft_lfc_level=updraft_lfc_level,
                pbl_level=pbl_level,
                p_cloud_levels_forced=p_cloud_levels_forced,
                cloud_moist_static_energy_forced=cloud_moist_static_energy_forced,
                environment_moist_static_energy_cloud_levels_forced=environment_moist_static_energy_cloud_levels_forced,
                dmoist_static_energydt=dmoist_static_energydt,
                convective_scale_velocity=convective_scale_velocity,
                cape_removal_time_scale=cape_removal_time_scale,
                arbitrary_numerical_parameter=arbitrary_numerical_parameter,
                f_dicycle_modified=f_dicycle_modified,
                cloud_workfunction_0=cloud_workfunction_0,
                cloud_workfunction_0_modified=cloud_workfunction_0_modified,
                cloud_workfunction_1=cloud_workfunction_1,
                cloud_workfunction_1_pbl=cloud_workfunction_1_pbl,
                xff_mid=xff_mid,
                CLOSURE_CHOICE=plume_dependent_constants.CLOSURE_CHOICE,
                plume=plume_dependent_constants.PLUME_INDEX,
            )

        if plume_dependent_constants.PLUME_INDEX == cumulus_parameterization_constants.SHALLOW:
            raise NotImplementedError(
                "Shallow plume options not implemented in LargeScaleForcing, but you"
                "should have been caught before getting here. Something is wrong with the config checker.,"
                "Beware, more untested paths may be executing without warning."
            )

        self._effective_precipitation(
            error_code=error_code,
            epsilon_forced=epsilon_forced,
            condensate_to_fall_forced=condensate_to_fall_forced,
            evaporate_in_downdraft_forced=evaporate_in_downdraft_forced,
            effective_condensate_to_fall_forced=effective_condensate_to_fall_forced,
            plume=plume_dependent_constants.PLUME_INDEX,
        )
