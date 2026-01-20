from ndsl import StencilFactory, QuantityFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.state import GF2020CumulusParameterizationState
from pyMoist.convection.GF_2020.cumulus_parameterization.locals import GF2020CumulusParameterizationLocals
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)
from ndsl.logging import ndsl_log
from ndsl.dsl.gt4py import computation, interval, FORWARD, K, function, PARALLEL, exp
from ndsl.dsl.typing import FloatField, FloatFieldIJ, Float, IntFieldIJ, Int, BoolFieldIJ
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import (
    FloatField_Plume,
    IntFieldIJ_Plume,
    FloatFieldIJ_Ensemble,
    FloatFieldIJ_Plume,
)

import pyMoist.constants as constants
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants


def cup_forcing_ens_3d(
    xf_ens: FloatFieldIJ_Ensemble,
    error_code: IntFieldIJ_Plume,
    error_code2: IntFieldIJ_Plume,
    error_code3: IntFieldIJ_Plume,
    plume: Int,
    cloud_workfunction_1: FloatFieldIJ,
    cloud_workfunction_0: FloatFieldIJ,
    cloud_workfunction_1_pbl: FloatFieldIJ,
    xff_ens3: FloatFieldIJ_Ensemble,
    updraft_lfc_level: IntFieldIJ_Plume,
    omega: FloatField,
    moisture_convergence: FloatFieldIJ,
    vapor_forced: FloatField,
    tau_ecmwf: FloatFieldIJ,
    ichoice: IntFieldIJ,
    cloud_workfunction_0_modified: FloatFieldIJ,
    arbitrary_numerical_parameter: FloatFieldIJ,
    ocean_fraction: FloatFieldIJ,
    precipitation_ensemble: FloatFieldIJ_Ensemble,
    f_dicycle_modified: FloatFieldIJ,
):
    from __externals__ import DTIME, ZERO_DIFF, DIURNAL_CYCLE

    with computation(FORWARD), interval(0, 1):
        ens_adj: FloatFieldIJ = 1.0
        member = 0
        while member < 16:
            xf_ens[0, 0][member] = 0.0
            member += 1
    with computation(FORWARD), interval(...):
        if error_code[0, 0][plume] == 0:
            xff0: FloatFieldIJ = (cloud_workfunction_1 - cloud_workfunction_0) / DTIME

            xff_ens3[0, 0][0] = max(0.0, (cloud_workfunction_1 - cloud_workfunction_0) / DTIME)

            xff_ens3[0, 0][1] = xff_ens3[0, 0][0]
            xff_ens3[0, 0][2] = xff_ens3[0, 0][0]
            xff_ens3[0, 0][15] = xff_ens3[0, 0][0]

            xomg: FloatFieldIJ = 0.0

            kk = 0

            xff_ens3[0, 0][3] = 0.0
            lev = max(0, updraft_lfc_level[0, 0][plume] - 1)
            while lev <= updraft_lfc_level[0, 0][plume] + 1:
                betajb = 1.0
                xomg = xomg - omega.at(K=lev) / constants.MAPL_GRAV / betajb
                kk = kk + 1
                lev += 1

            if kk > 0:
                xff_ens3[0, 0][3] = xomg / float(kk)

            xff_ens3[0, 0][3] = max(0.0, xff_ens3[0, 0][3])
            xff_ens3[0, 0][4] = xff_ens3[0, 0][3]
            xff_ens3[0, 0][5] = xff_ens3[0, 0][3]
            xff_ens3[0, 0][13] = xff_ens3[0, 0][3]

            moisture_convergence = (
                -omega.at(K=updraft_lfc_level[0, 0][plume])
                * vapor_forced.at(K=updraft_lfc_level[0, 0][plume])
                / constants.MAPL_GRAV
            )

            moisture_convergence = max(0.0, moisture_convergence)
            xff_ens3[0, 0][6] = moisture_convergence
            xff_ens3[0, 0][7] = xff_ens3[0, 0][6]
            xff_ens3[0, 0][8] = xff_ens3[0, 0][6]
            xff_ens3[0, 0][14] = xff_ens3[0, 0][6]

            xff_ens3[0, 0][9] = cloud_workfunction_1 / tau_ecmwf

            xff_ens3[0, 0][10] = xff_ens3[0, 0][9]
            xff_ens3[0, 0][11] = xff_ens3[0, 0][9]
            xff_ens3[0, 0][12] = xff_ens3[0, 0][9]

            if ichoice == 0:
                if xff0 < 0.0:
                    xff_ens3[0, 0][0] = 0.0
                    xff_ens3[0, 0][1] = 0.0
                    xff_ens3[0, 0][2] = 0.0
                    xff_ens3[0, 0][15] = 0.0

                    xff_ens3[0, 0][9] = 0.0
                    xff_ens3[0, 0][10] = 0.0
                    xff_ens3[0, 0][11] = 0.0
                    xff_ens3[0, 0][12] = 0.0

            xk: FloatFieldIJ = (
                cloud_workfunction_0_modified - (cloud_workfunction_1)
            ) / arbitrary_numerical_parameter

            if xk <= 0.0 and xk > -0.1 * arbitrary_numerical_parameter:
                xk = -0.1 * arbitrary_numerical_parameter

            if xk > 0.0 and xk < 1.0e-2:
                xk = 1.0e-2

            if ocean_fraction < 0.1:
                if error_code2[0, 0][plume] > 0.0 or error_code3[0, 0][plume] > 0:
                    member = 0
                    while member < 16:
                        xff_ens3[0, 0][member] = ens_adj * xff_ens3[0, 0][member]
                        member += 1

            if xk < 0.0:
                if xff_ens3[0, 0][0] > 0.0:
                    xf_ens[0, 0][0] = max(0.0, -xff_ens3[0, 0][0] / xk)

                if xff_ens3[0, 0][1] > 0.0:
                    xf_ens[0, 0][1] = max(0.0, -xff_ens3[0, 0][1] / xk)

                if xff_ens3[0, 0][2] > 0.0:
                    xf_ens[0, 0][2] = max(0.0, -xff_ens3[0, 0][2] / xk)

                if xff_ens3[0, 0][15] > 0.0:
                    xf_ens[0, 0][15] = max(0.0, -xff_ens3[0, 0][15] / xk)
            else:
                xff_ens3[0, 0][0] = 0.0
                xff_ens3[0, 0][1] = 0.0
                xff_ens3[0, 0][2] = 0.0
                xff_ens3[0, 0][15] = 0.0

            xf_ens[0, 0][3] = max(0.0, xff_ens3[0, 0][3])
            xf_ens[0, 0][4] = max(0.0, xff_ens3[0, 0][4])
            xf_ens[0, 0][5] = max(0.0, xff_ens3[0, 0][5])
            xf_ens[0, 0][13] = max(0.0, xff_ens3[0, 0][13])

            a1: FloatFieldIJ = max(1.0e-3, precipitation_ensemble[0, 0][6])
            xf_ens[0, 0][6] = max(0.0, xff_ens3[0, 0][6] / a1)

            a1 = max(1.0e-3, precipitation_ensemble[0, 0][7])
            xf_ens[0, 0][7] = max(0.0, xff_ens3[0, 0][7] / a1)

            a1 = max(1.0e-3, precipitation_ensemble[0, 0][8])
            xf_ens[0, 0][8] = max(0.0, xff_ens3[0, 0][8] / a1)

            a1 = max(1.0e-3, precipitation_ensemble[0, 0][14])
            xf_ens[0, 0][14] = max(0.0, xff_ens3[0, 0][14] / a1)

            if xk < 0.0:
                xf_ens[0, 0][9] = max(0.0, -xff_ens3[0, 0][9] / xk)
                xf_ens[0, 0][10] = max(0.0, -xff_ens3[0, 0][10] / xk)
                xf_ens[0, 0][11] = max(0.0, -xff_ens3[0, 0][11] / xk)
                xf_ens[0, 0][12] = max(0.0, -xff_ens3[0, 0][12] / xk)
            else:
                xf_ens[0, 0][9] = 0.0
                xf_ens[0, 0][10] = 0.0
                xf_ens[0, 0][11] = 0.0
                xf_ens[0, 0][12] = 0.0

            if ichoice >= 1:
                member = 0
                while member < 16:
                    xf_ens[0, 0][member] = xf_ens[0, 0][ichoice]
                    member += 1

            # Not sure if this stays inside conditional?
            if ZERO_DIFF == 0 and ichoice == 0:
                member = 0
                while member < 16:
                    xf_ens[0, 0][member] = (1.0 - ocean_fraction) * xf_ens[0, 0][9] + ocean_fraction * xf_ens[
                        0, 0
                    ][member]
                    member += 1

    with computation(FORWARD), interval(...):
        if DIURNAL_CYCLE == 1:
            f_dicycle_modified = 0.0
            if error_code[0, 0][plume] == 0:
                xff_dicycle = (cloud_workfunction_1 - cloud_workfunction_1_pbl) / tau_ecmwf
                if xk < 0:
                    f_dicycle_modified = max(0.0, -xff_dicycle / xk)
                f_dicycle_modified = xf_ens[0, 0][9] - f_dicycle_modified


def cup_forcing_ens_3d_mid(
    error_code: IntFieldIJ_Plume,
    plume: Int,
    k_x: FloatFieldIJ_Ensemble,
    cloud_workfunction_1: FloatFieldIJ,
    cloud_workfunction_0: FloatFieldIJ,
    cloud_workfunction_0_modified: FloatFieldIJ,
    cloud_workfunction_1_pbl: FloatFieldIJ,
    arbitrary_numerical_parameter: FloatFieldIJ,
    tau_ecmwf: FloatFieldIJ,
    f_dicycle_modified: FloatFieldIJ,
    xff_mid: FloatFieldIJ_Ensemble,
    updraft_origin_level: IntFieldIJ_Plume,
    pbl_level: IntFieldIJ,
    updraft_lfc_level: IntFieldIJ_Plume,
    moist_static_energy: FloatField,
    p_cloud_levels_forced: FloatField_Plume,
    cloud_moist_static_energy_forced: FloatField,
    env_moist_static_energy_cloud_levels_forced: FloatField,
    convective_scale_velocity: FloatFieldIJ,
    ichoice: IntFieldIJ,
):
    from __externals__ import DIURNAL_CYCLE, DTIME

    with computation(FORWARD), interval(0, 1):
        member = 0
        while member < 16:  # 16 should come from config file
            xff_mid[0, 0][member] = 0.0
            member += 1
        f_dicycle_modified: FloatFieldIJ = 0.0

    with computation(FORWARD), interval(...):
        if error_code[0, 0][plume] == 0:
            k_x[0, 0][0] = (
                cloud_workfunction_0_modified - cloud_workfunction_1
            ) / arbitrary_numerical_parameter

            if k_x[0, 0][0] <= 0.0 and k_x[0, 0][0] > -0.1 * arbitrary_numerical_parameter:
                k_x[0, 0][0] = -0.1 * arbitrary_numerical_parameter

            if k_x[0, 0][0] > 0.0 and k_x[0, 0][0] < 0.01:
                k_x[0, 0][0] = 0.01

            if DIURNAL_CYCLE == 1 or DIURNAL_CYCLE == 6 or DIURNAL_CYCLE == 0:
                xff_dicycle: FloatFieldIJ = cloud_workfunction_1_pbl / tau_ecmwf
                f_dicycle_modified = max(0.0, -xff_dicycle / k_x[0, 0][0])

            if DIURNAL_CYCLE == 5:
                xff_ens1: FloatFieldIJ = max(0.0, (cloud_workfunction_1 - cloud_workfunction_0) / DTIME)
                mf_ens1: FloatFieldIJ = 0.0

                if k_x[0, 0][0] < 0.0:
                    mf_ens1 = max(0.0, -xff_ens1 / k_x[0, 0][0])

                f_dicycle_modified = -(
                    max(0.0, -(cloud_workfunction_1_pbl - cloud_workfunction_0) / DTIME / k_x[0, 0][0])
                    - mf_ens1
                )

            if k_x[0, 0][0] < 0.0:
                xff_mid[0, 0][2] = max(0.0, -(cloud_workfunction_1 / tau_ecmwf) / k_x[0, 0][0])
                xff_mid[0, 0][3] = f_dicycle_modified

    with computation(FORWARD), interval(...):
        if error_code[0, 0][plume] == 0:
            if updraft_origin_level[0, 0][plume] < pbl_level + 1:
                blqe: FloatFieldIJ = 0.0

                lev = 0
                while lev >= 0 and lev <= updraft_lfc_level[0, 0][plume]:
                    blqe = (
                        blqe
                        + 100.0
                        * moist_static_energy.at(K=lev)
                        * (
                            p_cloud_levels_forced.at(K=lev, ddim=[plume])
                            - p_cloud_levels_forced.at(K=lev + 1, ddim=[plume])
                        )
                        / constants.MAPL_GRAV
                    )
                    lev += 1

                trash: FloatFieldIJ = max(
                    (
                        cloud_moist_static_energy_forced.at(K=updraft_lfc_level)
                        - env_moist_static_energy_cloud_levels_forced.at(K=updraft_lfc_level)
                    ),
                    1.0e1,
                )
                xff_mid[0, 0][1] = max(0.0, blqe / trash)

            xff_mid[0, 0][0] = 0.03 * convective_scale_velocity

    with computation(FORWARD), interval(0, 1):
        if DIURNAL_CYCLE > 0 and ichoice == 4:
            f_dicycle_modified = 0.0

    with computation(FORWARD), interval(...):
        if ichoice == 5:
            if DIURNAL_CYCLE > 0:
                if error_code[0, 0][plume] == 0:
                    xff_mid[0, 0][4] = 0.05 * xff_mid[0, 0][3]
                    f_dicycle_modified = 0.0

            elif DIURNAL_CYCLE == 0:
                if error_code[0, 0][plume] == 0:
                    xff_mid[0, 0][4] = 0.05 * xff_mid[0, 0][3]
                    f_dicycle_modified = 0.0


def update_condensate_to_fall(
    error_code: IntFieldIJ_Plume,
    plume: Int,
    effective_condensate_to_fall_forced: FloatField_Plume,
    condensate_to_fall_forced: FloatField_Plume,
    epsilon_forced: FloatFieldIJ_Plume,
    evaporate_in_downdraft_forced: FloatField_Plume,
):
    with computation(FORWARD), interval(0, -1):
        if error_code[0, 0][plume] == 0:
            effective_condensate_to_fall_forced[0, 0, 0][plume] = (
                condensate_to_fall_forced[0, 0, 0][plume]
                + epsilon_forced[0, 0][plume] * evaporate_in_downdraft_forced[0, 0, 0][plume]
            )


class CloudBaseMassFlux:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GF2020Config,
        cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
    ):
        self.config = config
        self.cumulus_parameterization_config = cumulus_parameterization_config

        self._cup_forcing_ens_3d = stencil_factory.from_dims_halo(
            func=cup_forcing_ens_3d,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "DIURNAL_CYCLE": cumulus_parameterization_config.DIURNAL_CYCLE,
                "DTIME": cumulus_parameterization_config.DTIME,
                "ZERO_DIFF": cumulus_parameterization_config.ZERO_DIFF,
                # "ENSEMBLE_MEMBERS": cumulus_parameterization_config.ENSEMBLE_MEMBERS,
            },
        )

        self._cup_forcing_ens_3d_mid = stencil_factory.from_dims_halo(
            func=cup_forcing_ens_3d_mid,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "DIURNAL_CYCLE": cumulus_parameterization_config.DIURNAL_CYCLE,
                "DTIME": cumulus_parameterization_config.DTIME,
                # "ENSEMBLE_MEMBERS": cumulus_parameterization_config.ENSEMBLE_MEMBERS,
            },
        )

        self._update_condensate_to_fall = stencil_factory.from_dims_halo(
            func=update_condensate_to_fall,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        if self.cumulus_parameterization_config.DIURNAL_CYCLE != 1:
            raise NotImplementedError("Expected DIURNAL_CYCLE == 1!!")

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        if plume_dependent_constants.PLUME_INDEX == 2:

            self._cup_forcing_ens_3d(
                xf_ens=locals.xf_ens,
                error_code=state.output.error_code,
                error_code2=state.output.error_code2,
                error_code3=state.output.error_code3,
                plume=plume_dependent_constants.PLUME_INDEX,
                cloud_workfunction_1=locals.cloud_workfunction_1,
                cloud_workfunction_0=locals.cloud_workfunction_0,
                cloud_workfunction_1_pbl=locals.cloud_work_function_1_pbl,
                xff_ens3=locals.xff_ens3,
                updraft_lfc_level=state.output.updraft_lfc_level,
                omega=state.input_output.omega,
                moisture_convergence=locals.moisture_convergence,
                vapor_forced=locals.vapor_forced,
                tau_ecmwf=locals.tau_ecmwf,
                ichoice=locals.ichoice,
                cloud_workfunction_0_modified=locals.cloud_workfunction_0_modified,
                arbitrary_numerical_parameter=locals.arbitrary_numerical_parameter,
                ocean_fraction=locals.ocean_fraction,
                precipitation_ensemble=locals.precipitation_ensemble,
                f_dicycle_modified=locals.f_dicycle_modified,
            )

        if plume_dependent_constants.PLUME_INDEX == 1:

            self._cup_forcing_ens_3d_mid(
                error_code=state.output.error_code,
                plume=plume_dependent_constants.PLUME_INDEX,
                k_x=locals.k_x,
                cloud_workfunction_1=locals.cloud_workfunction_1,
                cloud_workfunction_0=locals.cloud_workfunction_0,
                cloud_workfunction_0_modified=locals.cloud_workfunction_0_modified,
                cloud_workfunction_1_pbl=locals.cloud_work_function_1_pbl,
                arbitrary_numerical_parameter=locals.arbitrary_numerical_parameter,
                tau_ecmwf=locals.tau_ecmwf,
                f_dicycle_modified=locals.f_dicycle_modified,
                xff_mid=locals.xff_mid,
                updraft_origin_level=state.output.updraft_origin_level,
                pbl_level=state.input_output.pbl_level,
                updraft_lfc_level=state.output.updraft_lfc_level,
                moist_static_energy=locals.moist_static_energy,
                p_cloud_levels_forced=state.output.p_cloud_levels_forced,
                cloud_moist_static_energy_forced=locals.cloud_moist_static_energy_forced,
                env_moist_static_energy_cloud_levels_forced=locals.environment_moist_static_energy_cloud_levels_forced,
                convective_scale_velocity=state.input_output.convective_scale_velocity,
                ichoice=locals.ichoice,
            )

        if plume_dependent_constants.PLUME_INDEX == 0:
            raise NotImplementedError("Plume == 0 (shallow) not expected!!")

        self._update_condensate_to_fall(
            error_code=state.output.error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
            effective_condensate_to_fall_forced=state.output.effective_condensate_to_fall_forced,
            condensate_to_fall_forced=state.output.condensate_to_fall_forced,
            epsilon_forced=state.output.epsilon_forced,
            evaporate_in_downdraft_forced=state.output.evaporate_in_downdraft_forced,
        )
