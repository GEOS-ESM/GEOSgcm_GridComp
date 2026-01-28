from ndsl.dsl.gt4py import computation, interval, PARALLEL, FORWARD, K, exp, BACKWARD
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import (
    IntFieldIJ_Plume,
    FloatField_Tracers,
    FloatFieldIJ_Tracers,
    FloatField_Plume,
    FloatFieldIJ_Plume,
)
from ndsl.dsl.typing import Int, FloatField, FloatFieldIJ, Float
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl import StencilFactory, QuantityFactory, Quantity, Local
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)
import pyMoist.constants as constants
from pyMoist.convection.GF_2020.cumulus_parameterization.shared_stencils import tridiag
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
from pyMoist.convection.GF_2020.cumulus_parameterization.shared_functions import get_cloud_boundary_conditions


def get_incloud_sc_chem_up(
    error_code: IntFieldIJ_Plume,
    plume: Int,
    sc_up_chem: FloatField_Tracers,
    chemistry_tracers_cloud_levels: FloatField_Tracers,
    pw_up_chem: FloatField_Tracers,
    tot_pw_up_chem: FloatFieldIJ_Tracers,
    ocean_fraction: FloatFieldIJ,
    AVERAGE_LAYER_DEPTH: Float,
    updraft_origin_level: IntFieldIJ_Plume,
    po: FloatField,
    sc_b: FloatFieldIJ_Tracers,
    cloud_top_level: IntFieldIJ_Plume,
    normalized_massflux_updraft_forced: FloatField_Plume,
    mass_detrainment_updraft_forced: FloatField_Plume,
    mass_entrainment_updraft_forced: FloatField_Plume,
    chemistry_tracers: FloatField_Tracers,
    geopotential_height_cloud_levels: FloatField,
    vertical_velocity_3d: FloatField,
    CNV_Tracers_fscav: FloatFieldIJ_Tracers,
    p_cloud_levels_forced: FloatField_Plume,
):
    from __externals__ import k_end, BOUNDARY_CONDITION_METHOD, USE_TRACER_SCAVEN, NUMBER_OF_TRACERS

    with computation(FORWARD), interval(...):
        n = 0
        while n < NUMBER_OF_TRACERS:
            sc_up_chem[0, 0, 0][n] = chemistry_tracers_cloud_levels[0, 0, 0][n]
            pw_up_chem[0, 0, 0][n] = 0.0
            tot_pw_up_chem[0, 0][n] = 0.0
            n += 1

    with computation(PARALLEL), interval(...):
        # make garbage field so the get_cloud_boundary_conditions call does not break
        # this is never touched so long as compute_perturbation=False
        dummy_field_no_read = 0.0
        chemistry_tracers_cloud_levels_temp = 0.0

    # with computation(FORWARD), interval(...):
    #     if error_code[0, 0][plume] == 0:
    #         n = 0
    #         while n < NUMBER_OF_TRACERS:
    #             chemistry_tracers_cloud_levels_temp = chemistry_tracers_cloud_levels[0, 0, 0][n]

    #             sc_b[0, 0][n] = get_cloud_boundary_conditions(
    #                 field=chemistry_tracers_cloud_levels_temp,
    #                 scalar_perturbation=0,
    #                 p=po,
    #                 updraft_origin_level=updraft_origin_level[0, 0][plume],
    #                 ocean_fraction=ocean_fraction,
    #                 BOUNDARY_CONDITION_METHOD=BOUNDARY_CONDITION_METHOD,
    #                 AVERAGE_LAYER_DEPTH=AVERAGE_LAYER_DEPTH,
    #                 k_end=k_end,
    #                 compute_perturbation=False,
    #                 perturbation_field=dummy_field_no_read,
    #             )
    #             n += 1

    # with computation(FORWARD), interval(...):
    #     if error_code[0, 0][plume] == 0:
    #         n = 0
    #         while n < NUMBER_OF_TRACERS:
    #             if K <= updraft_origin_level[0, 0][plume]:
    #                 sc_up_chem[0, 0, 0][n] = sc_b[0, 0][n]
    #             n += 1

    with computation(FORWARD), interval(...):
        if error_code[0, 0][plume] == 0:
            if K >= updraft_origin_level[0, 0][plume] + 1 and K <= cloud_top_level[0, 0][plume] + 1:
                XZZ: FloatFieldIJ = normalized_massflux_updraft_forced.at(K=K - 1, ddim=[plume])
                XZD: FloatFieldIJ = 0.5 * mass_detrainment_updraft_forced.at(K=K - 1, ddim=[plume])
                XZE: FloatFieldIJ = mass_entrainment_updraft_forced.at(K=K - 1, ddim=[plume])
                denom: FloatFieldIJ = XZZ - XZD + XZE

                n = 0
                while n < NUMBER_OF_TRACERS:
                    if denom > 0.0:
                        sc_up_chem[0, 0, 0][n] = (
                            sc_up_chem.at(K=K - 1, ddim=[n]) * XZZ
                            - sc_up_chem.at(K=K - 1, ddim=[n]) * XZD
                            + chemistry_tracers.at(K=K - 1, ddim=[n]) * XZE
                        ) / denom
                    else:
                        sc_up_chem[0, 0, 0][n] = sc_up_chem.at(K=K - 1, ddim=[n])
                    n += 1

                if USE_TRACER_SCAVEN != 0 or plume != cumulus_parameterization_constants.SHALLOW:
                    dz = geopotential_height_cloud_levels - geopotential_height_cloud_levels.at(K=K - 1)

                    w_upd: FloatFieldIJ = vertical_velocity_3d

                    tracer = 0
                    while tracer < NUMBER_OF_TRACERS:
                        if CNV_Tracers_fscav[0, 0][tracer] > 1.0e-6:
                            if USE_TRACER_SCAVEN == 1:
                                pw_up_chem[0, 0, 0][tracer] = max(
                                    0.0,
                                    sc_up_chem[0, 0, 0][tracer]
                                    * (1.0 - exp(-CNV_Tracers_fscav[0, 0][tracer] * (dz / 1000.0))),
                                )

                            sc_up_chem[0, 0, 0][tracer] = (
                                sc_up_chem[0, 0, 0][tracer] - pw_up_chem[0, 0, 0][tracer]
                            )

                        tracer += 1

                dp = 100.0 * (p_cloud_levels_forced[0, 0, 0][plume] - p_cloud_levels_forced[0, 0, 1][plume])

                n = 0
                while n < NUMBER_OF_TRACERS:
                    tot_pw_up_chem[0, 0][n] = (
                        tot_pw_up_chem[0, 0][n] + pw_up_chem[0, 0, 0][n] * dp / constants.MAPL_GRAV
                    )
                    n += 1


def get_incloud_sc_chem_dd(
    error_code: IntFieldIJ_Plume,
    plume: Int,
    sc_dn: FloatField_Tracers,
    pw_dn: FloatField_Tracers,
    chemistry_tracers: FloatField_Tracers,
    tot_pw_dn_chem: FloatFieldIJ_Tracers,
    tot_pw_up_chem: FloatFieldIJ_Tracers,
    total_normalized_integrated_evaporate_forced: FloatFieldIJ_Plume,
    total_normalized_integrated_condensate_forced: FloatFieldIJ_Plume,
    downdraft_origin_level: IntFieldIJ_Plume,
    evaporate_in_downdraft_forced: FloatField_Plume,
    p_cloud_levels_forced: FloatField_Plume,
    chemistry_tracers_cloud_levels: FloatField_Tracers,
    normalized_massflux_downdraft_forced: FloatField_Plume,
    mass_detrainment_downdraft_forced: FloatField_Plume,
    mass_entrainment_downdraft_forced: FloatField_Plume,
):
    from __externals__ import NUMBER_OF_TRACERS, USE_TRACER_EVAP

    with computation(FORWARD), interval(...):
        n = 0
        while n < NUMBER_OF_TRACERS:
            sc_dn[0, 0, 0][n] = 0.0
            pw_dn[0, 0, 0][n] = 0.0
            tot_pw_dn_chem[0, 0][n] = 0.0
            n += 1

    with computation(FORWARD), interval(0, 1):
        if plume != cumulus_parameterization_constants.SHALLOW:
            if error_code[0, 0][plume] == 0:
                frac_evap: FloatFieldIJ = -total_normalized_integrated_evaporate_forced[0, 0][plume] / (
                    1.0e-16 + total_normalized_integrated_condensate_forced[0, 0][plume]
                )
                lev = downdraft_origin_level[0, 0][plume]

                pwdper: FloatFieldIJ = (
                    evaporate_in_downdraft_forced.at(K=lev, ddim=[plume])
                    / (1.0e-16 + total_normalized_integrated_evaporate_forced[0, 0][plume])
                    * frac_evap
                )

                if USE_TRACER_EVAP == 0:
                    pwdper = 0.0

                dp: FloatFieldIJ = 100.0 * (
                    p_cloud_levels_forced.at(K=lev, ddim=[plume])
                    - p_cloud_levels_forced.at(K=lev + 1, ddim=[plume])
                )

                n = 0
                while n < NUMBER_OF_TRACERS:
                    sc_dn[0, 0, lev][n] = chemistry_tracers_cloud_levels[0, 0, lev][n]
                    pw_dn[0, 0, lev][n] = -pwdper * tot_pw_up_chem[0, 0][n] * constants.MAPL_GRAV / dp
                    sc_dn[0, 0, lev][n] = sc_dn[0, 0, lev][n] - pw_dn[0, 0, lev][n]
                    tot_pw_dn_chem[0, 0][n] = (
                        tot_pw_dn_chem[0, 0][n] + pw_dn[0, 0, lev][n] * dp / constants.MAPL_GRAV
                    )
                    n += 1

    with computation(BACKWARD), interval(...):
        if plume != cumulus_parameterization_constants.SHALLOW:
            if error_code[0, 0][plume] == 0:
                if K <= downdraft_origin_level[0, 0][plume] - 1:
                    XZZ: FloatFieldIJ = normalized_massflux_downdraft_forced[0, 0, 1][plume]
                    XZD: FloatFieldIJ = 0.5 * mass_detrainment_downdraft_forced[0, 0, 0][plume]
                    XZE: FloatFieldIJ = mass_entrainment_downdraft_forced[0, 0, 0][plume]

                    denom: FloatFieldIJ = XZZ - XZD + XZE

                    n = 0
                    while n < NUMBER_OF_TRACERS:
                        if denom > 0.0:
                            sc_dn[0, 0, 0][n] = (
                                sc_dn.at(K=K + 1, ddim=[n]) * XZZ
                                - sc_dn.at(K=K + 1, ddim=[n]) * XZD
                                + chemistry_tracers[0, 0, 0][n] * XZE
                            ) / denom
                        else:
                            sc_dn[0, 0, 0][n] = sc_dn.at(K=K + 1, ddim=[n])
                        n += 1

                    if USE_TRACER_EVAP != 0:

                        dp = 100.0 * (
                            p_cloud_levels_forced[0, 0, 0][plume]
                            - p_cloud_levels_forced.at(K=K + 1, ddim=[plume])
                        )

                        pwdper = evaporate_in_downdraft_forced[0, 0, 0][plume] / (
                            1.0e-16 + total_normalized_integrated_evaporate_forced[0, 0][plume]
                        )

                        pwdper = pwdper * frac_evap

                        pwdper = min(1.0, max(pwdper, 0.0))

                        n = 0
                        while n < NUMBER_OF_TRACERS:
                            pw_dn[0, 0, 0][n] = -pwdper * tot_pw_up_chem[0, 0][n] * constants.MAPL_GRAV / dp

                            sc_dn[0, 0, 0][n] = sc_dn[0, 0, 0][n] - pw_dn[0, 0, 0][n]

                            tot_pw_dn_chem[0, 0][n] = (
                                tot_pw_dn_chem[0, 0][n] + pw_dn[0, 0, 0][n] * dp / constants.MAPL_GRAV
                            )

                            n += 1


def environment_cloud_levels_chemistry(
    error_code: IntFieldIJ_Plume,
    chemistry_tracers: FloatField_Tracers,
    chemistry_tracers_cloud_levels: FloatField_Tracers,
    plume: Int,
):
    from __externals__ import k_end, NUMBER_OF_TRACERS

    with computation(FORWARD), interval(0, 1):
        # set up internal constants
        # NOTE only cloud_level_option = 2 has been tested
        cloud_level_option = 2

    with computation(FORWARD), interval(1, -1):
        if error_code[0, 0][plume] == 0 and cloud_level_option == 1:
            tracer = 0
            while tracer < NUMBER_OF_TRACERS:
                chemistry_tracers_cloud_levels[0, 0, 0][tracer] = (
                    0.5 * chemistry_tracers[0, 0, -1][tracer] + chemistry_tracers[0, 0, 0][tracer]
                )
                tracer += 1

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0 and cloud_level_option == 1:
            tracer = 0
            while tracer < NUMBER_OF_TRACERS:
                chemistry_tracers_cloud_levels[0, 0, 0][tracer] = chemistry_tracers[0, 0, 0][tracer]
                tracer += 1

    with computation(FORWARD), interval(-1, None):
        if error_code[0, 0][plume] == 0 and cloud_level_option == 1:
            tracer = 0
            while tracer < NUMBER_OF_TRACERS:
                chemistry_tracers_cloud_levels[0, 0, 0][tracer] = chemistry_tracers[0, 0, 0][tracer]
                tracer += 1

    with computation(PARALLEL), interval(0, -1):
        if error_code[0, 0][plume] == 0 and cloud_level_option != 1:
            tracer = 0
            while tracer < NUMBER_OF_TRACERS:
                chemistry_tracers_cloud_levels[0, 0, 0][tracer] = chemistry_tracers[0, 0, 0][tracer]
                tracer += 1


def determine_vertical_transport1(
    error_code: IntFieldIJ_Plume,
    plume: Int,
    cloud_top_level: IntFieldIJ_Plume,
    environment_massflux: FloatField,
    p_cloud_levels_forced: FloatField_Plume,
    ddtr: FloatField_Tracers,
    chemistry_tracers: FloatField_Tracers,
    sc_up_chem: FloatField_Tracers,
    normalized_massflux_updraft_forced: FloatField_Plume,
    normalized_massflux_downdraft_forced: FloatField_Plume,
    epsilon_forced: FloatFieldIJ_Plume,
    sc_dn: FloatField_Tracers,
    out_chem: FloatField_Tracers,
    pw_dn: FloatField_Tracers,
    pw_up_chem: FloatField_Tracers,
    aa: FloatField,
    bb: FloatField,
    cc: FloatField,
):
    from __externals__ import (
        ALP1,
        USE_FLUX_FORM,
        DTIME,
        NUMBER_OF_TRACERS,
        USE_TRACER_EVAP,
        USE_TRACER_SCAVEN,
    )

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0:
            if USE_FLUX_FORM == 1 and ALP1 > 0.0:
                alp0 = 1.0 - ALP1
                if K <= cloud_top_level[0, 0][plume] + 1:
                    fp = 0.5 * (environment_massflux + abs(environment_massflux))
                    fm = 0.5 * (environment_massflux - abs(environment_massflux))

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0:
            if USE_FLUX_FORM == 1 and ALP1 > 0.0:
                if K <= cloud_top_level[0, 0][plume]:
                    dp = 100.0 * (
                        p_cloud_levels_forced[0, 0, 0][plume] - p_cloud_levels_forced[0, 0, 1][plume]
                    )
                    beta1 = DTIME * constants.MAPL_GRAV / dp
                    aa = ALP1 * beta1 * fm
                    bb = 1.0 + ALP1 * beta1 * (fp - fm[0, 0, 1])
                    cc = -ALP1 * beta1 * fp[0, 0, 1]

                    n = 0
                    while n < NUMBER_OF_TRACERS:
                        ddtr[0, 0, 0][n] = (
                            chemistry_tracers[0, 0, 0][n]
                            - (
                                normalized_massflux_updraft_forced[0, 0, 1][plume] * sc_up_chem[0, 0, 1][n]
                                - normalized_massflux_updraft_forced[0, 0, 0][plume] * sc_up_chem[0, 0, 0][n]
                            )
                            * beta1
                            + (
                                normalized_massflux_downdraft_forced[0, 0, 1][plume] * sc_dn[0, 0, 1][n]
                                - normalized_massflux_downdraft_forced[0, 0, 0][plume] * sc_dn[0, 0, 0][n]
                            )
                            * beta1
                            * epsilon_forced[0, 0][plume]
                        )
                        n += 1

                    if USE_TRACER_EVAP == 1 and plume != cumulus_parameterization_constants.SHALLOW:
                        n = 0
                        while n < NUMBER_OF_TRACERS:
                            out_chem[0, 0, 0][n] = (
                                out_chem[0, 0, 0][n]
                                - 0.5
                                * epsilon_forced[0, 0][plume]
                                * (
                                    normalized_massflux_downdraft_forced[0, 0, 0][plume] * pw_dn[0, 0, 0][n]
                                    + normalized_massflux_downdraft_forced[0, 0, 1][plume] * pw_dn[0, 0, 1][n]
                                )
                                * beta1
                            )
                            n += 1

                    if USE_TRACER_SCAVEN > 0 and plume != cumulus_parameterization_constants.SHALLOW:
                        n = 0
                        while n < NUMBER_OF_TRACERS:
                            out_chem[0, 0, 0][n] = (
                                out_chem[0, 0, 0][n]
                                - 0.5
                                * (
                                    normalized_massflux_updraft_forced[0, 0, 0][plume]
                                    * pw_up_chem[0, 0, 0][n]
                                    + normalized_massflux_updraft_forced[0, 0, 1][plume]
                                    * pw_up_chem[0, 0, 1][n]
                                )
                                * beta1
                            )
                            n += 1

                    n = 0
                    while n < NUMBER_OF_TRACERS:
                        ddtr[0, 0, 0][n] = (
                            ddtr[0, 0, 0][n]
                            + out_chem[0, 0, 0][n]
                            + alp0 * beta1 * (-fm * chemistry_tracers.at(K=max(0, K - 1), ddim=[n]))
                            + (fm[0, 0, 1] - fp) * chemistry_tracers[0, 0, 0][n]
                            + fp[0, 0, 1] * chemistry_tracers[0, 0, 1][n]
                        )
                        n += 1


# def determine_vertical_transport2(
#     error_code: IntFieldIJ_Plume,
#     plume: Int,
#     cloud_top_level: IntFieldIJ_Plume,
# ):
#     from __externals__ import (
#         ALP1,
#         USE_FLUX_FORM,
#         DTIME,
#         NUMBER_OF_TRACERS,
#     )

#     with computation(PARALLEL), interval(...):
#         if error_code[0, 0][plume] == 0:
#             if USE_FLUX_FORM == 1 and ALP1 > 0.0:
#                 if K <= cloud_top_level[0, 0][plume]:
#                     n = 0
#                     while n < NUMBER_OF_TRACERS:
#                         out_chem[0, 0, 0][n] = (ddtr[0, 0, 0][n] - chemistry_tracers[0, 0, 0][n]) / DTIME
#                         n += 1

# with computation(FORWARD), interval(...):
#     if error_code[0, 0][plume] == 0:
#         n = 0
#         while n < NUMBER_OF_TRACERS:
#             trash_[0, 0][n] = 0.0
#             trash2_[0, 0][n] = 0.0
#             evap_[0, 0][n] = 0.0
#             wetdep_[0, 0][n] = 0.0
#             residu_[0, 0][n] = 0.0

#             if K <= cloud_top_level[0, 0][plume]:
#                 dp = 100.0 * (
#                     p_cloud_levels_forced[0, 0, 0][plume] - p_cloud_levels_forced[0, 0, 1][plume]
#                 )
#                 evap = (
#                     -0.5
#                     * (
#                         normalized_massflux_downdraft_forced[0, 0, 0][plume] * pw_dn[0, 0, 0][n]
#                         + normalized_massflux_downdraft_forced[0, 0, 1][plume] * pw_dn[0, 0, 1][n]
#                     )
#                     * constants.MAPL_GRAV
#                     / dp
#                     * epsilon_forced[0, 0][plume]
#                 )
#                 wetdep = (
#                     0.5
#                     * (
#                         normalized_massflux_updraft_forced[0, 0, 0][plume] * pw_up_chem[0, 0, 0][n]
#                         + normalized_massflux_updraft_forced[0, 0, 1][plume] * pw_up_chem[0, 0, 1][n]
#                     )
#                     * constants.MAPL_GRAV
#                     / dp
#                 )

#                 evap_[0, 0][n] = evap_[0, 0][n] + evap * dp / constants.MAPL_GRAV
#                 wetdep_[0, 0][n] = wetdep_[0, 0][n] + wetdep * dp / constants.MAPL_GRAV
#                 residu_[0, 0][n] = residu_[0, 0][n] + (wetdep - evap) * dp / constants.MAPL_GRAV

#                 trash_[0, 0][n] = trash_[0, 0][n] + (out_chem[0, 0, 0][n]) * dp / constants.MAPL_GRAV

#                 trash2_[0, 0][n] = (
#                     trash2_[0, 0][n] + chemistry_tracers[0, 0, 0][n] * dp / constants.MAPL_GRAV
#                 )

#             if residu_[0, 0][n] < 0.0:
#                 beta1 = constants.MAPL_GRAV / (
#                     p_cloud_levels_forced.at(K=0, ddim=[plume])
#                     - p_cloud_levels_forced.at(K=cloud_top_level[0, 0][plume] + 1, ddim=[plume])
#                 )
#                 if K <= cloud_top_level[0, 0][plume]:
#                     out_chem[0, 0, 0][n] = out_chem[0, 0, 0][n] + residu_[0, 0][n] * beta1

#             n += 1


class ColdPoolParameterization:
    def __init__(self, cumulus_parameterization_config: GF2020CumulusParameterizationConfig):
        if cumulus_parameterization_config.CONVECTION_TRACER == 1:
            raise NotImplementedError(
                "The ColdPoolParameterization section has not been implemented. You should have been caught"
                "before getting here by the config checker. Beware, something likely failing in the config"
                "checker as well - you may be unknowingly calling other untested/unimplemented sections."
            )

    def __call__(self, *args, **kwds):
        pass


class AtmosphericComposition:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GF2020Config,
        cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
    ):
        # make configuration visible at runtime
        self.config = config
        self.cumulus_parameterization_config = cumulus_parameterization_config

        self._aa: Local = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._bb: Local = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._cc: Local = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._ddtr: Local = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        # construct stencils and functions
        self._environment_cloud_levels_chemistry = stencil_factory.from_dims_halo(
            func=environment_cloud_levels_chemistry,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"NUMBER_OF_TRACERS": config.NUMBER_OF_TRACERS},
        )

        self._get_incloud_sc_chem_up = stencil_factory.from_dims_halo(
            func=get_incloud_sc_chem_up,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "NUMBER_OF_TRACERS": config.NUMBER_OF_TRACERS,
                "BOUNDARY_CONDITION_METHOD": cumulus_parameterization_config.BOUNDARY_CONDITION_METHOD,
                "USE_TRACER_SCAVEN": cumulus_parameterization_config.USE_TRACER_SCAVEN,
            },
        )

        if self.cumulus_parameterization_config.USE_TRACER_SCAVEN != 1:
            raise NotImplementedError(
                "[NDSL] GF2020-->CumulusParameterization-->AtmosphericComposition called with an unimplemented path."
                "This should have been caught at initialization, but somehow you made it here."
                "Choose another option for USE_TRACER_SCAVEN or implement to continue."
            )

        self._determine_vertical_transport1 = stencil_factory.from_dims_halo(
            func=determine_vertical_transport1,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "NUMBER_OF_TRACERS": config.NUMBER_OF_TRACERS,
                "BOUNDARY_CONDITION_METHOD": cumulus_parameterization_config.BOUNDARY_CONDITION_METHOD,
                "USE_TRACER_SCAVEN": cumulus_parameterization_config.USE_TRACER_SCAVEN,
                "ALP1": cumulus_parameterization_config.ALP1,
                "USE_FLUX_FORM": cumulus_parameterization_config.USE_FLUX_FORM,
                "DTIME": cumulus_parameterization_config.DTIME,
                "USE_TRACER_EVAP": cumulus_parameterization_config.USE_TRACER_EVAP,
            },
        )

        if (
            self.cumulus_parameterization_config.USE_FLUX_FORM != 1
            or self.cumulus_parameterization_config.ALP1 != 1
        ):
            raise NotImplementedError(
                "[NDSL] GF2020-->CumulusParameterization-->AtmosphericComposition called with an unimplemented path."
                "This should have been caught at initialization, but somehow you made it here."
                "Choose another option for USE_FlUX_FORM and ALP1 or implement to continue."
            )

        self._get_incloud_sc_chem_dd = stencil_factory.from_dims_halo(
            func=get_incloud_sc_chem_dd,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "NUMBER_OF_TRACERS": config.NUMBER_OF_TRACERS,
                "USE_TRACER_EVAP": cumulus_parameterization_config.USE_TRACER_EVAP,
            },
        )

        self._tridiag = stencil_factory.from_dims_halo(
            func=tridiag,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        # self._determine_vertical_transport2 = stencil_factory.from_dims_halo(
        #     func=determine_vertical_transport2,
        #     compute_dims=[X_DIM, Y_DIM, Z_DIM],
        #     externals={
        #         "NUMBER_OF_TRACERS": config.NUMBER_OF_TRACERS,
        #         "DTIME": cumulus_parameterization_config.DTIME,
        #         "ALP1": cumulus_parameterization_config.ALP1,
        #         "USE_FLUX_FORM": cumulus_parameterization_config.USE_FLUX_FORM,
        #     },
        # )

    def __call__(
        self,
        error_code: Quantity,
        chemistry_tracers: Quantity,
        chemistry_tracers_cloud_levels: Quantity,
        plume_dependent_constants: GF2020PlumeDependentConstants,
        sc_up_chem: Quantity,
        pw_up_chem: Quantity,
        tot_pw_up_chem: Quantity,
        ocean_fraction: Quantity,
        AVERAGE_LAYER_DEPTH: Quantity,
        updraft_origin_level: Quantity,
        po: Quantity,
        sc_b: Quantity,
        cloud_top_level: Quantity,
        normalized_massflux_updraft_forced: Quantity,
        mass_detrainment_updraft_forced: Quantity,
        mass_entrainment_updraft_forced: Quantity,
        geopotential_height_cloud_levels: Quantity,
        vertical_velocity_3d: Quantity,
        CNV_Tracers_fscav: Quantity,
        CNV_Tracers_Vect_hcts: Quantity,
        p_cloud_levels_forced: Quantity,
        sc_dn: Quantity,
        pw_dn: Quantity,
        tot_pw_dn_chem: Quantity,
        total_normalized_integrated_evaporate_forced: Quantity,
        total_normalized_integrated_condensate_forced: Quantity,
        downdraft_origin_level: Quantity,
        evaporate_in_downdraft_forced: Quantity,
        normalized_massflux_downdraft_forced: Quantity,
        mass_detrainment_downdraft_forced: Quantity,
        mass_entrainment_downdraft_forced: Quantity,
        environment_massflux: Quantity,
        ddtr: Quantity,
        epsilon_forced: Quantity,
        out_chem: Quantity,
        trash_: Quantity,
        trash2_: Quantity,
        evap_: Quantity,
        wetdep_: Quantity,
        residu_: Quantity,
    ):

        self._environment_cloud_levels_chemistry(
            error_code=error_code,
            chemistry_tracers=chemistry_tracers,
            chemistry_tracers_cloud_levels=chemistry_tracers_cloud_levels,
            plume=plume_dependent_constants.PLUME_INDEX,
        )

        # NOTE Raise an error if CNV_Tracers_Vect_Hcts > 1.e-6
        # if CNV_Tracers_Vect_hcts.view[:].any() > float(1.0e-6):
        #     raise NotImplementedError(
        #         "[NDSL] GF2020-->CumulusParameterization-->AtmosphericComposition called with an unimplemented path."
        #         "This should have been caught at initialization, but somehow you made it here."
        #     )

        self._get_incloud_sc_chem_up(
            error_code=error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
            sc_up_chem=sc_up_chem,
            chemistry_tracers_cloud_levels=chemistry_tracers_cloud_levels,
            pw_up_chem=pw_up_chem,
            tot_pw_up_chem=tot_pw_up_chem,
            ocean_fraction=ocean_fraction,
            AVERAGE_LAYER_DEPTH=AVERAGE_LAYER_DEPTH,
            updraft_origin_level=updraft_origin_level,
            po=po,
            sc_b=sc_b,
            cloud_top_level=cloud_top_level,
            normalized_massflux_updraft_forced=normalized_massflux_updraft_forced,
            mass_detrainment_updraft_forced=mass_detrainment_updraft_forced,
            mass_entrainment_updraft_forced=mass_entrainment_updraft_forced,
            chemistry_tracers=chemistry_tracers,
            geopotential_height_cloud_levels=geopotential_height_cloud_levels,
            vertical_velocity_3d=vertical_velocity_3d,
            CNV_Tracers_fscav=CNV_Tracers_fscav,
            p_cloud_levels_forced=p_cloud_levels_forced,
        )

        self._get_incloud_sc_chem_dd(
            error_code=error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
            sc_dn=sc_dn,
            pw_dn=pw_dn,
            chemistry_tracers=chemistry_tracers,
            tot_pw_dn_chem=tot_pw_dn_chem,
            tot_pw_up_chem=tot_pw_up_chem,
            total_normalized_integrated_evaporate_forced=total_normalized_integrated_evaporate_forced,
            total_normalized_integrated_condensate_forced=total_normalized_integrated_condensate_forced,
            downdraft_origin_level=downdraft_origin_level,
            evaporate_in_downdraft_forced=evaporate_in_downdraft_forced,
            p_cloud_levels_forced=p_cloud_levels_forced,
            chemistry_tracers_cloud_levels=chemistry_tracers_cloud_levels,
            normalized_massflux_downdraft_forced=normalized_massflux_downdraft_forced,
            mass_detrainment_downdraft_forced=mass_detrainment_downdraft_forced,
            mass_entrainment_downdraft_forced=mass_entrainment_downdraft_forced,
        )

        self._determine_vertical_transport1(
            error_code=error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
            cloud_top_level=cloud_top_level,
            environment_massflux=environment_massflux,
            p_cloud_levels_forced=p_cloud_levels_forced,
            ddtr=self._ddtr,
            chemistry_tracers=chemistry_tracers,
            sc_up_chem=sc_up_chem,
            normalized_massflux_updraft_forced=normalized_massflux_updraft_forced,
            normalized_massflux_downdraft_forced=normalized_massflux_downdraft_forced,
            epsilon_forced=epsilon_forced,
            sc_dn=sc_dn,
            out_chem=out_chem,
            pw_dn=pw_dn,
            pw_up_chem=pw_up_chem,
            aa=self._aa,
            bb=self._bb,
            cc=self._cc,
        )

        self._tridiag(
            m=cloud_top_level,
            a=self._aa,
            b=self._bb,
            c=self._cc,
            f=self._ddtr,
            error_code=error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
        )

        # self._determine_vertical_transport2()
