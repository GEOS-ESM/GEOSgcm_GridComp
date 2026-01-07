from ndsl import StencilFactory, QuantityFactory, Quantity, Local
from ndsl.dsl.gt4py import computation, FORWARD, interval, PARALLEL, K
from ndsl.dsl.typing import FloatField, FloatFieldIJ, IntFieldIJ, Int
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import (
    FloatField_Plume,
    FloatFieldIJ_Plume,
    IntFieldIJ_Plume,
)
import pyMoist.constants as constants
from ndsl.stencils.column_operations import column_max
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from pyMoist.convection.GF_2020.cumulus_parameterization.shared_stencils import tridiag


def zero_tendencies(
    del_u_cloud_ensemble: FloatField,
    del_v_cloud_ensemble: FloatField,
    del_moist_static_energy_cloud_ensemble: FloatField,
    del_t_cloud_ensemble: FloatField,
    del_vapor_cloud_ensemble: FloatField,
    del_cloud_liquid_cloud_ensemble: FloatField,
    del_buoyancy_cloud_ensemble: FloatField,
    t_tendency_from_environmental_subsidence: FloatField,
    moist_static_energy_tendency_from_environmental_subsidence: FloatField,
    vapor_tendency_from_environmental_subsidence: FloatField,
):
    with computation(PARALLEL), interval(...):
        del_u_cloud_ensemble = 0.0
        del_v_cloud_ensemble = 0.0
        del_moist_static_energy_cloud_ensemble = 0.0
        del_t_cloud_ensemble = 0.0
        del_vapor_cloud_ensemble = 0.0
        del_cloud_liquid_cloud_ensemble = 0.0
        del_buoyancy_cloud_ensemble = 0.0
        t_tendency_from_environmental_subsidence = 0.0
        moist_static_energy_tendency_from_environmental_subsidence = 0.0
        vapor_tendency_from_environmental_subsidence = 0.0


def convective_transport_of_momentum(
    error_code: IntFieldIJ_Plume,
    cloud_top_level: IntFieldIJ_Plume,
    p_cloud_levels_forced: FloatField_Plume,
    normalized_massflux_updraft_forced: FloatField_Plume,
    normalized_massflux_downdraft_forced: FloatField_Plume,
    environment_massflux: FloatField,
    u: FloatField,
    v: FloatField,
    u_cloud_levels: FloatField,
    v_cloud_levels: FloatField,
    u_c: FloatField,
    v_c: FloatField,
    u_c_downdraft: FloatField,
    v_c_downdraft: FloatField,
    del_u_cloud_ensemble: FloatField,
    del_v_cloud_ensemble: FloatField,
    epsilon_forced: FloatFieldIJ_Plume,
    fp: FloatField,
    fm: FloatField,
    aa: FloatField,
    bb: FloatField,
    cc: FloatField,
    ddu: FloatField,
    ddv: FloatField,
    plume: Int,
):
    from __externals__ import VERTICAL_DISCRETIZATION_OPTION, ALP1, DTIME

    with computation(FORWARD), interval(0, 1):
        # prepare bounds for subsequent computation
        upper_bound: IntFieldIJ = cloud_top_level[0, 0][plume] + 1
        upper_bound_up_1: IntFieldIJ = upper_bound + 1

    with computation(FORWARD), interval(0, upper_bound):
        if (
            error_code[0, 0][plume] == 0 and VERTICAL_DISCRETIZATION_OPTION == 1 and ALP1 == 0.0
        ):  # fully time explicit
            dp = 100.0 * (p_cloud_levels_forced[0, 0, 0][plume] - p_cloud_levels_forced[0, 0, 1][plume])

            del_u_cloud_ensemble = (
                -(
                    normalized_massflux_updraft_forced[0, 0, 1][plume]
                    * (u_c[0, 0, 1] - u_cloud_levels[0, 0, 1])
                    - normalized_massflux_updraft_forced[0, 0, 0][plume] * (u_c - u_cloud_levels)
                )
                * constants.MAPL_GRAV
                / dp
                + (
                    normalized_massflux_downdraft_forced[0, 0, 1][plume]
                    * (u_c_downdraft[0, 0, 1] - u_cloud_levels[0, 0, 1])
                    - normalized_massflux_downdraft_forced[0, 0, 0][plume] * (u_c_downdraft - u_cloud_levels)
                )
                * constants.MAPL_GRAV
                / dp
                * epsilon_forced[0, 0][plume]
            )

            del_v_cloud_ensemble = (
                -(
                    normalized_massflux_updraft_forced[0, 0, 1][plume]
                    * (v_c_downdraft[0, 0, 1] - v_cloud_levels[0, 0, 1])
                    - normalized_massflux_updraft_forced[0, 0, 0][plume] * (v_c_downdraft - v_cloud_levels)
                )
                * constants.MAPL_GRAV
                / dp
                + (
                    normalized_massflux_downdraft_forced[0, 0, 1][plume]
                    * (v_c_downdraft[0, 0, 1] - v_cloud_levels[0, 0, 1])
                    - normalized_massflux_downdraft_forced[0, 0, 0][plume] * (v_c_downdraft - v_cloud_levels)
                )
                * constants.MAPL_GRAV
                / dp
                * epsilon_forced[0, 0][plume]
            )

    with computation(PARALLEL), interval(0, upper_bound_up_1):
        if (
            error_code[0, 0][plume] == 0 and VERTICAL_DISCRETIZATION_OPTION == 1 and ALP1 > 0.0
        ):  # time alp0*explict + alp1*implicit + upstream
            alp0 = 1.0 - ALP1
            fp = 0.5 * (environment_massflux + abs(environment_massflux))
            fm = 0.5 * (environment_massflux - abs(environment_massflux))

    with computation(FORWARD), interval(0, upper_bound):
        if error_code[0, 0][plume] == 0 and VERTICAL_DISCRETIZATION_OPTION == 1 and ALP1 > 0.0:
            dp = 100.0 * (p_cloud_levels_forced[0, 0, 0][plume] - p_cloud_levels_forced[0, 0, 1][plume])

            beta1 = DTIME * constants.MAPL_GRAV / dp
            aa = ALP1 * beta1 * fm
            bb = 1.0 + ALP1 * beta1 * (fp - fm[0, 0, 1])
            cc = -ALP1 * beta1 * fp[0, 0, 1]

            ddu = (
                u
                - (
                    normalized_massflux_updraft_forced[0, 0, 1][plume] * u_c[0, 0, 1]
                    - normalized_massflux_updraft_forced[0, 0, 0][plume] * u_c[0, 0, 0]
                )
                * beta1
                + (
                    normalized_massflux_downdraft_forced[0, 0, 1][plume] * u_c_downdraft[0, 0, 1]
                    - normalized_massflux_downdraft_forced[0, 0, 1][plume] * u_c_downdraft
                )
                * beta1
                * epsilon_forced[0, 0][plume]
            )

            _, max_index = column_max(u, 0, K - 1)
            ddu = ddu + alp0 * beta1 * (
                -fm * u.at(K=max_index) + (fm[0, 0, 1] - fp) * u + fp[0, 0, 1] * u[0, 0, 1]
            )

            ddv = (
                v
                - (
                    normalized_massflux_updraft_forced[0, 0, 1][plume] * v_c[0, 0, 1]
                    - normalized_massflux_updraft_forced[0, 0, 0][plume] * v_c
                )
                * beta1
                + (
                    normalized_massflux_downdraft_forced[0, 0, 1][plume] * v_c_downdraft[0, 0, 1]
                    - normalized_massflux_downdraft_forced[0, 0, 0][plume] * v_c_downdraft
                )
                * beta1
                * epsilon_forced[0, 0][plume]
            )

            _, max_index = column_max(u, 0, K - 1)
            ddv = ddv + alp0 * beta1 * (
                -fm * v.at(K=max_index) + (fm[0, 0, 1] - fp) * v + fp[0, 0, 1] * v[0, 0, 1]
            )


def update_after_tridiag(
    error_code: IntFieldIJ_Plume,
    cloud_top_level: IntFieldIJ_Plume,
    in_field: FloatField,
    out_field: FloatField,
    wind: FloatField,
    plume: Int,
):
    from __externals__ import DTIME

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0 and K <= cloud_top_level[0, 0][plume]:
            out_field = (in_field - wind) / DTIME


class VerticalDiscretization:
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

        # initalize local fields
        self._fp: Local = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._fm: Local = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._aa: Local = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._bb: Local = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._cc: Local = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._ddu: Local = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._ddv: Local = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        self._zero_tendencies = stencil_factory.from_dims_halo(
            func=zero_tendencies,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._convective_transport_of_momentum = stencil_factory.from_dims_halo(
            func=convective_transport_of_momentum,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "VERTICAL_DISCRETIZATION_OPTION": cumulus_parameterization_config.VERTICAL_DISCRETIZATION_OPTION,
                "ALP1": cumulus_parameterization_config.ALP1,
                "DTIME": cumulus_parameterization_config.DTIME,
            },
        )

        self._tridiag = stencil_factory.from_dims_halo(
            func=tridiag,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._update_after_tridiag = stencil_factory.from_dims_halo(
            func=update_after_tridiag,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"DTIME": cumulus_parameterization_config.DTIME},
        )

    def __call__(
        self,
        error_code: Quantity,
        cloud_top_level: Quantity,
        p_cloud_levels_forced: Quantity,
        normalized_massflux_updraft_forced: Quantity,
        normalized_massflux_downdraft_forced: Quantity,
        environment_massflux: Quantity,
        u: Quantity,
        v: Quantity,
        u_cloud_levels: Quantity,
        v_cloud_levels: Quantity,
        u_c: Quantity,
        v_c: Quantity,
        u_c_downdraft: Quantity,
        v_c_downdraft: Quantity,
        del_u_cloud_ensemble: Quantity,
        del_v_cloud_ensemble: Quantity,
        del_moist_static_energy_cloud_ensemble: Quantity,
        del_t_cloud_ensemble: Quantity,
        del_vapor_cloud_ensemble: Quantity,
        del_cloud_liquid_cloud_ensemble: Quantity,
        del_buoyancy_cloud_ensemble: Quantity,
        t_tendency_from_environmental_subsidence: Quantity,
        moist_static_energy_tendency_from_environmental_subsidence: Quantity,
        vapor_tendency_from_environmental_subsidence: Quantity,
        epsilon_forced: Quantity,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._zero_tendencies(
            del_u_cloud_ensemble=del_u_cloud_ensemble,
            del_v_cloud_ensemble=del_v_cloud_ensemble,
            del_moist_static_energy_cloud_ensemble=del_moist_static_energy_cloud_ensemble,
            del_t_cloud_ensemble=del_t_cloud_ensemble,
            del_vapor_cloud_ensemble=del_vapor_cloud_ensemble,
            del_cloud_liquid_cloud_ensemble=del_cloud_liquid_cloud_ensemble,
            del_buoyancy_cloud_ensemble=del_buoyancy_cloud_ensemble,
            t_tendency_from_environmental_subsidence=t_tendency_from_environmental_subsidence,
            moist_static_energy_tendency_from_environmental_subsidence=moist_static_energy_tendency_from_environmental_subsidence,
            vapor_tendency_from_environmental_subsidence=vapor_tendency_from_environmental_subsidence,
        )

        self._convective_transport_of_momentum(
            error_code=error_code,
            cloud_top_level=cloud_top_level,
            p_cloud_levels_forced=p_cloud_levels_forced,
            normalized_massflux_updraft_forced=normalized_massflux_updraft_forced,
            normalized_massflux_downdraft_forced=normalized_massflux_downdraft_forced,
            environment_massflux=environment_massflux,
            u=u,
            v=v,
            u_cloud_levels=u_cloud_levels,
            v_cloud_levels=v_cloud_levels,
            u_c=u_c,
            v_c=v_c,
            u_c_downdraft=u_c_downdraft,
            v_c_downdraft=v_c_downdraft,
            del_u_cloud_ensemble=del_u_cloud_ensemble,
            del_v_cloud_ensemble=del_v_cloud_ensemble,
            epsilon_forced=epsilon_forced,
            fp=self._fp,
            fm=self._fm,
            aa=self._aa,
            bb=self._bb,
            cc=self._cc,
            ddu=self._ddu,
            ddv=self._ddv,
            plume=plume_dependent_constants.PLUME_INDEX,
        )

        self._tridiag(
            m=cloud_top_level,
            a=self._aa,
            b=self._bb,
            c=self._cc,
            f=self._ddu,
            error_code=error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
        )

        self._update_after_tridiag(
            error_code=error_code,
            cloud_top_level=cloud_top_level,
            in_field=self._ddu,
            out_field=del_u_cloud_ensemble,
            wind=u,
            plume=plume_dependent_constants.PLUME_INDEX,
        )

        self._tridiag(
            m=cloud_top_level,
            a=self._aa,
            b=self._bb,
            c=self._cc,
            f=self._ddv,
            error_code=error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
        )

        self._update_after_tridiag(
            error_code=error_code,
            cloud_top_level=cloud_top_level,
            in_field=self._ddv,
            out_field=del_v_cloud_ensemble,
            wind=v,
            plume=plume_dependent_constants.PLUME_INDEX,
        )
