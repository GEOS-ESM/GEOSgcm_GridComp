from ndsl import StencilFactory, QuantityFactory, Local, Quantity
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import (
    GF2020CumulusParameterizationConfig,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.state import (
    GF2020CumulusParameterizationState,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.locals import (
    GF2020CumulusParameterizationLocals,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import IntFieldIJ_Plume, FloatField_Plume
from ndsl.dsl.typing import FloatField, FloatFieldIJ, Float, IntFieldIJ, Int
from ndsl.dsl.gt4py import computation, FORWARD, PARALLEL, interval, sqrt
from pyMoist.shared_generic_math import sigma
import pyMoist.constants as constants


def set_time_scales(
    error_code: IntFieldIJ_Plume,
    updraft_lfc_level: IntFieldIJ_Plume,
    cloud_top_level: IntFieldIJ_Plume,
    grid_length: FloatFieldIJ,
    ocean_fraction: FloatFieldIJ,
    topography_height_no_negative: FloatFieldIJ,
    geopotential_height_cloud_levels_forced: FloatField,
    u: FloatField,
    v: FloatField,
    vertical_velocity_2d: FloatFieldIJ,
    cape_removal_time_scale: FloatFieldIJ,
    cape_removal_time_scale_from_state: FloatFieldIJ,
    pbl_time_scale: FloatFieldIJ,
    pbl_time_scale_from_state: FloatFieldIJ,
    TAU_CAPE_REMOVAL: Float,
    plume: Int,
):
    from __externals__ import SGS_W_TIMESCALE, DTIME

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            if SGS_W_TIMESCALE == 0:
                # time-scale cape removal from Bechtold et al. 2008
                dz = geopotential_height_cloud_levels_forced.at(
                    K=cloud_top_level[0, 0][plume]
                ) - geopotential_height_cloud_levels_forced.at(K=updraft_lfc_level[0, 0][plume])
            else:
                # time-scale cape removal from Bechtold et al. 2008
                dz = geopotential_height_cloud_levels_forced.at(
                    K=cloud_top_level[0, 0][plume]
                ) - geopotential_height_cloud_levels_forced.at(K=updraft_lfc_level[0, 0][plume])
                cape_removal_time_scale = (
                    3600.0 * (sigma(grid_length))
                    + 10800.0 * (1.0 - sigma(grid_length))
                    + (dz / vertical_velocity_2d)
                )
                cape_removal_time_scale = max(DTIME, cape_removal_time_scale)

            if ocean_fraction > 0.99:  # over water
                umean = 2.0 + sqrt(
                    0.5
                    * (
                        u**2
                        + v**2
                        + u.at(K=updraft_lfc_level[0, 0][plume]) ** 2
                        + v.at(K=updraft_lfc_level[0, 0][plume]) ** 2
                    )
                )
                pbl_time_scale = (
                    geopotential_height_cloud_levels_forced.at(K=updraft_lfc_level[0, 0][plume])
                    - topography_height_no_negative
                ) / umean
            else:  # over land
                pbl_time_scale = (
                    geopotential_height_cloud_levels_forced.at(K=cloud_top_level[0, 0][plume])
                    - geopotential_height_cloud_levels_forced.at(K=updraft_lfc_level[0, 0][plume])
                ) / 3.0  # 3.0 m/s is estimated wmean

    with computation(FORWARD), interval(0, 1):
        cape_removal_time_scale_from_state = cape_removal_time_scale
        pbl_time_scale_from_state = pbl_time_scale


def cloud_workfunction_1_pbl(
    error_code: IntFieldIJ_Plume,
    pbl_level: IntFieldIJ,
    geopotential_height_cloud_levels_forced: FloatField,
    t_old: FloatField,
    t_new: FloatField,
    t_cloud_levels_forced: FloatField,
    vapor_old: FloatField,
    vapor_forced: FloatField,
    cloud_work_function_1_pbl: FloatFieldIJ,
    cloud_work_function_1_fa: FloatFieldIJ,
    plume: Int,
):
    from __externals__ import DTIME

    with computation(FORWARD), interval(0, 1):
        cloud_work_function_1_pbl = 0.0
        cloud_work_function_1_fa = 0.0

        top_bound: IntFieldIJ = pbl_level + 1

    with computation(FORWARD), interval(0, top_bound):
        if error_code[0, 0][plume] == 0:
            dz = (
                geopotential_height_cloud_levels_forced[0, 0, 1] - geopotential_height_cloud_levels_forced
            ) * constants.MAPL_GRAV
            da = dz * (t_new * (1.0 + 0.608 * vapor_forced) - t_old * (1.0 + 0.608 * vapor_old)) / DTIME

            cloud_work_function_1_pbl = cloud_work_function_1_pbl + da  # Units : J K / (kg seg)


def scale_cloud_workfunction_1_pbl(
    error_code: IntFieldIJ_Plume,
    pbl_time_scale: FloatFieldIJ,
    cloud_work_function_1_pbl: FloatFieldIJ,
    T_STAR: Float,
    plume: Int,
):
    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            cloud_work_function_1_pbl = cloud_work_function_1_pbl / T_STAR * pbl_time_scale


class DiurnalCycle:
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

        # construct stencils and functions
        self._set_time_scales = stencil_factory.from_dims_halo(
            func=set_time_scales,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "SGS_W_TIMESCALE": cumulus_parameterization_config.SGS_W_TIMESCALE,
                "DTIME": cumulus_parameterization_config.DTIME,
            },
        )

        self._cloud_workfunction_1_pbl = stencil_factory.from_dims_halo(
            func=cloud_workfunction_1_pbl,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "SGS_W_TIMESCALE": cumulus_parameterization_config.SGS_W_TIMESCALE,
                "DTIME": cumulus_parameterization_config.DTIME,
            },
        )

        self._scale_cloud_workfunction_1_pbl = stencil_factory.from_dims_halo(
            func=scale_cloud_workfunction_1_pbl,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        # initialize internal fields
        self._tau_ecmwf: Local = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        self._tau_bl: Local = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")

    def __call__(
        self,
        error_code: Quantity,
        updraft_lfc_level: Quantity,
        cloud_top_level: Quantity,
        pbl_level: Quantity,
        grid_length: Quantity,
        ocean_fraction: Quantity,
        topography_height_no_negative: Quantity,
        geopotential_height_cloud_levels_forced: Quantity,
        t_old: Quantity,
        t_new: Quantity,
        t_cloud_levels_forced: Quantity,
        vapor_old: Quantity,
        vapor_forced: Quantity,
        u: Quantity,
        v: Quantity,
        vertical_velocity_2d: Quantity,
        cape_removal_time_scale: Quantity,
        cape_removal_time_scale_from_state: Quantity,
        pbl_time_scale: Quantity,
        pbl_time_scale_from_state: Quantity,
        cloud_work_function_1_pbl: Quantity,
        cloud_work_function_1_fa: Quantity,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        # Bechtold et al 2008 time-scale of cape removal
        self._set_time_scales(
            error_code=error_code,
            updraft_lfc_level=updraft_lfc_level,
            cloud_top_level=cloud_top_level,
            grid_length=grid_length,
            ocean_fraction=ocean_fraction,
            topography_height_no_negative=topography_height_no_negative,
            geopotential_height_cloud_levels_forced=geopotential_height_cloud_levels_forced,
            u=u,
            v=v,
            vertical_velocity_2d=vertical_velocity_2d,
            cape_removal_time_scale=cape_removal_time_scale,
            cape_removal_time_scale_from_state=cape_removal_time_scale_from_state,
            pbl_time_scale=pbl_time_scale,
            pbl_time_scale_from_state=pbl_time_scale_from_state,
            TAU_CAPE_REMOVAL=plume_dependent_constants.TAU_CAPE_REMOVAL,
            plume=plume_dependent_constants.PLUME_INDEX,
        )

        if (
            self.cumulus_parameterization_config.DIURNAL_CYCLE == 1
            or self.cumulus_parameterization_config.DIURNAL_CYCLE == 6
        ) or (
            self.cumulus_parameterization_config.DIURNAL_CYCLE == 0
            and plume_dependent_constants.PLUME_INDEX == 1  # mid plume
        ):
            # calculate pcape from BL forcing only
            self._cloud_workfunction_1_pbl(
                error_code=error_code,
                pbl_level=pbl_level,
                geopotential_height_cloud_levels_forced=geopotential_height_cloud_levels_forced,
                t_old=t_old,
                t_new=t_new,
                t_cloud_levels_forced=t_cloud_levels_forced,
                vapor_old=vapor_old,
                vapor_forced=vapor_forced,
                cloud_work_function_1_pbl=cloud_work_function_1_pbl,
                cloud_work_function_1_fa=cloud_work_function_1_fa,
                plume=plume_dependent_constants.PLUME_INDEX,
            )

            if self.cumulus_parameterization_config.DIURNAL_CYCLE == 6:
                raise NotImplementedError(
                    "[NDSL] GF2020-->CumulusParameterization-->DiurnalCycle called with an unimplemented path."
                    "This should have been caught at initialization, but somehow you made it here."
                    "Choose another option for DIURNAL_CYCLE or implement to continue."
                )

            self._scale_cloud_workfunction_1_pbl(
                error_code=error_code,
                pbl_time_scale=pbl_time_scale,
                cloud_work_function_1_pbl=cloud_work_function_1_pbl,
                T_STAR=plume_dependent_constants.T_STAR,
                plume=plume_dependent_constants.PLUME_INDEX,
            )

        elif self.cumulus_parameterization_config.DIURNAL_CYCLE == 4:
            raise NotImplementedError(
                "[NDSL] GF2020-->CumulusParameterization-->DiurnalCycle called with an unimplemented path."
                "This should have been caught at initialization, but somehow you made it here."
                "Choose another option or for DIURNAL_CYCLE implement to continue."
            )

        if (
            self.cumulus_parameterization_config.DIURNAL_CYCLE == 5
            or self.cumulus_parameterization_config.DIURNAL_CYCLE == 2
        ):
            raise NotImplementedError(
                "[NDSL] GF2020-->CumulusParameterization-->DiurnalCycle called with an unimplemented path."
                "This should have been caught at initialization, but somehow you made it here."
                "Choose another option or for DIURNAL_CYCLE implement to continue."
            )
