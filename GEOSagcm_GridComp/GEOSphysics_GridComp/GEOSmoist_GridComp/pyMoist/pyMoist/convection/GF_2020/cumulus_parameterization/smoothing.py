from ndsl.dsl.typing import FloatField, FloatFieldIJ, Int
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import FloatField_Plume, IntFieldIJ_Plume
from ndsl.dsl.gt4py import computation, FORWARD, interval, PARALLEL, K


def smooth_tendencies(
    error_code: IntFieldIJ_Plume,
    cloud_top_level: IntFieldIJ_Plume,
    p_cloud_levels_forced: FloatField_Plume,
    del_moist_static_energy_cloud_ensemble: FloatField,
    del_vapor_cloud_ensemble: FloatField,
    del_cloud_liquid_cloud_ensemble: FloatField,
    del_u_cloud_ensemble: FloatField,
    del_v_cloud_ensemble: FloatField,
    plume: Int,
):
    from __externals__ import USE_SMOOTH_TENDENCIES

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0 and USE_SMOOTH_TENDENCIES >= 1:
            # initialize 2d temporaries
            internal_tendency_1_2d: FloatFieldIJ = 0.0
            internal_tendency_2_2d: FloatFieldIJ = 0.0
            internal_tendency_3_2d: FloatFieldIJ = 0.0
            internal_tendency_4_2d: FloatFieldIJ = 0.0
            internal_tendency_5_2d: FloatFieldIJ = 0.0
            rcount: FloatFieldIJ = 0.0
            dp: FloatFieldIJ = 0.0

    with computation(FORWARD), interval(...):
        # NOTE this entire block can be (and should be, for the sake of performance) rewritten with
        # dynamic intervals, but this cannot be done until the feature is more stable. will revisit later
        if error_code[0, 0][plume] == 0 and USE_SMOOTH_TENDENCIES >= 1 and K <= cloud_top_level[0, 0][plume]:
            rcount = 1.0e-8
            internal_tendency_1_2d = 0
            internal_tendency_2_2d = 0
            internal_tendency_3_2d = 0
            internal_tendency_4_2d = 0
            internal_tendency_5_2d = 0
            inner_loop = max(0, K - USE_SMOOTH_TENDENCIES)
            while inner_loop <= min(cloud_top_level[0, 0][plume], K + USE_SMOOTH_TENDENCIES):
                dp = p_cloud_levels_forced.at(K=inner_loop, ddim=[plume]) - p_cloud_levels_forced.at(
                    K=inner_loop + 1, ddim=[plume]
                )
                rcount = rcount + dp
                internal_tendency_1_2d = (
                    internal_tendency_1_2d + dp * del_moist_static_energy_cloud_ensemble.at(K=inner_loop)
                )
                internal_tendency_2_2d = internal_tendency_2_2d + dp * del_vapor_cloud_ensemble.at(
                    K=inner_loop
                )
                internal_tendency_3_2d = internal_tendency_3_2d + dp * del_cloud_liquid_cloud_ensemble.at(
                    K=inner_loop
                )
                internal_tendency_4_2d = internal_tendency_4_2d + dp * del_u_cloud_ensemble.at(K=inner_loop)
                internal_tendency_5_2d = internal_tendency_5_2d + dp * del_v_cloud_ensemble.at(K=inner_loop)
                inner_loop += 1

            internal_tendency_1_3d = internal_tendency_1_2d / rcount
            internal_tendency_2_3d = internal_tendency_2_2d / rcount
            internal_tendency_3_3d = internal_tendency_3_2d / rcount
            internal_tendency_4_3d = internal_tendency_4_2d / rcount
            internal_tendency_5_3d = internal_tendency_5_2d / rcount

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0 and USE_SMOOTH_TENDENCIES >= 1 and K <= cloud_top_level[0, 0][plume]:
            del_moist_static_energy_cloud_ensemble = internal_tendency_1_3d
            del_vapor_cloud_ensemble = internal_tendency_2_3d
            del_cloud_liquid_cloud_ensemble = internal_tendency_3_3d
            del_u_cloud_ensemble = internal_tendency_4_3d
            del_v_cloud_ensemble = internal_tendency_5_3d


class MakeSmoother:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass


class ApplySmoother:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass
