from ndsl.dsl.gt4py import BACKWARD, FORWARD, computation, interval, log
from ndsl.dsl.typing import BoolFieldIJ, Float, FloatField, FloatFieldIJ


def vertical_interpolation_interface(
    field: FloatField,
    interpolated_field: FloatFieldIJ,
    p_interface_mb: FloatField,
    target_pressure: Float,
):
    """
    Interpolate to a specific vertical level on interface.

    Only works for interface fields. Must be constructed using K_INTERFACE_DIM.

    Arguments:
        field (in): three dimensional field to be interpolated to a specific pressure
        interpolated_field (out): output two dimension field of interpolated values
        p_interface_mb (in): interface pressure in mb
        target_pressure (in): target pressure for interpolation in Pascals
    """

    with computation(FORWARD), interval(-1, None):
        pb: FloatFieldIJ = log(p_interface_mb * 100.0)

    with computation(BACKWARD), interval(0, -1):
        pt: FloatFieldIJ = log(p_interface_mb * 100.0)
        if log(target_pressure) > pt and log(target_pressure) <= pb:
            al = (pb - log(target_pressure)) / (pb - pt)
            interpolated_field = field * al + field[0, 0, 1] * (1.0 - al)
        pb = pt


def vertical_interpolation(
    field: FloatField,
    interpolated_field: FloatFieldIJ,
    p_interface_mb: FloatField,
    target_pressure: Float,
):
    """
    Interpolate to a specific vertical level.

    Only works for non-interface fields. Must be constructed using K_DIM.

    Arguments:
        field (in): three dimensional field to be interpolated to a specific pressure
        interpolated_field (out): output two dimension field of interpolated values
        p_interface_mb (in): interface pressure in mb
        target_pressure (in): target pressure for interpolation in Pascals
    """
    # mask tracks which points have been touched. check later on ensures that every point has been touched
    with computation(FORWARD), interval(0, 1):
        track_points: BoolFieldIJ = False

    # with computation(PARALLEL), interval(...):
    #     p = log(p_interface_mb * 100.0)

    with computation(FORWARD), interval(-1, None):
        pb: FloatFieldIJ = 0.5 * (log(p_interface_mb * 100.0) + log(p_interface_mb[0, 0, 1] * 100.0))

    with computation(BACKWARD), interval(1, None):
        pt: FloatFieldIJ = 0.5 * (log(p_interface_mb[0, 0, -1] * 100.0) + log(p_interface_mb * 100.0))
        if log(target_pressure) > pt and log(target_pressure) <= pb and not track_points:
            al = (pb - log(target_pressure)) / (pb - pt)
            interpolated_field = field[0, 0, -1] * al + field * (1.0 - al)
            track_points = True
        pb = pt

    with computation(FORWARD), interval(-1, None):
        pb2: FloatFieldIJ = 0.5 * (log(p_interface_mb * 100.0) + log(p_interface_mb[0, 0, 1] * 100.0))
        if log(target_pressure) > pb2 and log(target_pressure) <= log(p_interface_mb[0, 0, 1] * 100.0) and not track_points:
            interpolated_field = field
            track_points = True

        # ensure every point was actually touched
        if track_points == False:  # noqa
            interpolated_field = field
