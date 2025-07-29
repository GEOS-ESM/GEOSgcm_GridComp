from gt4py.cartesian.gtscript import THIS_K

from ndsl.dsl.gt4py import BACKWARD, FORWARD, PARALLEL, computation, interval, log
from ndsl.dsl.typing import BoolFieldIJ, Float, FloatField, FloatFieldIJ


def vertical_interpolation(
    field: FloatField,
    interpolated_field: FloatFieldIJ,
    p_interface_mb: FloatField,
    target_pressure: Float,
    pb: FloatFieldIJ,
    pt: FloatFieldIJ,
    boolean_2d_mask: BoolFieldIJ,
    interface: bool = False,
):
    """
    Interpolate to a specific vertical level.

    Only works for non-interface fields. Must be constructed using Z_INTERFACE_DIM.

    Arguments:
        field (in): three dimensional field to be interpolated to a specific pressure
        interpolated_field (out): output two dimension field of interpolated values
        p_interface_mb (in): interface pressure in mb
        target_pressure (in): target pressure for interpolation in Pascals
        pb (in): placeholder 2d quantity, can be removed onces 2d temporaries are available
        pt (in): placeholder 2d quantity, can be removed onces 2d temporaries are available
        boolean_2d_mask (in): boolean mask to track when each cell is modified
        interface (in): specifies if input 'field' is an interface (True) or non-interface (False) field
    """
    # from __externals__ import k_end

    # mask tracks which points have been touched. check later on ensures that every point has been touched
    with computation(FORWARD), interval(0, 1):
        boolean_2d_mask = False

    with computation(PARALLEL), interval(...):
        p = log(p_interface_mb * 100)

    with computation(FORWARD), interval(-1, None):
        if interface == True:  # noqa
            pb = p
    with computation(FORWARD), interval(-1, None):
        if interface == False:  # noqa
            pb = 0.5 * (p[0, 0, -1] + p)

    with computation(BACKWARD), interval(0, -1):
        # interval is (0, -1) instead of (...) becauase the stencil is built with Z_INTERFACE_DIM
        if interface == True:  # noqa
            pt = p
            if log(target_pressure) > pt and log(target_pressure) <= pb:
                al = (pb - log(target_pressure)) / (pb - pt)
                interpolated_field = field * al + field[0, 0, 1] * (1.0 - al)
                boolean_2d_mask = True
            pb = pt

    with computation(BACKWARD), interval(1, -1):
        if interface == False:  # noqa
            pt = 0.5 * (p[0, 0, -1] + p)
            if log(target_pressure) > pt and log(target_pressure) <= pb and boolean_2d_mask == False:
                al = (pb - log(target_pressure)) / (pb - pt)
                interpolated_field = field[0, 0, -1] * al + field * (1.0 - al)
                boolean_2d_mask = True
            pb = pt

    with computation(FORWARD), interval(-2, -1):
        if interface == False:  # noqa
            pt = 0.5 * (p + p[0, 0, -1])
            pb = 0.5 * (p + p[0, 0, 1])
            if log(target_pressure) > pb and log(target_pressure) <= p[0, 0, 1] and boolean_2d_mask == False:
                interpolated_field = field
                boolean_2d_mask = True

            # ensure every point was actually touched
            if boolean_2d_mask == False:  # noqa
                interpolated_field = field

    # reset masks and temporaries for later use
    with computation(FORWARD), interval(0, 1):
        boolean_2d_mask = False
        pb = 0
        pt = 0
