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
        if interface is True:
            pb = p
    with computation(FORWARD), interval(-1, None):
        if interface is False:
            pb = 0.5 * (p[0, 0, -1] + p)

    with computation(BACKWARD), interval(0, -1):
        if interface is True:
            pt = p.at(K=THIS_K)
            if log(target_pressure) > pt and log(target_pressure) <= pb:
                al = (pb - log(target_pressure)) / (pb - pt)
                interpolated_field = field.at(K=THIS_K) * al + field.at(K=THIS_K + 1) * (1.0 - al)
            pb = pt

    with computation(BACKWARD), interval(1, -1):
        if interface is False:
            pt = 0.5 * (p.at(K=THIS_K - 1) + p.at(K=THIS_K))
            if log(target_pressure) > pt and log(target_pressure) <= pb:
                al = (pb - log(target_pressure)) / (pb - pt)
                interpolated_field = field.at(K=THIS_K - 1) * al + field.at(K=THIS_K) * (1.0 - al)
                boolean_2d_mask = True
            pb = pt

    with computation(FORWARD), interval(-2, -1):
        if interface is False:
            pt = 0.5 * (p + p[0, 0, -1])
            pb = 0.5 * (p + p[0, 0, 1])
            if log(target_pressure) > pb and log(target_pressure) <= p[0, 0, 1]:
                interpolated_field = field[0, 0, -1]

            # ensure every point was actually touched
            if boolean_2d_mask is False:
                interpolated_field = p

    # reset masks and temporaries for later use
    with computation(FORWARD), interval(0, 1):
        boolean_2d_mask = False
        pb = 0
        pt = 0


# @function
# def cs_prof(

# )

# def cs_interpolator(
#     field: FloatField,
#     edge_height_above_surface: FloatField,
#     desired_height: Float,
#     interpolated_field: FloatFieldIJ,
#     minimum_value: Float,
#     # 2d temporaries
#     grid_ratio: FloatFieldIJ,

# ):
#     with computation(PARALLEL), interval(...):
#         dz = edge_height_above_surface - edge_height_above_surface[0,0,1]
#         field_copy = field


#     # start cs_prof
#     with computation(PARALLEL), interval(0, 1):
#         grid_radio = dz[0,0,1] - dz
#         bet = grid_radio*(grid_radio+0.5)
# q(i,1) = ( (grat+grat)*(grat+1.)*q2(i,1) + q2(i,2) ) / bet
# gam(i,1) = ( 1. + grat*(grat+1.5) ) / bet


#   do i=i1,i2
#          grat = delp(i,2) / delp(i,1)   ! grid ratio
#           bet = grat*(grat+0.5)
#        q(i,1) = ( (grat+grat)*(grat+1.)*q2(i,1) + q2(i,2) ) / bet
#      gam(i,1) = ( 1. + grat*(grat+1.5) ) / bet
#   enddo

#   do k=2,km
#      do i=i1,i2
#            d4(i) = delp(i,k-1) / delp(i,k)
#              bet =  2. + d4(i) + d4(i) - gam(i,k-1)
#           q(i,k) = ( 3.*(q2(i,k-1)+d4(i)*q2(i,k)) - q(i,k-1) )/bet
#         gam(i,k) = d4(i) / bet
#      enddo
#   enddo

#   do i=i1,i2
#          a_bot = 1. + d4(i)*(d4(i)+1.5)
#      q(i,km+1) = (2.*d4(i)*(d4(i)+1.)*q2(i,km)+q2(i,km-1)-a_bot*q(i,km))  &
#                / ( d4(i)*(d4(i)+0.5) - a_bot*gam(i,km) )
#   enddo

#   do k=km,1,-1
#      do i=i1,i2
#         q(i,k) = q(i,k) - gam(i,k)*q(i,k+1)
#      enddo
#   enddo

# ! Apply *large-scale* constraints
#   do i=i1,i2
#      q(i,2) = min( q(i,2), max(q2(i,1), q2(i,2)) )
#      q(i,2) = max( q(i,2), min(q2(i,1), q2(i,2)) )
#   enddo

#   do k=2,km
#      do i=i1,i2
#         gam(i,k) = q2(i,k) - q2(i,k-1)
#      enddo
#   enddo

# ! Interior:
#   do k=3,km-1
#      do i=i1,i2
#         if ( gam(i,k-1)*gam(i,k+1)>0. ) then
# ! Apply large-scale constraint to ALL fields if not local max/min
#              q(i,k) = min( q(i,k), max(q2(i,k-1),q2(i,k)) )
#              q(i,k) = max( q(i,k), min(q2(i,k-1),q2(i,k)) )
#         else
#           if ( gam(i,k-1) > 0. ) then
# ! There exists a local max
#                q(i,k) = max(q(i,k), min(q2(i,k-1),q2(i,k)))
#           else
# ! There exists a local min
#                q(i,k) = min(q(i,k), max(q2(i,k-1),q2(i,k)))
#                if ( iv==0 ) q(i,k) = max(0., q(i,k))
#           endif
#         endif
#      enddo
#   enddo

# ! Bottom:
#   do i=i1,i2
#      q(i,km) = min( q(i,km), max(q2(i,km-1), q2(i,km)) )
#      q(i,km) = max( q(i,km), min(q2(i,km-1), q2(i,km)) )
#   enddo
