import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
from ndsl.dsl.gt4py import function, log
from ndsl.stencils.column_operations import column_max
from pyMoist.shared.incloud_processes import ice_fraction


@function
def saturation_vapor_pressure(t):
    """
    Compute saturation vapor pressure

    Args:
        t (in): temperature in Kelvin
    """
    # CAUTION: This function has not been verified!!!
    t_celcius = t - 273.155
    if t_celcius < -20.0:
        toot = 273.16 / t
        toto = 1 / toot
        eilog = (
            -9.09718 * (toot - 1)
            - 3.56654 * (log(toot) / log(10.0))
            + 0.876793 * (1 - toto)
            + (log(6.1071) / log(10.0))
        )
        saturation_vapor_pressure = 10**eilog
    else:
        tsot = 373.16 / t
        ewlog = -7.90298 * (tsot - 1) + 5.02808 * (log(tsot) / log(10.0))
        ewlog2 = ewlog - 1.3816e-07 * (10 ** (11.344 * (1 - (1 / tsot))) - 1)
        ewlog3 = ewlog2 + 0.0081328 * (10 ** (-3.49149 * (tsot - 1)) - 1)
        ewlog4 = ewlog3 + (log(1013.246) / log(10.0))
        saturation_vapor_pressure = 10**ewlog4

    return saturation_vapor_pressure


@function
def updraft_origin_conditions_perturbation(value, perturbation, start_index, end_index):
    """
    Modify the value based by the maximum value of perturbation

    Args:
        value
        perturbation
        start_index: starting index for max value search
        end_index: ending index for max value search
    """
    max_value, _ = column_max(perturbation, start_index, end_index)
    return value + cumulus_parameterization_constants.CP * max_value


@function
def get_cloud_boundary_conditions(
    field,
    scalar_perturbation,
    p,
    updraft_origin_level,
    ocean_fraction,
    BOUNDARY_CONDITION_METHOD,
    AVERAGE_LAYER_DEPTH,
    k_end,
    compute_perturbation,
    perturbation_field,
):
    """
    Get the conditions of a field at the origin level of a updraft by taking an average across layers.

    Dimensions of the average:
        a) to pick the value at k22 level, instead of an average between
            (k22-order_aver, ..., k22-1, k22) set order_aver=kts
        b) to average between kts and k22 => set order_aver = k22

        order_aver = 4:
            BOUNDARY_CONDITION_METHOD == 0: average between updraft_origin_level,
                updraft_origin_level - 1, updraft_origin_level - 2
            BOUNDARY_CONDITION_METHOD == 1: average between updraft_origin_level + 1,
                updraft_origin_level, updraft_origin_level - 1

        order_aver = 0:
            average between updraft_origin_level, updraft_origin_level - 1, updraft_origin_level - 2

    Args:
        field: field from which value is extracted
        scalar_perturbation: scalar perturbation to be applied to the value after extraction
        p: pressure field
        updraft_origin_level: origin level of updraft (LCL, LFC, surface, etc.)
        ocean_fraction: fraction of surface in the grid cell which is ocean
        BOUNDARY_CONDITION_METHOD: switch which controls function
        AVERAGE_LAYER_DEPTH: average depth of a model layer
        k_end: number of model layers
        compute_perturbation: switch which controls if a full field perturbation is applied to the value
        perturbation_field: perturbation field which may be applied to the value (must supply dummy field if
            compute_perturbation is False)
    """
    # internal constant
    frac_ave_layer_ocean = 0.3

    if BOUNDARY_CONDITION_METHOD == 0:
        order_aver = 3
        local_order_aver = min(updraft_origin_level + 1, order_aver)

        source_parcel_value = 0.0
        count = 0
        while count < local_order_aver:
            source_parcel_value = source_parcel_value + field.at(K=updraft_origin_level - count)
            count += 1

        source_parcel_value = source_parcel_value / local_order_aver

    elif BOUNDARY_CONDITION_METHOD == 1:
        effec_frac = (1.0 - ocean_fraction) + ocean_fraction * frac_ave_layer_ocean
        average_layer = AVERAGE_LAYER_DEPTH * effec_frac

        start_index = 0
        level = 0
        p_at_origin = p.at(K=updraft_origin_level) + 0.5 * average_layer
        while level <= k_end - 1:
            new = abs(p.at(K=level) - p_at_origin)
            old = abs(p.at(K=start_index) - p_at_origin)
            if new < old:
                start_index = level
            level += 1

        end_index = 0
        level = 0
        p_at_origin = p.at(K=updraft_origin_level) - 0.5 * average_layer
        while level <= k_end - 1:
            new = abs(p.at(K=level) - p_at_origin)
            old = abs(p.at(K=end_index) - p_at_origin)
            if new < old:
                end_index = level
            level += 1

        if start_index >= end_index:
            source_parcel_value = field.at(K=updraft_origin_level)
            dp_layer = 1.0
            level_c = start_index

        else:
            dp_layer = 1.0e-06
            source_parcel_value = 0.0
            level_c = 0
            stop_computation = False
            level = start_index
            while level <= k_end - 1 and stop_computation == False:
                dp = -(p.at(K=level + 1) - p.at(K=level))
                if dp_layer + dp <= average_layer:
                    dp_layer = dp_layer + dp
                    source_parcel_value = source_parcel_value + field.at(K=level) * dp
                    level += 1
                else:
                    dp = average_layer - dp_layer
                    dp_layer = dp_layer + dp
                    source_parcel_value = source_parcel_value + field.at(K=level) * dp
                    stop_computation = True
                    level += 1

            source_parcel_value = source_parcel_value / dp_layer
            level_c = max(start_index, level)

        # this perturbation is only used for moist static energy
        if compute_perturbation == True:
            source_parcel_value = updraft_origin_conditions_perturbation(
                value=source_parcel_value,
                perturbation=perturbation_field,
                start_index=start_index,
                end_index=level_c,
            )

    return source_parcel_value + scalar_perturbation


@function
def compute_dewpoint(p, vapor):
    rr = vapor + 1e-8
    es = p * rr / (0.622 + rr)
    esln = log(es)
    return (35.86 * esln - 4947.2325) / (esln - 23.6837)


@function
def liquid_fraction(
    t,
    convection_fraction,
    surface_type,
    FRAC_MODIS,
):

    if FRAC_MODIS == 1:
        liquid_fraction = 1.0 - ice_fraction(t, convection_fraction, surface_type)
    else:
        liquid_fraction = min(
            1.0,
            (
                max(0.0, (t - cumulus_parameterization_constants.T_ICE))
                / (cumulus_parameterization_constants.T_0 - cumulus_parameterization_constants.T_ICE)
            )
            ** 2,
        )

    return liquid_fraction


# @function
# def get_delmix(
#     ocean_fraction,
#     sub_cloud_level,
#     p_forced,
#     in_var,
#     out_var,
# ):
#     internal_var = out_var.at(K=0)
#     x1 = 0.0
#     x2 = 0.0
#     level = 0
#     while level <= sub_cloud_level:
#         dp = p_forced[0,0,1] - p_forced
#         x1 = x1 + dp*in_var
#         x2 = x2 + dp
#         level += 1

#     delta = abs(internal_var-x1/(x2 + 1.e-12))

#     level = 0
#     while level <= sub_cloud_level:
#         out_var =
