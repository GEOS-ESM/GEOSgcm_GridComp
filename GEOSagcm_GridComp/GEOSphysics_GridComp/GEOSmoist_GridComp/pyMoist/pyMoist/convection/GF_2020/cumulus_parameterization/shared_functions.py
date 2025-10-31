from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField, IntField, Float, Int
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
import pyMoist.constants as constants
from gt4py.cartesian.gtscript import (
    PARALLEL,
    computation,
    interval,
    int32,
    log,
    exp,
)
from ndsl.dsl.gt4py import function


@function
def saturation_vapor_pressure(t: Float):
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
def column_max(field, start_index, end_index):
    max_index = 0
    level = start_index
    while level <= end_index:
        new = field.at(K=level)
        old = field.at(K=max_index)
        if new > old:
            max_index = level
        level += 1

    return field.at(K=max_index), max_index


@function
def column_min(field, start_index, end_index):
    min_index = 0
    level = start_index
    while level <= end_index:
        new = field.at(K=level)
        old = field.at(K=min_index)
        if new < old:
            min_index = level
        level += 1

    return field.at(K=min_index), min_index


@function
def cloud_boundary_conditions_perturbation(average, perturbation, start_index, end_index):
    max_value, _ = column_max(perturbation, start_index, end_index)
    return average + cumulus_parameterization_constants.CP * max_value


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
    # internal constant
    frac_ave_layer_ocean = 0.3

    if BOUNDARY_CONDITION_METHOD == 0:
        order_aver = 3
        local_order_aver = min(updraft_origin_level, order_aver)

        source_parcel_value = 0.0
        level = 0
        while level < local_order_aver:
            source_parcel_value = source_parcel_value + field.at(K=updraft_origin_level - level + 1)

        source_parcel_value = source_parcel_value / local_order_aver

    elif BOUNDARY_CONDITION_METHOD == 1:
        effec_frac = (1.0 - ocean_fraction) + ocean_fraction * frac_ave_layer_ocean
        x_ave_layer = AVERAGE_LAYER_DEPTH * effec_frac

        start_index = 0
        level = 0
        while level <= k_end - 1:
            new = abs(p.at(K=level) - (p.at(K=updraft_origin_level) + 0.5 * x_ave_layer))
            old = abs(p.at(K=start_index) - (p.at(K=updraft_origin_level) + 0.5 * x_ave_layer))
            if new < old:
                start_index = level
            level += 1

        end_index = 0
        level = 0
        while level <= k_end - 1:
            new = abs(p.at(K=level) - (p.at(K=updraft_origin_level) - 0.5 * x_ave_layer))
            old = abs(p.at(K=end_index) - (p.at(K=updraft_origin_level) - 0.5 * x_ave_layer))
            if new < old:
                end_index = level
            level += 1

        start_index = min(k_end - 1, max(start_index, 0))
        end_index = min(k_end - 1, max(end_index, 0))

        if start_index >= end_index:
            source_parcel_value = field.at(K=updraft_origin_level)
            dp_layer = 0.0
            level_c = start_index

        else:
            dp_layer = 1.0e-06
            source_parcel_value = 0.0
            level_c = 0
            stop_computation = False
            level = start_index
            while level <= k_end - 1 and stop_computation == False:
                dp = -(p[0, 0, 1] - p)
                if dp_layer + dp <= x_ave_layer:
                    dp_layer = dp_layer + dp
                    source_parcel_value = source_parcel_value + field * dp

                else:
                    dp = x_ave_layer - dp_layer
                    dp_layer = dp_layer + dp
                    source_parcel_value = source_parcel_value + field * dp
                    stop_computation = True

        source_parcel_value = source_parcel_value / dp_layer
        level_c = max(start_index, level)

        # this perturbation is only used for moist static energy
        if compute_perturbation == True:
            source_parcel_value = cloud_boundary_conditions_perturbation(
                average=source_parcel_value,
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
