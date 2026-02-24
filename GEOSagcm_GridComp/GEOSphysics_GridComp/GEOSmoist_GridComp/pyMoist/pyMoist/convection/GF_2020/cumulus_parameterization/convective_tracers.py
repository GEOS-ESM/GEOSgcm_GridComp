from ndsl.dsl.gt4py import computation, interval, PARALLEL, FORWARD, K, exp, BACKWARD
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import (
    IntFieldIJ_Plume,
    FloatField_ConvectionTracers,
    FloatFieldIJ_ConvectionTracers,
    FloatField_Plume,
    FloatFieldIJ_Plume,
    FloatField_ConvectionTracers_Plume,
)
from pyMoist.field_types import (
    ConvectionTracerMetaDataTable_Float,
    ConvectionTracerMetaDataTable_Bool,
    ConvectionTracerMetaDataTable_x4,
)
from ndsl.dsl.typing import IntFieldIJ, Int, FloatField, FloatFieldIJ, Float
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl import StencilFactory, QuantityFactory, Quantity, Local, NDSLRuntime
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)
import pyMoist.constants as constants
from pyMoist.convection.GF_2020.cumulus_parameterization.shared_stencils import tridiag
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
from pyMoist.convection.GF_2020.cumulus_parameterization.shared_functions import get_cloud_boundary_conditions
from pyMoist.convection_tracers import ConvectionTracers


def environment_cloud_levels_chemistry(
    error_code: IntFieldIJ_Plume,
    chemistry_tracers: FloatField_ConvectionTracers,
    chemistry_tracers_cloud_levels: FloatField_ConvectionTracers,
    plume: Int,
):
    """Set chemistry tracers at cloud levels. Behavior is determined by CLOUD_LEVEL_OPTION external.
    Only options 1 and 2 have been implemented, only option 2 has been tested.

    Args:
        error_code (IntFieldIJ_Plume)
        chemistry_tracers (FloatField_ConvectionTracers)
        chemistry_tracers_cloud_levels (FloatField_ConvectionTracers)
        plume (Int)
    """
    from __externals__ import k_end, NUMBER_OF_TRACERS, CLOUD_LEVEL_OPTION

    with computation(FORWARD), interval(1, -1):
        if error_code[0, 0][plume] == 0 and CLOUD_LEVEL_OPTION == 1:
            tracer = 0
            while tracer < NUMBER_OF_TRACERS:
                chemistry_tracers_cloud_levels[0, 0, 0][tracer] = (
                    0.5 * chemistry_tracers[0, 0, -1][tracer] + chemistry_tracers[0, 0, 0][tracer]
                )
                tracer += 1

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0 and CLOUD_LEVEL_OPTION == 1:
            tracer = 0
            while tracer < NUMBER_OF_TRACERS:
                chemistry_tracers_cloud_levels[0, 0, 0][tracer] = chemistry_tracers[0, 0, 0][tracer]
                tracer += 1

    with computation(FORWARD), interval(-1, None):
        if error_code[0, 0][plume] == 0 and CLOUD_LEVEL_OPTION == 1:
            tracer = 0
            while tracer < NUMBER_OF_TRACERS:
                chemistry_tracers_cloud_levels[0, 0, 0][tracer] = chemistry_tracers[0, 0, 0][tracer]
                tracer += 1

    with computation(PARALLEL), interval(0, -1):
        if error_code[0, 0][plume] == 0 and CLOUD_LEVEL_OPTION != 1:
            tracer = 0
            while tracer < NUMBER_OF_TRACERS:
                chemistry_tracers_cloud_levels[0, 0, 0][tracer] = chemistry_tracers[0, 0, 0][tracer]
                tracer += 1


def updraft_chemistry(
    error_code: IntFieldIJ_Plume,
    updraft_origin_level: IntFieldIJ_Plume,
    cloud_top_level: IntFieldIJ_Plume,
    ocean_fraction: FloatFieldIJ,
    p_forced: FloatField,
    p_cloud_levels_forced: FloatField_Plume,
    geopotential_height_cloud_levels: FloatField,
    vertical_velocity_3d: FloatField,
    mass_entrainment_updraft_forced: FloatField_Plume,
    mass_detrainment_updraft_forced: FloatField_Plume,
    normalized_massflux_updraft_forced: FloatField_Plume,
    chemistry_tracers: FloatField_ConvectionTracers,
    chemistry_tracers_sc_updraft: FloatField_ConvectionTracers,
    chemistry_tracers_cloud_levels: FloatField_ConvectionTracers,
    chemistry_tracers_pw_updraft: FloatField_ConvectionTracers,
    chemistry_tracers_total_pw_updraft: FloatFieldIJ_ConvectionTracers,
    convection_tracers_vect_hcts: ConvectionTracerMetaDataTable_x4,
    convection_tracers_fscav: ConvectionTracerMetaDataTable_Float,
    convection_tracers_use_gcc_washout: ConvectionTracerMetaDataTable_Bool,
    tracer_cloud_boundary: FloatFieldIJ_ConvectionTracers,
    AVERAGE_LAYER_DEPTH: Float,
    plume: Int,
):
    """Compute chemistry effects in the updraft.

    Args:
        error_code (IntFieldIJ_Plume)
        updraft_origin_level (IntFieldIJ_Plume)
        cloud_top_level (IntFieldIJ_Plume)
        ocean_fraction (FloatFieldIJ)
        p_forced (FloatField)
        p_cloud_levels_forced (FloatField_Plume)
        geopotential_height_cloud_levels (FloatField)
        vertical_velocity_3d (FloatField)
        mass_entrainment_updraft_forced (FloatField_Plume)
        mass_detrainment_updraft_forced (FloatField_Plume)
        normalized_massflux_updraft_forced (FloatField_Plume)
        chemistry_tracers (FloatField_ConvectionTracers)
        chemistry_tracers_sc_updraft (FloatField_ConvectionTracers)
        chemistry_tracers_cloud_levels (FloatField_ConvectionTracers)
        chemistry_tracers_pw_updraft (FloatField_ConvectionTracers)
        chemistry_tracers_total_pw_updraft (FloatFieldIJ_ConvectionTracers)
        convection_tracers_vect_hcts (ConvectionTracerMetaDataTable_x4)
        convection_tracers_fscav (ConvectionTracerMetaDataTable_Float)
        convection_tracers_use_gcc_washout (ConvectionTracerMetaDataTable_Bool)
        tracer_cloud_boundary (FloatFieldIJ_ConvectionTracers)
        AVERAGE_LAYER_DEPTH (Float)
        plume (Int)
    """
    from __externals__ import k_end, BOUNDARY_CONDITION_METHOD, USE_TRACER_SCAVENGE, NUMBER_OF_TRACERS

    with computation(FORWARD), interval(...):
        tracer = 0
        while tracer < NUMBER_OF_TRACERS:
            # initialization
            chemistry_tracers_sc_updraft[0, 0, 0][tracer] = chemistry_tracers_cloud_levels[0, 0, 0][tracer]
            chemistry_tracers_pw_updraft[0, 0, 0][tracer] = 0.0
            chemistry_tracers_total_pw_updraft[0, 0][tracer] = 0.0
            tracer += 1

    with computation(PARALLEL), interval(...):
        # make garbage field so the get_cloud_boundary_conditions call does not break
        # this is never touched so long as compute_perturbation=False
        dummy_field_no_read = 0.0
        chemistry_tracers_cloud_levels_3d = 0.0

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            # NOTE this double loop is disgusting. need a better solution for it
            tracer = 0
            while tracer < NUMBER_OF_TRACERS:
                level = 0
                while level <= k_end:
                    chemistry_tracers_cloud_levels_3d[0, 0, level] = chemistry_tracers_cloud_levels[
                        0, 0, level
                    ][tracer]
                    level += 1
                tracer_cloud_boundary[0, 0][tracer] = get_cloud_boundary_conditions(
                    field=chemistry_tracers_cloud_levels_3d,
                    scalar_perturbation=0,
                    p=p_forced,
                    updraft_origin_level=updraft_origin_level[0, 0][plume],
                    ocean_fraction=ocean_fraction,
                    BOUNDARY_CONDITION_METHOD=BOUNDARY_CONDITION_METHOD,
                    AVERAGE_LAYER_DEPTH=AVERAGE_LAYER_DEPTH,
                    k_end=k_end,
                    compute_perturbation=False,
                    perturbation_field=dummy_field_no_read,
                )
                tracer += 1

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0:
            tracer = 0
            while tracer < NUMBER_OF_TRACERS:
                if K <= updraft_origin_level[0, 0][plume]:
                    chemistry_tracers_sc_updraft[0, 0, 0][tracer] = tracer_cloud_boundary[0, 0][tracer]
                tracer += 1

    with computation(FORWARD), interval(1, None):
        if (
            error_code[0, 0][plume] == 0
            and K >= updraft_origin_level[0, 0][plume] + 1
            and K <= cloud_top_level[0, 0][plume] + 1
        ):
            # entrainment, detrainment, mass flux
            massflux_internal = normalized_massflux_updraft_forced[0, 0, -1][plume]
            entrainment_internal = mass_entrainment_updraft_forced[0, 0, -1][plume]
            detrainment_internal = 0.5 * mass_detrainment_updraft_forced[0, 0, -1][plume]
            denom = massflux_internal - detrainment_internal + entrainment_internal

            # transport + mixing
            tracer = 0
            while tracer < NUMBER_OF_TRACERS:
                if denom > 0.0:
                    chemistry_tracers_sc_updraft[0, 0, 0][tracer] = (
                        chemistry_tracers_sc_updraft[0, 0, -1][tracer] * massflux_internal
                        - chemistry_tracers_sc_updraft[0, 0, -1][tracer] * detrainment_internal
                        + chemistry_tracers[0, 0, -1][tracer] * entrainment_internal
                    ) / denom
                else:
                    chemistry_tracers_sc_updraft[0, 0, 0][tracer] = chemistry_tracers_sc_updraft[0, 0, -1][
                        tracer
                    ]
                tracer += 1

            # scavenging section - skip if USE_TRACER_SCAVENGE = 0 or on shallow plume
            if not (USE_TRACER_SCAVENGE == 0 or plume == cumulus_parameterization_constants.SHALLOW):
                dz = geopotential_height_cloud_levels - geopotential_height_cloud_levels[0, 0, -1]

                # in-cloud vert velocity for scavenging formulation 2
                updraft_vertical_velosity_internal = vertical_velocity_3d

                tracer = 0
                while tracer < NUMBER_OF_TRACERS:
                    is_gcc = convection_tracers_use_gcc_washout.A[tracer]
                    # aerosol scavenging
                    if convection_tracers_fscav.A[tracer] > 1.0e-6:
                        # formulation 1 as in GOCART with RAS conv_par
                        if USE_TRACER_SCAVENGE == 1:
                            chemistry_tracers_pw_updraft[0, 0, 0][tracer] = max(
                                0.0,
                                chemistry_tracers_sc_updraft[0, 0, 0][tracer]
                                * (1.0 - exp(-convection_tracers_fscav.A[tracer] * (dz / 1000.0))),
                            )

                        # formulation 2 as in GOCART
                        elif USE_TRACER_SCAVENGE == 2:
                            option_not_implemented = True

                        # formulation 3 - orignal GF conv_par
                        elif USE_TRACER_SCAVENGE == 3:
                            option_not_implemented = True

                        # (in cloud) total mixing ratio in gas and aqueous phases
                        chemistry_tracers_sc_updraft[0, 0, 0][tracer] = (
                            chemistry_tracers_sc_updraft[0, 0, 0][tracer]
                            - chemistry_tracers_pw_updraft[0, 0, 0][tracer]
                        )

                    # tracer gas phase scavenging
                    elif convection_tracers_vect_hcts.A[tracer, 1] > 1.0e-6:
                        option_not_implemented = True
                        if is_gcc == True:
                            option_not_implemented = True
                        if USE_TRACER_SCAVENGE == 3:
                            option_not_implemented = True
                        else:
                            if is_gcc == True:
                                option_not_implemented = True
                            else:
                                option_not_implemented = True

                    tracer += 1

                dp = 100.0 * (p_cloud_levels_forced[0, 0, 0][plume] - p_cloud_levels_forced[0, 0, 1][plume])
                tracer = 0
                while tracer < NUMBER_OF_TRACERS:
                    chemistry_tracers_total_pw_updraft[0, 0][tracer] = (
                        chemistry_tracers_total_pw_updraft[0, 0][tracer]
                        + chemistry_tracers_pw_updraft[0, 0, 0][tracer] * dp / constants.MAPL_GRAV
                    )
                    tracer += 1


def downdraft_chemistry(
    error_code: IntFieldIJ_Plume,
    downdraft_origin_level: IntFieldIJ_Plume,
    p_cloud_levels_forced: FloatField_Plume,
    evaporate_in_downdraft_forced: FloatField_Plume,
    total_normalized_integrated_evaporate_forced: FloatFieldIJ,
    total_normalized_integrated_condensate_forced: FloatFieldIJ_Plume,
    normalized_massflux_downdraft_forced: FloatField_Plume,
    mass_entrainment_downdraft_forced: FloatField_Plume,
    mass_detrainment_downdraft_forced: FloatField_Plume,
    chemistry_tracers: FloatField_ConvectionTracers,
    chemistry_tracers_cloud_levels: FloatField_ConvectionTracers,
    chemistry_tracers_sc_downdraft: FloatField_ConvectionTracers,
    chemistry_tracers_pw_downdraft: FloatField_ConvectionTracers,
    chemistry_tracers_total_pw_downdraft: FloatFieldIJ_ConvectionTracers,
    chemistry_tracers_total_pw_updraft: FloatFieldIJ_ConvectionTracers,
    plume: Int,
):
    """Compute chemistry effects in the downdraft.

    Args:
        error_code (IntFieldIJ_Plume)
        downdraft_origin_level (IntFieldIJ_Plume)
        p_cloud_levels_forced (FloatField_Plume)
        evaporate_in_downdraft_forced (FloatField_Plume)
        total_normalized_integrated_evaporate_forced (FloatFieldIJ)
        total_normalized_integrated_condensate_forced (FloatFieldIJ_Plume)
        normalized_massflux_downdraft_forced (FloatField_Plume)
        mass_entrainment_downdraft_forced (FloatField_Plume)
        mass_detrainment_downdraft_forced (FloatField_Plume)
        chemistry_tracers (FloatField_ConvectionTracers)
        chemistry_tracers_cloud_levels (FloatField_ConvectionTracers)
        chemistry_tracers_sc_downdraft (FloatField_ConvectionTracers)
        chemistry_tracers_pw_downdraft (FloatField_ConvectionTracers)
        chemistry_tracers_total_pw_downdraft (FloatFieldIJ_ConvectionTracers)
        chemistry_tracers_total_pw_updraft (FloatFieldIJ_ConvectionTracers)
        plume (Int)
    """
    from __externals__ import NUMBER_OF_TRACERS, USE_TRACER_EVAPORATION

    with computation(FORWARD), interval(...):
        tracer = 0
        while tracer < NUMBER_OF_TRACERS:
            chemistry_tracers_sc_downdraft[0, 0, 0][tracer] = 0.0
            chemistry_tracers_pw_downdraft[0, 0, 0][tracer] = 0.0
            chemistry_tracers_total_pw_downdraft[0, 0][tracer] = 0.0
            tracer += 1

    with computation(FORWARD), interval(0, 1):
        if plume != cumulus_parameterization_constants.SHALLOW and error_code[0, 0][plume] == 0:
            # fration of the total rain that was evaporated
            evaporation_fraction: FloatFieldIJ = -total_normalized_integrated_evaporate_forced / (
                1.0e-16 + total_normalized_integrated_condensate_forced[0, 0][plume]
            )

            # scalar concentration in-cloud - downdraft

            # at downdraft_origin_level
            level = downdraft_origin_level[0, 0][plume]
            if USE_TRACER_EVAPORATION == 0:
                internal_precip: FloatFieldIJ = 0.0
            else:
                internal_precip: FloatFieldIJ = (
                    evaporate_in_downdraft_forced.at(K=level, ddim=[plume])
                    / (1.0e-16 + total_normalized_integrated_evaporate_forced)
                    * evaporation_fraction
                )

            dp = 100.0 * (
                p_cloud_levels_forced.at(K=level, ddim=[plume])
                - p_cloud_levels_forced.at(K=level + 1, ddim=[plume])
            )

            # downdrafts will be initiated with a mixture of 50% environmental and in-cloud concentrations
            tracer = 0
            while tracer < NUMBER_OF_TRACERS:
                chemistry_tracers_sc_downdraft[0, 0, level][tracer] = chemistry_tracers_cloud_levels[
                    0, 0, level
                ][tracer]
                chemistry_tracers_pw_downdraft[0, 0, level][tracer] = (
                    -internal_precip
                    * chemistry_tracers_total_pw_updraft[0, 0][tracer]
                    * constants.MAPL_GRAV
                    / dp
                )
                chemistry_tracers_sc_downdraft[0, 0, level][tracer] = (
                    chemistry_tracers_sc_downdraft[0, 0, level][tracer]
                    - chemistry_tracers_pw_downdraft[0, 0, level][tracer]
                )
                chemistry_tracers_total_pw_downdraft[0, 0][tracer] = (
                    chemistry_tracers_total_pw_downdraft[0, 0][tracer]
                    + chemistry_tracers_pw_downdraft[0, 0, level][tracer] * dp / constants.MAPL_GRAV
                )
                tracer += 1

    # calculate downdraft mass terms
    with computation(BACKWARD), interval(0, -1):
        if (
            plume != cumulus_parameterization_constants.SHALLOW
            and error_code[0, 0][plume] == 0
            and K <= downdraft_origin_level[0, 0][plume] - 1
        ):
            massflux_internal = normalized_massflux_downdraft_forced[0, 0, 1][plume]
            entrainment_internal = mass_entrainment_downdraft_forced[0, 0, 0][plume]
            detrainment_internal = 0.5 * mass_detrainment_downdraft_forced[0, 0, 0][plume]
            denom = massflux_internal - detrainment_internal + entrainment_internal

            # transport + mixing
            tracer = 0
            while tracer < NUMBER_OF_TRACERS:
                if denom > 0.0:
                    chemistry_tracers_sc_downdraft[0, 0, 0][tracer] = (
                        chemistry_tracers_sc_downdraft[0, 0, 1][tracer] * massflux_internal
                        - chemistry_tracers_sc_downdraft[0, 0, 1][tracer] * detrainment_internal
                        + chemistry_tracers[0, 0, 0][tracer] * entrainment_internal
                    ) / denom
                else:
                    chemistry_tracers_sc_downdraft[0, 0, 0][tracer] = chemistry_tracers_sc_downdraft.at(
                        K=K + 1, ddim=[tracer]
                    )
                tracer += 1

            # evaporation term
            if USE_TRACER_EVAPORATION != 0:
                dp = 100.0 * (p_cloud_levels_forced[0, 0, 0][plume] - p_cloud_levels_forced[0, 0, 1][plume])

                # fraction of evaporated precip per layer
                internal_precip = evaporate_in_downdraft_forced[0, 0, 0][plume] / (
                    1.0e-16 + total_normalized_integrated_evaporate_forced
                )

                # fraction of the total precip that was actually evaporated at layer k
                internal_precip = internal_precip * evaporation_fraction

                # sanity check
                internal_precip = min(1.0, max(internal_precip, 0.0))

                tracer = 0
                while tracer < NUMBER_OF_TRACERS:
                    # amount evaporated by the downdraft from the precipitation
                    chemistry_tracers_pw_downdraft[0, 0, 0][tracer] = (
                        -internal_precip
                        * chemistry_tracers_total_pw_updraft[0, 0][tracer]
                        * constants.MAPL_GRAV
                        / dp
                    )

                    # final tracer in the downdraft
                    chemistry_tracers_sc_downdraft[0, 0, 0][tracer] = (
                        chemistry_tracers_sc_downdraft[0, 0, 0][tracer]
                        - chemistry_tracers_pw_downdraft[0, 0, 0][tracer]
                    )

                    # total evaporated tracer
                    chemistry_tracers_total_pw_downdraft[0, 0][tracer] = (
                        chemistry_tracers_total_pw_downdraft[0, 0][tracer]
                        + chemistry_tracers_pw_downdraft[0, 0, 0][tracer] * dp / constants.MAPL_GRAV
                    )

                    tracer += 1


def vertical_transport_part_1(
    error_code: IntFieldIJ_Plume,
    cloud_top_level: IntFieldIJ_Plume,
    p_cloud_levels_forced: FloatField_Plume,
    environment_massflux: FloatField,
    normalized_massflux_updraft_forced: FloatField_Plume,
    normalized_massflux_downdraft_forced: FloatField_Plume,
    epsilon_forced: FloatFieldIJ_Plume,
    chemistry_tracers: FloatField_ConvectionTracers,
    chemistry_tracers_output: FloatField_ConvectionTracers_Plume,
    chemistry_tracers_cloud_levels: FloatField_ConvectionTracers,
    chemistry_tracers_sc_updraft: FloatField_ConvectionTracers,
    chemistry_tracers_sc_downdraft: FloatField_ConvectionTracers,
    chemistry_tracers_pw_updraft: FloatField_ConvectionTracers,
    chemistry_tracers_pw_downdraft: FloatField_ConvectionTracers,
    dd_tracers: FloatField_ConvectionTracers,
    aa: FloatField,
    bb: FloatField,
    cc: FloatField,
    plume: Int,
):
    """Advect tracers up/down the column - part 1/2 (before tridiag).

    Args:
        error_code (IntFieldIJ_Plume)
        cloud_top_level (IntFieldIJ_Plume)
        p_cloud_levels_forced (FloatField_Plume)
        environment_massflux (FloatField)
        normalized_massflux_updraft_forced (FloatField_Plume)
        normalized_massflux_downdraft_forced (FloatField_Plume)
        epsilon_forced (FloatFieldIJ_Plume)
        chemistry_tracers (FloatField_ConvectionTracers)
        chemistry_tracers_output (FloatField_ConvectionTracers_Plume)
        chemistry_tracers_cloud_levels (FloatField_ConvectionTracers)
        chemistry_tracers_sc_updraft (FloatField_ConvectionTracers)
        chemistry_tracers_sc_downdraft (FloatField_ConvectionTracers)
        chemistry_tracers_pw_updraft (FloatField_ConvectionTracers)
        chemistry_tracers_pw_downdraft (FloatField_ConvectionTracers)
        dd_tracers (FloatField_ConvectionTracers)
        aa (FloatField)
        bb (FloatField)
        cc (FloatField)
        plume (Int)
    """
    from __externals__ import (
        ALP1,
        DTIME,
        NUMBER_OF_TRACERS,
        USE_TRACER_EVAPORATION,
        USE_TRACER_SCAVENGE,
        USE_FLUX_FORM,
        USE_FCT,
    )

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            alp0: FloatFieldIJ = 0.0

    # flux form + source/sink terms + time explicit + FCT
    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0 and USE_FLUX_FORM == 1 and ALP1 == 0:
            option_not_implemented = True

    # flux form + source/sink terms + time explicit/implicit + upstream
    with computation(FORWARD), interval(...):
        if error_code[0, 0][plume] == 0 and USE_FLUX_FORM == 1 and ALP1 > 0.0:
            alp0 = 1.0 - ALP1
            if K <= cloud_top_level[0, 0][plume] + 1:
                fp = 0.5 * (environment_massflux + abs(environment_massflux))
                fm = 0.5 * (environment_massflux - abs(environment_massflux))

    with computation(FORWARD), interval(...):
        if (
            error_code[0, 0][plume] == 0
            and USE_FLUX_FORM == 1
            and ALP1 > 0.0
            and K <= cloud_top_level[0, 0][plume]
        ):
            dp = 100.0 * (p_cloud_levels_forced[0, 0, 0][plume] - p_cloud_levels_forced[0, 0, 1][plume])
            beta = DTIME * constants.MAPL_GRAV / dp
            aa = ALP1 * beta * fm
            bb = 1.0 + ALP1 * beta * (fp - fm[0, 0, 1])
            cc = -ALP1 * beta * fp[0, 0, 1]

            tracer = 0
            while tracer < NUMBER_OF_TRACERS:
                dd_tracers[0, 0, 0][tracer] = (
                    chemistry_tracers[0, 0, 0][tracer]
                    - (
                        normalized_massflux_updraft_forced[0, 0, 1][plume]
                        * chemistry_tracers_sc_updraft[0, 0, 1][tracer]
                        - normalized_massflux_updraft_forced[0, 0, 0][plume]
                        * chemistry_tracers_sc_updraft[0, 0, 0][tracer]
                    )
                    * beta
                    + (
                        (
                            normalized_massflux_downdraft_forced[0, 0, 1][plume]
                            * chemistry_tracers_sc_downdraft[0, 0, 1][tracer]
                            - normalized_massflux_downdraft_forced[0, 0, 0][plume]
                            * chemistry_tracers_sc_downdraft[0, 0, 0][tracer]
                        )
                        * beta
                        * epsilon_forced[0, 0][plume]
                    )
                )
                tracer += 1

            # include evaporation
            if USE_TRACER_EVAPORATION == 1 and plume != cumulus_parameterization_constants.SHALLOW:
                tracer = 0
                while tracer < NUMBER_OF_TRACERS:
                    # evaporated ( chemistry_tracer_pw_downdraft < 0 => E_dn > 0)
                    chemistry_tracers_output[0, 0, 0][plume, tracer] = (
                        chemistry_tracers_output[0, 0, 0][plume, tracer]
                        - 0.5
                        * epsilon_forced[0, 0][plume]
                        * (
                            normalized_massflux_downdraft_forced[0, 0, 0][plume]
                            * chemistry_tracers_pw_downdraft[0, 0, 0][tracer]
                            + normalized_massflux_downdraft_forced[0, 0, 1][plume]
                            * chemistry_tracers_pw_downdraft[0, 0, 1][tracer]
                        )
                        * beta
                    )
                    tracer += 1

            # include scavenging
            if USE_TRACER_SCAVENGE > 0 and plume != cumulus_parameterization_constants.SHALLOW:
                tracer = 0
                while tracer < NUMBER_OF_TRACERS:
                    # incorporated in rainfall (<0)
                    chemistry_tracers_output[0, 0, 0][plume, tracer] = (
                        chemistry_tracers_output[0, 0, 0][plume, tracer]
                        - 0.5
                        * (
                            normalized_massflux_updraft_forced[0, 0, 0][plume]
                            * chemistry_tracers_pw_updraft[0, 0, 0][tracer]
                            + normalized_massflux_updraft_forced[0, 0, 1][plume]
                            * chemistry_tracers_pw_updraft[0, 0, 1][tracer]
                        )
                        * beta
                    )
                    tracer += 1

            tracer = 0
            while tracer < NUMBER_OF_TRACERS:
                dd_tracers[0, 0, 0][tracer] = (
                    dd_tracers[0, 0, 0][tracer]
                    + chemistry_tracers_output[0, 0, 0][plume, tracer]
                    + alp0
                    * beta
                    * (
                        -fm * chemistry_tracers.at(K=max(0, K - 1), ddim=[tracer])
                        + (fm[0, 0, 1] - fp) * chemistry_tracers[0, 0, 0][tracer]
                        + fp[0, 0, 1] * chemistry_tracers[0, 0, 1][tracer]
                    )
                )
                tracer += 1

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0 and (USE_FLUX_FORM == 2 or USE_FLUX_FORM == 3):
            option_not_implemented = True


def update_after_tridiag(
    error_code: IntFieldIJ_Plume,
    cloud_top_level: IntFieldIJ_Plume,
    dd: FloatField_ConvectionTracers,
    chemistry_tracers: FloatField_ConvectionTracers,
    chemistry_tracers_output: FloatField_ConvectionTracers_Plume,
    tracer: Int,
    plume: Int,
):
    from __externals__ import DTIME

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0 and K <= cloud_top_level[0, 0][plume]:
            chemistry_tracers_output[0, 0, 0][plume, tracer] = (
                dd[0, 0, 0][tracer] - chemistry_tracers[0, 0, 0][tracer]
            ) / DTIME


def vertical_transport_part_2(
    error_code: IntFieldIJ_Plume,
    cloud_top_level: IntFieldIJ_Plume,
    p_cloud_levels_forced: FloatField_Plume,
    normalized_massflux_updraft_forced: FloatField_Plume,
    normalized_massflux_downdraft_forced: FloatField_Plume,
    epsilon_forced: FloatFieldIJ_Plume,
    chemistry_tracers: FloatField_ConvectionTracers,
    chemistry_tracers_output: FloatField_ConvectionTracers_Plume,
    chemistry_tracers_pw_updraft: FloatField_ConvectionTracers,
    chemistry_tracers_pw_downdraft: FloatField_ConvectionTracers,
    dd_tracers: FloatField_ConvectionTracers,
    residual: FloatFieldIJ_ConvectionTracers,
    plume: Int,
):
    """Advect tracers up/down the column - part 2/2 (after tridiag).

    Args:
        error_code (IntFieldIJ_Plume)
        cloud_top_level (IntFieldIJ_Plume)
        p_cloud_levels_forced (FloatField_Plume)
        normalized_massflux_updraft_forced (FloatField_Plume)
        normalized_massflux_downdraft_forced (FloatField_Plume)
        epsilon_forced (FloatFieldIJ_Plume)
        chemistry_tracers (FloatField_ConvectionTracers)
        chemistry_tracers_output (FloatField_ConvectionTracers_Plume)
        chemistry_tracers_pw_updraft (FloatField_ConvectionTracers)
        chemistry_tracers_pw_downdraft (FloatField_ConvectionTracers)
        dd_tracers (FloatField_ConvectionTracers)
        residual (FloatFieldIJ_ConvectionTracers)
        plume (Int)
    """
    from __externals__ import NUMBER_OF_TRACERS

    with computation(FORWARD), interval(0, 1):
        tracer = 0
        while tracer < NUMBER_OF_TRACERS:
            residual[0, 0][tracer] = 0.0
            tracer += 1

    with computation(FORWARD), interval(0, -1):
        if error_code[0, 0][plume] == 0 and K <= cloud_top_level[0, 0][plume]:
            tracer = 0
            while tracer < NUMBER_OF_TRACERS:
                dp = 100.0 * (p_cloud_levels_forced[0, 0, 0][plume] - p_cloud_levels_forced[0, 0, 1][plume])
                evap = (
                    -0.5
                    * (
                        normalized_massflux_downdraft_forced[0, 0, 0][plume]
                        * chemistry_tracers_pw_downdraft[0, 0, 0][tracer]
                        + normalized_massflux_downdraft_forced[0, 0, 1][plume]
                        * chemistry_tracers_pw_downdraft[0, 0, 1][tracer]
                    )
                    * constants.MAPL_GRAV
                    / dp
                    * epsilon_forced[0, 0][plume]
                )
                wetdep = (
                    0.5
                    * (
                        normalized_massflux_updraft_forced[0, 0, 0][plume]
                        * chemistry_tracers_pw_updraft[0, 0, 0][tracer]
                        + normalized_massflux_updraft_forced[0, 0, 1][plume]
                        * chemistry_tracers_pw_updraft[0, 0, 1][tracer]
                    )
                    * constants.MAPL_GRAV
                    / dp
                )

                residual[0, 0][tracer] = residual[0, 0][tracer] + (wetdep - evap) * dp / constants.MAPL_GRAV

                tracer += 1

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0:
            tracer = 0
            while tracer < NUMBER_OF_TRACERS:
                if residual[0, 0][tracer] < 0:
                    beta = constants.MAPL_GRAV / (
                        p_cloud_levels_forced.at(K=0, ddim=[plume])
                        - p_cloud_levels_forced.at(K=cloud_top_level[0, 0][plume] + 1, ddim=[plume])
                    )
                    chemistry_tracers_output[0, 0, 0][plume, tracer] = (
                        chemistry_tracers_output[0, 0, 0][plume, tracer] + residual[0, 0][tracer] * beta
                    )
                tracer += 1


class ColdPoolParameterization:
    def __init__(self, config: GF2020Config):
        if config.CONVECTION_TRACER == 1:
            raise NotImplementedError(
                "[NDSL] GF2020-->CumulusParameterization-->ColdPoolParameterization: The"
                "ColdPoolParameterization section has not been implemented. You should have been caught"
                "before getting here by the config checker. Beware, something likely failing in the config"
                "checker as well - you may be unknowingly calling other untested/unimplemented sections."
            )

    def __call__(self, *args, **kwds):
        pass


class AtmosphericComposition(NDSLRuntime):
    """Consider the effects of and on chemistry/tracers by convection."""

    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GF2020Config,
        cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
    ):
        # init NDSLRuntime
        super().__init__(stencil_factory)

        # make configuration visible at runtime
        self.config = config
        self.cumulus_parameterization_config = cumulus_parameterization_config

        quantity_factory.add_data_dimensions({"convection_tracers": config.NUMBER_OF_TRACERS})
        self._aa: Local = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._bb: Local = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._cc: Local = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._dd_tracers: Local = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM, "convection_tracers"], "n/a")
        self._tracer_cloud_boundary: Local = quantity_factory.zeros(
            [X_DIM, Y_DIM, "convection_tracers"], "n/a"
        )
        self._residual: Local = quantity_factory.zeros([X_DIM, Y_DIM, "convection_tracers"], "n/a")

        # construct stencils and functions
        self._environment_cloud_levels_chemistry = stencil_factory.from_dims_halo(
            func=environment_cloud_levels_chemistry,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"NUMBER_OF_TRACERS": config.NUMBER_OF_TRACERS, "CLOUD_LEVEL_OPTION": 2},
        )

        self._updraft_chemistry = stencil_factory.from_dims_halo(
            func=updraft_chemistry,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "NUMBER_OF_TRACERS": config.NUMBER_OF_TRACERS,
                "BOUNDARY_CONDITION_METHOD": cumulus_parameterization_config.BOUNDARY_CONDITION_METHOD,
                "USE_TRACER_SCAVENGE": cumulus_parameterization_config.USE_TRACER_SCAVENGE,
            },
        )

        self._downdraft_chemistry = stencil_factory.from_dims_halo(
            func=downdraft_chemistry,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "NUMBER_OF_TRACERS": config.NUMBER_OF_TRACERS,
                "USE_TRACER_EVAPORATION": cumulus_parameterization_config.USE_TRACER_EVAPORATION,
            },
        )

        self._vertical_transport_part_1 = stencil_factory.from_dims_halo(
            func=vertical_transport_part_1,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "ALP1": cumulus_parameterization_config.ALP1,
                "DTIME": cumulus_parameterization_config.DTIME,
                "NUMBER_OF_TRACERS": config.NUMBER_OF_TRACERS,
                "USE_TRACER_EVAPORATION": cumulus_parameterization_config.USE_TRACER_EVAPORATION,
                "USE_TRACER_SCAVENGE": cumulus_parameterization_config.USE_TRACER_SCAVENGE,
                "USE_FLUX_FORM": cumulus_parameterization_config.USE_FLUX_FORM,
                "USE_FCT": cumulus_parameterization_config.USE_FCT,
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

        self._vertical_transport_part_2 = stencil_factory.from_dims_halo(
            func=vertical_transport_part_2,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"NUMBER_OF_TRACERS": config.NUMBER_OF_TRACERS},
        )

    def __call__(
        self,
        error_code: Quantity,
        cloud_top_level: Quantity,
        updraft_origin_level: Quantity,
        downdraft_origin_level: Quantity,
        ocean_fraction: Quantity,
        p_forced: Quantity,
        p_cloud_levels_forced: Quantity,
        geopotential_height_cloud_levels: Quantity,
        environment_massflux: Quantity,
        normalized_massflux_updraft_forced: Quantity,
        normalized_massflux_downdraft_forced: Quantity,
        mass_entrainment_updraft_forced: Quantity,
        mass_detrainment_updraft_forced: Quantity,
        mass_entrainment_downdraft_forced: Quantity,
        mass_detrainment_downdraft_forced: Quantity,
        vertical_velocity_3d: Quantity,
        total_normalized_integrated_condensate_forced: Quantity,
        total_normalized_integrated_evaporate_forced: Quantity,
        evaporate_in_downdraft_forced: Quantity,
        epsilon_forced: Quantity,
        chemistry_tracers: Quantity,
        chemistry_tracers_output: Quantity,
        chemistry_tracers_cloud_levels: Quantity,
        chemistry_tracers_sc_updraft: Quantity,
        chemistry_tracers_sc_downdraft: Quantity,
        chemistry_tracers_pw_updraft: Quantity,
        chemistry_tracers_pw_downdraft: Quantity,
        chemistry_tracers_total_pw_updraft: Quantity,
        chemistry_tracers_total_pw_downdraft: Quantity,
        convection_tracers: ConvectionTracers,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        """Apply the effects of convection to the convection/chemistry tracers.

        Args:
            error_code (Quantity)
            cloud_top_level (Quantity)
            updraft_origin_level (Quantity)
            downdraft_origin_level (Quantity)
            ocean_fraction (Quantity)
            p_forced (Quantity)
            p_cloud_levels_forced (Quantity)
            geopotential_height_cloud_levels (Quantity)
            environment_massflux (Quantity)
            normalized_massflux_updraft_forced (Quantity)
            normalized_massflux_downdraft_forced (Quantity)
            mass_entrainment_updraft_forced (Quantity)
            mass_detrainment_updraft_forced (Quantity)
            mass_entrainment_downdraft_forced (Quantity)
            mass_detrainment_downdraft_forced (Quantity)
            vertical_velocity_3d (Quantity)
            total_normalized_integrated_condensate_forced (Quantity)
            total_normalized_integrated_evaporate_forced (Quantity)
            evaporate_in_downdraft_forced (Quantity)
            epsilon_forced (Quantity)
            chemistry_tracers (Quantity)
            chemistry_tracers_output (Quantity)
            chemistry_tracers_cloud_levels (Quantity)
            chemistry_tracers_sc_updraft (Quantity)
            chemistry_tracers_sc_downdraft (Quantity)
            chemistry_tracers_pw_updraft (Quantity)
            chemistry_tracers_pw_downdraft (Quantity)
            chemistry_tracers_total_pw_updraft (Quantity)
            chemistry_tracers_total_pw_downdraft (Quantity)
            convection_tracers (ConvectionTracers): Collection of tracers from the rest of the model which
                will be updated within convection. These may come from a variaty of sources, and need to be
                collected into the expected ConvectionTracers data type before being passed down.
            plume_dependent_constants (GF2020PlumeDependentConstants)

        Raises:
            NotImplementedError: _description_
        """

        # 1) get mass mixing ratios at the cloud levels
        self._environment_cloud_levels_chemistry(
            error_code=error_code,
            chemistry_tracers=chemistry_tracers,
            chemistry_tracers_cloud_levels=chemistry_tracers_cloud_levels,
            plume=plume_dependent_constants.PLUME_INDEX,
        )

        if (convection_tracers.vect_hcts.field[:, 0] > Float(1.0e-6)).any():
            raise NotImplementedError(
                "[NDSL] GF2020-->CumulusParameterization-->AtmosphericComposition called with "
                "convection_tracers.vect_hcts > 1.0e-6. This condition requires an unimplemented option in "
                "updraft_chemistry. Please implement, then remove this error to continue."
            )

        # 2) determine in-cloud tracer mixing ratios
        # a) updraft chemistry
        self._updraft_chemistry(
            error_code=error_code,
            updraft_origin_level=updraft_origin_level,
            cloud_top_level=cloud_top_level,
            ocean_fraction=ocean_fraction,
            p_forced=p_forced,
            p_cloud_levels_forced=p_cloud_levels_forced,
            geopotential_height_cloud_levels=geopotential_height_cloud_levels,
            vertical_velocity_3d=vertical_velocity_3d,
            mass_entrainment_updraft_forced=mass_entrainment_updraft_forced,
            mass_detrainment_updraft_forced=mass_detrainment_updraft_forced,
            normalized_massflux_updraft_forced=normalized_massflux_updraft_forced,
            chemistry_tracers=chemistry_tracers,
            chemistry_tracers_sc_updraft=chemistry_tracers_sc_updraft,
            chemistry_tracers_cloud_levels=chemistry_tracers_cloud_levels,
            chemistry_tracers_pw_updraft=chemistry_tracers_pw_updraft,
            chemistry_tracers_total_pw_updraft=chemistry_tracers_total_pw_updraft,
            convection_tracers_vect_hcts=convection_tracers.vect_hcts,
            convection_tracers_fscav=convection_tracers.fscav,
            convection_tracers_use_gcc_washout=convection_tracers.use_gcc_washout,
            tracer_cloud_boundary=self._tracer_cloud_boundary,
            AVERAGE_LAYER_DEPTH=plume_dependent_constants.AVERAGE_LAYER_DEPTH,
            plume=plume_dependent_constants.PLUME_INDEX,
        )

        # b) downdraft chemistry
        self._downdraft_chemistry(
            error_code=error_code,
            downdraft_origin_level=downdraft_origin_level,
            p_cloud_levels_forced=p_cloud_levels_forced,
            evaporate_in_downdraft_forced=evaporate_in_downdraft_forced,
            total_normalized_integrated_evaporate_forced=total_normalized_integrated_evaporate_forced,
            total_normalized_integrated_condensate_forced=total_normalized_integrated_condensate_forced,
            normalized_massflux_downdraft_forced=normalized_massflux_downdraft_forced,
            mass_entrainment_downdraft_forced=mass_entrainment_downdraft_forced,
            mass_detrainment_downdraft_forced=mass_detrainment_downdraft_forced,
            chemistry_tracers=chemistry_tracers,
            chemistry_tracers_cloud_levels=chemistry_tracers_cloud_levels,
            chemistry_tracers_sc_downdraft=chemistry_tracers_sc_downdraft,
            chemistry_tracers_pw_downdraft=chemistry_tracers_pw_downdraft,
            chemistry_tracers_total_pw_downdraft=chemistry_tracers_total_pw_downdraft,
            chemistry_tracers_total_pw_updraft=chemistry_tracers_total_pw_updraft,
            plume=plume_dependent_constants.PLUME_INDEX,
        )

        # 3) determine the vertical transport including mixing, scavenging and evaporation
        self._vertical_transport_part_1(
            error_code=error_code,
            cloud_top_level=cloud_top_level,
            p_cloud_levels_forced=p_cloud_levels_forced,
            environment_massflux=environment_massflux,
            normalized_massflux_updraft_forced=normalized_massflux_updraft_forced,
            normalized_massflux_downdraft_forced=normalized_massflux_downdraft_forced,
            epsilon_forced=epsilon_forced,
            chemistry_tracers=chemistry_tracers,
            chemistry_tracers_output=chemistry_tracers_output,
            chemistry_tracers_cloud_levels=chemistry_tracers_cloud_levels,
            chemistry_tracers_sc_updraft=chemistry_tracers_sc_updraft,
            chemistry_tracers_sc_downdraft=chemistry_tracers_sc_downdraft,
            chemistry_tracers_pw_updraft=chemistry_tracers_pw_updraft,
            chemistry_tracers_pw_downdraft=chemistry_tracers_pw_downdraft,
            dd_tracers=self._dd_tracers,
            aa=self._aa,
            bb=self._bb,
            cc=self._cc,
            plume=plume_dependent_constants.PLUME_INDEX,
        )

        if (
            self.cumulus_parameterization_config.USE_FLUX_FORM == 1
            and self.cumulus_parameterization_config.ALP1 > 0.0
        ):
            for tracer in range(self.config.NUMBER_OF_TRACERS):
                # NOTE, this code is an extention of the vertical transport section for option
                # which executes when (USE_FLUX_FORM == 1 and ALP1 > 0.0)
                self._tridiag(
                    m=cloud_top_level,
                    a=self._aa,
                    b=self._bb,
                    c=self._cc,
                    f=self._dd_tracers.data[:, :, :, tracer],
                    error_code=error_code,
                    plume=plume_dependent_constants.PLUME_INDEX,
                )

                self._update_after_tridiag(
                    error_code=error_code,
                    cloud_top_level=cloud_top_level,
                    dd=self._dd_tracers,
                    chemistry_tracers=chemistry_tracers,
                    chemistry_tracers_output=chemistry_tracers_output,
                    tracer=Int(tracer),
                    plume=plume_dependent_constants.PLUME_INDEX,
                )

        self._vertical_transport_part_2(
            error_code=error_code,
            cloud_top_level=cloud_top_level,
            p_cloud_levels_forced=p_cloud_levels_forced,
            normalized_massflux_updraft_forced=normalized_massflux_updraft_forced,
            normalized_massflux_downdraft_forced=normalized_massflux_downdraft_forced,
            epsilon_forced=epsilon_forced,
            chemistry_tracers=chemistry_tracers,
            chemistry_tracers_output=chemistry_tracers_output,
            chemistry_tracers_pw_updraft=chemistry_tracers_pw_updraft,
            chemistry_tracers_pw_downdraft=chemistry_tracers_pw_downdraft,
            dd_tracers=self._dd_tracers,
            residual=self._residual,
            plume=plume_dependent_constants.PLUME_INDEX,
        )
