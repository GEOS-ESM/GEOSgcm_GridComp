import pyMoist.constants as constants
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
from ndsl import NDSLRuntime, QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.gt4py import BACKWARD, FORWARD, PARALLEL, K, abs, computation, function, interval, max, min
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ, Int, IntFieldIJ
from ndsl.stencils.column_operations import column_min
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import (
    FloatField_ConvectionTracers,
    FloatField_ConvectionTracers_Plume,
    FloatField_Plume,
    FloatFieldIJ_Plume,
    IntFieldIJ_Plume,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.state import GF2020CumulusParameterizationState
from pyMoist.convection.GF_2020.locals import GF2020Locals
from pyMoist.convection.GF_2020.state import GF2020State
from pyMoist.convection_tracers import ConvectionTracers
from pyMoist.saturation_tables import (
    GlobalTable_saturation_tables,
    saturation_specific_humidity,
    saturation_specific_humidity_liquid_surface,
)
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.shared_incloud_processes import ice_fraction


def copy_from_cumulus_parameterization_state(
    pbl_time_scale_from_cumulus_parameterization: FloatFieldIJ,
    pbl_time_scale: FloatFieldIJ,
    cape_removal_time_scale_from_cumulus_parameterization: FloatFieldIJ,
    cape_removal_time_scale: FloatFieldIJ,
    cloud_workfunction_0_from_cumulus_parameterization: FloatFieldIJ,
    cloud_workfunction_0: FloatFieldIJ,
    cloud_workfunction_1: FloatFieldIJ,
    cloud_workfunction_1_from_cumulus_parameterization: FloatFieldIJ,
):
    """Copy fields from the cumulus parameterization state to the overarching model state.

    Args:
        pbl_time_scale_from_cumulus_parameterization (FloatFieldIJ)
        pbl_time_scale (FloatFieldIJ)
        cape_removal_time_scale_from_cumulus_parameterization (FloatFieldIJ)
        cape_removal_time_scale (FloatFieldIJ)
        cloud_workfunction_0_from_cumulus_parameterization (FloatFieldIJ)
        cloud_workfunction_0 (FloatFieldIJ)
        cloud_workfunction_1 (FloatFieldIJ)
        cloud_workfunction_1_from_cumulus_parameterization (FloatFieldIJ)
    """
    with computation(FORWARD), interval(0, 1):
        pbl_time_scale = pbl_time_scale_from_cumulus_parameterization
        cape_removal_time_scale = cape_removal_time_scale_from_cumulus_parameterization

        cloud_workfunction_0 = cloud_workfunction_0_from_cumulus_parameterization
        cloud_workfunction_1 = cloud_workfunction_1_from_cumulus_parameterization


def flag_computed_plumes_and_columns(
    error_code: IntFieldIJ_Plume,
    do_this_column: IntFieldIJ,
):
    """Flag which plumes were computed, and which columns made it through the entire scheme
    (error code = 0) for one or more columns

    Args:
        error_code (IntFieldIJ_Plume): contains error codes from cumulus parameterization core
        do_this_column (IntFieldIJ): mask which shows "successful" columns
    """
    from __externals__ import ENABLE_DEEP, ENABLE_MID, ENABLE_SHALLOW

    with computation(FORWARD), interval(0, 1):
        if ENABLE_SHALLOW == 0:
            error_code[0, 0][cumulus_parameterization_constants.SHALLOW] = -99
        if ENABLE_MID == 0:
            error_code[0, 0][cumulus_parameterization_constants.MID] = -99
        if ENABLE_DEEP == 0:
            error_code[0, 0][cumulus_parameterization_constants.DEEP] = -99

    # generate a mask (which has been initialized to zero in GF2020Setup) which flags only columns
    # which made it through the entire cumulus parameterization scheme at one or more plumes
    with computation(FORWARD), interval(0, 1):
        plume = 0
        while plume < cumulus_parameterization_constants.NUMBER_OF_PLUMES:
            if error_code[0, 0][plume] == 0:
                do_this_column = 1
            plume += 1


def check_vapor_mixing_ratio(
    vapor_current: FloatField,
    dvapordt: FloatField_Plume,
    t_tendency_from_vapor: FloatField,
    fix_out_vapor: FloatFieldIJ,
    do_this_column: IntFieldIJ,
):
    """Ensure that output water vapor mixing ratio is a reasonable value

    Args:
        vapor_current (FloatField): vapor mixing ratio before GF2020CumulusParameterization call
        dvapordt (FloatField_Plume): vapor tendency output from GF2020CumulusParameterization call
        t_tendency_from_vapor (FloatField): temperature tendency from water vapor
        fix_out_vapor (FloatFieldIJ): modification to "fix" water vapor mixing ratio
        do_this_column (IntFieldIJ): mask which shows "successful" columns
    """
    from __externals__ import DT_MOIST, k_end

    with computation(PARALLEL), interval(0, -1):
        if do_this_column != 0:
            t_tendency_from_vapor = (
                dvapordt[0, 0, 0][cumulus_parameterization_constants.SHALLOW]
                + dvapordt[0, 0, 0][cumulus_parameterization_constants.MID]
                + dvapordt[0, 0, 0][cumulus_parameterization_constants.DEEP]
            )

            distance = vapor_current + t_tendency_from_vapor * DT_MOIST

    with computation(FORWARD), interval(0, 1):
        # prefill
        fix_out_vapor = 1.0

        if do_this_column != 0:
            min_value, min_index = column_min(distance, 0, k_end - 1)
            if min_value < 0.0:
                if abs(t_tendency_from_vapor.at(K=min_index) * DT_MOIST) < constants.FLOAT_TINY:
                    fix_out_vapor = 0.999999
                else:
                    fix_out_vapor = (
                        cumulus_parameterization_constants.smaller_qv - vapor_current.at(K=min_index)
                    ) / (t_tendency_from_vapor.at(K=min_index) * DT_MOIST)
                fix_out_vapor = max(0.0, min(fix_out_vapor, 1.0))


def feedback(
    fix_out_vapor: FloatFieldIJ,
    precip: FloatFieldIJ,
    precip_from_cumulus_parameterization: FloatFieldIJ_Plume,
    evaporation_sublimation_tendency: FloatField,
    evaporation_sublimation_tendency_from_cumulus_parameterization: FloatField,
    convective_precip_flux: FloatField,
    convective_precip_flux_from_cumulus_parameterization: FloatField,
    dtdt: FloatField,
    dtdt_from_cumulus_parameterization: FloatField_Plume,
    dvapordt: FloatField,
    dvapordt_from_cumulus_parameterization: FloatField_Plume,
    dcloudicedt: FloatField,
    dcloudicedt_from_cumulus_parameterization: FloatField_Plume,
    dudt: FloatField,
    dudt_from_cumulus_parameterization: FloatField_Plume,
    dvdt: FloatField,
    dvdt_from_cumulus_parameterization: FloatField_Plume,
    dlarge_scale_icedt: FloatField,
    dlarge_scale_icedt_from_cumulus_parameterization: FloatField_Plume,
    dconvective_icedt: FloatField,
    dconvective_icedt_from_cumulus_parameterization: FloatField_Plume,
    dlarge_scale_liquiddt: FloatField,
    dlarge_scale_liquiddt_from_cumulus_parameterization: FloatField_Plume,
    dconvective_liquiddt: FloatField,
    dconvective_liquiddt_from_cumulus_parameterization: FloatField_Plume,
    dlarge_scale_cloud_fractiondt: FloatField,
    dlarge_scale_cloud_fractiondt_from_cumulus_parameterization: FloatField_Plume,
    dconvective_cloud_fractiondt: FloatField,
    dconvective_cloud_fractiondt_from_cumulus_parameterization: FloatField_Plume,
    dbuoyancydt_from_cumulus_parameterization: FloatField_Plume,
    dbuoyancydt: FloatField,
    do_this_column: IntFieldIJ,
):
    """Feedback output from the cumulus parameterization core to the local state.
    This cannot be fed straight into the model as the local state and model state have opposing k-axis
    directions, and further calculations may be performed on one or more of these outputs before this flip
    and final exchange is performed.

    These could be condensed, but have been kept separate for readibility and clarity.

    Args:
        fix_out_vapor (FloatFieldIJ)
        evaporation_sublimation_tendency (FloatField)
        evaporation_sublimation_tendency_from_cumulus_parameterization (FloatField)
        convective_precip_flux (FloatField)
        convective_precip_flux_from_cumulus_parameterization (FloatField_Plume)
        dtdt (FloatField)
        dtdt_from_cumulus_parameterization (FloatField_Plume)
        dvapordt (FloatField)
        dvapordt_from_cumulus_parameterization (FloatField_Plume)
        dcloudicedt (FloatField)
        dcloudicedt_from_cumulus_parameterization (FloatField_Plume)
        dudt (FloatField)
        dudt_from_cumulus_parameterization (FloatField_Plume)
        dvdt (FloatField)
        dvdt_from_cumulus_parameterization (FloatField_Plume)
        dlarge_scale_icedt (FloatField)
        dlarge_scale_icedt_from_cumulus_parameterization (FloatField_Plume)
        dconvective_icedt (FloatField)
        dconvective_icedt_from_cumulus_parameterization (FloatField_Plume)
        dlarge_scale_liquiddt (FloatField)
        dlarge_scale_liquiddt_from_cumulus_parameterization (FloatField_Plume)
        dconvective_liquiddt (FloatField)
        dconvective_liquiddt_from_cumulus_parameterization (FloatField_Plume)
        dlarge_scale_cloud_fractiondt (FloatField)
        dlarge_scale_cloud_fractiondt_from_cumulus_parameterization (FloatField_Plume)
        dconvective_cloud_fractiondt (FloatField)
        dconvective_cloud_fractiondt_from_cumulus_parameterization (FloatField_Plume)
        dbuoyancydt_from_cumulus_parameterization (FloatField_Plume)
        dbuoyancydt (FloatField)
        do_this_column (IntFieldIJ)
    """
    from __externals__ import APPLY_SUBSIDENCE_MICROPHYSICS, CONVECTION_TRACER, USE_MOMENTUM_TRANSPORT, k_end

    with computation(FORWARD), interval(0, 1):
        if do_this_column != 0:
            precip = (
                precip_from_cumulus_parameterization[0, 0][cumulus_parameterization_constants.SHALLOW]
                + precip_from_cumulus_parameterization[0, 0][cumulus_parameterization_constants.MID]
                + precip_from_cumulus_parameterization[0, 0][cumulus_parameterization_constants.DEEP]
            ) * fix_out_vapor

    with computation(PARALLEL), interval(...):
        # combining effects of shallow + mid + deep convection
        if do_this_column != 0:
            # feedback the tendencies from convection
            dtdt = (
                dtdt_from_cumulus_parameterization[0, 0, 0][cumulus_parameterization_constants.SHALLOW]
                + dtdt_from_cumulus_parameterization[0, 0, 0][cumulus_parameterization_constants.MID]
                + dtdt_from_cumulus_parameterization[0, 0, 0][cumulus_parameterization_constants.DEEP]
            ) * fix_out_vapor
            dvapordt = (
                dvapordt_from_cumulus_parameterization[0, 0, 0][cumulus_parameterization_constants.SHALLOW]
                + dvapordt_from_cumulus_parameterization[0, 0, 0][cumulus_parameterization_constants.MID]
                + dvapordt_from_cumulus_parameterization[0, 0, 0][cumulus_parameterization_constants.DEEP]
            ) * fix_out_vapor
            dcloudicedt = (
                dcloudicedt_from_cumulus_parameterization[0, 0, 0][cumulus_parameterization_constants.SHALLOW]
                + dcloudicedt_from_cumulus_parameterization[0, 0, 0][cumulus_parameterization_constants.MID]
                + dcloudicedt_from_cumulus_parameterization[0, 0, 0][cumulus_parameterization_constants.DEEP]
            ) * fix_out_vapor
            evaporation_sublimation_tendency = (
                evaporation_sublimation_tendency_from_cumulus_parameterization * fix_out_vapor
            )  # already contains deep and mid amounts.

            # precip flux is only computed for deep plume
            convective_precip_flux = (
                convective_precip_flux_from_cumulus_parameterization * fix_out_vapor
            )  # ice/liquid precip flux of the deep plume

            if USE_MOMENTUM_TRANSPORT > 0:
                dudt = (
                    dudt_from_cumulus_parameterization[0, 0, 0][cumulus_parameterization_constants.SHALLOW]
                    + dudt_from_cumulus_parameterization[0, 0, 0][cumulus_parameterization_constants.MID]
                    + dudt_from_cumulus_parameterization[0, 0, 0][cumulus_parameterization_constants.DEEP]
                ) * fix_out_vapor
                dvdt = (
                    dvdt_from_cumulus_parameterization[0, 0, 0][cumulus_parameterization_constants.SHALLOW]
                    + dvdt_from_cumulus_parameterization[0, 0, 0][cumulus_parameterization_constants.MID]
                    + dvdt_from_cumulus_parameterization[0, 0, 0][cumulus_parameterization_constants.DEEP]
                ) * fix_out_vapor

            if APPLY_SUBSIDENCE_MICROPHYSICS == 1:
                dlarge_scale_icedt = (
                    dlarge_scale_icedt_from_cumulus_parameterization[0, 0, 0][
                        cumulus_parameterization_constants.SHALLOW
                    ]
                    + dlarge_scale_icedt_from_cumulus_parameterization[0, 0, 0][
                        cumulus_parameterization_constants.MID
                    ]
                    + dlarge_scale_icedt_from_cumulus_parameterization[0, 0, 0][
                        cumulus_parameterization_constants.DEEP
                    ]
                ) * fix_out_vapor
                dconvective_icedt = (
                    dconvective_icedt_from_cumulus_parameterization[0, 0, 0][
                        cumulus_parameterization_constants.SHALLOW
                    ]
                    + dconvective_icedt_from_cumulus_parameterization[0, 0, 0][
                        cumulus_parameterization_constants.MID
                    ]
                    + dconvective_icedt_from_cumulus_parameterization[0, 0, 0][
                        cumulus_parameterization_constants.DEEP
                    ]
                ) * fix_out_vapor
                dlarge_scale_liquiddt = (
                    dlarge_scale_liquiddt_from_cumulus_parameterization[0, 0, 0][
                        cumulus_parameterization_constants.SHALLOW
                    ]
                    + dlarge_scale_liquiddt_from_cumulus_parameterization[0, 0, 0][
                        cumulus_parameterization_constants.MID
                    ]
                    + dlarge_scale_liquiddt_from_cumulus_parameterization[0, 0, 0][
                        cumulus_parameterization_constants.DEEP
                    ]
                ) * fix_out_vapor
                dconvective_liquiddt = (
                    dconvective_liquiddt_from_cumulus_parameterization[0, 0, 0][
                        cumulus_parameterization_constants.SHALLOW
                    ]
                    + dconvective_liquiddt_from_cumulus_parameterization[0, 0, 0][
                        cumulus_parameterization_constants.MID
                    ]
                    + dconvective_liquiddt_from_cumulus_parameterization[0, 0, 0][
                        cumulus_parameterization_constants.DEEP
                    ]
                ) * fix_out_vapor
                dlarge_scale_cloud_fractiondt = (
                    dlarge_scale_cloud_fractiondt_from_cumulus_parameterization[0, 0, 0][
                        cumulus_parameterization_constants.SHALLOW
                    ]
                    + dlarge_scale_cloud_fractiondt_from_cumulus_parameterization[0, 0, 0][
                        cumulus_parameterization_constants.MID
                    ]
                    + dlarge_scale_cloud_fractiondt_from_cumulus_parameterization[0, 0, 0][
                        cumulus_parameterization_constants.DEEP
                    ]
                ) * fix_out_vapor
                dconvective_cloud_fractiondt = (
                    dconvective_cloud_fractiondt_from_cumulus_parameterization[0, 0, 0][
                        cumulus_parameterization_constants.SHALLOW
                    ]
                    + dconvective_cloud_fractiondt_from_cumulus_parameterization[0, 0, 0][
                        cumulus_parameterization_constants.MID
                    ]
                    + dconvective_cloud_fractiondt_from_cumulus_parameterization[0, 0, 0][
                        cumulus_parameterization_constants.DEEP
                    ]
                ) * fix_out_vapor

    with computation(PARALLEL), interval(...):
        if do_this_column != 0 and CONVECTION_TRACER == 1:
            dbuoyancydt = (
                dbuoyancydt_from_cumulus_parameterization[0, 0, 0][cumulus_parameterization_constants.SHALLOW]
                + dbuoyancydt_from_cumulus_parameterization[0, 0, 0][cumulus_parameterization_constants.MID]
                + dbuoyancydt_from_cumulus_parameterization[0, 0, 0][cumulus_parameterization_constants.DEEP]
            ) * fix_out_vapor


def feedback_tracers(
    tracer: Int,
    fix_out_vapor: FloatFieldIJ,
    dconvection_tracersdt: FloatField_ConvectionTracers,
    dconvection_tracersdt_from_cumulus_parameterization: FloatField_ConvectionTracers_Plume,
    chemistry_tracers_from_cumulus_parameterization: FloatField_ConvectionTracers,
    do_this_column: IntFieldIJ,
):
    """See documentation for "feedback". This is a separate stencil to accomodate the tracer data dimension.

    Args:
        tracer (Int)
        fix_out_vapor (FloatFieldIJ)
        dconvection_tracersdt (FloatField_ConvectionTracers)
        dconvection_tracersdt_from_cumulus_parameterization (FloatField_ConvectionTracers)
        chemistry_tracers_from_cumulus_parameterization (FloatField_ConvectionTracers)
        do_this_column (IntFieldIJ)
    """
    from __externals__ import DT_MOIST, USE_TRACER_TRANSPORT, k_end

    with computation(PARALLEL), interval(...):
        if do_this_column != 0 and USE_TRACER_TRANSPORT == 1:
            dconvection_tracersdt[0, 0, 0][tracer] = (
                dconvection_tracersdt_from_cumulus_parameterization[0, 0, 0][
                    cumulus_parameterization_constants.SHALLOW, tracer
                ]
                + dconvection_tracersdt_from_cumulus_parameterization[0, 0, 0][
                    cumulus_parameterization_constants.MID, tracer
                ]
                + dconvection_tracersdt_from_cumulus_parameterization[0, 0, 0][
                    cumulus_parameterization_constants.DEEP, tracer
                ]
            ) * fix_out_vapor

    with computation(PARALLEL), interval(0, -1):
        if do_this_column != 0 and USE_TRACER_TRANSPORT == 1:
            # constrain positivity for tracers
            distance = (
                chemistry_tracers_from_cumulus_parameterization[0, 0, 0][tracer]
                + dconvection_tracersdt[0, 0, 0][tracer] * DT_MOIST
            )

    with computation(FORWARD), interval(0, 1):
        # ensure temporary is initialized properly
        fix_tracers: FloatFieldIJ = 0.0

        if do_this_column != 0 and USE_TRACER_TRANSPORT == 1:
            # fixer for mass of tracer
            min_value, min_index = column_min(distance, 0, k_end - 1)
            if min_value < 0.0:
                if (
                    abs(dconvection_tracersdt.at(K=min_index, ddim=[tracer]) * DT_MOIST)
                    < constants.FLOAT_TINY
                ):
                    fix_tracers = 0.999999
                else:
                    fix_tracers = (
                        constants.FLOAT_TINY
                        - chemistry_tracers_from_cumulus_parameterization.at(K=min_index, ddim=[tracer])
                    ) / (dconvection_tracersdt.at(K=min_index, ddim=[tracer]) * DT_MOIST)
                if fix_tracers > 1.0 or fix_tracers < 0.0:
                    fix_tracers = 0.0

    with computation(PARALLEL), interval(0, -1):
        if do_this_column != 0 and USE_TRACER_TRANSPORT == 1:
            # apply fixer
            dconvection_tracersdt[0, 0, 0][tracer] = fix_tracers * dconvection_tracersdt[0, 0, 0][tracer]


def cloud_workfunction_output(
    precip: FloatFieldIJ_Plume,
    fix_out_vapor: FloatFieldIJ,
    cloud_workfunction_2: FloatFieldIJ,
    cloud_workfunction_3: FloatFieldIJ,
):
    """Fill cloud workfunction 2 and 3 with data from the cumulus parameterization core

    Args:
        precip (FloatFieldIJ_Plume)
        fix_out_vapor (FloatFieldIJ)
        cloud_workfunction_2 (FloatFieldIJ)
        cloud_workfunction_3 (FloatFieldIJ)
    """
    with computation(FORWARD), interval(0, 1):
        cloud_workfunction_2 = precip[0, 0][cumulus_parameterization_constants.MID] * fix_out_vapor
        cloud_workfunction_3 = precip[0, 0][cumulus_parameterization_constants.DEEP] * fix_out_vapor


def prefill_entrainment(
    lateral_entrainment_rate_shallow: FloatField,
    lateral_entrainment_rate_mid: FloatField,
    lateral_entrainment_rate_deep: FloatField,
):
    with computation(PARALLEL), interval(...):
        lateral_entrainment_rate_shallow = constants.MAPL_UNDEF
        lateral_entrainment_rate_mid = constants.MAPL_UNDEF
        lateral_entrainment_rate_deep = constants.MAPL_UNDEF


def feed_3d_model(
    do_this_column: IntFieldIJ,
    precip: FloatFieldIJ,
    convective_precipitation_GF: FloatFieldIJ,
    evaporation_sublimation_tendency_from_cumulus_parameterization: FloatField,
    evaporation_sublimation_tendency: FloatField,
    convective_precip_flux_from_cumulus_parameterization: FloatField,
    convective_precip_flux: FloatField,
    dconvection_tracersdt: FloatField_ConvectionTracers,
    convection_tracers: FloatField_ConvectionTracers,
):
    """Update the model state with feedback from the plume-independent cumulus parameterization core fields.

    Args:
        do_this_column (IntFieldIJ)
        precip (FloatFieldIJ)
        convective_precipitation_GF (FloatFieldIJ)
        evaporation_sublimation_tendency_from_cumulus_parameterization (FloatField)
        evaporation_sublimation_tendency (FloatField)
        convective_precip_flux_from_cumulus_parameterization (FloatField)
        convective_precip_flux (FloatField)
        dconvection_tracersdt (FloatField_ConvectionTracers)
        convection_tracers (FloatField_ConvectionTracers)
    """
    from __externals__ import DT_MOIST, FEED_3D_MODEL, ITEST, USE_TRACER_TRANSPORT, k_end

    with computation(FORWARD), interval(0, 1):
        if FEED_3D_MODEL == True and do_this_column != 0:
            # convective precip rate: mm/s = kg m-2 s-1
            if ITEST == 0:
                convective_precipitation_GF = 0.0
            else:
                convective_precipitation_GF = precip

    with computation(PARALLEL), interval(...):
        if FEED_3D_MODEL == True and do_this_column != 0:
            # sublimation/evaporation tendencies (kg/kg/s)
            evaporation_sublimation_tendency = (
                evaporation_sublimation_tendency_from_cumulus_parameterization.at(K=k_end - K)
            )
            # preciptation fluxes (kg/kg/s)
            convective_precip_flux = convective_precip_flux_from_cumulus_parameterization.at(K=k_end - K)

            if USE_TRACER_TRANSPORT == 1:
                # update tracer mass mixing ratios
                tracer = 0
                while tracer < constants.NUMBER_OF_TRACERS:
                    convection_tracers[0, 0, 0][tracer] = convection_tracers[0, 0, 0][
                        tracer
                    ] + DT_MOIST * dconvection_tracersdt.at(K=k_end - K, ddim=[tracer])

                    # final check for negative tracer mass mixing ratio
                    convection_tracers[0, 0, 0][tracer] = max(
                        convection_tracers[0, 0, 0][tracer], constants.FLOAT_TINY
                    )
                    tracer += 1


def feed_3d_model_from_plumes(
    plume: Int,
    do_this_column: IntFieldIJ,
    dz: FloatField,
    air_density: FloatField,
    p_flipped: FloatField,
    vapor_flipped: FloatField,
    dcloudicedt: FloatField,
    lateral_entrainment_rate_shallow: FloatField,
    lateral_entrainment_rate_mid: FloatField,
    lateral_entrainment_rate_deep: FloatField,
    entrainment_rate_from_cumulus_parameterization: FloatField_Plume,
    error_code_from_cumulus_parameterization: IntFieldIJ_Plume,
    cloud_top_level_from_cumulus_parameterization: IntFieldIJ_Plume,
    t_updraft_from_cumulus_parameterization: FloatField_Plume,
    mass_detrainment_updraft_forced_from_cumulus_parameterization: FloatField_Plume,
    mass_entrainment_updraft_forced_from_cumulus_parameterization: FloatField_Plume,
    normalized_massflux_updraft_forced_from_cumulus_parameterization: FloatField_Plume,
    cloud_liquid_after_rain_forced_from_cumulus_parameterization: FloatField_Plume,
    condensate_to_fall_forced_from_cumulus_parameterization: FloatField_Plume,
    epsilon_forced_from_cumulus_parameterization: FloatFieldIJ_Plume,
    evaporate_in_downdraft_forced_from_cumulus_parameterization: FloatField_Plume,
    total_water_flux_deep_convection: FloatField,
    convective_condensate_source: FloatField,
    mass_flux_cloud_base: FloatField,
    mass_flux_deep_updraft_detrained: FloatField,
    mass_flux_deep_updraft_interface: FloatField,
    entrainment_parameter: FloatField,
    updraft_vertical_velocity: FloatField,
    convective_condensate_grid_mean: FloatField,
    convective_precipitation_RAS: FloatField,
    updraft_areal_fraction: FloatField,
    esw: GlobalTable_saturation_tables,
    estlqu: Float,
):
    """Update the model state with feedback from the cumulus parameterization core.

    Args:
        plume (Float)
        do_this_column (IntFieldIJ)
        dz (FloatField)
        air_density (FloatField)
        p_flipped (FloatField)
        vapor_flipped (FloatField)
        dcloudicedt (FloatField)
        lateral_entrainment_rate_shallow (FloatField)
        lateral_entrainment_rate_mid (FloatField)
        lateral_entrainment_rate_deep (FloatField)
        error_code_from_cumulus_parameterization (FloatFieldIJ_Plume)
        cloud_top_level_from_cumulus_parameterization (FloatFieldIJ_Plume)
        lateral_entrainment_rate_from_cumulus_parameterization (FloatField_Plume)
        t_updraft_from_cumulus_parameterization (FloatField_Plume)
        mass_detrainment_updraft_forced_from_cumulus_parameterization (FloatField_Plume)
        mass_entrainment_updraft_forced_from_cumulus_parameterization (FloatField_Plume)
        normalized_massflux_updraft_forced_from_cumulus_parameterization (FloatField_Plume)
        cloud_liquid_after_rain_forced_from_cumulus_parameterization (FloatField_Plume)
        condensate_to_fall_forced_from_cumulus_parameterization (FloatField_Plume)
        epsilon_forced_from_cumulus_parameterization (FloatFieldIJ_Plume)
        evaporate_in_downdraft_forced_from_cumulus_parameterization (FloatField_Plume)
        total_water_flux_deep_convection (FloatField)
        convective_condensate_source (FloatField)
        mass_flux_cloud_base (FloatField)
        mass_flux_deep_updraft_detrained (FloatField)
        mass_flux_deep_updraft_interface (FloatField)
        entrainment_parameter (FloatField)
        updraft_vertical_velocity (FloatField)
        convective_condensate_grid_mean (FloatField)
        convective_precipitation_RAS (FloatField)
        esw (GlobalTable_saturation_tables)
        estlqu (Float)
    """
    from __externals__ import DT_MOIST, FEED_3D_MODEL, k_end

    with computation(FORWARD), interval(...):
        if (
            FEED_3D_MODEL == True
            and K >= k_end - cloud_top_level_from_cumulus_parameterization[0, 0][plume] - 1
            and error_code_from_cumulus_parameterization[0, 0][plume] == 0
        ):
            # deep convective total water flux
            # assumes .033 fractional area
            saturation_specific_humidity_updraft, _ = saturation_specific_humidity_liquid_surface(
                esw=esw,
                lqu=estlqu,
                t=t_updraft_from_cumulus_parameterization.at(K=k_end - K, ddim=[plume]),
                p=p_flipped.at(K=k_end - K),
                pressure_correction=True,
            )

            saturation_specific_humidity_updraft = (
                saturation_specific_humidity_updraft
                + cloud_liquid_after_rain_forced_from_cumulus_parameterization.at(K=k_end - K, ddim=[plume])
                / 0.033
            )
            total_water_flux_deep_convection[0, 0, 1] = total_water_flux_deep_convection[
                0, 0, 1
            ] + normalized_massflux_updraft_forced_from_cumulus_parameterization.at(
                K=k_end - K, ddim=[plume]
            ) * (
                saturation_specific_humidity_updraft - vapor_flipped.at(K=k_end - K)
            )

    with computation(PARALLEL), interval(...):
        if (
            FEED_3D_MODEL == True
            and K >= k_end - cloud_top_level_from_cumulus_parameterization[0, 0][plume] - 1
            and error_code_from_cumulus_parameterization[0, 0][plume] == 0
        ):
            if plume == cumulus_parameterization_constants.SHALLOW:
                # export entrainment rates used by GF
                lateral_entrainment_rate_shallow = entrainment_rate_from_cumulus_parameterization.at(
                    K=k_end - K, ddim=[plume]
                )
            if plume == cumulus_parameterization_constants.MID:
                # export entrainment rates used by GF
                lateral_entrainment_rate_mid = entrainment_rate_from_cumulus_parameterization.at(
                    K=k_end - K, ddim=[plume]
                )
            if plume == cumulus_parameterization_constants.DEEP:
                # export entrainment rates used by GF
                lateral_entrainment_rate_deep = entrainment_rate_from_cumulus_parameterization.at(
                    K=k_end - K, ddim=[plume]
                )
            # special treatment for convective_condensate_source
            # units = 'kg m-2 s-1',
            # dcloudicedt contains contributions from all plumes, so no need to accumulate across levels
            convective_condensate_source = dcloudicedt.at(K=k_end - K) * dz * air_density

            # detraining_mass_flux
            # units = 'kg m-2 s-1'
            mass_flux_deep_updraft_detrained = (
                mass_flux_deep_updraft_detrained
                + mass_detrainment_updraft_forced_from_cumulus_parameterization.at(K=k_end - K, ddim=[plume])
            )

            # cloud_base_mass_flux
            # units = 'kg m-2 s-1'
            mass_flux_cloud_base = (
                mass_flux_cloud_base
                + normalized_massflux_updraft_forced_from_cumulus_parameterization.at(
                    K=k_end - K, ddim=[plume]
                )
            )

            if (
                normalized_massflux_updraft_forced_from_cumulus_parameterization.at(K=k_end - K, ddim=[plume])
                > 1.0e-6
            ):
                # entrainment parameter
                # units ='m-1',
                entrainment_parameter = entrainment_parameter + (
                    mass_entrainment_updraft_forced_from_cumulus_parameterization.at(
                        K=k_end - K, ddim=[plume]
                    )
                    / (
                        dz
                        * normalized_massflux_updraft_forced_from_cumulus_parameterization.at(
                            K=k_end - K, ddim=[plume]
                        )
                    )
                )

                #     # updraft_vertical_velocity
                #     # units = 'hPa s-1',
                updraft_vertical_velocity = -0.2  # hPa/s =>  4 m/s

            # convective_condensate_grid_mean
            # units ='kg kg-1'
            convective_condensate_grid_mean = (
                convective_condensate_grid_mean
                + cloud_liquid_after_rain_forced_from_cumulus_parameterization.at(K=k_end - K, ddim=[plume])
            )

            #  not using progno-cloud to calculate the precip from the convective column
            #  if CNV_PRC3 will be send to progno-cloud, set CNPCPRATE = zero
            # 'convective_precipitation_from_GF',UNITS     = 'kg m-2 s-1',
            #  JAN/17/2017 : the units above are wrong. The correct are kg[precip water]/kg[air]
            convective_precipitation_RAS = convective_precipitation_RAS + (
                condensate_to_fall_forced_from_cumulus_parameterization.at(K=k_end - K, ddim=[plume])
                + epsilon_forced_from_cumulus_parameterization[0, 0][plume]
                * evaporate_in_downdraft_forced_from_cumulus_parameterization.at(K=k_end - K, ddim=[plume])
            ) * DT_MOIST / (dz * air_density)

            # updraft_area_fraction
            if (
                normalized_massflux_updraft_forced_from_cumulus_parameterization.at(K=k_end - K, ddim=[plume])
                > 1.0e-6
            ):
                updraft_areal_fraction = 0.033

    with computation(BACKWARD), interval(...):
        # this must be done in a separate computation because the offset write is incompatable with PARALLEL
        if FEED_3D_MODEL == True:
            if (
                K >= k_end - cloud_top_level_from_cumulus_parameterization[0, 0][plume] - 1
                and error_code_from_cumulus_parameterization[0, 0][plume] == 0
            ):
                # convective mass flux - only updraft
                # units = 'kg m-2 s-1'
                mass_flux_deep_updraft_interface[0, 0, 1] = mass_flux_deep_updraft_interface[
                    0, 0, 1
                ] + normalized_massflux_updraft_forced_from_cumulus_parameterization.at(
                    K=k_end - K, ddim=[plume]
                )


def update_convection_tracer(
    land_fraction: FloatFieldIJ,
    dbuoyancydt: FloatField,
    convection_tracer: FloatField,
):
    """Update the convection tracer field with output from the cumulus parameterizaiton core.

    Args:
        land_fraction (FloatFieldIJ)
        dbuoyancydt (FloatField)
        convection_tracer (FloatField)
    """
    from __externals__ import CONVECTION_TRACER, DT_MOIST, k_end

    # cold pool/"convection tracer"
    with computation(PARALLEL), interval(...):
        if CONVECTION_TRACER == 1:
            cold_pool_timescale = land_fraction * (6.0 / 3600.0) + (1 - land_fraction) * (6.0 / 3600.0)

            # sink term (exp decay 1h)
            sink = DT_MOIST * abs(convection_tracer) / cold_pool_timescale
            # source term
            # downdraft detrainment of buoyancy [ J/kg s^{-1}]
            # negative sign => source for updraft lifting
            source = DT_MOIST * min(0.0, dbuoyancydt.at(K=k_end - K))

            # 'continuity' equation = ADV + SRC - SINK
            convection_tracer = convection_tracer + source - sink


def update_outputs(
    p_flipped: FloatField,
    error_code_from_cumulus_parameterization: IntFieldIJ_Plume,
    cloud_top_level_from_cumulus_parameterization: IntFieldIJ_Plume,
    cloud_base_mass_flux_modified_from_cumulus_parameterization: FloatFieldIJ_Plume,
    normalized_massflux_updraft_forced_from_cumulus_parameterization: FloatField_Plume,
    normalized_massflux_downdraft_forced_from_cumulus_parameterization: FloatField_Plume,
    epsilon_forced_from_cumulus_parameterization: FloatFieldIJ_Plume,
    scale_dependence_factor_from_cumulus_parameterizaiton: FloatFieldIJ_Plume,
    pressure_shallow_convective_cloud_top: FloatFieldIJ,
    pressure_mid_convective_cloud_top: FloatFieldIJ,
    pressure_deep_convective_cloud_top: FloatFieldIJ,
    mass_flux_cloud_base_shallow: FloatFieldIJ,
    mass_flux_cloud_base_mid: FloatFieldIJ,
    mass_flux_cloud_base_deep: FloatFieldIJ,
    mass_flux_shallow: FloatField,
    mass_flux_mid: FloatField,
    mass_flux_deep_updraft: FloatField,
    mass_flux_deep_downdraft: FloatField,
    sigma_mid: FloatFieldIJ,
    sigma_deep: FloatFieldIJ,
):
    """Update various outputs in preparation for a push back to the model state.

    Args:
        p_flipped (FloatField)
        error_code_from_cumulus_parameterization (IntFieldIJ_Plume)
        cloud_top_level_from_cumulus_parameterization (IntFieldIJ_Plume)
        cloud_base_mass_flux_modified_from_cumulus_parameterization (FloatFieldIJ_Plume)
        normalized_massflux_updraft_forced_from_cumulus_parameterization (FloatField_Plume)
        normalized_massflux_downdraft_forced_from_cumulus_parameterization (FloatField_Plume)
        epsilon_forced_from_cumulus_parameterization (FloatFieldIJ_Plume)
        scale_dependence_factor_from_cumulus_parameterizaiton (FloatFieldIJ)
        pressure_shallow_convective_cloud_top (FloatFieldIJ)
        pressure_mid_convective_cloud_top (FloatFieldIJ)
        pressure_deep_convective_cloud_top (FloatFieldIJ)
        mass_flux_cloud_base_shallow (FloatFieldIJ)
        mass_flux_cloud_base_mid (FloatFieldIJ)
        mass_flux_cloud_base_deep (FloatFieldIJ)
        mass_flux_shallow (FloatField)
        mass_flux_mid (FloatField)
        mass_flux_deep_updraft (FloatField)
        mass_flux_deep_downdraft (FloatField)
        sigma_mid (FloatFieldIJ)
        sigma_deep (FloatFieldIJ)
    """
    from __externals__ import ENABLE_DEEP, ENABLE_MID, ENABLE_SHALLOW, k_end

    with computation(FORWARD), interval(0, 1):
        pressure_shallow_convective_cloud_top = constants.MAPL_UNDEF
        pressure_mid_convective_cloud_top = constants.MAPL_UNDEF
        pressure_deep_convective_cloud_top = constants.MAPL_UNDEF

    with computation(FORWARD), interval(0, 1):
        if (
            ENABLE_SHALLOW == True
            and error_code_from_cumulus_parameterization[0, 0][cumulus_parameterization_constants.SHALLOW]
            == 0
        ):
            pressure_shallow_convective_cloud_top = p_flipped.at(
                K=cloud_top_level_from_cumulus_parameterization[0, 0][
                    cumulus_parameterization_constants.SHALLOW
                ]
            )
            mass_flux_cloud_base_shallow = cloud_base_mass_flux_modified_from_cumulus_parameterization[0, 0][
                cumulus_parameterization_constants.SHALLOW
            ]

        if (
            ENABLE_MID == True
            and error_code_from_cumulus_parameterization[0, 0][cumulus_parameterization_constants.MID] == 0
        ):
            pressure_mid_convective_cloud_top = p_flipped.at(
                K=cloud_top_level_from_cumulus_parameterization[0, 0][cumulus_parameterization_constants.MID]
            )
            mass_flux_cloud_base_mid = cloud_base_mass_flux_modified_from_cumulus_parameterization[0, 0][
                cumulus_parameterization_constants.MID
            ]
            sigma_mid = scale_dependence_factor_from_cumulus_parameterizaiton[0, 0][
                cumulus_parameterization_constants.MID
            ]

        if (
            ENABLE_DEEP == True
            and error_code_from_cumulus_parameterization[0, 0][cumulus_parameterization_constants.DEEP] == 0
        ):
            pressure_deep_convective_cloud_top = p_flipped.at(
                K=cloud_top_level_from_cumulus_parameterization[0, 0][cumulus_parameterization_constants.DEEP]
            )
            mass_flux_cloud_base_deep = cloud_base_mass_flux_modified_from_cumulus_parameterization[0, 0][
                cumulus_parameterization_constants.DEEP
            ]
            sigma_deep = scale_dependence_factor_from_cumulus_parameterizaiton[0, 0][
                cumulus_parameterization_constants.DEEP
            ]

    with computation(PARALLEL), interval(...):
        if (
            ENABLE_SHALLOW == True
            and error_code_from_cumulus_parameterization[0, 0][cumulus_parameterization_constants.SHALLOW]
            == 0
        ):
            mass_flux_shallow = normalized_massflux_updraft_forced_from_cumulus_parameterization.at(
                K=k_end - K, ddim=[cumulus_parameterization_constants.SHALLOW]
            )

        if (
            ENABLE_MID == True
            and error_code_from_cumulus_parameterization[0, 0][cumulus_parameterization_constants.MID] == 0
        ):
            mass_flux_mid = normalized_massflux_updraft_forced_from_cumulus_parameterization.at(
                K=k_end - K, ddim=[cumulus_parameterization_constants.MID]
            )

        if (
            ENABLE_DEEP == True
            and error_code_from_cumulus_parameterization[0, 0][cumulus_parameterization_constants.DEEP] == 0
        ):
            mass_flux_deep_updraft = normalized_massflux_updraft_forced_from_cumulus_parameterization.at(
                K=k_end - K, ddim=[cumulus_parameterization_constants.DEEP]
            )
            mass_flux_deep_downdraft = (
                normalized_massflux_downdraft_forced_from_cumulus_parameterization.at(
                    K=k_end - K, ddim=[cumulus_parameterization_constants.DEEP]
                )
                * epsilon_forced_from_cumulus_parameterization[0, 0][cumulus_parameterization_constants.DEEP]
            )

    with computation(FORWARD), interval(0, 1):
        # for output purposes, error_code=0 (convection scheme runs without error) will be changed to 1
        # and all other internal errors are forced to zero
        plume = 0
        while plume < cumulus_parameterization_constants.NUMBER_OF_PLUMES:
            if error_code_from_cumulus_parameterization[0, 0][plume] == 0:
                error_code_from_cumulus_parameterization[0, 0][plume] = 1
            if error_code_from_cumulus_parameterization[0, 0][plume] > 1:
                error_code_from_cumulus_parameterization[0, 0][plume] = 0
            plume += 1


def update_tendencies(
    dtdt: FloatField,
    dvapordt: FloatField,
    dudt: FloatField,
    dvdt: FloatField,
    dtdt_deep_convection: FloatField,
    dvapordt_deep_convection: FloatField,
    dudt_deep_convection: FloatField,
    dvdt_deep_convection: FloatField,
):
    """Push tendency values computed in previous stencils (in the finalize class) back to the model state.

    Args:
        dtdt (FloatField)
        dvapordt (FloatField)
        dudt (FloatField)
        dvdt (FloatField)
        dtdt_deep_convection (FloatField)
        dvapordt_deep_convection (FloatField)
        dudt_deep_convection (FloatField)
        dvdt_deep_convection (FloatField)
    """
    from __externals__ import k_end

    with computation(PARALLEL), interval(...):
        # tendencies
        dtdt_deep_convection = dtdt.at(K=k_end - K)
        dvapordt_deep_convection = dvapordt.at(K=k_end - K)
        dudt_deep_convection = dudt.at(K=k_end - K)
        dvdt_deep_convection = dvdt.at(K=k_end - K)


def update_convection_codes(
    error_code_from_cumulus_parameterization: IntFieldIJ_Plume,
    convection_code_shallow: FloatFieldIJ,
    convection_code_mid: FloatFieldIJ,
    convection_code_deep: FloatFieldIJ,
):
    with computation(FORWARD), interval(0, 1):
        # error codes
        convection_code_shallow = error_code_from_cumulus_parameterization[0, 0][
            cumulus_parameterization_constants.SHALLOW
        ]
        convection_code_mid = error_code_from_cumulus_parameterization[0, 0][
            cumulus_parameterization_constants.MID
        ]
        convection_code_deep = error_code_from_cumulus_parameterization[0, 0][
            cumulus_parameterization_constants.DEEP
        ]


def update_state_with_tendencies(
    convection_fraction: FloatFieldIJ,
    surface_type: FloatFieldIJ,
    u: FloatField,
    v: FloatField,
    vapor: FloatField,
    t: FloatField,
    p: FloatField,
    p_kappa: FloatField,
    mass: FloatField,
    mass_flux_deep_updraft_detrained: FloatField,
    mass_flux_deep_updraft_interface: FloatField,
    total_cumulative_mass_flux_interface: FloatField,
    total_detraining_mass_flux: FloatField,
    dudt_deep_convection: FloatField,
    dvdt_deep_convection: FloatField,
    dvapordt_deep_convection: FloatField,
    dtdt_deep_convection: FloatField,
    dliquiddt_deep_convection: FloatField,
    dicedt_deep_convection: FloatField,
    dcloudfractiondt_deep_convection: FloatField,
    convective_condensate_source: FloatField,
    evaporation_sublimation_tendency: FloatField,
    convective_precip_flux: FloatField,
    sublimation_of_convective_precipitation: FloatField,
    evaporation_of_convective_precipitation: FloatField,
    ice_fraction_in_convective_tower: FloatField,
    ice_precip_flux_interface: FloatField,
    liquid_precip_flux_interface: FloatField,
    convective_liquid: FloatField,
    convective_ice: FloatField,
    convective_cloud_fraction: FloatField,
    convective_rainwater_source: FloatField,
    convective_precipitation_RAS: FloatField,
    ese: GlobalTable_saturation_tables,
    esx: GlobalTable_saturation_tables,
):
    """Update the model state (excluding the few fields which have already been updated in earlier
    stencils) with the output from the cumulus parameterization core.

    Containts a call to saturation_specific_humidity, which is techincally a port of the GEOS_QSAT function.
    In fortran

    Args:
        convection_fraction (FloatFieldIJ)
        surface_type (FloatFieldIJ)
        u (FloatField)
        v (FloatField)
        vapor (FloatField)
        t (FloatField)
        p (FloatField)
        p_kappa (FloatField)
        mass (FloatField)
        mass_flux_deep_updraft_detrained (FloatField)
        mass_flux_deep_updraft_interface (FloatField)
        total_cumulative_mass_flux_interface (FloatField)
        total_detraining_mass_flux (FloatField)
        dudt_deep_convection (FloatField)
        dvdt_deep_convection (FloatField)
        dvapordt_deep_convection (FloatField)
        dtdt_deep_convection (FloatField)
        dliquiddt_deep_convection (FloatField)
        dicedt_deep_convection (FloatField)
        dcloudfractiondt_deep_convection (FloatField)
        convective_condensate_source (FloatField)
        evaporation_sublimation_tendency (FloatField)
        convective_precip_flux (FloatField)
        sublimation_of_convective_precipitation (FloatField)
        evaporation_of_convective_precipitation (FloatField)
        ice_fraction_in_convective_tower (FloatField)
        ice_precip_flux_interface (FloatField)
        liquid_precip_flux_interface (FloatField)
        convective_liquid (FloatField)
        convective_ice (FloatField)
        convective_cloud_fraction (FloatField)
        convective_rainwater_source (FloatField)
        convective_precipitation_RAS (FloatField)
        ese (GlobalTable_saturation_tables)
        esx (GlobalTable_saturation_tables)
    """
    from __externals__ import DT_MOIST, FIX_CONVECTIVE_CLOUD, SCLM_DEEP

    with computation(PARALLEL), interval(...):
        u = u + dudt_deep_convection * DT_MOIST
        v = v + dvdt_deep_convection * DT_MOIST
        vapor = vapor + dvapordt_deep_convection * DT_MOIST
        t = t + dtdt_deep_convection * DT_MOIST

        # update deep cumulus liquid/ice/cloud fraction tendencies
        fraction_ice = ice_fraction(t, convection_fraction, surface_type)
        condensate_per_mass = convective_condensate_source / mass
        dliquiddt_deep_convection = (1.0 - fraction_ice) * condensate_per_mass
        dicedt_deep_convection = fraction_ice * condensate_per_mass
        dcloudfractiondt_deep_convection = mass_flux_deep_updraft_detrained * SCLM_DEEP / mass

        # sublimation/evaporation tendencies (kg/kg/s)
        sublimation_of_convective_precipitation = evaporation_sublimation_tendency * fraction_ice
        evaporation_of_convective_precipitation = evaporation_sublimation_tendency * (1.0 - fraction_ice)

    with computation(FORWARD), interval(...):
        # preciptation fluxes (kg/kg/s)
        ice_precip_flux_interface[0, 0, 1] = convective_precip_flux * fraction_ice
        liquid_precip_flux_interface[0, 0, 1] = convective_precip_flux * (1.0 - fraction_ice)

    with computation(PARALLEL), interval(...):
        # add liquid/ice/cloud fraction tendencies
        convective_liquid = convective_liquid + dliquiddt_deep_convection * DT_MOIST
        convective_ice = convective_ice + dicedt_deep_convection * DT_MOIST
        convective_cloud_fraction = max(
            min(convective_cloud_fraction + dcloudfractiondt_deep_convection * DT_MOIST, 1.0), 0.0
        )

        ice_fraction_in_convective_tower = fraction_ice

        # fix convective cloud fraction
        if FIX_CONVECTIVE_CLOUD == True:
            saturation_humidity, _ = saturation_specific_humidity(t, p, ese, esx)

            if convective_cloud_fraction < 1.0:
                modification = (vapor - saturation_humidity * convective_cloud_fraction) / (
                    1.0 - convective_cloud_fraction
                )
            min_saturation_humidity = 0.001
            if (
                modification - min_saturation_humidity * saturation_humidity
            ) < 0.0 and convective_cloud_fraction > 0.0:
                convective_cloud_fraction = (vapor - min_saturation_humidity * saturation_humidity) / (
                    saturation_humidity * (1.0 - min_saturation_humidity)
                )
            # if a suitable environment relative humidity cannot be made then destroy anvil
            if convective_cloud_fraction < 0.0:
                convective_cloud_fraction = 0.0
                dliquiddt_deep_convection = dliquiddt_deep_convection - (convective_liquid) / DT_MOIST
                dicedt_deep_convection = dicedt_deep_convection - (convective_ice) / DT_MOIST
                dvapordt_deep_convection = (
                    dvapordt_deep_convection + (convective_liquid + convective_ice) / DT_MOIST
                )
                vapor = vapor + (convective_liquid + convective_ice)
                modification = (
                    constants.MAPL_ALHL * convective_liquid + constants.MAPL_ALHS * convective_ice
                ) / constants.MAPL_CP
                dtdt_deep_convection = dtdt_deep_convection - modification / DT_MOIST
                t = t - modification
                convective_liquid = 0.0
                convective_ice = 0.0

        total_cumulative_mass_flux_interface = (
            total_cumulative_mass_flux_interface + mass_flux_deep_updraft_interface
        )
        total_detraining_mass_flux = total_detraining_mass_flux + mass_flux_deep_updraft_detrained
        convective_rainwater_source = convective_precipitation_RAS / DT_MOIST


class GF2020Finalize(NDSLRuntime):
    """This class performs the entire finalization sequence for the GF2020 convection parameterization scheme

    In the source Fortran codee, this code is split across three subroutines nested as follows:

    - GF_Run
        - GF2020_INTERFACE
            - GF2020_DRV beginning at the end (outside) of the of the plume loop (after the call to CUP_gf)

    This python implementation simplifies this structure by bringing all setup calculations to the same level.
    An effort has been made to reduce duplicate/unnecessary locals where possible, but some have been retained
    for the sake of readibility.
    """

    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GF2020Config,
        cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
        saturation_tables: SaturationVaporPressureTable,
    ):
        super().__init__(stencil_factory)

        # make status of plumes visible at runtime
        self._plume_status = [
            cumulus_parameterization_config.ENABLE_SHALLOW,
            cumulus_parameterization_config.ENABLE_MID,
            cumulus_parameterization_config.ENABLE_DEEP,
        ]

        # make saturation table data visible at runtime
        # NOTE: this is an orchestration workaround. Direct call to
        #   `self.tables.X` fails closure capture for
        #   argument reconstruction at call time
        self._ese = saturation_tables.ese
        self._esw = saturation_tables.esw
        self._esx = saturation_tables.esx
        self._estfrz = saturation_tables.frz
        self._estlqu = saturation_tables.lqu

        # construct stencils
        self._copy_from_cumulus_parameterization_state = stencil_factory.from_dims_halo(
            func=copy_from_cumulus_parameterization_state,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._flag_computed_plumes_and_columns = stencil_factory.from_dims_halo(
            func=flag_computed_plumes_and_columns,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "ENABLE_SHALLOW": cumulus_parameterization_config.ENABLE_SHALLOW,
                "ENABLE_MID": cumulus_parameterization_config.ENABLE_MID,
                "ENABLE_DEEP": cumulus_parameterization_config.ENABLE_DEEP,
            },
        )

        self._check_vapor_mixing_ratio = stencil_factory.from_dims_halo(
            func=check_vapor_mixing_ratio,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"DT_MOIST": config.DT_MOIST},
        )

        self._feedback = stencil_factory.from_dims_halo(
            func=feedback,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "USE_MOMENTUM_TRANSPORT": config.USE_MOMENTUM_TRANSPORT,
                "APPLY_SUBSIDENCE_MICROPHYSICS": config.APPLY_SUBSIDENCE_MICROPHYSICS,
                "CONVECTION_TRACER": config.CONVECTION_TRACER,
            },
        )

        self._feedback_tracers = stencil_factory.from_dims_halo(
            func=feedback_tracers,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "USE_TRACER_TRANSPORT": config.USE_TRACER_TRANSPORT,
                "DT_MOIST": config.DT_MOIST,
            },
        )

        self._cloud_workfunction_output = stencil_factory.from_dims_halo(
            func=cloud_workfunction_output,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._prefill_entrainment = stencil_factory.from_dims_halo(
            func=prefill_entrainment,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._feed_3d_model = stencil_factory.from_dims_halo(
            func=feed_3d_model,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "FEED_3D_MODEL": cumulus_parameterization_config.FEED_3D_MODEL,
                "ITEST": cumulus_parameterization_config.ITEST,
                "USE_TRACER_TRANSPORT": config.USE_TRACER_TRANSPORT,
                "DT_MOIST": config.DT_MOIST,
            },
        )

        self._feed_3d_model_from_plumes = stencil_factory.from_dims_halo(
            func=feed_3d_model_from_plumes,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "FEED_3D_MODEL": cumulus_parameterization_config.FEED_3D_MODEL,
                "ITEST": cumulus_parameterization_config.ITEST,
                "USE_TRACER_TRANSPORT": config.USE_TRACER_TRANSPORT,
                "DT_MOIST": config.DT_MOIST,
            },
        )

        self._update_convection_tracer = stencil_factory.from_dims_halo(
            func=update_convection_tracer,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"CONVECTION_TRACER": config.CONVECTION_TRACER, "DT_MOIST": config.DT_MOIST},
        )

        self._update_outputs = stencil_factory.from_dims_halo(
            func=update_outputs,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "ENABLE_SHALLOW": cumulus_parameterization_config.ENABLE_SHALLOW,
                "ENABLE_MID": cumulus_parameterization_config.ENABLE_MID,
                "ENABLE_DEEP": cumulus_parameterization_config.ENABLE_DEEP,
            },
        )

        self._update_tendencies = stencil_factory.from_dims_halo(
            func=update_tendencies,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._update_convection_codes = stencil_factory.from_dims_halo(
            func=update_convection_codes,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._update_state_with_tendencies = stencil_factory.from_dims_halo(
            func=update_state_with_tendencies,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "SCLM_DEEP": config.SCLM_DEEP,
                "DT_MOIST": config.DT_MOIST,
                "FIX_CONVECTIVE_CLOUD": config.FIX_CONVECTIVE_CLOUD,
            },
        )

    def __call__(
        self,
        state: GF2020State,
        locals: GF2020Locals,
        cumulus_parameterization_state: GF2020CumulusParameterizationState,
        convection_tracers: ConvectionTracers,
    ):
        """Finish the GF2020 Cumulus Parameterization scheme

        Args:
            state (GF2020State): NDSL State containing all model fields required for GF2020.
            locals (GF2020Locals): NDSL LocalState containing all locals for GF2020.
            cumulus_parameterization_state (GF2020CumulusParameterizationState): NDSL State containing all
                fields required for the CumulusParameterization.
            convection_tracers (ConvectionTracers): Collection of tracers from the rest of the model which
                will be updated within convection. These may come from a variaty of sources, and need to be
                collected into the expected ConvectionTracers data type before being passed down.
        """

        self._copy_from_cumulus_parameterization_state(
            pbl_time_scale_from_cumulus_parameterization=cumulus_parameterization_state.output.pbl_time_scale,
            pbl_time_scale=state.pbl_time_scale,
            cape_removal_time_scale_from_cumulus_parameterization=cumulus_parameterization_state.output.cape_removal_time_scale,
            cape_removal_time_scale=state.cape_removal_time_scale,
            cloud_workfunction_0_from_cumulus_parameterization=cumulus_parameterization_state.output.cloud_workfunction_0,
            cloud_workfunction_0=state.cloud_workfunction_0,
            cloud_workfunction_1_from_cumulus_parameterization=cumulus_parameterization_state.output.cloud_workfunction_1,
            cloud_workfunction_1=state.cloud_workfunction_1,
        )

        self._flag_computed_plumes_and_columns(
            error_code=cumulus_parameterization_state.output.error_code,
            do_this_column=locals.do_this_column,
        )

        self._check_vapor_mixing_ratio(
            do_this_column=locals.do_this_column,
            vapor_current=locals.flipped_copy.vapor_current,
            dvapordt=cumulus_parameterization_state.output.dvapordt,
            t_tendency_from_vapor=locals.t_tendency_from_vapor,
            fix_out_vapor=locals.fix_out_vapor,
        )

        self._feedback(
            fix_out_vapor=locals.fix_out_vapor,
            precip=locals.precip,
            precip_from_cumulus_parameterization=cumulus_parameterization_state.output.precip,
            evaporation_sublimation_tendency=locals.evaporation_sublimation_tendency,
            evaporation_sublimation_tendency_from_cumulus_parameterization=cumulus_parameterization_state.output.evaporation_sublimation_tendency,
            convective_precip_flux=locals.convective_precip_flux,
            convective_precip_flux_from_cumulus_parameterization=cumulus_parameterization_state.output.convective_precip_flux,
            dtdt=locals.dtdt,
            dtdt_from_cumulus_parameterization=cumulus_parameterization_state.output.dtdt,
            dvapordt=locals.dvapordt,
            dvapordt_from_cumulus_parameterization=cumulus_parameterization_state.output.dvapordt,
            dcloudicedt=locals.dcloudicedt,
            dcloudicedt_from_cumulus_parameterization=cumulus_parameterization_state.output.dcloudicedt,
            dudt=locals.dudt,
            dudt_from_cumulus_parameterization=cumulus_parameterization_state.output.dudt,
            dvdt=locals.dvdt,
            dvdt_from_cumulus_parameterization=cumulus_parameterization_state.output.dvdt,
            dlarge_scale_icedt=locals.dlarge_scale_icedt,
            dlarge_scale_icedt_from_cumulus_parameterization=cumulus_parameterization_state.output.dlargescaleicedt,
            dconvective_icedt=locals.dconvective_icedt,
            dconvective_icedt_from_cumulus_parameterization=cumulus_parameterization_state.output.dconvectiveicedt,
            dlarge_scale_liquiddt=locals.dlarge_scale_liquiddt,
            dlarge_scale_liquiddt_from_cumulus_parameterization=cumulus_parameterization_state.output.dlargescaleliquiddt,
            dconvective_liquiddt=locals.dconvective_liquiddt,
            dconvective_liquiddt_from_cumulus_parameterization=cumulus_parameterization_state.output.dconvectiveliquiddt,
            dlarge_scale_cloud_fractiondt=locals.dlarge_scale_cloud_fractiondt,
            dlarge_scale_cloud_fractiondt_from_cumulus_parameterization=cumulus_parameterization_state.output.dlargescalecloudfractiondt,
            dconvective_cloud_fractiondt=locals.dconvective_cloud_fractiondt,
            dconvective_cloud_fractiondt_from_cumulus_parameterization=cumulus_parameterization_state.output.dconvectivecloudfractiondt,
            dbuoyancydt_from_cumulus_parameterization=cumulus_parameterization_state.output.dbuoyancydt,
            dbuoyancydt=locals.dbuoyancydt,
            do_this_column=locals.do_this_column,
        )

        for tracer in range(constants.NUMBER_OF_TRACERS):
            self._feedback_tracers(
                tracer=Int(tracer),
                fix_out_vapor=locals.fix_out_vapor,
                dconvection_tracersdt=locals.dconvection_tracersdt,
                dconvection_tracersdt_from_cumulus_parameterization=cumulus_parameterization_state.input_output.chemistry_tracers_output,
                chemistry_tracers_from_cumulus_parameterization=cumulus_parameterization_state.input_output.chemistry_tracers,
                do_this_column=locals.do_this_column,
            )

        self._cloud_workfunction_output(
            precip=cumulus_parameterization_state.output.precip,
            fix_out_vapor=locals.fix_out_vapor,
            cloud_workfunction_2=state.cloud_workfunction_2,
            cloud_workfunction_3=state.cloud_workfunction_3,
        )

        self._prefill_entrainment(
            lateral_entrainment_rate_shallow=state.lateral_entrainment_rate_shallow,
            lateral_entrainment_rate_mid=state.lateral_entrainment_rate_mid,
            lateral_entrainment_rate_deep=state.lateral_entrainment_rate_deep,
        )

        self._feed_3d_model(
            do_this_column=locals.do_this_column,
            precip=locals.precip,
            convective_precipitation_GF=state.convective_precipitation_GF,
            evaporation_sublimation_tendency_from_cumulus_parameterization=cumulus_parameterization_state.output.evaporation_sublimation_tendency,
            evaporation_sublimation_tendency=locals.evaporation_sublimation_tendency,
            convective_precip_flux_from_cumulus_parameterization=cumulus_parameterization_state.output.convective_precip_flux,
            convective_precip_flux=locals.convective_precip_flux,
            dconvection_tracersdt=locals.dconvection_tracersdt,
            convection_tracers=convection_tracers.tracers,
        )

        for plume in range(cumulus_parameterization_constants.NUMBER_OF_PLUMES):
            # Only call the function if the current plume is enabled self._plume_status is fixed to
            # [shallow, mid, deep] order. This conditional ensures the correct is checked in plume_status
            # and the correct plume is written in the stencil even if the overarching data order changes
            if plume == cumulus_parameterization_constants.SHALLOW:
                index = 0
            elif plume == cumulus_parameterization_constants.MID:
                index = 1
            elif plume == cumulus_parameterization_constants.DEEP:
                index = 2
            if self._plume_status[index] == 1:
                self._feed_3d_model_from_plumes(
                    plume=Int(index),
                    do_this_column=locals.do_this_column,
                    dz=locals.derived_state.dz,
                    air_density=locals.derived_state.air_density,
                    p_flipped=locals.flipped_copy.p,
                    vapor_flipped=locals.flipped_copy.vapor,
                    dcloudicedt=locals.dcloudicedt,
                    lateral_entrainment_rate_shallow=state.lateral_entrainment_rate_shallow,
                    lateral_entrainment_rate_mid=state.lateral_entrainment_rate_mid,
                    lateral_entrainment_rate_deep=state.lateral_entrainment_rate_deep,
                    entrainment_rate_from_cumulus_parameterization=cumulus_parameterization_state.output.entrainment_rate,
                    error_code_from_cumulus_parameterization=cumulus_parameterization_state.output.error_code,
                    cloud_top_level_from_cumulus_parameterization=cumulus_parameterization_state.output.cloud_top_level,
                    t_updraft_from_cumulus_parameterization=cumulus_parameterization_state.output.t_updraft,
                    mass_detrainment_updraft_forced_from_cumulus_parameterization=cumulus_parameterization_state.output.mass_detrainment_updraft_forced,
                    mass_entrainment_updraft_forced_from_cumulus_parameterization=cumulus_parameterization_state.output.mass_entrainment_updraft_forced,
                    normalized_massflux_updraft_forced_from_cumulus_parameterization=cumulus_parameterization_state.output.normalized_massflux_updraft_forced,
                    cloud_liquid_after_rain_forced_from_cumulus_parameterization=cumulus_parameterization_state.output.cloud_liquid_after_rain_forced,
                    condensate_to_fall_forced_from_cumulus_parameterization=cumulus_parameterization_state.output.condensate_to_fall_forced,
                    epsilon_forced_from_cumulus_parameterization=cumulus_parameterization_state.output.epsilon_forced,
                    evaporate_in_downdraft_forced_from_cumulus_parameterization=cumulus_parameterization_state.output.evaporate_in_downdraft_forced,
                    total_water_flux_deep_convection=state.total_water_flux_deep_convection_interface,
                    convective_condensate_source=state.convective_condensate_source,
                    mass_flux_cloud_base=state.mass_flux_cloud_base,
                    mass_flux_deep_updraft_detrained=state.mass_flux_deep_updraft_detrained,
                    mass_flux_deep_updraft_interface=state.mass_flux_deep_updraft_interface,
                    entrainment_parameter=state.entrainment_parameter,
                    updraft_vertical_velocity=state.updraft_vertical_velocity,
                    convective_condensate_grid_mean=state.convective_condensate_grid_mean,
                    convective_precipitation_RAS=state.convective_precipitation_RAS,
                    updraft_areal_fraction=state.updraft_areal_fraction,
                    esw=self._esw,
                    estlqu=self._estlqu,
                )

        self._update_convection_tracer(
            land_fraction=state.land_fraction,
            dbuoyancydt=locals.dbuoyancydt,
            convection_tracer=state.convection_tracer,
        )

        self._update_outputs(
            p_flipped=locals.flipped_copy.p,
            error_code_from_cumulus_parameterization=cumulus_parameterization_state.output.error_code,
            cloud_top_level_from_cumulus_parameterization=cumulus_parameterization_state.output.cloud_top_level,
            cloud_base_mass_flux_modified_from_cumulus_parameterization=cumulus_parameterization_state.output.cloud_base_mass_flux_modified,
            normalized_massflux_updraft_forced_from_cumulus_parameterization=cumulus_parameterization_state.output.normalized_massflux_updraft_forced,
            normalized_massflux_downdraft_forced_from_cumulus_parameterization=cumulus_parameterization_state.output.normalized_massflux_downdraft_forced,
            epsilon_forced_from_cumulus_parameterization=cumulus_parameterization_state.output.epsilon_forced,
            scale_dependence_factor_from_cumulus_parameterizaiton=cumulus_parameterization_state.output.scale_dependence_factor,
            pressure_shallow_convective_cloud_top=state.pressure_shallow_convective_cloud_top,
            pressure_mid_convective_cloud_top=state.pressure_mid_convective_cloud_top,
            pressure_deep_convective_cloud_top=state.pressure_deep_convective_cloud_top,
            mass_flux_cloud_base_shallow=state.mass_flux_cloud_base_shallow,
            mass_flux_cloud_base_mid=state.mass_flux_cloud_base_mid,
            mass_flux_cloud_base_deep=state.mass_flux_cloud_base_deep,
            mass_flux_shallow=state.mass_flux_shallow,
            mass_flux_mid=state.mass_flux_mid,
            mass_flux_deep_updraft=state.mass_flux_deep_updraft,
            mass_flux_deep_downdraft=state.mass_flux_deep_downdraft,
            sigma_mid=state.sigma_mid,
            sigma_deep=state.sigma_deep,
        )

        self._update_tendencies(
            dtdt=locals.dtdt,
            dvapordt=locals.dvapordt,
            dudt=locals.dudt,
            dvdt=locals.dvdt,
            dtdt_deep_convection=state.dtdt_deep_convection,
            dvapordt_deep_convection=state.dvapordt_deep_convection,
            dudt_deep_convection=state.dudt_deep_convection,
            dvdt_deep_convection=state.dvdt_deep_convection,
        )

        self._update_convection_codes(
            error_code_from_cumulus_parameterization=cumulus_parameterization_state.output.error_code,
            convection_code_shallow=state.convection_code_shallow,
            convection_code_mid=state.convection_code_mid,
            convection_code_deep=state.convection_code_deep,
        )

        self._update_state_with_tendencies(
            convection_fraction=state.convection_fraction,
            surface_type=state.surface_type,
            u=state.u,
            v=state.v,
            vapor=state.vapor,
            t=state.t,
            p=locals.derived_state.p,
            p_kappa=locals.derived_state.p_kappa,
            mass=locals.derived_state.mass,
            mass_flux_deep_updraft_detrained=state.mass_flux_deep_updraft_detrained,
            mass_flux_deep_updraft_interface=state.mass_flux_deep_updraft_interface,
            total_cumulative_mass_flux_interface=state.total_cumulative_mass_flux_interface,
            total_detraining_mass_flux=state.total_detraining_mass_flux,
            dudt_deep_convection=state.dudt_deep_convection,
            dvdt_deep_convection=state.dvdt_deep_convection,
            dvapordt_deep_convection=state.dvapordt_deep_convection,
            dtdt_deep_convection=state.dtdt_deep_convection,
            dliquiddt_deep_convection=state.dliquiddt_deep_convection,
            dicedt_deep_convection=state.dicedt_deep_convection,
            dcloudfractiondt_deep_convection=state.dcloudfractiondt_deep_convection,
            convective_condensate_source=state.convective_condensate_source,
            evaporation_sublimation_tendency=locals.evaporation_sublimation_tendency,
            convective_precip_flux=locals.convective_precip_flux,
            sublimation_of_convective_precipitation=state.sublimation_of_convective_precipitation,
            evaporation_of_convective_precipitation=state.evaporation_of_convective_precipitation,
            ice_fraction_in_convective_tower=state.ice_fraction_in_convective_tower,
            ice_precip_flux_interface=state.ice_precip_flux_interface,
            liquid_precip_flux_interface=state.liquid_precip_flux_interface,
            convective_liquid=state.convective_liquid,
            convective_ice=state.convective_ice,
            convective_cloud_fraction=state.convective_cloud_fraction,
            convective_rainwater_source=state.convective_rainwater_source,
            convective_precipitation_RAS=state.convective_precipitation_RAS,
            ese=self._ese,
            esx=self._esx,
        )
