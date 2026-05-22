import dataclasses

from ndsl import Local, LocalState, NDSLRuntime, QuantityFactory, StencilFactory
from ndsl.constants import I_DIM, J_DIM, K_DIM, K_INTERFACE_DIM
from ndsl.dsl.gt4py import BACKWARD, FORWARD, PARALLEL, computation, exp, function, interval, log, max, sqrt
from ndsl.dsl.typing import Bool, BoolFieldIJ, Float, FloatField, FloatFieldIJ
from ndsl.stencils import set_IJ_mask_value, set_value, set_value_2D

from pyMoist.microphysics.GFDL_1M.config import GFDL1MConfig
from pyMoist.microphysics.GFDL_1M.driver.config_constants import GFDL1MDriverConfigDependentConstants
from pyMoist.microphysics.GFDL_1M.driver.constants import constants
from pyMoist.microphysics.GFDL_1M.driver.sat_tables import GFDL_driver_tables, GlobalTable_driver_qsat
from pyMoist.microphysics.GFDL_1M.driver.stencils import implicit_fall, wqs2


@function
def revap_racc(
    t: Float,
    mixing_ratio_vapor: Float,
    mixing_ratio_liquid: Float,
    mixing_ratio_rain: Float,
    mixing_ratio_ice: Float,
    mixing_ratio_snow: Float,
    mixing_ratio_graupel: Float,
    cloud_fraction: Float,
    density: Float,
    density_factor: Float,
    rh_limited: Float,
    table1: Float,
    table2: Float,
    table3: Float,
    table4: Float,
    des1: Float,
    des2: Float,
    des3: Float,
    des4: Float,
    c_air: Float,
    c_vap: Float,
    cracw: Float,
    crevp_0: Float,
    crevp_1: Float,
    crevp_2: Float,
    crevp_3: Float,
    crevp_4: Float,
    d0_vap: Float,
    lv00: Float,
    tau_revp: Float,
    dts: Float,
):
    """Evaporate rain

    reference Fortran: gfdl_cloud_microphys.F90: subroutine revap_racc

    Args:
        t (Float)
        mixing_ratio_vapor (Float)
        mixing_ratio_liquid (Float)
        mixing_ratio_rain (Float)
        mixing_ratio_ice (Float)
        mixing_ratio_snow (Float)
        mixing_ratio_graupel (Float)
        cloud_fraction (Float)
        density (Float)
        density_factor (Float)
        rh_limited (Float)
        table1 (Float)
        table2 (Float)
        table3 (Float)
        table4 (Float)
        des1 (Float)
        des2 (Float)
        des3 (Float)
        des4 (Float)
        c_air (Float)
        c_vap (Float)
        cracw (Float)
        crevp_0 (Float)
        crevp_1 (Float)
        crevp_2 (Float)
        crevp_3 (Float)
        crevp_4 (Float)
        d0_vap (Float)
        lv00 (Float)
        tau_revp (Float)
        dts (Float)

    Returns:
        Float: t
        Float: mixing_ratio_vapor
        Float: mixing_ratio_liquid
        Float: mixing_ratio_rain
        Float: mixing_ratio_ice
        Float: mixing_ratio_snow
        Float: mixing_ratio_graupel
        Float: cloud_fraction
        Float: revap
    """
    revap = 0.0

    if t > constants.T_WFR and mixing_ratio_rain > constants.QPMIN:
        # area and timescale efficiency on revap
        fac_revp = 1.0 - exp(-(0.5 * dts) / tau_revp)

        # -----------------------------------------------------------------------
        # define heat capacity and latent heat coefficient
        # -----------------------------------------------------------------------

        lhl = lv00 + d0_vap * t
        q_liq = mixing_ratio_liquid + mixing_ratio_rain
        q_sol = mixing_ratio_ice + mixing_ratio_snow + mixing_ratio_graupel
        cvm = c_air + mixing_ratio_vapor * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
        lcpk = lhl / cvm

        tin = t - lcpk * mixing_ratio_liquid  # presence of clouds suppresses the rain evap
        qpz = mixing_ratio_vapor + mixing_ratio_liquid
        qsat, dqsdt = wqs2(tin, density, table2, des2)
        dqh = max(mixing_ratio_liquid, rh_limited * max(qpz, constants.QCMIN))
        dqh = min(dqh, 0.2 * qpz)  # new limiter
        dqv = qsat - mixing_ratio_vapor  # use this to prevent super - sat the gird box
        q_minus = qpz - dqh
        q_plus = qpz + dqh

        # -----------------------------------------------------------------------
        # qsat must be > q_minus to activate evaporation
        # qsat must be < q_plus to activate accretion
        # -----------------------------------------------------------------------

        # -----------------------------------------------------------------------
        # rain evaporation
        # -----------------------------------------------------------------------

        if dqv > constants.QVMIN and qsat > q_minus:
            if qsat > q_plus:
                dq = qsat - qpz
            else:
                # -----------------------------------------------------------------------
                # q_minus < qsat < q_plus
                # dq == dqh if qsat == q_minus
                # -----------------------------------------------------------------------
                dq = 0.25 * (q_minus - qsat) ** 2 / dqh
            qden = mixing_ratio_rain * density
            t2 = tin * tin
            evap = crevp_0 * t2 * dq * (crevp_1 * sqrt(qden) + crevp_2 * exp(0.725 * log(qden))) / (crevp_3 * t2 + crevp_4 * qsat * density)
            evap = min(mixing_ratio_rain, min(0.5 * dts * fac_revp * evap, dqv / (1.0 + lcpk * dqsdt)))
            mixing_ratio_rain = mixing_ratio_rain - evap
            mixing_ratio_vapor = mixing_ratio_vapor + evap
            q_liq = q_liq - evap
            cvm = c_air + mixing_ratio_vapor * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
            t = t - evap * lhl / cvm
            revap = evap / (0.5 * dts)

        # -----------------------------------------------------------------------
        # accretion: pracc
        # -----------------------------------------------------------------------

        if mixing_ratio_rain > constants.QPMIN and mixing_ratio_liquid > constants.QCMIN and qsat < q_minus:
            sink = 0.5 * dts * density_factor * cracw * exp(0.95 * log(mixing_ratio_rain * density))
            sink = sink / (1.0 + sink) * mixing_ratio_liquid
            mixing_ratio_liquid = mixing_ratio_liquid - sink
            mixing_ratio_rain = mixing_ratio_rain + sink

    return (
        t,
        mixing_ratio_vapor,
        mixing_ratio_rain,
        mixing_ratio_liquid,
        mixing_ratio_ice,
        mixing_ratio_snow,
        mixing_ratio_graupel,
        cloud_fraction,
        revap,
    )


def warm_rain_step_1(
    dp: FloatField,
    dz: FloatField,
    t: FloatField,
    mixing_ratio_vapor: FloatField,
    mixing_ratio_liquid: FloatField,
    mixing_ratio_rain: FloatField,
    mixing_ratio_ice: FloatField,
    mixing_ratio_snow: FloatField,
    mixing_ratio_graupel: FloatField,
    cloud_fraction: FloatField,
    ccn: FloatField,
    density: FloatField,
    density_factor: FloatField,
    c_praut: FloatField,
    terminal_speed_rain: FloatField,
    driver_evaporation: FloatField,
    driver_liquid_precip_flux: FloatField,
    rh_limited: FloatField,
    estimated_inversion_strength: FloatFieldIJ,
    one_minus_sigma: FloatFieldIJ,
    z_interface: FloatField,
    dmass: FloatField,
    precip_fall: BoolFieldIJ,
    table1: GlobalTable_driver_qsat,
    table2: GlobalTable_driver_qsat,
    table3: GlobalTable_driver_qsat,
    table4: GlobalTable_driver_qsat,
    des1: GlobalTable_driver_qsat,
    des2: GlobalTable_driver_qsat,
    des3: GlobalTable_driver_qsat,
    des4: GlobalTable_driver_qsat,
):
    """Warm rain cloud microphysics: evaporation, accretion

    first half of the timestep

    Args:
        dp (FloatField)
        dz (FloatField)
        t (FloatField)
        mixing_ratio_vapor (FloatField)
        mixing_ratio_liquid (FloatField)
        mixing_ratio_rain (FloatField)
        mixing_ratio_ice (FloatField)
        mixing_ratio_snow (FloatField)
        mixing_ratio_graupel (FloatField)
        cloud_fraction (FloatField)
        ccn (FloatField)
        density (FloatField)
        density_factor (FloatField)
        c_praut (FloatField)
        terminal_speed_rain (FloatField)
        driver_evaporation (FloatField)
        driver_liquid_precip_flux (FloatField)
        rh_limited (FloatField)
        estimated_inversion_strength (FloatFieldIJ)
        one_minus_sigma (FloatFieldIJ)
        z_interface (FloatField)
        dmass (FloatField)
        precip_fall (BoolFieldIJ)
        table1 (GlobalTable_driver_qsat)
        table2 (GlobalTable_driver_qsat)
        table3 (GlobalTable_driver_qsat)
        table4 (GlobalTable_driver_qsat)
        des1 (GlobalTable_driver_qsat)
        des2 (GlobalTable_driver_qsat)
        des3 (GlobalTable_driver_qsat)
        des4 (GlobalTable_driver_qsat)
    """
    from __externals__ import (
        c_air,
        c_vap,
        const_vr,
        cracw,
        crevp_0,
        crevp_1,
        crevp_2,
        crevp_3,
        crevp_4,
        d0_vap,
        do_qa,
        do_sedi_w,
        dts,
        irain_f,
        lv00,
        ql0_max,
        rthreshs,
        rthreshu,
        tau_revp,
        vr_fac,
        vr_max,
        z_slope_liq,
    )

    with computation(PARALLEL), interval(...):
        # reset to zero, clearing data from previous call
        driver_liquid_precip_flux = 0.0

    with computation(FORWARD), interval(0, 1):
        # ensure mask is clear of previous value
        precip_fall = False

    # reference Fortran: gfdl_cloud_microphys.F90: subroutine check_column
    # determine if any precip falls in the column
    # if it falls anywhere in the column, the entire column becomes true
    # initialized to 0 (false), potentially changed to 1 (true)
    with computation(FORWARD), interval(...):
        if mixing_ratio_rain > constants.QPMIN:
            precip_fall = True
    # end reference Fortran: gfdl_cloud_microphys.F90: subroutine check_column

    # -----------------------------------------------------------------------
    # auto - conversion
    # assuming linear subgrid vertical distribution of cloud water
    # following lin et al. 1994, mwr
    # -----------------------------------------------------------------------

    with computation(PARALLEL), interval(...):
        # Use In-Cloud condensates
        if do_qa == False:  # noqa
            cloud_fraction_limited = max(cloud_fraction, constants.QCMIN)
        else:
            cloud_fraction_limited = 1.0
        mixing_ratio_liquid = mixing_ratio_liquid / cloud_fraction_limited
        mixing_ratio_ice = mixing_ratio_ice / cloud_fraction_limited

        fac_rc = min(1.0, estimated_inversion_strength / 15.0) ** 2  # Estimated inversion strength determine stable regime
        fac_rc = constants.RC * (rthreshs * fac_rc + rthreshu * (1.0 - fac_rc)) ** 3
        # NOTE: the multiplication "constants.RC * (result of parenthetical)" produces different results
        # in Fortran and Python, despite constants.RC and (result of parenthetical) being identical.
        # This creates errors later throughout warm_rain and the larger driver.

    with computation(PARALLEL), interval(...):
        if irain_f != 0:
            # -----------------------------------------------------------------------
            # no subgrid variability
            # -----------------------------------------------------------------------
            if cloud_fraction_limited > one_minus_sigma:
                if t > constants.T_WFR:
                    qc = fac_rc * ccn / density
                    dq = mixing_ratio_liquid - qc
                    if dq > 0.0:
                        sink = min(
                            dq,
                            dts * c_praut * density * exp(constants.SO3 * log(mixing_ratio_liquid)),
                        )
                        sink = min(ql0_max, min(mixing_ratio_liquid, max(0.0, sink)))
                        mixing_ratio_liquid = mixing_ratio_liquid - sink
                        mixing_ratio_rain = mixing_ratio_rain + sink * cloud_fraction_limited

    with computation(FORWARD), interval(1, None):
        if irain_f == 0:
            # -----------------------------------------------------------------------
            # with subgrid variability
            # -----------------------------------------------------------------------

            # begin reference Fortran: gfdl_cloud_microphys.F90: subroutine linear_prof
            # definition of vertical subgrid variability
            # used for cloud ice and cloud water autoconversion
            # qi -- > ql & ql -- > qr
            # edges: qe == qbar + / - dm
            if z_slope_liq == True:  # noqa
                dql = 0.5 * (mixing_ratio_liquid - mixing_ratio_liquid[0, 0, -1])

    # -----------------------------------------------------------------------
    # use twice the strength of the positive definiteness limiter (lin et al 1994)
    # -----------------------------------------------------------------------

    with computation(FORWARD), interval(1, -1):
        if irain_f == 0:
            if z_slope_liq == True:  # noqa
                dl = 0.5 * min(abs(dql + dql[0, 0, 1]), 0.5 * mixing_ratio_liquid)
                if dql * dql[0, 0, 1] <= 0.0:
                    if dql > 0.0:  # local max
                        dl = min(dl, min(dql, -dql[0, 0, 1]))
                    else:
                        dl = 0.0

    with computation(FORWARD), interval(0, 1):
        if irain_f == 0:
            if z_slope_liq == True:  # noqa
                dl = 0

    with computation(FORWARD), interval(-1, None):
        if irain_f == 0:
            if z_slope_liq == True:  # noqa
                dl = 0

    # -----------------------------------------------------------------------
    # impose a presumed background horizontal variability
    # that is proportional to the value itself
    # -----------------------------------------------------------------------
    with computation(PARALLEL), interval(...):
        if irain_f == 0:
            if z_slope_liq == True:  # noqa
                dl = max(dl, max(constants.QVMIN, rh_limited * mixing_ratio_liquid))
            if z_slope_liq == False:  # noqa
                dl = max(constants.QVMIN, rh_limited * mixing_ratio_liquid)

            # end reference Fortran: gfdl_cloud_microphys.F90: subroutine linear_prof

    with computation(PARALLEL), interval(...):
        if irain_f == 0:
            if cloud_fraction_limited > one_minus_sigma:
                if t > constants.T_WFR + constants.DT_FR:
                    dl = min(max(constants.QCMIN, dl), 0.5 * mixing_ratio_liquid)
                    # --------------------------------------------------------------------
                    # as in klein's gfdl am2 stratiform scheme (with subgrid variations)
                    # --------------------------------------------------------------------
                    qc = fac_rc * ccn / density
                    dq = 0.5 * (mixing_ratio_liquid + dl - qc)
                    # --------------------------------------------------------------------
                    # dq = dl if qc == q_minus = ql - dl
                    # dq = 0 if qc == q_plus = ql + dl
                    # --------------------------------------------------------------------
                    if dq > 0.0:  # q_plus > qc
                        # --------------------------------------------------------------------
                        # revised continuous form: linearly decays
                        # (with subgrid dl) to zero at qc == ql + dl
                        # --------------------------------------------------------------------
                        sink = min(1.0, dq / dl) * dts * c_praut * density * exp(constants.SO3 * log(mixing_ratio_liquid))
                        sink = min(ql0_max, min(mixing_ratio_liquid, max(0.0, sink)))
                        mixing_ratio_liquid = mixing_ratio_liquid - sink
                        mixing_ratio_rain = mixing_ratio_rain + sink * cloud_fraction_limited

        # Revert In-Cloud condensate
        mixing_ratio_liquid = mixing_ratio_liquid * cloud_fraction_limited
        mixing_ratio_ice = mixing_ratio_ice * cloud_fraction_limited

        # -----------------------------------------------------------------------
        # fall speed of rain
        # -----------------------------------------------------------------------

        if precip_fall == False:  # noqa
            terminal_speed_rain = constants.VF_MIN
        elif const_vr == True:  # noqa
            terminal_speed_rain = vr_fac  # ifs_2016: 4.0
        else:
            qden = mixing_ratio_rain * density
            if mixing_ratio_rain < constants.THR:
                terminal_speed_rain = constants.VR_MIN
            else:
                terminal_speed_rain = vr_fac * constants.VCONR * sqrt(min(10.0, constants.SFCRHO / density)) * exp(0.2 * log(qden / constants.NORMR))
                terminal_speed_rain = min(vr_max, max(constants.VR_MIN, terminal_speed_rain))

    with computation(FORWARD), interval(-1, None):
        z_interface[0, 0, 1] = constants.ZS

    with computation(BACKWARD), interval(...):
        z_interface = z_interface[0, 0, 1] - dz  # dz < 0

    # -----------------------------------------------------------------------
    # evaporation and accretion of rain for the first 1 / 2 time step
    # -----------------------------------------------------------------------
    with computation(PARALLEL), interval(...):
        (
            t,
            mixing_ratio_vapor,
            mixing_ratio_rain,
            mixing_ratio_liquid,
            mixing_ratio_ice,
            mixing_ratio_snow,
            mixing_ratio_graupel,
            cloud_fraction,
            revap,
        ) = revap_racc(
            t,
            mixing_ratio_vapor,
            mixing_ratio_liquid,
            mixing_ratio_rain,
            mixing_ratio_ice,
            mixing_ratio_snow,
            mixing_ratio_graupel,
            cloud_fraction,
            density,
            density_factor,
            rh_limited,
            table1,
            table2,
            table3,
            table4,
            des1,
            des2,
            des3,
            des4,
            c_air,
            c_vap,
            cracw,
            crevp_0,
            crevp_1,
            crevp_2,
            crevp_3,
            crevp_4,
            d0_vap,
            lv00,
            tau_revp,
            dts,
        )

        driver_evaporation = revap

    with computation(PARALLEL), interval(...):
        if do_sedi_w == True:  # noqa
            dmass = dp * (1.0 + mixing_ratio_vapor + mixing_ratio_liquid + mixing_ratio_rain + mixing_ratio_ice + mixing_ratio_snow + mixing_ratio_graupel)


def warm_rain_step_2(
    t: FloatField,
    mixing_ratio_vapor: FloatField,
    mixing_ratio_liquid: FloatField,
    mixing_ratio_rain: FloatField,
    mixing_ratio_ice: FloatField,
    mixing_ratio_snow: FloatField,
    mixing_ratio_graupel: FloatField,
    cloud_fraction: FloatField,
    density: FloatField,
    density_factor: FloatField,
    terminal_speed_rain: FloatField,
    driver_evaporation: FloatField,
    driver_liquid_precip_flux: FloatField,
    w: FloatField,
    rh_limited: FloatField,
    dmass: FloatField,
    table1: GlobalTable_driver_qsat,
    table2: GlobalTable_driver_qsat,
    table3: GlobalTable_driver_qsat,
    table4: GlobalTable_driver_qsat,
    des1: GlobalTable_driver_qsat,
    des2: GlobalTable_driver_qsat,
    des3: GlobalTable_driver_qsat,
    des4: GlobalTable_driver_qsat,
):
    """Warm rain cloud microphysics: evaporation, accretion

    second half of the timestep

    Args:
        t (FloatField)
        mixing_ratio_vapor (FloatField)
        mixing_ratio_liquid (FloatField)
        mixing_ratio_rain (FloatField)
        mixing_ratio_ice (FloatField)
        mixing_ratio_snow (FloatField)
        mixing_ratio_graupel (FloatField)
        cloud_fraction (FloatField)
        density (FloatField)
        density_factor (FloatField)
        terminal_speed_rain (FloatField)
        driver_evaporation (FloatField)
        driver_liquid_precip_flux (FloatField)
        w (FloatField)
        rh_limited (FloatField)
        dmass (FloatField)
        table1 (GlobalTable_driver_qsat)
        table2 (GlobalTable_driver_qsat)
        table3 (GlobalTable_driver_qsat)
        table4 (GlobalTable_driver_qsat)
        des1 (GlobalTable_driver_qsat)
        des2 (GlobalTable_driver_qsat)
        des3 (GlobalTable_driver_qsat)
        des4 (GlobalTable_driver_qsat)
    """
    from __externals__ import c_air, c_vap, cracw, crevp_0, crevp_1, crevp_2, crevp_3, crevp_4, d0_vap, do_sedi_w, dts, lv00, tau_revp

    # -----------------------------------------------------------------------
    # vertical velocity transportation during sedimentation
    # -----------------------------------------------------------------------
    with computation(FORWARD), interval(0, 1):
        if do_sedi_w == True:  # noqa
            w = (dmass * w + driver_liquid_precip_flux * terminal_speed_rain) / (dmass - driver_liquid_precip_flux)

    with computation(FORWARD), interval(1, None):
        if do_sedi_w == True:  # noqa
            w = (dmass * w - driver_liquid_precip_flux[0, 0, -1] * terminal_speed_rain[0, 0, -1] + driver_liquid_precip_flux * terminal_speed_rain) / (
                dmass + driver_liquid_precip_flux[0, 0, -1] - driver_liquid_precip_flux
            )

    # -----------------------------------------------------------------------
    # evaporation and accretion of rain for the remaing 1 / 2 time step
    # -----------------------------------------------------------------------
    with computation(PARALLEL), interval(...):
        (
            t,
            mixing_ratio_vapor,
            mixing_ratio_rain,
            mixing_ratio_liquid,
            mixing_ratio_ice,
            mixing_ratio_snow,
            mixing_ratio_graupel,
            cloud_fraction,
            revap,
        ) = revap_racc(
            t,
            mixing_ratio_vapor,
            mixing_ratio_liquid,
            mixing_ratio_rain,
            mixing_ratio_ice,
            mixing_ratio_snow,
            mixing_ratio_graupel,
            cloud_fraction,
            density,
            density_factor,
            rh_limited,
            table1,
            table2,
            table3,
            table4,
            des1,
            des2,
            des3,
            des4,
            c_air,
            c_vap,
            cracw,
            crevp_0,
            crevp_1,
            crevp_2,
            crevp_3,
            crevp_4,
            d0_vap,
            lv00,
            tau_revp,
            dts,
        )

        driver_evaporation = driver_evaporation + revap


def update_outputs(
    driver_liquid_precip_flux: FloatField,
    driver_ice_precip_flux: FloatField,
    driver_rain: FloatFieldIJ,
    driver_evaporation: FloatField,
    evaporation: FloatField,
    liquid_precip_flux: FloatField,
    ice_precip_flux: FloatField,
    mass: FloatField,
    rain: FloatFieldIJ,
):
    """Ensure that information is pushed back to the rest of the model

    Args:
        driver_liquid_precip_flux (FloatField)
        driver_ice_precip_flux (FloatField)
        driver_rain (FloatFieldIJ)
        driver_evaporation (FloatField)
        evaporation (FloatField)
        liquid_precip_flux (FloatField)
        ice_precip_flux (FloatField)
        mass (FloatField)
        rain (FloatFieldIJ)
    """
    with computation(PARALLEL), interval(...):
        evaporation = evaporation + driver_evaporation
        mass = mass + driver_liquid_precip_flux + driver_ice_precip_flux

    with computation(FORWARD), interval(...):
        liquid_precip_flux[0, 0, 1] = liquid_precip_flux[0, 0, 1] + driver_liquid_precip_flux
        ice_precip_flux[0, 0, 1] = ice_precip_flux[0, 0, 1] + driver_ice_precip_flux

    with computation(FORWARD), interval(0, 1):
        rain = rain + driver_rain


@dataclasses.dataclass
class Locals(LocalState):
    dmass: Local = dataclasses.field(
        metadata={
            "name": "dm",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    unused_m1: Local = dataclasses.field(
        metadata={
            "name": "unused_m1",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    z_interface: Local = dataclasses.field(
        metadata={
            "name": "z_interface",
            "dims": [I_DIM, J_DIM, K_INTERFACE_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    precip_fall: Local = dataclasses.field(
        metadata={
            "name": "precip_fall",
            "dims": [I_DIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Bool,
        }
    )


class GFDL1MWarmRain(NDSLRuntime):
    """
    Warm rain cloud microphysics: evaporation, accretion
    """

    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GFDL1MConfig,
        config_dependent_constants: GFDL1MDriverConfigDependentConstants,
        saturation_tables: GFDL_driver_tables,
    ):
        """Initialize the warm rain module

        Args:
            stencil_factory (StencilFactory)
            quantity_factory (QuantityFactory)
            config (GFDL1MConfig)
            config_dependent_constants (GFDL1MDriverConfigDependentConstants)
            saturation_tables (GFDL_driver_tables)
        """
        # initialize NDSLRuntime
        super().__init__(stencil_factory)

        # make saturation tables and config available at runtime
        self.config = config
        self.saturation_tables = saturation_tables

        # initialize locals
        self._locals = Locals.make_locals(quantity_factory)

        # construct stencils
        self._warm_rain_step_1 = stencil_factory.from_dims_halo(
            func=warm_rain_step_1,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "dts": config_dependent_constants.DTS,
                "do_qa": config.DO_QA,
                "rthreshs": config.RTHRESHS,
                "rthreshu": config.RTHRESHU,
                "irain_f": config.IRAIN_F,
                "ql0_max": config.QL0_MAX,
                "z_slope_liq": config.Z_SLOPE_LIQ,
                "vr_fac": config.VR_FAC,
                "const_vr": config.CONST_VR,
                "vr_max": config.VR_MAX,
                "tau_revp": config.TAU_REVP,
                "lv00": config_dependent_constants.LV00,
                "d0_vap": config_dependent_constants.D0_VAP,
                "c_air": config_dependent_constants.C_AIR,
                "c_vap": config_dependent_constants.C_VAP,
                "crevp_0": config_dependent_constants.CREVP_0,
                "crevp_1": config_dependent_constants.CREVP_1,
                "crevp_2": config_dependent_constants.CREVP_2,
                "crevp_3": config_dependent_constants.CREVP_3,
                "crevp_4": config_dependent_constants.CREVP_4,
                "cracw": config_dependent_constants.CRACW,
                "do_sedi_w": config.DO_SEDI_W,
                "use_ppm": config.USE_PPM,
            },
        )

        self._implicit_fall = stencil_factory.from_dims_halo(
            func=implicit_fall,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "dts": config_dependent_constants.DTS,
                "use_ppm": config.USE_PPM,
            },
        )

        self._warm_rain_step_2 = stencil_factory.from_dims_halo(
            func=warm_rain_step_2,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "dts": config_dependent_constants.DTS,
                "do_qa": config.DO_QA,
                "rthreshs": config.RTHRESHS,
                "rthreshu": config.RTHRESHU,
                "irain_f": config.IRAIN_F,
                "ql0_max": config.QL0_MAX,
                "z_slope_liq": config.Z_SLOPE_LIQ,
                "vr_fac": config.VR_FAC,
                "const_vr": config.CONST_VR,
                "vr_max": config.VR_MAX,
                "tau_revp": config.TAU_REVP,
                "lv00": config_dependent_constants.LV00,
                "d0_vap": config_dependent_constants.D0_VAP,
                "c_air": config_dependent_constants.C_AIR,
                "c_vap": config_dependent_constants.C_VAP,
                "crevp_0": config_dependent_constants.CREVP_0,
                "crevp_1": config_dependent_constants.CREVP_1,
                "crevp_2": config_dependent_constants.CREVP_2,
                "crevp_3": config_dependent_constants.CREVP_3,
                "crevp_4": config_dependent_constants.CREVP_4,
                "cracw": config_dependent_constants.CRACW,
                "do_sedi_w": config.DO_SEDI_W,
                "use_ppm": config.USE_PPM,
            },
        )

        self._update_outputs = stencil_factory.from_dims_halo(
            func=update_outputs,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )
        self._set_value_IJ = stencil_factory.from_dims_halo(
            func=set_value_2D,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )
        self._set_value = stencil_factory.from_dims_halo(
            func=set_value,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )
        self._set_value_K_interface = stencil_factory.from_dims_halo(
            func=set_value,
            compute_dims=[I_DIM, J_DIM, K_INTERFACE_DIM],
        )
        self._set_IJ_mask = stencil_factory.from_dims_halo(
            func=set_IJ_mask_value,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )

        # Dev NOTE: this is an orchestration workaround. Direct call to
        #           `self.saturation_tables.X` fails closure capture for
        #           argument reconstruction at call time
        self._table1 = (self.saturation_tables.table1,)
        self._table2 = (self.saturation_tables.table2,)
        self._table3 = (self.saturation_tables.table3,)
        self._table4 = (self.saturation_tables.table4,)
        self._des1 = (self.saturation_tables.des1,)
        self._des2 = (self.saturation_tables.des2,)
        self._des3 = (self.saturation_tables.des3,)
        self._des4 = (self.saturation_tables.des4,)

    def __call__(
        self,
        t: FloatField,
        dp: FloatField,
        dz: FloatField,
        w: FloatField,
        mixing_ratio_vapor: FloatField,
        mixing_ratio_liquid: FloatField,
        mixing_ratio_rain: FloatField,
        mixing_ratio_ice: FloatField,
        mixing_ratio_snow: FloatField,
        mixing_ratio_graupel: FloatField,
        cloud_fraction: FloatField,
        ccn: FloatField,
        density: FloatField,
        density_factor: FloatField,
        c_praut: FloatField,
        terminal_speed_rain: FloatField,
        rh_limited: FloatField,
        estimated_inversion_strength: FloatFieldIJ,
        one_minus_sigma: FloatFieldIJ,
        mass: FloatField,
        rain: FloatFieldIJ,
        driver_rain: FloatFieldIJ,
        ice_precip_flux: FloatField,
        driver_ice_precip_flux: FloatField,
        liquid_precip_flux: FloatField,
        driver_liquid_precip_flux: FloatField,
        evaporation: FloatField,
        driver_evaporation: FloatField,
    ):
        """
        Warm rain cloud microphysics: evaporation, accretion

        Args:
            t (FloatField): (inout) temperature (K)
            dp (FloatField): (in) change in pressure between model levels (Pa)
            dz (FloatField): (in) change in height between model levels (m)
            w (FloatField): (inout) vertical motion (m/s)
            mixing_ratio_vapor (FloatField): (inout) water vapor mixing ratio (kg/kg)
            mixing_ratio_liquid (FloatField): (inout) liquid water mixing ratio (kg/kg)
            mixing_ratio_rain (FloatField): (inout) rain mixing ratio (kg/kg)
            mixing_ratio_ice (FloatField): (inout) ice mixing ratio (kg/kg)
            mixing_ratio_snow (FloatField): (inout) snow mixing ratio (kg/kg)
            mixing_ratio_graupel (FloatField): (inout) graupel mixing ratio (kg/kg)
            cloud_fraction (FloatField): (inout) cloud fraction
            ccn (FloatField): (in) cloud condensation nuclei
            density (FloatField): (in) density of the grid cell (kg m^-3)
            density_factor (FloatField): (in) details unknown
            c_praut (FloatField): (in) details unknown
            terminal_speed_rain (FloatField): (inout) terminal speed of rain (m/s)
            rh_limited (FloatField): (in) relative humidity with limits imposed
            estimated_inversion_strength (FloatFieldIJ): (in) estimated inversion strength (K)
            one_minus_sigma (FloatFieldIJ): (out) details unknown
            mass (FloatField): (inout) mass of grid cell
            rain (FloatFieldIJ): (out) model at-large rain precipitation (kg m^-2 s^-1)
            driver_rain (FloatFieldIJ): (out) in-driver rain precipitation (kg m^-2 s^-1)
            ice_precip_flux (FloatField): (out) model at-large non-anvil large scale ice precip flux
            (kg m^-2 s^-1)
            driver_ice_precip_flux (FloatField): (out) in-driver non-anvil large scale ice precip flux
            (kg m^-2 s^-1)
            liquid_precip_flux (FloatField): (out) model at-large non-anvil large scare liquid precip flux
            (kg m^-2 s^-1)
            driver_liquid_precip_flux (FloatField): (out) in-driver non-anvil large scare liquid precip flux
            (kg m^-2 s^-1)
            evaporation (FloatField): (out) model at-large non-anvil large scale evaporation (kg kg-1 s-1)
            driver_evaporation (FloatField): (out) in-driver non-anvil large scale evaporation (kg kg-1 s-1)
        """
        # reset locals
        self._set_value(self._locals.dmass, Float(0))
        self._set_value(self._locals.unused_m1, Float(0))
        self._set_value_K_interface(self._locals.z_interface, Float(0))
        self._set_IJ_mask(self._locals.precip_fall, False)

        self._warm_rain_step_1(
            dp=dp,
            dz=dz,
            t=t,
            mixing_ratio_vapor=mixing_ratio_vapor,
            mixing_ratio_liquid=mixing_ratio_liquid,
            mixing_ratio_rain=mixing_ratio_rain,
            mixing_ratio_ice=mixing_ratio_ice,
            mixing_ratio_snow=mixing_ratio_snow,
            mixing_ratio_graupel=mixing_ratio_graupel,
            cloud_fraction=cloud_fraction,
            ccn=ccn,
            density=density,
            density_factor=density_factor,
            c_praut=c_praut,
            terminal_speed_rain=terminal_speed_rain,
            driver_evaporation=driver_evaporation,
            driver_liquid_precip_flux=driver_liquid_precip_flux,
            rh_limited=rh_limited,
            estimated_inversion_strength=estimated_inversion_strength,
            one_minus_sigma=one_minus_sigma,
            z_interface=self._locals.z_interface,
            dmass=self._locals.dmass,
            precip_fall=self._locals.precip_fall,
            table1=self._table1,
            table2=self._table2,
            table3=self._table3,
            table4=self._table4,
            des1=self._des1,
            des2=self._des2,
            des3=self._des3,
            des4=self._des4,
        )

        if not self.config.USE_PPM:
            self._implicit_fall(
                mixing_ratio=mixing_ratio_rain,
                terminal_speed=terminal_speed_rain,
                z_interface=self._locals.z_interface,
                dp=dp,
                mass=self._locals.unused_m1,
                precip_flux=driver_liquid_precip_flux,
                precip=driver_rain,
                precip_fall=self._locals.precip_fall,
            )

        self._warm_rain_step_2(
            t=t,
            mixing_ratio_vapor=mixing_ratio_vapor,
            mixing_ratio_liquid=mixing_ratio_liquid,
            mixing_ratio_rain=mixing_ratio_rain,
            mixing_ratio_ice=mixing_ratio_ice,
            mixing_ratio_snow=mixing_ratio_snow,
            mixing_ratio_graupel=mixing_ratio_graupel,
            cloud_fraction=cloud_fraction,
            density=density,
            density_factor=density_factor,
            terminal_speed_rain=terminal_speed_rain,
            driver_evaporation=driver_evaporation,
            driver_liquid_precip_flux=driver_liquid_precip_flux,
            w=w,
            rh_limited=rh_limited,
            dmass=self._locals.dmass,
            table1=self._table1,
            table2=self._table2,
            table3=self._table3,
            table4=self._table4,
            des1=self._des1,
            des2=self._des2,
            des3=self._des3,
            des4=self._des4,
        )

        self._update_outputs(
            driver_liquid_precip_flux=driver_liquid_precip_flux,
            driver_ice_precip_flux=driver_ice_precip_flux,
            driver_rain=driver_rain,
            driver_evaporation=driver_evaporation,
            evaporation=evaporation,
            liquid_precip_flux=liquid_precip_flux,
            ice_precip_flux=ice_precip_flux,
            mass=mass,
            rain=rain,
        )
