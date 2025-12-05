import dataclasses

from ndsl import NDSLRuntime, Quantity, QuantityFactory, State, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.gt4py import BACKWARD, FORWARD, PARALLEL, computation, exp, function, interval, log, max, sqrt
from ndsl.dsl.typing import Bool, BoolFieldIJ, Float, FloatField, FloatFieldIJ
from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.GFDL_1M.driver.config_constants import GFDL1MDriverConfigDependentConstants
from pyMoist.GFDL_1M.driver.constants import constants
from pyMoist.GFDL_1M.driver.sat_tables import GFDL_driver_tables, GlobalTable_driver_qsat
from pyMoist.GFDL_1M.driver.stencils import implicit_fall, wqs2


@function
def revap_racc(
    t1: Float,
    qv1: Float,
    ql1: Float,
    qr1: Float,
    qi1: Float,
    qs1: Float,
    qg1: Float,
    qa1: Float,
    den1: Float,
    denfac: Float,
    rh_limited: Float,
    table1: GlobalTable_driver_qsat,
    table2: GlobalTable_driver_qsat,
    table3: GlobalTable_driver_qsat,
    table4: GlobalTable_driver_qsat,
    des1: GlobalTable_driver_qsat,
    des2: GlobalTable_driver_qsat,
    des3: GlobalTable_driver_qsat,
    des4: GlobalTable_driver_qsat,
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
    """
    Evaporate rain

    reference Fortran: gfdl_cloud_microphys.F90: subroutine revap_racc
    """

    revap = 0.0

    if t1 > constants.T_WFR and qr1 > constants.QPMIN:
        # area and timescale efficiency on revap
        fac_revp = 1.0 - exp(-(0.5 * dts) / tau_revp)

        # -----------------------------------------------------------------------
        # define heat capacity and latent heat coefficient
        # -----------------------------------------------------------------------

        lhl = lv00 + d0_vap * t1
        q_liq = ql1 + qr1
        q_sol = qi1 + qs1 + qg1
        cvm = c_air + qv1 * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
        lcpk = lhl / cvm

        tin = t1 - lcpk * ql1  # presence of clouds suppresses the rain evap
        qpz = qv1 + ql1
        qsat, dqsdt = wqs2(tin, den1, table2, des2)
        dqh = max(ql1, rh_limited * max(qpz, constants.QCMIN))
        dqh = min(dqh, 0.2 * qpz)  # new limiter
        dqv = qsat - qv1  # use this to prevent super - sat the gird box
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
            qden = qr1 * den1
            t2 = tin * tin
            evap = (
                crevp_0
                * t2
                * dq
                * (crevp_1 * sqrt(qden) + crevp_2 * exp(0.725 * log(qden)))
                / (crevp_3 * t2 + crevp_4 * qsat * den1)
            )
            evap = min(qr1, min(0.5 * dts * fac_revp * evap, dqv / (1.0 + lcpk * dqsdt)))
            qr1 = qr1 - evap
            qv1 = qv1 + evap
            q_liq = q_liq - evap
            cvm = c_air + qv1 * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
            t1 = t1 - evap * lhl / cvm
            revap = evap / (0.5 * dts)

        # -----------------------------------------------------------------------
        # accretion: pracc
        # -----------------------------------------------------------------------

        if qr1 > constants.QPMIN and ql1 > constants.QCMIN and qsat < q_minus:
            sink = 0.5 * dts * denfac * cracw * exp(0.95 * log(qr1 * den1))
            sink = sink / (1.0 + sink) * ql1
            ql1 = ql1 - sink
            qr1 = qr1 + sink

    return (
        t1,
        qv1,
        qr1,
        ql1,
        qi1,
        qs1,
        qg1,
        qa1,
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
    """
    Warm rain cloud microphysics: evaporation, accretion

    first half of the timestep
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
        # ensure mask is clear of any previous values
        driver_liquid_precip_flux = 0

    with computation(FORWARD), interval(0, 1):
        # ensure mask is clear of previous value
        precip_fall = False

    # reference Fortran: gfdl_cloud_microphys.F90: subroutine check_column
    # determine if any precip falls in the column
    # if it falls anywhere in the column, the entire column becomes true
    # initalized to 0 (false), potentially changed to 1 (true)
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

        fac_rc = (
            min(1.0, estimated_inversion_strength / 15.0) ** 2
        )  # Estimated inversion strength determine stable regime
        fac_rc = constants.RC * (rthreshs * fac_rc + rthreshu * (1.0 - fac_rc)) ** 3
        # NOTE: the multiplication "constants.RC * (result of parenthetical)" produces different results
        # in Fortran and Python, despite constants.RC and (result of parenthetical) being identical.
        # This creates errors later throughout warm_rain and the larger driver.

    with computation(PARALLEL), interval(...):
        if irain_f != 0:
            # -----------------------------------------------------------------------
            # no subgrid varaibility
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
                        sink = (
                            min(1.0, dq / dl)
                            * dts
                            * c_praut
                            * density
                            * exp(constants.SO3 * log(mixing_ratio_liquid))
                        )
                        sink = min(ql0_max, min(mixing_ratio_liquid, max(0.0, sink)))
                        mixing_ratio_liquid = mixing_ratio_liquid - sink
                        mixing_ratio_rain = mixing_ratio_rain + sink * cloud_fraction_limited

        # Revert In-Cloud condensate
        mixing_ratio_liquid = mixing_ratio_liquid * cloud_fraction_limited
        mixing_ratio_ice = mixing_ratio_ice * cloud_fraction_limited

        # -----------------------------------------------------------------------
        # fall speed of rain
        # -----------------------------------------------------------------------

        if precip_fall == False:
            terminal_speed_rain = constants.VF_MIN
        elif const_vr == True:  # noqa
            terminal_speed_rain = vr_fac  # ifs_2016: 4.0
        else:
            qden = mixing_ratio_rain * density
            if mixing_ratio_rain < constants.THR:
                terminal_speed_rain = constants.VR_MIN
            else:
                terminal_speed_rain = (
                    vr_fac
                    * constants.VCONR
                    * sqrt(min(10.0, constants.SFCRHO / density))
                    * exp(0.2 * log(qden / constants.NORMR))
                )
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
            dmass = dp * (
                1.0
                + mixing_ratio_vapor
                + mixing_ratio_liquid
                + mixing_ratio_rain
                + mixing_ratio_ice
                + mixing_ratio_snow
                + mixing_ratio_graupel
            )


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
    """
    Warm rain cloud microphysics: evaporation, accretion

    first half of the timestep
    """
    from __externals__ import (
        c_air,
        c_vap,
        cracw,
        crevp_0,
        crevp_1,
        crevp_2,
        crevp_3,
        crevp_4,
        d0_vap,
        do_sedi_w,
        dts,
        lv00,
        tau_revp,
    )

    # -----------------------------------------------------------------------
    # vertical velocity transportation during sedimentation
    # -----------------------------------------------------------------------
    with computation(FORWARD), interval(0, 1):
        if do_sedi_w == True:  # noqa
            w = (dmass * w + driver_liquid_precip_flux * terminal_speed_rain) / (
                dmass - driver_liquid_precip_flux
            )

    with computation(FORWARD), interval(1, None):
        if do_sedi_w == True:  # noqa
            w = (
                dmass * w
                - driver_liquid_precip_flux[0, 0, -1] * terminal_speed_rain[0, 0, -1]
                + driver_liquid_precip_flux * terminal_speed_rain
            ) / (dmass + driver_liquid_precip_flux[0, 0, -1] - driver_liquid_precip_flux)

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
    """
    Ensure that information is pushed back to the rest of the model
    """
    with computation(PARALLEL), interval(...):
        evaporation = evaporation + driver_evaporation
        liquid_precip_flux = liquid_precip_flux + driver_liquid_precip_flux
        ice_precip_flux = ice_precip_flux + driver_ice_precip_flux
        mass = mass + driver_liquid_precip_flux + driver_ice_precip_flux

    with computation(FORWARD), interval(0, 1):
        rain = rain + driver_rain


@dataclasses.dataclass
class Locals(State):
    dmass: Quantity = dataclasses.field(
        metadata={
            "name": "dm",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    unused_m1: Quantity = dataclasses.field(
        metadata={
            "name": "unused_m1",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    z_interface: Quantity = dataclasses.field(
        metadata={
            "name": "z_interface",
            "dims": [X_DIM, Y_DIM, Z_INTERFACE_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    precip_fall: Quantity = dataclasses.field(
        metadata={
            "name": "precip_fall",
            "dims": [X_DIM, Y_DIM],
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
        # initalize NDSLRuntime
        super().__init__(stencil_factory)

        # make saturation tables and config available at runtime
        self.config = config
        self.saturation_tables = saturation_tables

        # initalize temporaries
        self._locals = Locals.zeros(quantity_factory)

        # construct stencils
        self._warm_rain_step_1 = stencil_factory.from_dims_halo(
            func=warm_rain_step_1,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
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
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "dts": config_dependent_constants.DTS,
                "use_ppm": config.USE_PPM,
            },
        )

        self._warm_rain_step_2 = stencil_factory.from_dims_halo(
            func=warm_rain_step_2,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
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
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

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
            t (inout): temperature (K)
            dp (in): change in pressure between model levels (Pa)
            dz (in): change in height between model levels (m)
            w (inout): vertical motion (m/s)
            mixing_ratio_vapor (inout): water vapor mixing ratio (kg/kg)
            mixing_ratio_liquid (inout): liquid water mixing ratio (kg/kg)
            mixing_ratio_rain (inout): rain mixing ratio (kg/kg)
            mixing_ratio_ice (inout): ice mixing ratio (kg/kg)
            mixing_ratio_snow (inout): snow mixing ratio (kg/kg)
            mixing_ratio_graupel (inout): graupel mixing ratio (kg/kg)
            cloud_fraction (inout): cloud fraction
            ccn (in): cloud condensation nuclei
            density (in): density of the grid cell (kg m^-3)
            density_factor (in): details unknown
            c_praut (in): details unknown
            terminal_speed_rain (inout): terminal speed of rain (m/s)
            rh_limited (in): relative humidity with limits imposed
            estimated_inversion_strength (in): estimated inversion strength (K)
            one_minus_sigma (out): details unknown
            mass (inout): mass of grid cell
            rain (out): model at-large rain precipitation (kg m^-2 s^-1)
            driver_rain (): in-driver rain precipitation (kg m^-2 s^-1)
            ice_precip_flux (out): model at-large non-anvil large scale ice precip flux (kg m^-2 s^-1)
            driver_ice_precip_flux (out): in-driver non-anvil large scale ice precip flux (kg m^-2 s^-1)
            liquid_precip_flux (out): model at-large non-anvil large scare liquid precip flux (kg m^-2 s^-1)
            driver_liquid_precip_flux (out): in-driver non-anvil large scare liquid precip flux (kg m^-2 s^-1)
            evaporation (out): model at-large non-anvil large scale evaporation (kg kg-1 s-1)
            driver_evaporation (out): in-driver non-anvil large scale evaporation (kg kg-1 s-1)
        """
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
            table1=self.saturation_tables.table1,
            table2=self.saturation_tables.table2,
            table3=self.saturation_tables.table3,
            table4=self.saturation_tables.table4,
            des1=self.saturation_tables.des1,
            des2=self.saturation_tables.des2,
            des3=self.saturation_tables.des3,
            des4=self.saturation_tables.des4,
        )
        if self.config.USE_PPM is False:
            # NOTE: somehow errors pop up in rain1 and m1_rain within implicit fall, despite all of the
            # imputs being correct and implicit_fall verifying at three separate calls
            # within the terminal_fall module. May be a similar issue to the warm_rain_part_1 error
            # (different result despite inputs being identical, possible registry issue??).
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
            table1=self.saturation_tables.table1,
            table2=self.saturation_tables.table2,
            table3=self.saturation_tables.table3,
            table4=self.saturation_tables.table4,
            des1=self.saturation_tables.des1,
            des2=self.saturation_tables.des2,
            des3=self.saturation_tables.des3,
            des4=self.saturation_tables.des4,
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
