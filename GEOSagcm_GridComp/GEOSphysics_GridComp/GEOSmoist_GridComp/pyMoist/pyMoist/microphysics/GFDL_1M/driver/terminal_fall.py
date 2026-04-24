import dataclasses

from ndsl import Local, LocalState, NDSLRuntime, Quantity, QuantityFactory, StencilFactory
from ndsl.constants import I_DIM, J_DIM, K_DIM, K_INTERFACE_DIM
from ndsl.dsl.gt4py import BACKWARD, FORWARD, PARALLEL, computation, exp, function, interval
from ndsl.dsl.typing import Bool, BoolFieldIJ, Float, FloatField, FloatFieldIJ
from ndsl.stencils import set_IJ_mask_value, set_value, set_value_2D

from pyMoist.microphysics.GFDL_1M.config import GFDL1MConfig
from pyMoist.microphysics.GFDL_1M.driver.config_constants import GFDL1MDriverConfigDependentConstants
from pyMoist.microphysics.GFDL_1M.driver.constants import constants
from pyMoist.microphysics.GFDL_1M.driver.stencils import implicit_fall


def check_precip_get_zt(
    mixing_ratio: FloatField,
    terminal_speed: FloatField,
    z_interface: FloatField,
    z_interface_modified: FloatField,
    internal_precip: FloatFieldIJ,
    precip_fall: BoolFieldIJ,
):
    """Check if a critical mixing ratio is met for precipitation to occur and compute the updated
    grid edge height.

    Args:
        mixing_ratio (FloatField)
        terminal_speed (FloatField)
        z_interface (FloatField)
        z_interface_modified (FloatField)
        internal_precip (FloatFieldIJ)
        precip_fall (BoolFieldIJ)
    """
    from __externals__ import dts

    # melting of falling snow into rain
    with computation(FORWARD), interval(...):
        # determine if any precip falls in the column
        # if it falls anywhere in the column, the entire column becomes true
        # precip_fall initialized to 0 (false), potentially changed to 1 (true)
        if mixing_ratio > constants.QPMIN:
            precip_fall = True

    with computation(FORWARD), interval(0, 1):
        if precip_fall == False:  # noqa
            internal_precip = 0

    with computation(FORWARD), interval(1, None):
        if precip_fall == True:  # noqa
            z_interface_modified = z_interface - dts * (terminal_speed[0, 0, -1] + terminal_speed) / 2.0

    with computation(FORWARD), interval(-1, None):
        if precip_fall == True:  # noqa
            z_interface_modified[0, 0, 1] = constants.ZS - dts * terminal_speed

    with computation(FORWARD), interval(...):
        if precip_fall == True:  # noqa
            if z_interface_modified[0, 0, 1] >= z_interface_modified:
                z_interface_modified[0, 0, 1] = z_interface_modified - constants.DZ_MIN


def update_dmass(
    dmass: FloatField,
    dp: FloatField,
    mixing_ratio_vapor: FloatField,
    mixing_ratio_liquid: FloatField,
    mixing_ratio_rain: FloatField,
    mixing_ratio_ice: FloatField,
    mixing_ratio_snow: FloatField,
    mixing_ratio_graupel: FloatField,
    precip_fall: BoolFieldIJ,
):
    """Compute the change in mass based on mixing ratios

    Args:
        dmass (FloatField)
        dp (FloatField)
        mixing_ratio_vapor (FloatField)
        mixing_ratio_liquid (FloatField)
        mixing_ratio_rain (FloatField)
        mixing_ratio_ice (FloatField)
        mixing_ratio_snow (FloatField)
        mixing_ratio_graupel (FloatField)
        precip_fall (BoolFieldIJ)
    """
    from __externals__ import do_sedi_w

    with computation(PARALLEL), interval(...):
        if precip_fall == True:  # noqa
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


def update_w(
    w: FloatField,
    dmass: FloatField,
    mass: FloatField,
    terminal_speed: FloatField,
    precip_fall: BoolFieldIJ,
):
    """Update vertical motion

    Args:
        w (FloatField)
        dmass (FloatField)
        mass (FloatField)
        terminal_speed (FloatField)
        precip_fall (BoolFieldIJ)
    """
    from __externals__ import do_sedi_w

    with computation(FORWARD), interval(0, 1):
        if precip_fall == True:  # noqa
            if do_sedi_w:
                w = (dmass * w + mass * terminal_speed) / (dmass - mass)

    with computation(FORWARD), interval(1, None):
        if precip_fall == True:  # noqa
            if do_sedi_w:
                w = (dmass * w - mass[0, 0, -1] * terminal_speed[0, 0, -1] + mass * terminal_speed) / (
                    dmass + mass[0, 0, -1] - mass
                )


def reset(
    mass: FloatField,
    precip_fall: BoolFieldIJ,
):
    """Reset masks

    Args:
        mass (FloatField)
        precip_fall (BoolFieldIJ)
    """
    # reset masks and temporaries
    with computation(FORWARD), interval(0, 1):
        precip_fall = False

    with computation(PARALLEL), interval(...):
        # must be reset to zero everywhere in case there is no precipitate
        # in a column and no calculations are performed
        mass = 0.0


@function
def prefall_melting(
    t: Float,
    vapor: Float,
    liquid: Float,
    rain: Float,
    graupel: Float,
    snow: Float,
    ice: Float,
    ice_precip_flux: Float,
    is_frozen: Float,
    c_air: Float,
    c_vap: Float,
    d0_vap: Float,
    dts: Float,
    lv00: Float,
    ql_mlt: Float,
    tau_imlt: Float,
):
    """Melt excess cloud ice before any precipitation occurs.

    reference Fortran: gfdl_cloud_microphys.F90: subroutine terminal_fall

    Args:
        t (Float): _description_
        vapor (Float): _description_
        liquid (Float): _description_
        rain (Float): _description_
        graupel (Float): _description_
        snow (Float): _description_
        ice (Float): _description_
        ice_precip_flux (Float): _description_
        is_frozen (bool): _description_
        c_air (Float): _description_
        c_vap (Float): _description_
        d0_vap (Float): _description_
        dts (Float): _description_
        lv00 (Float): _description_
        ql_mlt (Float): _description_
        tau_imlt (Float): _description_

    Returns:
        Float: t
        Float: liquid
        Float: rain
        Float: ice
        Float: cvm
        Float: lhi
        Float: icpk
        Float: ice_precip_flux
    """

    fac_imlt = 1.0 - exp(-dts / tau_imlt)

    # -----------------------------------------------------------------------
    # define heat capacity and latent heat coefficient
    # -----------------------------------------------------------------------

    ice_precip_flux = 0.0
    lhl = lv00 + d0_vap * t
    lhi = constants.LI00 + constants.DC_ICE * t
    q_liq = liquid + rain
    q_sol = ice + snow + graupel
    cvm = c_air + vapor * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
    lcpk = lhl / cvm
    icpk = lhi / cvm

    # -----------------------------------------------------------------------
    # melting of cloud_ice (before fall) :
    # -----------------------------------------------------------------------

    tc = t - constants.TICE
    if is_frozen == False and (ice > constants.QCMIN and tc > 0.0):  # noqa
        sink = min(ice, fac_imlt * tc / icpk)
        if ql_mlt - liquid > 0:
            ans = ql_mlt - liquid
        else:
            ans = 0
        tmp = min(sink, ans)
        liquid = liquid + tmp
        rain = rain + sink - tmp
        ice = ice - sink
        q_liq = q_liq + sink
        q_sol = q_sol - sink
        cvm = c_air + vapor * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
        t = t - sink * lhi / cvm

    return t, liquid, rain, ice, cvm, lhi, icpk, ice_precip_flux


def setup(
    t: FloatField,
    mixing_ratio_vapor: FloatField,
    mixing_ratio_liquid: FloatField,
    mixing_ratio_rain: FloatField,
    mixing_ratio_graupel: FloatField,
    mixing_ratio_snow: FloatField,
    mixing_ratio_ice: FloatField,
    dz: FloatField,
    ice_precip_flux: FloatField,
    z_interface: FloatField,
    z_interface_modified: FloatField,
    lhi: FloatField,
    icpk: FloatField,
    cvm: FloatField,
    internal_rain: FloatFieldIJ,
    internal_graupel: FloatFieldIJ,
    internal_snow: FloatFieldIJ,
    internal_ice: FloatFieldIJ,
):
    """Perform initial calculations of the terminal fall module. Get extra parameters
    and compute prefall melting.

    Args:
        t (FloatField): temperature (Kelvin)
        mixing_ratio_vapor (FloatField): water vapor mixing ratio (kg/kg)
        mixing_ratio_liquid (FloatField): liquid water mixing ratio (kg/kg)
        mixing_ratio_rain (FloatField): rain mixing ratio (kg/kg)
        mixing_ratio_graupel (FloatField): graupel mixing ratio (kg/kg)
        mixing_ratio_snow (FloatField): snow mixing ratio (kg/kg)
        mixing_ratio_ice (FloatField): ice mixing ratio (kg/kg)
        dz (FloatField): change in height between model layers (m)
        ice_precip_flux (FloatField): ice precipitation flux (kg m^-2 s^-1)
        z_interface (FloatField): grid center height (m)
        z_interface_modified (FloatField): grid edge height (m)
        lhi (FloatField): _description_
        icpk (FloatField): _description_
        cvm (FloatField): _description_
        internal_rain (FloatFieldIJ): precipitable rain from within the terminal fall
            module (kg m^-2 s^-1)
        internal_graupel (FloatFieldIJ): precipitable graupel from within the terminal fall
            module (kg m^-2 s^-1)
        internal_snow (FloatFieldIJ): precipitable snow from within the terminal fall
            module (kg m^-2 s^-1)
        internal_ice (FloatFieldIJ): precipitable ice from within the terminal fall
            module (kg m^-2 s^-1)
    """
    from __externals__ import c_air, c_vap, d0_vap, dts, lv00, ql_mlt, tau_imlt

    # determine frozen levels
    # later operations will only be executed if frozen/melted
    # True = frozen, False = melted
    with computation(PARALLEL), interval(...):
        if t <= constants.TICE:
            is_frozen = True
        else:
            is_frozen = False

    # we only want the melting layer farthest from the surface
    with computation(FORWARD), interval(1, None):
        if is_frozen[0, 0, -1] == False and is_frozen == True:  # type: ignore[index] # noqa
            is_frozen = False

    # force surface to "melt" for later calculations
    with computation(PARALLEL), interval(-1, None):
        is_frozen = False

    with computation(PARALLEL), interval(...):
        (
            t,
            mixing_ratio_liquid,
            mixing_ratio_rain,
            mixing_ratio_ice,
            cvm,
            lhi,
            icpk,
            ice_precip_flux,
        ) = prefall_melting(
            t,
            mixing_ratio_vapor,
            mixing_ratio_liquid,
            mixing_ratio_rain,
            mixing_ratio_graupel,
            mixing_ratio_snow,
            mixing_ratio_ice,
            ice_precip_flux,
            is_frozen,
            c_air,
            c_vap,
            d0_vap,
            dts,
            lv00,
            ql_mlt,
            tau_imlt,
        )

    # if timestep is too small turn off melting (only calculate at surface)
    with computation(PARALLEL), interval(0, -1):
        if dts < 300.0:
            is_frozen = True

    with computation(BACKWARD), interval(...):
        z_interface = z_interface[0, 0, 1] - dz  # dz < 0

    with computation(FORWARD), interval(0, 1):
        z_interface_modified = z_interface

    # -----------------------------------------------------------------------
    # update capacity heat and latent heat coefficient
    # -----------------------------------------------------------------------

    with computation(PARALLEL), interval(...):
        if is_frozen == False:  # noqa
            lhi = constants.LI00 + constants.DC_ICE * t
            icpk = lhi / cvm

    with computation(FORWARD), interval(-1, None):
        lhi = constants.LI00 + constants.DC_ICE * t
        icpk = lhi / cvm

    # zero local precipitation values
    with computation(FORWARD), interval(0, 1):
        internal_rain = 0
        internal_graupel = 0
        internal_snow = 0
        internal_ice = 0


def update_outputs(
    rain: FloatFieldIJ,
    graupel: FloatFieldIJ,
    snow: FloatFieldIJ,
    ice: FloatFieldIJ,
    internal_rain: FloatFieldIJ,
    internal_graupel: FloatFieldIJ,
    internal_snow: FloatFieldIJ,
    internal_ice: FloatFieldIJ,
):
    """Ensure information for all precipitates is pushed back to the rest of the model

    Args:
        rain (FloatFieldIJ): total precipitated rain from microphysics (kg m^-2 s^-1)
        graupel (FloatFieldIJ): total precipitated graupel from microphysics (kg m^-2 s^-1)
        snow (FloatFieldIJ): total precipitated snow from microphysics (kg m^-2 s^-1)
        ice (FloatFieldIJ): total precipitated ice from microphysics (kg m^-2 s^-1)
        internal_rain (FloatFieldIJ): precipitated rain from the terminal fall module (kg m^-2 s^-1)
        internal_graupel (FloatFieldIJ): precipitated graupel from the terminal fall module (kg m^-2 s^-1)
        internal_snow (FloatFieldIJ): precipitated snow from the terminal fall module (kg m^-2 s^-1)
        internal_ice (FloatFieldIJ): precipitated ice from the terminal fall module (kg m^-2 s^-1)
    """
    with computation(FORWARD), interval(0, 1):
        rain = rain + internal_rain  # from melted snow & ice that reached the ground
        graupel = graupel + internal_graupel
        snow = snow + internal_snow
        ice = ice + internal_ice


@dataclasses.dataclass
class Locals(LocalState):
    lhi: Local = dataclasses.field(
        metadata={
            "name": "lhi",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    icpk: Local = dataclasses.field(
        metadata={
            "name": "icpk",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    cvm: Local = dataclasses.field(
        metadata={
            "name": "cvm",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    mass: Local = dataclasses.field(
        metadata={
            "name": "m1",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dmass: Local = dataclasses.field(
        metadata={
            "name": "dm",
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
    z_interface_modified: Local = dataclasses.field(
        metadata={
            "name": "z_interface_modified",
            "dims": [I_DIM, J_DIM, K_INTERFACE_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    rain: Local = dataclasses.field(
        metadata={
            "name": "rain",
            "dims": [I_DIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    graupel: Local = dataclasses.field(
        metadata={
            "name": "graupel",
            "dims": [I_DIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    snow: Local = dataclasses.field(
        metadata={
            "name": "snow",
            "dims": [I_DIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    ice: Local = dataclasses.field(
        metadata={
            "name": "ice",
            "dims": [I_DIM, J_DIM],
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


class GFDL1MTerminalFall(NDSLRuntime):
    """
    Calculate terminal fall speed, accounting for melting of ice, snow, and graupel during fall
    """

    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GFDL1MConfig,
        config_dependent_constants: GFDL1MDriverConfigDependentConstants,
    ):
        """Initialize the TerminalFall module

        Args:
            stencil_factory (StencilFactory)
            quantity_factory (QuantityFactory)
            config (GFDL1MConfig)
            config_dependent_constants (GFDL1MDriverConfigDependentConstants)
        """
        # initialize NDSLRuntime
        super().__init__(stencil_factory)

        self.config = config

        # initialize temporaries
        self._locals: Locals = Locals.make_locals(quantity_factory)

        # construct stencils
        self._setup = stencil_factory.from_dims_halo(
            func=setup,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "dts": config_dependent_constants.DTS,
                "tau_imlt": config.TAU_IMLT,
                "ql_mlt": config.QL_MLT,
                "c_air": config_dependent_constants.C_AIR,
                "c_vap": config_dependent_constants.C_VAP,
                "d0_vap": config_dependent_constants.D0_VAP,
                "lv00": config_dependent_constants.LV00,
            },
        )

        self._check_precip_get_zt = stencil_factory.from_dims_halo(
            func=check_precip_get_zt,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "dts": config_dependent_constants.DTS,
            },
        )

        self._update_dmass = stencil_factory.from_dims_halo(
            func=update_dmass,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "do_sedi_w": config.DO_SEDI_W,
            },
        )

        self._implicit_fall = stencil_factory.from_dims_halo(
            func=implicit_fall,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "dts": config_dependent_constants.DTS,
            },
        )

        self._update_w = stencil_factory.from_dims_halo(
            func=update_w,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "do_sedi_w": config.DO_SEDI_W,
            },
        )

        self._reset = stencil_factory.from_dims_halo(
            func=reset,
            compute_dims=[I_DIM, J_DIM, K_DIM],
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

    def __call__(
        self,
        t: Quantity,
        w: Quantity,
        mixing_ratio_vapor: Quantity,
        mixing_ratio_liquid: Quantity,
        mixing_ratio_rain: Quantity,
        mixing_ratio_graupel: Quantity,
        mixing_ratio_snow: Quantity,
        mixing_ratio_ice: Quantity,
        dz: Quantity,
        dp: Quantity,
        terminal_velocity_graupel: Quantity,
        terminal_velocity_snow: Quantity,
        terminal_velocity_ice: Quantity,
        rain: Quantity,
        graupel: Quantity,
        snow: Quantity,
        ice: Quantity,
        ice_precip_flux: Quantity,
    ):
        """
        Calculate terminal fall speed, accounting for melting of ice, snow, and graupel during fall

        Args:
            t (Quantity): (inout) temperature (K)
            w (Quantity): (inout) vertical motion (m/s)
            mixing_ratio_vapor (Quantity): (inout) water vapor mixing ratio (kg/kg)
            mixing_ratio_liquid (Quantity): (inout) liquid water mixing ratio (kg/kg)
            mixing_ratio_rain (Quantity): (inout) rain mixing ratio (kg/kg)
            mixing_ratio_graupel (Quantity): (inout) graupel mixing ratio (kg/kg)
            mixing_ratio_snow (Quantity): (inout) snow mixing ratio (kg/kg)
            mixing_ratio_ice (Quantity): (inout) ice mixing ratio (kg/kg)
            dz (Quantity): (in) change in height between model layers (m)
            dp (Quantity): (in) change in pressure between model layers (Pa)
            terminal_velocity_graupel (Quantity): (in) terminal fall speed of graupel
            terminal_velocity_snow (Quantity): (in) terminal fall speed of snow
            terminal_velocity_ice (Quantity): (in) terminal fall speed of ice
            rain (Quantity): (out) rain precipitation (kg m^-2 s^-1)
            graupel (Quantity): (out) graupel precipitation (kg m^-2 s^-1)
            snow (Quantity): (out) snow precipitation (kg m^-2 s^-1)
            ice (Quantity): (out) ice precipitation (kg m^-2 s^-1)
            ice_precip_flux (Quantity): (out) ice precipitation flux (kg m^-2 s^-1)
        """

        # Reset locals
        self._set_value(self._locals.lhi, Float(0))
        self._set_value(self._locals.icpk, Float(0))
        self._set_value(self._locals.cvm, Float(0))
        self._set_value(self._locals.mass, Float(0))
        self._set_value(self._locals.dmass, Float(0))
        self._set_value_K_interface(self._locals.z_interface, Float(0))
        self._set_value_K_interface(self._locals.z_interface_modified, Float(0))
        self._set_value_IJ(self._locals.rain, Float(0))
        self._set_value_IJ(self._locals.graupel, Float(0))
        self._set_value_IJ(self._locals.snow, Float(0))
        self._set_value_IJ(self._locals.ice, Float(0))
        self._set_IJ_mask(self._locals.precip_fall, False)

        self._setup(
            t=t,
            mixing_ratio_vapor=mixing_ratio_vapor,
            mixing_ratio_liquid=mixing_ratio_liquid,
            mixing_ratio_rain=mixing_ratio_rain,
            mixing_ratio_graupel=mixing_ratio_graupel,
            mixing_ratio_snow=mixing_ratio_snow,
            mixing_ratio_ice=mixing_ratio_ice,
            dz=dz,
            ice_precip_flux=ice_precip_flux,
            z_interface=self._locals.z_interface,
            z_interface_modified=self._locals.z_interface_modified,
            lhi=self._locals.lhi,
            icpk=self._locals.icpk,
            cvm=self._locals.cvm,
            internal_rain=self._locals.rain,
            internal_graupel=self._locals.graupel,
            internal_snow=self._locals.snow,
            internal_ice=self._locals.ice,
        )

        # melting of falling cloud ice into rain
        if self.config.VI_FAC < 1.0e-5:
            self._set_value_IJ(self._locals.ice, 0)
        else:
            self._check_precip_get_zt(
                mixing_ratio=mixing_ratio_ice,
                terminal_speed=terminal_velocity_ice,
                z_interface=self._locals.z_interface,
                z_interface_modified=self._locals.z_interface_modified,
                internal_precip=self._locals.ice,
                precip_fall=self._locals.precip_fall,
            )

            self._update_dmass(
                dmass=self._locals.dmass,
                dp=dp,
                mixing_ratio_vapor=mixing_ratio_vapor,
                mixing_ratio_liquid=mixing_ratio_liquid,
                mixing_ratio_rain=mixing_ratio_rain,
                mixing_ratio_ice=mixing_ratio_ice,
                mixing_ratio_snow=mixing_ratio_snow,
                mixing_ratio_graupel=mixing_ratio_graupel,
                precip_fall=self._locals.precip_fall,
            )

            if not self.config.USE_PPM:
                self._implicit_fall(
                    mixing_ratio=mixing_ratio_ice,
                    terminal_speed=terminal_velocity_ice,
                    z_interface=self._locals.z_interface,
                    dp=dp,
                    mass=self._locals.mass,
                    precip_flux=ice_precip_flux,
                    precip=self._locals.ice,
                    precip_fall=self._locals.precip_fall,
                )

            self._update_w(
                w=w,
                dmass=self._locals.dmass,
                mass=self._locals.mass,
                terminal_speed=terminal_velocity_ice,
                precip_fall=self._locals.precip_fall,
            )

            self._reset(
                mass=self._locals.mass,
                precip_fall=self._locals.precip_fall,
            )

        # melting of falling snow into rain
        self._check_precip_get_zt(
            mixing_ratio=mixing_ratio_snow,
            terminal_speed=terminal_velocity_snow,
            z_interface=self._locals.z_interface,
            z_interface_modified=self._locals.z_interface_modified,
            internal_precip=self._locals.snow,
            precip_fall=self._locals.precip_fall,
        )

        self._update_dmass(
            dmass=self._locals.dmass,
            dp=dp,
            mixing_ratio_vapor=mixing_ratio_vapor,
            mixing_ratio_liquid=mixing_ratio_liquid,
            mixing_ratio_rain=mixing_ratio_rain,
            mixing_ratio_ice=mixing_ratio_ice,
            mixing_ratio_snow=mixing_ratio_snow,
            mixing_ratio_graupel=mixing_ratio_graupel,
            precip_fall=self._locals.precip_fall,
        )

        if not self.config.USE_PPM:
            self._implicit_fall(
                mixing_ratio=mixing_ratio_snow,
                terminal_speed=terminal_velocity_snow,
                z_interface=self._locals.z_interface,
                dp=dp,
                mass=self._locals.mass,
                precip_flux=ice_precip_flux,
                precip=self._locals.snow,
                precip_fall=self._locals.precip_fall,
            )

        self._update_w(
            w=w,
            dmass=self._locals.dmass,
            mass=self._locals.mass,
            terminal_speed=terminal_velocity_snow,
            precip_fall=self._locals.precip_fall,
        )

        self._reset(
            mass=self._locals.mass,
            precip_fall=self._locals.precip_fall,
        )

        # melting of falling graupel into rain
        self._check_precip_get_zt(
            mixing_ratio=mixing_ratio_graupel,
            terminal_speed=terminal_velocity_graupel,
            z_interface=self._locals.z_interface,
            z_interface_modified=self._locals.z_interface_modified,
            internal_precip=self._locals.graupel,
            precip_fall=self._locals.precip_fall,
        )

        self._update_dmass(
            dmass=self._locals.dmass,
            dp=dp,
            mixing_ratio_vapor=mixing_ratio_vapor,
            mixing_ratio_liquid=mixing_ratio_liquid,
            mixing_ratio_rain=mixing_ratio_rain,
            mixing_ratio_ice=mixing_ratio_ice,
            mixing_ratio_snow=mixing_ratio_snow,
            mixing_ratio_graupel=mixing_ratio_graupel,
            precip_fall=self._locals.precip_fall,
        )

        if not self.config.USE_PPM:
            self._implicit_fall(
                mixing_ratio=mixing_ratio_graupel,
                terminal_speed=terminal_velocity_graupel,
                z_interface=self._locals.z_interface,
                dp=dp,
                mass=self._locals.mass,
                precip_flux=ice_precip_flux,
                precip=self._locals.graupel,
                precip_fall=self._locals.precip_fall,
            )

        self._update_w(
            w=w,
            dmass=self._locals.dmass,
            mass=self._locals.mass,
            terminal_speed=terminal_velocity_graupel,
            precip_fall=self._locals.precip_fall,
        )

        self._reset(
            mass=self._locals.mass,
            precip_fall=self._locals.precip_fall,
        )

        # ensure information for all precipitates is pushed back to the rest of the model
        self._update_outputs(
            rain=rain,
            graupel=graupel,
            snow=snow,
            ice=ice,
            internal_rain=self._locals.rain,
            internal_graupel=self._locals.graupel,
            internal_snow=self._locals.snow,
            internal_ice=self._locals.ice,
        )
