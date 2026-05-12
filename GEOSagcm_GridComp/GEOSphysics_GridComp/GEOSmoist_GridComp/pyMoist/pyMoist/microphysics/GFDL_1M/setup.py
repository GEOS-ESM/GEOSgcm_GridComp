import dataclasses

from ndsl import Local, LocalState, NDSLRuntime, Quantity, QuantityFactory, StencilFactory
from ndsl.constants import I_DIM, J_DIM, K_DIM, K_INTERFACE_DIM
from ndsl.dsl.gt4py import BACKWARD, FORWARD, PARALLEL, K, computation, function, interval, log
from ndsl.dsl.typing import BoolFieldIJ, Float, FloatField, FloatFieldIJ, IntFieldIJ

from pyMoist.constants import MAPL_ALHL, MAPL_CP, MAPL_CPDRY, MAPL_CPVAP, MAPL_GRAV, MAPL_KAPPA, MAPL_P00, MAPL_RGAS, MAPL_RVAP
from pyMoist.microphysics.GFDL_1M.config import GFDL1MConfig
from pyMoist.microphysics.GFDL_1M.shared_stencils import prepare_tendencies
from pyMoist.saturation_tables import GlobalTable_saturation_tables, SaturationVaporPressureTable, saturation_specific_humidity
from pyMoist.shared.interpolations import vertical_interpolation


def set_unused_to_zero(
    shallow_convective_precipitation: FloatFieldIJ,
    deep_convective_precipitation: FloatFieldIJ,
    anvil_precipitation: FloatFieldIJ,
    shallow_convective_snow: FloatFieldIJ,
    deep_convective_snow: FloatFieldIJ,
    anvil_snow: FloatFieldIJ,
):
    """Set unused fields to zero. These fields are read in by the fortran,
    immediately set to zero, and never touched again.

    Args:
        shallow_convective_precipitation (FloatFieldIJ)
        deep_convective_precipitation (FloatFieldIJ)
        anvil_precipitation (FloatFieldIJ)
        shallow_convective_snow (FloatFieldIJ)
        deep_convective_snow (FloatFieldIJ)
        anvil_snow (FloatFieldIJ)
    """
    with computation(FORWARD), interval(0, 1):
        shallow_convective_precipitation = 0.0
        deep_convective_precipitation = 0.0
        anvil_precipitation = 0.0
        shallow_convective_snow = 0.0
        deep_convective_snow = 0.0
        anvil_snow = 0.0


def calculate_derived_states(
    p_interface: FloatField,
    p_interface_mb: FloatField,
    p_mb: FloatField,
    geopotential_height_interface: FloatField,
    edge_height_above_surface: FloatField,
    layer_height_above_surface: FloatField,
    layer_thickness: FloatField,
    layer_thickness_negative: FloatField,
    dp: FloatField,
    mass: FloatField,
    mass_inverse: FloatField,
    t: FloatField,
    esx: GlobalTable_saturation_tables,
    sat: FloatField,
    dsat: FloatField,
    u: FloatField,
    u_unmodified: FloatField,
    v: FloatField,
    v_unmodified: FloatField,
    th: FloatField,
):
    """
    Computes derived state fields required for the rest of the GFDL single moment
    microphysics module.

    Args:
        p_interface (FloatField)
        p_interface_mb (FloatField)
        p_mb (FloatField)
        geopotential_height_interface (FloatField)
        edge_height_above_surface (FloatField)
        layer_height_above_surface (FloatField)
        layer_thickness (FloatField)
        layer_thickness_negative (FloatField)
        dp (FloatField)
        mass (FloatField)
        mass_inverse (FloatField)
        t (FloatField)
        esx (GlobalTable_saturation_tables)
        sat (FloatField)
        dsat (FloatField)
        u (FloatField)
        u_unmodified (FloatField)
        v (FloatField)
        v_unmodified (FloatField)
        th (FloatField)
    """
    from __externals__ import k_end

    with computation(PARALLEL), interval(...):
        p_interface_mb = p_interface * 0.01
        edge_height_above_surface = geopotential_height_interface - geopotential_height_interface.at(K=k_end)
    with computation(PARALLEL), interval(0, -1):
        p_mb = 0.5 * (p_interface_mb + p_interface_mb[0, 0, 1])
        layer_height_above_surface = 0.5 * (edge_height_above_surface + edge_height_above_surface[0, 0, 1])
        layer_thickness = edge_height_above_surface - edge_height_above_surface[0, 0, 1]
        layer_thickness_negative = -1.0 * layer_thickness
        dp = p_interface[0, 0, 1] - p_interface
        mass = dp / MAPL_GRAV
        mass_inverse = 1 / mass
        sat, dsat = saturation_specific_humidity(t=t, p=p_mb * 100.0, esx=esx)
        u_unmodified = u
        v_unmodified = v
        th = (100.0 * p_mb / MAPL_P00) ** (MAPL_KAPPA)
        th = t / th


@function
def find_t_lcl(
    t: Float,
    rh: Float,
):
    """
    Computes the LCL temperature

    Arguments:
        t (Float): temperature at surface (K)
        rh (Float): relative humidity at surface

    Returns:
        tlcl: LCL temperature
    """
    term1 = 1.0 / (t - 55.0)
    term2 = log(max(0.1, rh) / 100.0) / 2840.0
    denom = term1 - term2
    tlcl = (1.0 / denom) + 55.0
    return tlcl


def find_lcl_level(
    t: FloatField,
    p_mb: FloatField,
    vapor: FloatField,
    esx: GlobalTable_saturation_tables,
    lcl_level: IntFieldIJ,
):
    """
    Find the level of the lifted condensation level (LCL).

    Arguments:
        t (FloatField): (in) Atmospheric temperature (K)
        p_mb (FloatField): (in) pressure (mb)
        vapor (FloatField): (in) water vapor mixing radio (kg/kg)
        esx (GlobalTable_saturation_tables): (in) saturation vapor pressure table, details unknown
        lcl_level (IntFieldIJ): (out) LCL level
    """
    from __externals__ import k_end

    # set up mask to stop computation
    with computation(FORWARD), interval(0, 1):
        found_level: BoolFieldIJ = False

    # get LCL pressure
    with computation(PARALLEL), interval(-1, None):
        qsat, _ = saturation_specific_humidity(t=t, p=p_mb * 100.0, esx=esx)
        rhsfc = 100.0 * vapor / qsat
        tlcl = find_t_lcl(t=t, rh=rhsfc)
        rm = (1.0 - vapor) * MAPL_RGAS + vapor * MAPL_RVAP
        cpm = (1.0 - vapor) * MAPL_CPDRY + vapor * MAPL_CPVAP
        plcl = p_mb * ((tlcl / t) ** (cpm / rm))

    # find nearest level <= LCL pressure
    with computation(BACKWARD), interval(...):
        if found_level == False:  # noqa
            lcl_level = K
        if p_mb <= plcl.at(K=k_end):
            found_level = True


def update_lcl_height(
    layer_height_above_surface: FloatField,
    lcl_level: IntFieldIJ,
    lcl_height: FloatFieldIJ,
):
    """Update LCL height

    Args:
        layer_height_above_surface (FloatField)
        lcl_level (IntFieldIJ)
        lcl_height (FloatFieldIJ)
    """
    with computation(FORWARD), interval(0, 1):
        lcl_height = layer_height_above_surface.at(K=lcl_level)


def compute_estimated_inversion_strength(
    t: FloatField,
    th: FloatField,
    layer_height_above_surface: FloatField,
    t700: FloatFieldIJ,
    th700: FloatFieldIJ,
    z700: FloatFieldIJ,
    lcl_level: IntFieldIJ,
    esx: GlobalTable_saturation_tables,
    lower_tropospheric_stability: FloatFieldIJ,
    estimated_inversion_strength: FloatFieldIJ,
):
    """
    Find estimated inversion strength. Returns Estimated Inversion
    Strength (K) according to Wood and Betherton, J.Climate, 2006.
    Based on Fortran code written by Donifan Barahona.

    Args:
        t (FloatField)
        th (FloatField)
        layer_height_above_surface (FloatField)
        t700 (FloatFieldIJ)
        th700 (FloatFieldIJ)
        z700 (FloatFieldIJ)
        lcl_level (IntFieldIJ)
        esx (GlobalTable_saturation_tables)
        lower_tropospheric_stability (FloatFieldIJ)
        estimated_inversion_strength (FloatFieldIJ)
    """
    with computation(FORWARD), interval(0, 1):
        from __externals__ import k_end

        lower_tropospheric_stability = th700 - th.at(K=k_end)
        lcl_height = layer_height_above_surface.at(K=lcl_level - 1)

        # Simplified single adiabat eq4 of https://doi.org/10.1175/JCLI3988.1
        t850 = 0.5 * (t.at(K=k_end) + t700)
        qs850, _ = saturation_specific_humidity(t=t850, p=100.0 * 850.0, esx=esx)
        gamma850 = (1.0 + (MAPL_ALHL * qs850 / (MAPL_RGAS * t850))) / (1.0 + (MAPL_ALHL * MAPL_ALHL * qs850 / (MAPL_CP * MAPL_RVAP * t850 * t850)))
        gamma850 = MAPL_GRAV / MAPL_CP * (1.0 - gamma850)
        estimated_inversion_strength = lower_tropospheric_stability - gamma850 * (z700 - lcl_height)


def update_precipitation(
    mixing_ratio: FloatField,
    shallow_convection_values: FloatField,
):
    """Update precipitate mixing ratio

    Args:
        mixing_ratio (FloatField)
        shallow_convection_values (FloatField)
    """
    from __externals__ import DT_MOIST

    with computation(PARALLEL), interval(...):
        mixing_ratio = mixing_ratio + shallow_convection_values * DT_MOIST


@dataclasses.dataclass
class GFDL1MSetupLocals(LocalState):
    th: Local = dataclasses.field(
        metadata={
            "name": "th",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    t700: Local = dataclasses.field(
        metadata={
            "name": "t700",
            "dims": [I_DIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    th700: Local = dataclasses.field(
        metadata={
            "name": "th700",
            "dims": [I_DIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    z700: Local = dataclasses.field(
        metadata={
            "name": "z700",
            "dims": [I_DIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )


class GFDL1MSetup(NDSLRuntime):
    """
    Perform the following functions to setup GFDL Single Moment microphysics:

    prepare_tendencies: preloads macrophysics tendencies for post-phase_change calculations
    calculate_derived_states: computes fields required for the module but not provided by the module
    find_k_lcl: identifies the LCL level
    update_z_lcl (conditional): computes the geometric height of the LCL and returns it to the model
    vertical_interpolation: interpolates various fields to the desired geometric height
    find_eis: computes the estimated inversion strength
    update_precipitation (conditional): updates precipitation (rain and snow) using shallow convection values
    """

    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GFDL1MConfig,
        saturation_tables: SaturationVaporPressureTable,
    ):
        """Initialize the GFDL1M microphysics setup class

        Args:
            stencil_factory (StencilFactory)
            quantity_factory (QuantityFactory)
            config (GFDL1MConfig)
            saturation_tables (SaturationVaporPressureTable)
        """
        # init NDSLRuntime
        super().__init__(stencil_factory)

        # make configuration and saturation tables visible at runtime
        self.config = config
        self.saturation_tables = saturation_tables

        # initialize locals
        self._locals = GFDL1MSetupLocals.make_locals(quantity_factory)

        # construct stencils
        self._set_unused_to_zero = stencil_factory.from_dims_halo(
            func=set_unused_to_zero,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )

        self._prepare_tendencies = stencil_factory.from_dims_halo(
            func=prepare_tendencies,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )

        self._calculate_derived_states = stencil_factory.from_dims_halo(
            func=calculate_derived_states,
            compute_dims=[I_DIM, J_DIM, K_INTERFACE_DIM],
        )

        self._find_lcl_level = stencil_factory.from_dims_halo(
            func=find_lcl_level,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )

        self._update_lcl_height = stencil_factory.from_dims_halo(
            func=update_lcl_height,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )

        self._vertical_interpolation = stencil_factory.from_dims_halo(
            func=vertical_interpolation,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )

        self._compute_estimated_inversion_strength = stencil_factory.from_dims_halo(
            func=compute_estimated_inversion_strength,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )

        self._update_precipitation = stencil_factory.from_dims_halo(
            func=update_precipitation,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "DT_MOIST": config.DT_MOIST,
            },
        )

        # Dev NOTE: this is an orchestration workaround. Direct call to
        #           `self.saturation_tables.X` fails closure capture for
        #           argument reconstruction at call time
        self._esx = self.saturation_tables.esx

    def __call__(
        self,
        p_interface: Quantity,
        z_interface: Quantity,
        u: Quantity,
        v: Quantity,
        t: Quantity,
        lcl_height: Quantity,
        lower_tropospheric_stability: Quantity,
        estimated_inversion_strength: Quantity,
        mixing_ratio_vapor: Quantity,
        mixing_ratio_rain: Quantity,
        mixing_ratio_snow: Quantity,
        mixing_ratio_graupel: Quantity,
        mixing_ratio_convective_liquid: Quantity,
        mixing_ratio_convective_ice: Quantity,
        mixing_ratio_large_scale_liquid: Quantity,
        mixing_ratio_large_scale_ice: Quantity,
        cloud_fraction_convective: Quantity,
        cloud_fraction_large_scale: Quantity,
        shallow_convection_rain: Quantity,
        shallow_convection_snow: Quantity,
        dudt_macro: Quantity,
        dvdt_macro: Quantity,
        dtdt_macro: Quantity,
        dvapordt_macro: Quantity,
        dliquiddt_macro: Quantity,
        dicedt_macro: Quantity,
        dcloud_fractiondt_macro: Quantity,
        draindt_macro: Quantity,
        dsnowdt_macro: Quantity,
        dgraupeldt_macro: Quantity,
        shallow_convective_precipitation: Quantity,
        deep_convective_precipitation: Quantity,
        anvil_precipitation: Quantity,
        shallow_convective_snow: Quantity,
        deep_convective_snow: Quantity,
        anvil_snow: Quantity,
        local_p_mb: Quantity,
        local_p_interface_mb: Quantity,
        local_edge_height_above_surface: Quantity,
        local_layer_height_above_surface: Quantity,
        local_layer_thickness: Quantity,
        local_layer_thickness_negative: Quantity,
        local_dp: Quantity,
        local_mass: Quantity,
        local_mass_inverse: Quantity,
        local_saturation_specific_humidity: Quantity,
        local_dsaturation_specific_humidity: Quantity,
        local_u_unmodified: Quantity,
        local_v_unmodified: Quantity,
        local_lcl_level: Quantity,
    ):
        """Setup the GFDL1M microphysics module

        Args:
            p_interface (Quantity)
            z_interface (Quantity)
            u (Quantity)
            v (Quantity)
            t (Quantity)
            lcl_height (Quantity)
            lower_tropospheric_stability (Quantity)
            estimated_inversion_strength (Quantity)
            mixing_ratio_vapor (Quantity)
            mixing_ratio_rain (Quantity)
            mixing_ratio_snow (Quantity)
            mixing_ratio_graupel (Quantity)
            mixing_ratio_convective_liquid (Quantity)
            mixing_ratio_convective_ice (Quantity)
            mixing_ratio_large_scale_liquid (Quantity)
            mixing_ratio_large_scale_ice (Quantity)
            cloud_fraction_convective (Quantity)
            cloud_fraction_large_scale (Quantity)
            shallow_convection_rain (Quantity)
            shallow_convection_snow (Quantity)
            dudt_macro (Quantity)
            dvdt_macro (Quantity)
            dtdt_macro (Quantity)
            dvapordt_macro (Quantity)
            dliquiddt_macro (Quantity)
            dicedt_macro (Quantity)
            dcloud_fractiondt_macro (Quantity)
            draindt_macro (Quantity)
            dsnowdt_macro (Quantity)
            dgraupeldt_macro (Quantity)
            shallow_convective_precipitation (Quantity)
            deep_convective_precipitation (Quantity)
            anvil_precipitation (Quantity)
            shallow_convective_snow (Quantity)
            deep_convective_snow (Quantity)
            anvil_snow (Quantity)
            local_p_mb (Quantity)
            local_p_interface_mb (Quantity)
            local_edge_height_above_surface (Quantity)
            local_layer_height_above_surface (Quantity)
            local_layer_thickness (Quantity)
            local_layer_thickness_negative (Quantity)
            local_dp (Quantity)
            local_mass (Quantity)
            local_mass_inverse (Quantity)
            local_saturation_specific_humidity (Quantity)
            local_dsaturation_specific_humidity (Quantity)
            local_u_unmodified (Quantity)
            local_v_unmodified (Quantity)
            local_lcl_level (Quantity)
        """
        # set unused fields to zero
        self._set_unused_to_zero(
            shallow_convective_precipitation=shallow_convective_precipitation,
            deep_convective_precipitation=deep_convective_precipitation,
            anvil_precipitation=anvil_precipitation,
            shallow_convective_snow=shallow_convective_snow,
            deep_convective_snow=deep_convective_snow,
            anvil_snow=anvil_snow,
        )

        # prepare macrophysics tendencies
        self._prepare_tendencies(
            u=u,
            v=v,
            t=t,
            vapor=mixing_ratio_vapor,
            rain=mixing_ratio_rain,
            snow=mixing_ratio_snow,
            graupel=mixing_ratio_graupel,
            convective_liquid=mixing_ratio_convective_liquid,
            convective_ice=mixing_ratio_convective_ice,
            large_scale_liquid=mixing_ratio_large_scale_liquid,
            large_scale_ice=mixing_ratio_large_scale_ice,
            convective_cloud_fraction=cloud_fraction_convective,
            large_scale_cloud_fraction=cloud_fraction_large_scale,
            du_dt=dudt_macro,
            dv_dt=dvdt_macro,
            dt_dt=dtdt_macro,
            dvapor_dt=dvapordt_macro,
            dliquid_dt=dliquiddt_macro,
            dice_dt=dicedt_macro,
            dcloud_fraction_dt=dcloud_fractiondt_macro,
            drain_dt=draindt_macro,
            dsnow_dt=dsnowdt_macro,
            dgraupel_dt=dgraupeldt_macro,
        )
        self._calculate_derived_states(
            p_interface=p_interface,
            p_interface_mb=local_p_interface_mb,
            p_mb=local_p_mb,
            geopotential_height_interface=z_interface,
            edge_height_above_surface=local_edge_height_above_surface,
            layer_height_above_surface=local_layer_height_above_surface,
            layer_thickness=local_layer_thickness,
            layer_thickness_negative=local_layer_thickness_negative,
            dp=local_dp,
            mass=local_mass,
            mass_inverse=local_mass_inverse,
            t=t,
            esx=self.saturation_tables.esx,
            sat=local_saturation_specific_humidity,
            dsat=local_dsaturation_specific_humidity,
            u=u,
            u_unmodified=local_u_unmodified,
            v=v,
            v_unmodified=local_v_unmodified,
            th=self._locals.th,
        )

        self._find_lcl_level(
            t=t,
            p_mb=local_p_mb,
            vapor=mixing_ratio_vapor,
            esx=self._esx,
            lcl_level=local_lcl_level,
        )

        if lcl_height is not None:
            self._update_lcl_height(
                layer_height_above_surface=local_layer_height_above_surface,
                lcl_level=local_lcl_level,
                lcl_height=lcl_height,
            )

        self._vertical_interpolation(
            field=self._locals.th,
            interpolated_field=self._locals.th700,
            p_interface_mb=local_p_interface_mb,
            target_pressure=Float(70000.0),
        )

        self._vertical_interpolation(
            field=t,
            interpolated_field=self._locals.t700,
            p_interface_mb=local_p_interface_mb,
            target_pressure=Float(70000.0),
        )

        self._vertical_interpolation(
            field=local_layer_height_above_surface,
            interpolated_field=self._locals.z700,
            p_interface_mb=local_p_interface_mb,
            target_pressure=Float(70000.0),
        )

        self._compute_estimated_inversion_strength(
            t=t,
            th=self._locals.th,
            layer_height_above_surface=local_layer_height_above_surface,
            t700=self._locals.t700,
            th700=self._locals.th700,
            z700=self._locals.z700,
            lcl_level=local_lcl_level,
            esx=self._esx,
            lower_tropospheric_stability=lower_tropospheric_stability,
            estimated_inversion_strength=estimated_inversion_strength,
        )

        if shallow_convection_rain is not None:
            self._update_precipitation(mixing_ratio_rain, shallow_convection_rain)

        if shallow_convection_snow is not None:
            self._update_precipitation(mixing_ratio_snow, shallow_convection_snow)
