from ndsl import StencilFactory, QuantityFactory, orchestrate, State, Quantity, NDSLRuntime
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.gt4py import BACKWARD, FORWARD, PARALLEL, K, computation, function, interval, log
from ndsl.dsl.typing import BoolFieldIJ, Float, FloatField, FloatFieldIJ, IntFieldIJ
from pyMoist.constants import (
    MAPL_ALHL,
    MAPL_CP,
    MAPL_CPDRY,
    MAPL_CPVAP,
    MAPL_GRAV,
    MAPL_KAPPA,
    MAPL_P00,
    MAPL_RGAS,
    MAPL_RVAP,
)
from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.interpolations import vertical_interpolation
from pyMoist.saturation_tables import (
    GlobalTable_saturation_tables,
    SaturationVaporPressureTable,
    saturation_specific_humidity,
)
from pyMoist.GFDL_1M.state import GFDL1MState
from pyMoist.GFDL_1M.locals import GFDL1MLocals
import dataclasses


def calculate_derived_states(
    p_interface: FloatField,
    p_interface_mb: FloatField,
    p_mb: FloatField,
    geopotential_height_interface: FloatField,
    edge_height_above_surface: FloatField,
    layer_height_above_surface: FloatField,
    layer_thickness: FloatField,
    layer_thinkness_negative: FloatField,
    dp: FloatField,
    mass: FloatField,
    mass_inverse: FloatField,
    t: FloatField,
    ese: GlobalTable_saturation_tables,
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

    This stencil MUST be built using Z_INTERFACE_DIM to function properly.
    """
    from __externals__ import k_end

    with computation(PARALLEL), interval(...):
        p_interface_mb = p_interface * 0.01
        edge_height_above_surface = geopotential_height_interface - geopotential_height_interface.at(K=k_end)
    with computation(PARALLEL), interval(0, -1):
        p_mb = 0.5 * (p_interface_mb + p_interface_mb[0, 0, 1])
        layer_height_above_surface = 0.5 * (edge_height_above_surface + edge_height_above_surface[0, 0, 1])
        layer_thickness = edge_height_above_surface - edge_height_above_surface[0, 0, 1]
        layer_thinkness_negative = -1.0 * layer_thickness
        dp = p_interface[0, 0, 1] - p_interface
        mass = dp / MAPL_GRAV
        mass_inverse = 1 / mass
        sat, dsat = saturation_specific_humidity(t=t, p=p_mb * 100, ese=ese, esx=esx)
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
        t: temperature at surface (K)
        rh: relative humidity at surface

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
    ese: GlobalTable_saturation_tables,
    esx: GlobalTable_saturation_tables,
    lcl_level: IntFieldIJ,
):
    """
    Find the level of the lifted condensation level (LCL).

    Arguments:
        t (in): Atmospheric temperature (K)
        p_mb (in): pressure (mb)
        vapor (in): water vapor mixing radio (kg/kg)
        ese (in): saturation vapor pressure table, details unknown
        esx (in): saturation vapor pressure table, details unknown
        k_lcl (out): LCL level
    """
    from __externals__ import k_end

    # set up mask to stop computation
    with computation(FORWARD), interval(0, 1):
        found_level: BoolFieldIJ = False

    # get LCL pressure
    with computation(PARALLEL), interval(-1, None):
        qsat, _ = saturation_specific_humidity(t=t, p=p_mb * 100, ese=ese, esx=esx)
        rhsfc = 100 * vapor / qsat
        tlcl = find_t_lcl(t=t, rh=rhsfc)
        rm = (1 - vapor) * MAPL_RGAS + vapor * MAPL_RVAP
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
    ese: GlobalTable_saturation_tables,
    esx: GlobalTable_saturation_tables,
    lower_tropospheric_stability: FloatFieldIJ,
    estimated_inversion_strength: FloatFieldIJ,
):
    """
    Find estimated inversion strength. Returns Estimated Inversion
    Strength (K) according to Wood and Betherton, J.Climate, 2006.
    Based on Fortran code written by Donifan Barahona.
    """
    with computation(FORWARD), interval(-1, None):
        lower_tropospheric_stability = th700 - th
        lcl_height = layer_height_above_surface.at(K=lcl_level - 1)

        # Simplified single adiabat eq4 of https://doi.org/10.1175/JCLI3988.1
        t850 = 0.5 * (t + t700)
        qs850, _ = saturation_specific_humidity(t=t850, p=100 * 850, ese=ese, esx=esx)
        gamma850 = (1.0 + (MAPL_ALHL * qs850 / (MAPL_RGAS * t850))) / (
            1.0 + (MAPL_ALHL * MAPL_ALHL * qs850 / (MAPL_CP * MAPL_RVAP * t850 * t850))
        )
        gamma850 = MAPL_GRAV / MAPL_CP * (1.0 - gamma850)
        estimated_inversion_strength = lower_tropospheric_stability - gamma850 * (z700 - lcl_height)


def update_precipitaiton(
    mixing_ratio: FloatField,
    shallow_convection_values: FloatField,
):
    from __externals__ import DT_MOIST

    with computation(PARALLEL), interval(...):
        mixing_ratio = mixing_ratio + shallow_convection_values * DT_MOIST


@dataclasses.dataclass
class GFDL1MSetupLocals(State):
    th: Quantity = dataclasses.field(
        metadata={
            "name": "th",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    t700: Quantity = dataclasses.field(
        metadata={
            "name": "t700",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    th700: Quantity = dataclasses.field(
        metadata={
            "name": "th700",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    z700: Quantity = dataclasses.field(
        metadata={
            "name": "z700",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )


class Setup(NDSLRuntime):
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
        prepare_tendencies,
    ):
        # make configuration and saturation tables visible at runtime
        self.config = config
        self.saturation_tables = saturation_tables

        # make the pre-built stencil visible at runtime
        self._prepare_tendencies = prepare_tendencies

        # initalize locals
        self._locals = GFDL1MSetupLocals.zeros(quantity_factory)

        # NOTE disabled orchestration because the data map has changed
        # orchestrate(
        #     obj=self,
        #     config=stencil_factory.config.dace_config,
        #     dace_compiletime_args=[
        #         "mixing_ratios",
        #         "cloud_fractions",
        #         "masks",
        #         "outputs",
        #         "temporaries",
        #     ],
        # )

        # construct stencils
        self._calculate_derived_states = stencil_factory.from_dims_halo(
            func=calculate_derived_states,
            compute_dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM],
        )

        self._find_lcl_level = stencil_factory.from_dims_halo(
            func=find_lcl_level,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._update_lcl_height = stencil_factory.from_dims_halo(
            func=update_lcl_height,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._vertical_interpolation = stencil_factory.from_dims_halo(
            func=vertical_interpolation,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._compute_estimated_inversion_strength = stencil_factory.from_dims_halo(
            func=compute_estimated_inversion_strength,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._update_precipitaiton = stencil_factory.from_dims_halo(
            func=update_precipitaiton,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "DT_MOIST": config.DT_MOIST,
            },
        )

        # Dev NOTE: this is an orchestration workaround. Direct call to
        #           `self.saturation_tables.X` fails closure capture for
        #           argument reconstruction at call time
        self._ese = self.saturation_tables.ese
        self._esx = self.saturation_tables.esx

    def __call__(
        self,
        state: GFDL1MState,
        locals: GFDL1MLocals,
    ):
        # prepare macrophysics tendencies
        self._prepare_tendencies(
            u=state.u,
            v=state.v,
            t=state.t,
            vapor=state.mixing_ratio.vapor,
            rain=state.mixing_ratio.rain,
            snow=state.mixing_ratio.snow,
            graupel=state.mixing_ratio.graupel,
            convective_liquid=state.mixing_ratio.convective_liquid,
            convective_ice=state.mixing_ratio.convective_ice,
            large_scale_liquid=state.mixing_ratio.large_scale_liquid,
            large_scale_ice=state.mixing_ratio.large_scale_ice,
            convective_cloud_fraction=state.cloud_fraction.convective,
            large_scale_cloud_fraction=state.cloud_fraction.large_scale,
            du_dt=state.tendencies.dudt_macro,
            dv_dt=state.tendencies.dvdt_macro,
            dt_dt=state.tendencies.dtdt_macro,
            dvapor_dt=state.tendencies.dvapordt_macro,
            dliquid_dt=state.tendencies.dliquiddt_macro,
            dice_dt=state.tendencies.dicedt_macro,
            dcloud_fraction_dt=state.tendencies.dcloud_fractiondt_macro,
            drain_dt=state.tendencies.draindt_macro,
            dsnow_dt=state.tendencies.dsnowdt_macro,
            dgraupel_dt=state.tendencies.dgraupeldt_macro,
        )

        self._calculate_derived_states(
            p_interface=state.p_interface,
            p_interface_mb=locals.p_interface_mb,
            p_mb=locals.p_mb,
            geopotential_height_interface=state.z_interface,
            edge_height_above_surface=locals.edge_height_above_surface,
            layer_height_above_surface=locals.layer_height_above_surface,
            layer_thickness=locals.layer_thickness,
            layer_thinkness_negative=locals.layer_thickness_negative,
            dp=locals.dp,
            mass=locals.mass,
            mass_inverse=locals.mass_inverse,
            t=state.t,
            ese=self._ese,
            esx=self.saturation_tables.esx,
            sat=locals.saturation_specific_humidity,
            dsat=locals.dsaturation_specific_humidity,
            u=state.u,
            u_unmodified=locals.u_unmodified,
            v=state.v,
            v_unmodified=locals.v_unmodified,
            th=self._locals.th,
        )

        self._find_lcl_level(
            t=state.t,
            p_mb=locals.p_mb,
            vapor=state.mixing_ratio.vapor,
            ese=self._ese,
            esx=self._esx,
            lcl_level=locals.lcl_level,
        )

        # NOTE need a new way to resolve this now that it is a state field. it will never be none
        if state.lcl_height is not None:
            self._update_lcl_height(
                layer_height_above_surface=locals.layer_height_above_surface,
                lcl_level=locals.lcl_level,
                lcl_height=state.lcl_height,
            )

        self._vertical_interpolation(
            field=self._locals.th,
            interpolated_field=self._locals.th700,
            p_interface_mb=locals.p_interface_mb,
            target_pressure=Float(70000.0),
        )

        self._vertical_interpolation(
            field=state.t,
            interpolated_field=self._locals.t700,
            p_interface_mb=locals.p_interface_mb,
            target_pressure=Float(70000.0),
        )

        self._vertical_interpolation(
            field=locals.layer_height_above_surface,
            interpolated_field=self._locals.z700,
            p_interface_mb=locals.p_interface_mb,
            target_pressure=Float(70000.0),
        )

        self._compute_estimated_inversion_strength(
            t=state.t,
            th=self._locals.th,
            layer_height_above_surface=locals.layer_height_above_surface,
            t700=self._locals.t700,
            th700=self._locals.th700,
            z700=self._locals.z700,
            lcl_level=locals.lcl_level,
            ese=self._ese,
            esx=self._esx,
            lower_tropospheric_stability=state.lower_tropospheric_stability,
            estimated_inversion_strength=state.estimated_inversion_strength,
        )

        # NOTE need a new way to resolve this now that it is a state field. it will never be none
        if state.shallow_convection_rain is not None:
            self._update_precipitaiton(state.mixing_ratio.rain, state.shallow_convection_rain)

        # NOTE need a new way to resolve this now that it is a state field. it will never be none
        if state.shallow_convection_snow is not None:
            self._update_precipitaiton(state.mixing_ratio.snow, state.shallow_convection_snow)
