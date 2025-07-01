from gt4py.cartesian.gtscript import THIS_K

from ndsl import StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.gt4py import BACKWARD, FORWARD, PARALLEL, computation, function, interval, log
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
from pyMoist.field_types import GlobalTable_saturaion_tables
from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.GFDL_1M.masks import Masks
from pyMoist.GFDL_1M.outputs import Outputs
from pyMoist.GFDL_1M.state import CloudFractions, MixingRatios
from pyMoist.GFDL_1M.temporaries import Temporaries
from pyMoist.interpolations import vertical_interpolation
from pyMoist.saturation_tables.qsat_functions import saturation_specific_humidity
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable


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
    t: FloatField,
    ese: GlobalTable_saturaion_tables,
    esx: GlobalTable_saturaion_tables,
    qsat: FloatField,
    dqsat: FloatField,
    u: FloatField,
    u_unmodified: FloatField,
    v: FloatField,
    v_unmodified: FloatField,
    temporary_3d: FloatField,
    th: FloatField,
):
    """
    Computes derived state fields required for the rest of the GFDL single moment
    microphysics module.

    This stencil MUST be built using Z_INTERFACE_DIM to funciton properly.
    """
    from __externals__ import k_end

    with computation(PARALLEL), interval(...):
        p_interface_mb = p_interface * 0.01
        edge_height_above_surface = geopotential_height_interface - geopotential_height_interface.at(K=k_end)
    with computation(FORWARD), interval(0, -1):
        p_mb = 0.5 * (p_interface_mb + p_interface_mb[0, 0, 1])
        layer_height_above_surface = 0.5 * (edge_height_above_surface + edge_height_above_surface[0, 0, 1])
        layer_thickness = edge_height_above_surface - edge_height_above_surface[0, 0, 1]
        layer_thinkness_negative = -1.0 * layer_thickness
        dp = p_interface[0, 0, 1] - p_interface
        mass = dp / MAPL_GRAV
        qsat, dqsat = saturation_specific_humidity(t=t, p=p_mb * 100, ese=ese, esx=esx)
        u_unmodified = u
        v_unmodified = v
        temporary_3d = (100.0 * p_mb / MAPL_P00) ** (MAPL_KAPPA)
        th = t / temporary_3d

    # reset temporary field for later uses
    with computation(PARALLEL), interval(...):
        temporary_3d = 0


@function
def find_tlcl(
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


def find_k_lcl(
    t: FloatField,
    p_mb: FloatField,
    vapor: FloatField,
    ese: GlobalTable_saturaion_tables,
    esx: GlobalTable_saturaion_tables,
    found_level: BoolFieldIJ,
    k_lcl: IntFieldIJ,
):
    """
    Find the level of the lifted condensation level (LCL).

    Arguments:
        t (in): Atmospheric temperature (K)
        p_mb (in): pressure (mb)
        vapor (in): water vapor mixing radio (kg/kg)
        ese (in): saturation vapor pressure table, details unknown
        esx (in): saturation vapor pressure table, details unknown
        found_level (in): boolean mask used to stop calculation after LCL is identified
        k_lcl (out): LCL level
    """
    from __externals__ import k_end

    # get LCL pressure
    with computation(FORWARD), interval(-1, None):
        qsat, _ = saturation_specific_humidity(t=t, p=p_mb * 100, ese=ese, esx=esx)
        rhsfc = 100 * vapor / qsat
        tlcl = find_tlcl(t=t, rh=rhsfc)
        rm = (1 - vapor) * MAPL_RGAS + vapor * MAPL_RVAP
        cpm = (1.0 - vapor) * MAPL_CPDRY + vapor * MAPL_CPVAP
        plcl = p_mb * ((tlcl / t) ** (cpm / rm))

    # find nearest level <= LCL pressure
    with computation(BACKWARD), interval(...):
        if found_level == False:  # noqa
            k_lcl = THIS_K
        if p_mb <= plcl.at(K=k_end):
            found_level = True

    # Reset mask for future use
    with computation(FORWARD), interval(0, 1):
        found_level = False


def update_z_lcl(
    layer_height_above_surface: FloatField,
    k_lcl: IntFieldIJ,
    z_lcl: FloatFieldIJ,
):
    with computation(FORWARD), interval(0, 1):
        z_lcl = layer_height_above_surface.at(K=k_lcl)


def find_eis(
    t: FloatField,
    th: FloatField,
    layer_height_above_surface: FloatField,
    t700: FloatFieldIJ,
    th700: FloatFieldIJ,
    z700: FloatFieldIJ,
    k_lcl: IntFieldIJ,
    ese: GlobalTable_saturaion_tables,
    esx: GlobalTable_saturaion_tables,
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
        zlcl = layer_height_above_surface.at(K=k_lcl - 1)

        # Simplified single adiabat eq4 of https://doi.org/10.1175/JCLI3988.1
        t850 = 0.5 * (t + t700)
        qs850, _ = saturation_specific_humidity(t=t850, p=100 * 850, ese=ese, esx=esx)
        gamma850 = (1.0 + (MAPL_ALHL * qs850 / (MAPL_RGAS * t850))) / (
            1.0 + (MAPL_ALHL * MAPL_ALHL * qs850 / (MAPL_CP * MAPL_RVAP * t850 * t850))
        )
        gamma850 = MAPL_GRAV / MAPL_CP * (1.0 - gamma850)
        estimated_inversion_strength = lower_tropospheric_stability - gamma850 * (z700 - zlcl)


def update_precipitaiton(
    mixing_ratio: FloatField,
    shallow_convection_values: FloatField,
):
    from __externals__ import DT_MOIST

    with computation(PARALLEL), interval(...):
        mixing_ratio = mixing_ratio + shallow_convection_values * DT_MOIST


class Setup:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        GFDL_1M_config: GFDL1MConfig,
        saturation_tables: SaturationVaporPressureTable,
        prepare_tendencies,
    ):
        self.GFDL_1M_config = GFDL_1M_config
        self.saturation_tables = saturation_tables
        self.prepare_tendencies = prepare_tendencies

        # construct stencils
        self.calculate_derived_states = stencil_factory.from_dims_halo(
            func=calculate_derived_states,
            compute_dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM],
        )

        self.find_k_lcl = stencil_factory.from_dims_halo(
            func=find_k_lcl,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self.update_z_lcl = stencil_factory.from_dims_halo(
            func=update_z_lcl,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self.vertical_interpolation = stencil_factory.from_dims_halo(
            func=vertical_interpolation,
            compute_dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM],
        )

        self.find_eis = stencil_factory.from_dims_halo(
            func=find_eis,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self.update_precipitaiton = stencil_factory.from_dims_halo(
            func=update_precipitaiton,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "DT_MOIST": GFDL_1M_config.DT_MOIST,
            },
        )

    def __call__(
        self,
        geopotential_height_interface: FloatField,
        p_interface: FloatField,
        t: FloatField,
        u: FloatField,
        v: FloatField,
        shallow_convective_rain: FloatField,
        shallow_convective_snow: FloatField,
        mixing_ratios: MixingRatios,
        cloud_fractions: CloudFractions,
        masks: Masks,
        outputs: Outputs,
        temporaries: Temporaries,
    ):

        # prepare macrophysics tendencies
        self.prepare_tendencies(
            u=u,
            v=v,
            t=t,
            vapor=mixing_ratios.vapor,
            rain=mixing_ratios.rain,
            snow=mixing_ratios.snow,
            graupel=mixing_ratios.graupel,
            convective_liquid=mixing_ratios.convective_liquid,
            convective_ice=mixing_ratios.convective_ice,
            large_scale_liquid=mixing_ratios.large_scale_liquid,
            large_scale_ice=mixing_ratios.large_scale_ice,
            convective_cloud_fraction=cloud_fractions.convective,
            large_scale_cloud_fraction=cloud_fractions.large_scale,
            du_dt=outputs.du_dt_macro,
            dv_dt=outputs.dv_dt_macro,
            dt_dt=outputs.dt_dt_macro,
            dvapor_dt=outputs.dvapor_dt_macro,
            dliquid_dt=outputs.dliquid_dt_macro,
            dice_dt=outputs.dice_dt_macro,
            dcloud_fraction_dt=outputs.dcloud_fraction_dt_macro,
            drain_dt=outputs.drain_dt_macro,
            dsnow_dt=outputs.dsnow_dt_macro,
            dgraupel_dt=outputs.dgraupel_dt_macro,
        )

        self.calculate_derived_states(
            p_interface=p_interface,
            p_interface_mb=temporaries.p_interface_mb,
            p_mb=temporaries.p_mb,
            geopotential_height_interface=geopotential_height_interface,
            edge_height_above_surface=temporaries.edge_height_above_surface,
            layer_height_above_surface=temporaries.layer_height_above_surface,
            layer_thickness=temporaries.layer_thickness,
            layer_thinkness_negative=temporaries.layer_thickness_negative,
            dp=temporaries.dp,
            mass=temporaries.mass,
            t=t,
            ese=self.saturation_tables.ese,
            esx=self.saturation_tables.esx,
            qsat=temporaries.qsat,
            dqsat=temporaries.dqsat,
            u=u,
            u_unmodified=temporaries.u_unmodified,
            v=v,
            v_unmodified=temporaries.v_unmodified,
            temporary_3d=temporaries.temporary_3d,
            th=temporaries.th,
        )

        self.find_k_lcl(
            t=t,
            p_mb=temporaries.p_mb,
            vapor=mixing_ratios.vapor,
            ese=self.saturation_tables.ese,
            esx=self.saturation_tables.esx,
            found_level=masks.boolean_2d_mask,
            k_lcl=temporaries.k_lcl,
        )

        if outputs.z_lcl is not None:
            self.update_z_lcl(
                temporaries.layer_height_above_surface,
                temporaries.k_lcl,
                outputs.z_lcl,
            )

        self.vertical_interpolation(
            field=temporaries.th,
            interpolated_field=temporaries.th700,
            p_interface_mb=temporaries.p_interface_mb,
            target_pressure=Float(70000.0),
            pb=temporaries.temporary_2d_1,
            pt=temporaries.temporary_2d_2,
            boolean_2d_mask=masks.boolean_2d_mask,
        )

        self.vertical_interpolation(
            field=t,
            interpolated_field=temporaries.t700,
            p_interface_mb=temporaries.p_interface_mb,
            target_pressure=Float(70000.0),
            pb=temporaries.temporary_2d_1,
            pt=temporaries.temporary_2d_2,
            boolean_2d_mask=masks.boolean_2d_mask,
        )
        self.vertical_interpolation(
            field=temporaries.layer_height_above_surface,
            interpolated_field=temporaries.z700,
            p_interface_mb=temporaries.p_interface_mb,
            target_pressure=Float(70000.0),
            pb=temporaries.temporary_2d_1,
            pt=temporaries.temporary_2d_2,
            boolean_2d_mask=masks.boolean_2d_mask,
        )

        self.find_eis(
            t=t,
            th=temporaries.th,
            layer_height_above_surface=temporaries.layer_height_above_surface,
            t700=temporaries.t700,
            th700=temporaries.th700,
            z700=temporaries.z700,
            k_lcl=temporaries.k_lcl,
            ese=self.saturation_tables.ese,
            esx=self.saturation_tables.esx,
            lower_tropospheric_stability=outputs.lower_tropospheric_stability,
            estimated_inversion_strength=outputs.estimated_inversion_strength,
        )

        if shallow_convective_rain is not None:
            self.update_precipitaiton(mixing_ratios.rain, shallow_convective_rain)

        if shallow_convective_snow is not None:
            self.update_precipitaiton(mixing_ratios.snow, shallow_convective_snow)
