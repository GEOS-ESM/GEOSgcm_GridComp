from ndsl.dsl.gt4py import PARALLEL, computation, exp, float64, function, interval, sqrt
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ, Int

import pyMoist.constants as constants
from pyMoist.saturation_tables import (
    GlobalTable_saturation_tables,
    saturation_specific_humidity,
    saturation_specific_humidity_frozen_surface,
    saturation_specific_humidity_liquid_surface,
)
from pyMoist.shared.incloud_processes import ice_fraction


@function
def pdffrac(
    pdfshape: Int,
    qtmean: Float,
    sigmaqt1: Float,
    sigmaqt2: Float,
    qstar: Float,
):
    """Determine cloud fraction

    Args:
        pdfshape (Int)
        qtmean (Float)
        sigmaqt1 (Float)
        sigmaqt2 (Float)
        qstar (Float)

    Returns:
        Float: cloud_fraction
    """
    if pdfshape == 1:
        if (qtmean + sigmaqt1) < qstar:
            cloud_fraction = 0.0
        else:
            if sigmaqt1 > 0.0:
                cloud_fraction = min(qtmean + sigmaqt1 - qstar, 2.0 * sigmaqt1) / (2.0 * sigmaqt1)
            else:
                cloud_fraction = 1.0

    # Above code only executes when pdfshape = 1. Fortran code exists for pdfshape = 2

    return cloud_fraction


@function
def pdfcondensate(
    pdfshape: Int,
    qtmean: Float,
    sigmaqt1: Float,
    sigmaqt2: Float,
    qstar: Float,
):
    """Quantify new condensate

    Args:
        pdfshape (Int)
        qtmean (Float)
        sigmaqt1 (Float)
        sigmaqt2 (Float)
        qstar (Float)

    Returns:
        Float: condensate
    """
    qtmean_64: float64 = qtmean
    sigmaqt1_64: float64 = sigmaqt1
    sigmaqt2_64: float64 = sigmaqt2
    qstar_64: float64 = qstar

    condensate: float64 = 0.0
    if pdfshape == 1:
        if (qtmean_64 + sigmaqt1_64) < qstar_64:
            condensate = float64(0.0)
        elif qstar_64 > (qtmean_64 - sigmaqt1_64):
            if sigmaqt1_64 > 0.0:
                condensate = (min(qtmean_64 + sigmaqt1_64 - qstar_64, float64(2.0) * sigmaqt1_64) ** 2) / (float64(4.0) * sigmaqt1_64)
            else:
                condensate = qtmean_64 - qstar_64
        else:
            condensate = qtmean_64 - qstar_64

    # Above code only executes when pdfshape = 1. Fortran code exists for pdfshape = 2

    return condensate


@function
def bergeron_partition(
    DT_MOIST: Float,
    p_mb: Float,
    t: Float,
    mixing_ratio_vapor: Float,
    mixing_ratio_large_scale_ice: Float,
    mixing_ratio_convective_ice: Float,
    mixing_ratio_large_scale_liquid: Float,
    mixing_ratio_convective_liquid: Float,
    concentration_ice: Float,
    dq_all: Float,
    convection_fraction: Float,
    surface_type: Float,
    ese: GlobalTable_saturation_tables,
    esw: GlobalTable_saturation_tables,
    frz: Float,
    lqu: Float,
):
    """Partition the new condensates

    Args:
        DT_MOIST (Float)
        p_mb (Float)
        t (Float)
        mixing_ratio_vapor (Float)
        mixing_ratio_large_scale_ice (Float)
        mixing_ratio_convective_ice (Float)
        mixing_ratio_large_scale_liquid (Float)
        mixing_ratio_convective_liquid (Float)
        concentration_ice (Float)
        dq_all (Float)
        convection_fraction (Float)
        surface_type (Float)
        ese (GlobalTable_saturation_tables)
        esw (GlobalTable_saturation_tables)
        frz (Float)
        lqu (Float)

    Returns:
        Float: fraction_ice
        Float: dq_all
    """
    internal_mixing_ratio_ice = mixing_ratio_large_scale_ice + mixing_ratio_convective_ice  # necessary because NI is for convective and large scale
    internal_mixing_ratio_liquid = mixing_ratio_large_scale_liquid + mixing_ratio_convective_liquid
    internal_total_precipitate = internal_mixing_ratio_ice + internal_mixing_ratio_liquid
    if internal_total_precipitate > 0.0:
        f_qa = (mixing_ratio_convective_ice + mixing_ratio_large_scale_ice) / internal_total_precipitate
    else:
        f_qa = 0.0
    ni_x = (1.0 - f_qa) * concentration_ice

    dq_all = dq_all / DT_MOIST
    t_c = t - constants.MAPL_TICE

    # Completely glaciated cloud:
    if t >= constants.iT_ICE_MAX:  # liquid cloud
        fraction_ice = 0.0
    elif t <= constants.iT_ICE_ALL:  # ice cloud
        fraction_ice = 1.0
    else:  # mixed phase cloud
        fraction_ice = 0.0
        if mixing_ratio_large_scale_ice <= 0.0:
            fraction_ice = ice_fraction(t, convection_fraction, surface_type)
        else:
            qv_inc = mixing_ratio_vapor
            qs_liq, _ = saturation_specific_humidity_liquid_surface(esw, lqu, t, p_mb * 100.0)
            qs_ice, dq_si = saturation_specific_humidity_frozen_surface(ese, frz, t, p_mb * 100.0)
            qv_inc = min(qv_inc, qs_liq)  # limit to below water saturation

            # Calculate deposition onto preexisting ice

            diff = (0.211 * 1013.25 / (p_mb + 0.1)) * (((t + 0.1) / constants.MAPL_TICE) ** 1.94) * 1e-4  # From Seinfeld and Pandis 2006
            den_air = p_mb * 100.0 / constants.MAPL_RDRY / t
            den_ice = 1000.0 * (0.9167 - 1.75e-4 * t_c - 5.0e-7 * t_c * t_c)  # From PK 97
            # lh_corr = 1.0 + dq_si * constants.MAPL_LATENT_HEAT_SUBLIMATION / (
            #     constants.MAPL_RDRY / constants.MAPL_KAPPA
            # )  # must be ice deposition
            lh_corr = 1.0 + dq_si * constants.MAPL_ALHS / constants.MAPL_CP  # must be ice deposition

            if ni_x > 1.0 and mixing_ratio_large_scale_ice > 1.0e-10:
                dc = max(
                    (mixing_ratio_large_scale_ice / (ni_x * den_ice * constants.MAPL_PI)) ** 0.333,
                    20.0e-6,
                )  # Assume monodisperse size distribution
            else:
                dc = 20.0e-6

            teff = ni_x * den_air * 2.0 * constants.MAPL_PI * diff * dc / lh_corr  # 1/Dep time scale

            dep = 0.0
            if teff > 0.0 and mixing_ratio_large_scale_ice > 1.0e-14:
                aux = max(min(DT_MOIST * teff, 20.0), 0.0)
                dep = (qv_inc - qs_ice) * (1.0 - exp(-aux)) / DT_MOIST
            dep = max(dep, -mixing_ratio_large_scale_ice / DT_MOIST)  # only existing ice can be sublimated

            d_mixing_ratio_ice = 0.0
            d_mixing_ratio_liquid = 0.0
            fraction_ice = 0.0
            # Partition DQALL accounting for Bergeron-Findensen process
            if dq_all >= 0:  # net condensation. Note: do not allow bergeron with QLCN
                if dep > 0.0:
                    d_mixing_ratio_ice = min(dep, dq_all + mixing_ratio_large_scale_liquid / DT_MOIST)
                    d_mixing_ratio_liquid = dq_all - d_mixing_ratio_ice
                else:
                    d_mixing_ratio_liquid = dq_all  # could happen because the PDF allows
                    # condensation in subsaturated conditions
                    d_mixing_ratio_ice = 0.0
            if dq_all < 0.0:  # net evaporation. Water evaporates first regardless of DEP
                d_mixing_ratio_liquid = max(dq_all, -mixing_ratio_large_scale_liquid / DT_MOIST)
                d_mixing_ratio_ice = max(dq_all - d_mixing_ratio_liquid, -mixing_ratio_large_scale_ice / DT_MOIST)
            if dq_all != 0.0:
                fraction_ice = max(min(d_mixing_ratio_ice / dq_all, 1.0), 0.0)

    return fraction_ice, dq_all


def hydrostatic_pdf(
    alpha: FloatField,
    convection_fraction: FloatFieldIJ,
    surface_type: FloatFieldIJ,
    p_mb: FloatField,
    mixing_ratio_vapor: FloatField,
    mixing_ratio_large_scale_liquid: FloatField,
    mixing_ratio_convective_liquid: FloatField,
    mixing_ratio_large_scale_ice: FloatField,
    mixing_ratio_convective_ice: FloatField,
    t: FloatField,
    large_scale_cloud_fraction: FloatField,
    convective_cloud_fraction: FloatField,
    ice_concentration: FloatField,
    relative_humidity: FloatField,
    pdf_iters: FloatField,
    ese: GlobalTable_saturation_tables,
    esw: GlobalTable_saturation_tables,
    esx: GlobalTable_saturation_tables,
    estfrz: Float,
    estlqu: Float,
):
    """
    Hydrostatic PDF for use in the phase change section of GFDL_1M microphysics

    Arguments:
        alpha (FloatField): (in) details unknown
        convection_fraction (FloatFieldIJ): (in) grid cell convection fraction
        surface_type (FloatFieldIJ): (in) surface type
        p_mb (FloatField): (in) pressure (millibars)
        mixing_ratio_vapor (FloatField): (inout) water vapor mixing ratio (kg/kg)
        mixing_ratio_large_scale_liquid (FloatField): (inout) large scale liquid mass fraction (unitless)
        mixing_ratio_convective_liquid (FloatField): (inout) convective liquid mass fraction (unitless)
        mixing_ratio_large_scale_ice (FloatField): (inout) large scale ice mass fraction (unitless)
        mixing_ratio_convective_ice (FloatField): (inout) convective ice mass fraction (unitless)
        t (FloatField): (inout) temperature
        large_scale_cloud_fraction (FloatField): (inout) large scale cloud fraction
        convective_cloud_fraction (FloatField): (inout) convective cloud fraction
        ice_concentration (FloatField): (in) ice concentration
        relative_humidity (FloatField): (inout) relative humidity after pdf
        pdf_iters (FloatField): (out) number of iterations executed
        ese (GlobalTable_saturation_tables): saturation table
        esw (GlobalTable_saturation_tables): saturation table
        esx (GlobalTable_saturation_tables): saturation table
        estfrz (Float): saturation table value at -40 C
        estlqu (Float): saturation table value at 0 C
    """
    from __externals__ import DT_MOIST, FLOAT_TINY, PDF_SHAPE, USE_BERGERON

    with computation(PARALLEL), interval(...):
        if convective_cloud_fraction < 1.0:
            inv_clcn = 1.0 / (1.0 - convective_cloud_fraction)
        else:
            inv_clcn = 0.0
        if convective_cloud_fraction > FLOAT_TINY:
            qa_x = (mixing_ratio_convective_liquid + mixing_ratio_convective_ice) / convective_cloud_fraction
        else:
            qa_x = 0.0
        cf_n = large_scale_cloud_fraction * inv_clcn
        qc_n = (mixing_ratio_large_scale_liquid + mixing_ratio_large_scale_ice) * inv_clcn
        qc_i = mixing_ratio_large_scale_ice * inv_clcn
        t_n = t
        qs_x, _ = saturation_specific_humidity(t=t, p=p_mb * 100.0, ese=ese, esx=esx)
        qv_n = (mixing_ratio_vapor - qs_x * convective_cloud_fraction) * inv_clcn

        qt = qc_n + qv_n  # Total LS water after microphysics

        count = 1
        while count <= 20:
            qv_p = qv_n
            qc_p = qc_n
            cf_p = cf_n
            t_p = t_n
            qs_n, dqs = saturation_specific_humidity(t=t_n, p=p_mb * 100.0, ese=ese, esx=esx)

            if PDF_SHAPE < 3:  # 1 = top-hat 2 = triangulat
                sigmaqt1 = alpha * qs_n
                sigmaqt2 = alpha * qs_n
            elif PDF_SHAPE == 4:  # lognormal (sigma is dimensionless)
                sigmaqt1 = max(alpha / sqrt(3.0), 0.001)

            if PDF_SHAPE < 5:
                cf_n = pdffrac(PDF_SHAPE, qt, sigmaqt1, sigmaqt2, qs_n)
                qc_n = pdfcondensate(PDF_SHAPE, qt, sigmaqt1, sigmaqt2, qs_n)

            if USE_BERGERON:
                dq_all = qc_n - qc_p
                Nfac = 100.0 * p_mb * constants.R_AIR / t_n  # density times conversion factor
                NIv = ice_concentration / Nfac
                f_qi, dq_all = bergeron_partition(
                    DT_MOIST,
                    p_mb,
                    t_n,
                    qt,
                    mixing_ratio_large_scale_ice,
                    mixing_ratio_convective_ice,
                    mixing_ratio_large_scale_liquid,
                    mixing_ratio_convective_liquid,
                    NIv,
                    dq_all,
                    convection_fraction,
                    surface_type,
                    ese,
                    esw,
                    estfrz,
                    estlqu,
                )
            else:
                f_qi = ice_fraction(t_n, convection_fraction, surface_type)

            latent_heat_factor = (1.0 - f_qi) * constants.ALHLBCP + f_qi * constants.ALHSBCP
            if PDF_SHAPE == 1:
                qc_n = qc_p + (qc_n - qc_p) / (1.0 - (cf_n * (alpha - 1.0) - (qc_n / qs_n)) * dqs * latent_heat_factor)
            elif PDF_SHAPE == 5:
                qc_n = qc_p + 0.5 * (qc_n - qc_p)

            if convective_cloud_fraction > 0:
                qa_o = qa_x
            else:
                qa_o = 0

            qv_n = qv_p - (qc_n - qc_p)
            t_n = (
                t_p
                + (1.0 - f_qi) * constants.ALHLBCP * ((qc_n - qc_p) * (1.0 - convective_cloud_fraction) + (qa_o - qa_x) * convective_cloud_fraction)
                + f_qi * constants.ALHSBCP * ((qc_n - qc_p) * (1.0 - convective_cloud_fraction) + (qa_o - qa_x) * convective_cloud_fraction)
            )

            # NOTE Differences in constants between fortran and python cause
            # calculation of t_n to be slightly off (at order 10^-5), causing
            # the following conditional to occasionsally trigger incorrectly,
            # leading to occasional errors in various fields at the end of the
            # PDF stencil.
            pdf_iters = count
            if abs(t_n - t_p) < 0.00001:
                count = 21  # break out of loop
            else:
                count = count + 1

        if convective_cloud_fraction < 1.0:
            large_scale_cloud_fraction = cf_n * (1.0 - convective_cloud_fraction)
            qc_n = qc_n * (1.0 - convective_cloud_fraction)
            qa_o = qa_o * convective_cloud_fraction
        else:
            # Special case CLCN=1, i.e., box filled with anvil.
            # - Note: no guarantee QV_box > QS_box
            large_scale_cloud_fraction = 0.0  # Remove any LS cloud
            qa_o = (
                mixing_ratio_convective_liquid + mixing_ratio_convective_ice + mixing_ratio_large_scale_liquid + mixing_ratio_large_scale_ice
            )  # Add all LS condensate to anvil type
            qc_n = 0.0  # Remove same from new LS
            qt = qa_o + mixing_ratio_vapor  # Update total water
            # Now set anvil condensate to any excess of total water
            # over QSx (saturation value at top)
            qa_o = max(qt - qs_x, 0.0)

        # Now take {\em New} condensate and partition into ice and liquid

        # large scale
        qc_x = qc_n - (mixing_ratio_large_scale_liquid + mixing_ratio_large_scale_ice)
        if qc_x < 0.0:  # net evaporation
            d_qlls = max(qc_x, -mixing_ratio_large_scale_liquid)  # Water evaporates first
            d_qils = max(qc_x - d_qlls, -mixing_ratio_large_scale_ice)  # Then sublimation
        else:
            d_qlls = (1.0 - f_qi) * qc_x
            d_qils = f_qi * qc_x

        # convective
        qa_x = qa_o - (mixing_ratio_convective_liquid + mixing_ratio_convective_ice)
        if qa_x < 0.0:  # net evaporation
            d_qlcn = max(qa_x, -mixing_ratio_convective_liquid)  # Water evaporates first
            d_qicn = max(qa_x - d_qlcn, -mixing_ratio_convective_ice)  # Then sublimation
        else:
            d_qlcn = (1.0 - f_qi) * qa_x
            d_qicn = f_qi * qa_x

        # Clean-up cloud if fractions are too small
        if convective_cloud_fraction < 1.0e-5:
            d_qicn = -mixing_ratio_convective_ice
            d_qlcn = -mixing_ratio_convective_liquid
        if large_scale_cloud_fraction < 1.0e-5:
            d_qils = -mixing_ratio_large_scale_ice
            d_qlls = -mixing_ratio_large_scale_liquid

        mixing_ratio_convective_ice = mixing_ratio_convective_ice + d_qicn
        mixing_ratio_convective_liquid = mixing_ratio_convective_liquid + d_qlcn
        mixing_ratio_large_scale_ice = mixing_ratio_large_scale_ice + d_qils
        mixing_ratio_large_scale_liquid = mixing_ratio_large_scale_liquid + d_qlls
        mixing_ratio_vapor = mixing_ratio_vapor - (d_qicn + d_qils + d_qlcn + d_qlls)
        t = (
            t
            + constants.MAPL_LATENT_HEAT_VAPORIZATION / constants.MAPL_CPDRY * (d_qicn + d_qils + d_qlcn + d_qlls)
            + constants.MAPL_LATENT_HEAT_FUSION / constants.MAPL_CPDRY * (d_qicn + d_qils)
        )

        # We need to take care of situations where QS moves past QA
        # during QSAT iteration. This should be only when QA/AF is small
        # to begin with. Effect is to make QAo negative. So, we
        # "evaporate" offending QA's
        # We get rid of anvil fraction also, although strictly
        # speaking, PDF-wise, we should not do this.

        if qa_o <= 0.0:
            mixing_ratio_vapor = mixing_ratio_vapor + mixing_ratio_convective_ice + mixing_ratio_convective_liquid
            t = (
                t
                - constants.MAPL_LATENT_HEAT_SUBLIMATION / constants.MAPL_CPDRY * mixing_ratio_convective_ice
                - constants.MAPL_LATENT_HEAT_VAPORIZATION / constants.MAPL_CPDRY * mixing_ratio_convective_liquid
            )
            mixing_ratio_convective_ice = 0.0
            mixing_ratio_convective_liquid = 0.0
            convective_cloud_fraction = 0.0

        denom, _ = saturation_specific_humidity(t=t, p=p_mb * 100.0, ese=ese, esx=esx)
        relative_humidity = mixing_ratio_vapor / denom
