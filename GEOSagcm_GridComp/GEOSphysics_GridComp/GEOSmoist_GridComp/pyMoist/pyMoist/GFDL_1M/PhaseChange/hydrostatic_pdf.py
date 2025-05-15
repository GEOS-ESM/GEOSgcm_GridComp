import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import (
    PARALLEL,
    computation,
    exp,
    interval,
    sqrt,
)

import pyMoist.constants as constants
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ, Int
from pyMoist.field_types import FloatField_VaporSaturationTable
from pyMoist.saturation_tables.qsat_functions import (
    QSat_Float,
    QSat_Float_Ice,
    QSat_Float_Liquid,
)
from pyMoist.shared_incloud_processes import ice_fraction


@gtscript.function
def pdffrac(
    pdfshape: Int,
    qtmean: Float,
    sigmaqt1: Float,
    sigmaqt2: Float,
    qstar: Float,
    clfrac: Float,
):
    if pdfshape == 1:
        if (qtmean + sigmaqt1) < qstar:
            clfrac = 0.0
        else:
            if sigmaqt1 > 0.0:
                clfrac = min((qtmean + sigmaqt1 - qstar), 2.0 * sigmaqt1) / (
                    2.0 * sigmaqt1
                )
            else:
                clfrac = 1.0

    # Above code only executes when pdfshape = 1. Fortran code exists for pdfshape = 2

    return clfrac


@gtscript.function
def pdfcondensate(
    pdfshape: Int,
    qtmean: Float,
    sigmaqt1: Float,
    sigmaqt2: Float,
    qstar: Float,
    condensate: Float,
):
    if pdfshape == 1:
        if (qtmean + sigmaqt1) < qstar:
            condensate = 0.0
        elif qstar > (qtmean - sigmaqt1):
            if sigmaqt1 > 0.0:
                condensate = (min(qtmean + sigmaqt1 - qstar, 2.0 * sigmaqt1) ** 2) / (
                    4.0 * sigmaqt1
                )
            else:
                condensate = qtmean - qstar
        else:
            condensate = qtmean - qstar

    # Above code only executes when pdfshape = 1. Fortran code exists for pdfshape = 2

    return condensate


@gtscript.function
def bergeron_partition(
    DTIME: Float,
    PL: Float,
    TE: Float,
    Q: Float,
    QILS: Float,
    QICN: Float,
    QLLS: Float,
    QLCN: Float,
    NI: Float,
    DQALL: Float,
    CNV_FRC: Float,
    SRF_TYPE: Float,
    ese: FloatField_VaporSaturationTable,
    esw: FloatField_VaporSaturationTable,
    estfrz: Float,
    estlqu: Float,
    count: Float,
):
    QI = QILS + QICN  # neccesary because NI is for convective and large scale
    QL = QLLS + QLCN
    QTOT = QI + QLLS
    if QTOT > 0.0:
        FQA = (QICN + QILS) / QTOT
    else:
        FQA = 0.0
    NIX = (1.0 - FQA) * NI

    DQALL = DQALL / DTIME
    TC = TE - constants.MAPL_TICE

    # Completelely glaciated cloud:
    if TE >= constants.iT_ICE_MAX:  # liquid cloud
        fQI = 0.0
    elif TE <= constants.iT_ICE_ALL:  # ice cloud
        fQI = 1.0
    else:  # mixed phase cloud
        fQI = 0.0
        if QILS <= 0.0:
            fQI = ice_fraction(TE, CNV_FRC, SRF_TYPE)
        else:
            QVINC = Q
            QSLIQ, _ = QSat_Float_Liquid(esw, estlqu, TE, PL * 100.0)
            QSICE, DQSI = QSat_Float_Ice(ese, estfrz, TE, PL * 100.0, compute_dq=True)
            QVINC = min(QVINC, QSLIQ)  # limit to below water saturation
            # Calculate deposition onto preexisting ice

            DIFF = (
                (0.211 * 1013.25 / (PL + 0.1))
                * (((TE + 0.1) / constants.MAPL_TICE) ** 1.94)
                * 1e-4
            )  # From Seinfeld and Pandis 2006
            DENAIR = PL * 100.0 / constants.MAPL_RDRY / TE
            DENICE = 1000.0 * (0.9167 - 1.75e-4 * TC - 5.0e-7 * TC * TC)  # From PK 97
            LHcorr = 1.0 + DQSI * constants.MAPL_LATENT_HEAT_SUBLIMATION / (
                constants.MAPL_RDRY / constants.MAPL_KAPPA
            )  # must be ice deposition

            if NIX > 1.0 and QILS > 1.0e-10:
                DC = max(
                    (QILS / (NIX * DENICE * constants.MAPL_PI)) ** 0.333, 20.0e-6
                )  # Assumme monodisperse size dsitribution
            else:
                DC = 20.0e-6

            TEFF = (
                NIX * DENAIR * 2.0 * constants.MAPL_PI * DIFF * DC / LHcorr
            )  # 1/Dep time scale

            DEP = 0.0
            if TEFF > 0.0 and QILS > 1.0e-14:
                AUX = max(min(DTIME * TEFF, 20.0), 0.0)
                DEP = (QVINC - QSICE) * (1.0 - exp(-AUX)) / DTIME
            DEP = max(DEP, -QILS / DTIME)  # only existing ice can be sublimated

            DQI = 0.0
            DQL = 0.0
            fQI = 0.0
            # Partition DQALL accounting for Bergeron-Findensen process
            if DQALL >= 0:  # net condensation. Note: do not allow bergeron with QLCN
                if DEP > 0.0:
                    DQI = min(DEP, DQALL + QLLS / DTIME)
                    DQL = DQALL - DQI
                else:
                    DQL = DQALL  # could happen because the PDF allows
                    # condensation in subsaturated conditions
                    DQI = 0.0
            if (
                DQALL < 0.0
            ):  # net evaporation. Water evaporates first regaardless of DEP
                DQL = max(DQALL, -QLLS / DTIME)
                DQI = max(DQALL - DQL, -QILS / DTIME)
            if DQALL != 0.0:
                fQI = max(min(DQI / DQALL, 1.0), 0.0)

    return fQI, DQALL


def hydrostatic_pdf(
    alpha: FloatField,
    cnv_frc: FloatFieldIJ,
    srf_type: FloatFieldIJ,
    p_mb: FloatField,
    q: FloatField,
    qlls: FloatField,
    qlcn: FloatField,
    qils: FloatField,
    qicn: FloatField,
    t: FloatField,
    clls: FloatField,
    clcn: FloatField,
    nacti: FloatField,
    rhx: FloatField,
    ese: FloatField_VaporSaturationTable,
    esw: FloatField_VaporSaturationTable,
    esx: FloatField_VaporSaturationTable,
    estfrz: Float,
    estlqu: Float,
):
    from __externals__ import DT_MOIST, PDF_SHAPE, USE_BERGERON

    # Reference Fortran: Process_Library.F90: subroutine hystpdf
    # with PDFSHAPE = 1, USE_BERGERON = True, and SC_ICE = False
    with computation(PARALLEL), interval(...):
        if clcn < 1.0:
            inv_clcn = 1.0 / (1.0 - clcn)
        else:
            inv_clcn = 0.0
        if (
            clcn > 0.0
        ):  # need to make sure this is equivilant to fortran CLCN > tiny(0.0)
            qa_x = (qlcn + qicn) / clcn
        else:
            qa_x = 0.0
        cf_n = clls * inv_clcn
        qc_n = (qlls + qils) * inv_clcn
        qc_i = qils * inv_clcn
        t_n = t
        qs_x, _ = QSat_Float(ese, esx, t, p_mb)
        qv_n = (q - qs_x * clcn) * inv_clcn

        qt = qc_n + qv_n  # Total LS water after microphysics

        count = 1
        while count <= 20:
            qv_p = qv_n
            qc_p = qc_n
            cf_p = cf_n
            t_p = t_n
            qs_n, dqs = QSat_Float(ese, esx, t_n, p_mb, DQSAT_trigger=True)

            if PDF_SHAPE < 3:  # 1 = top-hat 2 = triangulat
                sigmaqt1 = alpha * qs_n
                sigmaqt2 = alpha * qs_n
            elif PDF_SHAPE == 4:  # lognormal (sigma is dimensionless)
                sigmaqt1 = max(alpha / sqrt(3.0), 0.001)

            if PDF_SHAPE < 5:
                cf_n = pdffrac(PDF_SHAPE, qt, sigmaqt1, sigmaqt2, qs_n, cf_n)
                qc_n = pdfcondensate(PDF_SHAPE, qt, sigmaqt1, sigmaqt2, qs_n, qc_n)

            if USE_BERGERON:
                dq_all = qc_n - qc_p
                Nfac = (
                    100.0 * p_mb * constants.R_AIR / t_n
                )  # density times conversion factor
                NIv = nacti / Nfac
                f_qi, dq_all = bergeron_partition(
                    DT_MOIST,
                    p_mb,
                    t_n,
                    qt,
                    qils,
                    qicn,
                    qlls,
                    qlcn,
                    NIv,
                    dq_all,
                    cnv_frc,
                    srf_type,
                    ese,
                    esw,
                    estfrz,
                    estlqu,
                    count,
                )
            else:
                f_qi = ice_fraction(t_n, cnv_frc, srf_type)

            latent_heat_factor = (1.0 - f_qi) * (
                constants.MAPL_LATENT_HEAT_VAPORIZATION / constants.MAPL_CP
            ) + f_qi * (constants.MAPL_LATENT_HEAT_SUBLIMATION / constants.MAPL_CP)
            if PDF_SHAPE == 1:
                qc_n = qc_p + (qc_n - qc_p) / (
                    1.0
                    - (cf_n * (alpha - 1.0) - (qc_n / qs_n)) * dqs * latent_heat_factor
                )
            elif PDF_SHAPE == 5:
                qc_n = qc_p + 0.5 * (qc_n - qc_p)

            if clcn > 0:
                qa_o = qa_x
            else:
                qa_o = 0

            qv_n = qv_p - (qc_n - qc_p)
            t_n = (
                t_p
                + (1.0 - f_qi)
                * constants.MAPL_LATENT_HEAT_VAPORIZATION
                / constants.MAPL_CP
                * ((qc_n - qc_p) * (1.0 - clcn) + (qa_o - qa_x) * clcn)
                + f_qi
                * constants.MAPL_LATENT_HEAT_SUBLIMATION
                / constants.MAPL_CP
                * ((qc_n - qc_p) * (1.0 - clcn) + (qa_o - qa_x) * clcn)
            )

            # NOTE Differences in constants between fortran and python cause
            # calculation of t_n to be slightly off (at order 10^-5), causing
            # the following conditional to occasionsally trigger incorrectly,
            # leading to occasional errors in various fields at the end of the
            # PDF stencil.
            if abs(t_n - t_p) < 0.00001:
                count = 21  # break out of loop
            else:
                count = count + 1

        if clcn < 1.0:
            clls = cf_n * (1.0 - clcn)
            qc_n = qc_n * (1.0 - clcn)
            qa_o = qa_o * clcn
        else:
            # Special case CLCN=1, i.e., box filled with anvil.
            # - Note: no guarantee QV_box > QS_box
            clls = 0.0  # Remove any LS cloud
            qa_o = qlcn + qicn + qlls + qils  # Add all LS condensate to anvil type
            qc_n = 0.0  # Remove same from new LS
            qt = qa_o + q  # Update total water
            # Now set anvil condensate to any excess of total water
            # over QSx (saturation value at top)
            qa_o = max(qt - qs_x, 0.0)

        # Now take {\em New} condensate and partition into ice and liquid

        # large scale
        qc_xCx = qc_n - (qlls + qils)
        if qc_xCx < 0.0:  # net evaporation
            d_qlls = max(qc_xCx, -qlls)  # Water evaporates first
            d_qils = max(qc_xCx - d_qlls, -qils)  # Then sublimation
        else:
            d_qlls = (1.0 - f_qi) * qc_xCx
            d_qils = f_qi * qc_xCx

        # convective
        qa_x = qa_o - (qlcn + qicn)
        if qa_x < 0.0:  # net evaporation
            d_qlcn = max(qa_x, -qlcn)  # Water evaporates first
            d_qicn = max(qa_x - d_qlcn, -qicn)  # Then sublimation
        else:
            d_qlcn = (1.0 - f_qi) * qa_x
            d_qicn = f_qi * qa_x

        # Clean-up cloud if fractions are too small
        if clcn < 1.0e-5:
            d_qicn = -qicn
            d_qlcn = -qlcn
        if clls < 1.0e-5:
            d_qils = -qils
            d_qlls = -qlls

        qicn = qicn + d_qicn
        qlcn = qlcn + d_qlcn
        qils = qils + d_qils
        qlls = qlls + d_qlls
        q = q - (d_qicn + d_qils + d_qlcn + d_qlls)
        t = (
            t
            + constants.MAPL_LATENT_HEAT_VAPORIZATION
            / constants.MAPL_CPDRY
            * (d_qicn + d_qils + d_qlcn + d_qlls)
            + constants.MAPL_LATENT_HEAT_FUSION
            / constants.MAPL_CPDRY
            * (d_qicn + d_qils)
        )

        # We need to take care of situations where QS moves past QA
        # during QSAT iteration. This should be only when QA/AF is small
        # to begin with. Effect is to make QAo negative. So, we
        # "evaporate" offending QA's
        # We get rid of anvil fraction also, although strictly
        # speaking, PDF-wise, we should not do this.

        if qa_o <= 0.0:
            q = q + qicn + qlcn
            t = (
                t
                - constants.MAPL_LATENT_HEAT_SUBLIMATION / constants.MAPL_CPDRY * qicn
                - constants.MAPL_LATENT_HEAT_VAPORIZATION / constants.MAPL_CPDRY * qlcn
            )
            qicn = 0.0
            qlcn = 0.0
            clcn = 0.0

        denom, _ = QSat_Float(ese, esx, t, p_mb)
        rhx = q / denom
