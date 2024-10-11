"""Functions and stencils in this file make up the core of the evap/subl/pdf loop,
called before the GFDL_1M driver in the Fortran. Functions in this file are
unique to GFDL_1M and need not be visible to the rest pyMoist."""

import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import (
    PARALLEL,
    atan,
    computation,
    exp,
    interval,
    sqrt,
    tan,
)

import pyMoist.pyMoist_constants as constants
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ, Int
from pyMoist.field_types import FloatField_VaporSaturationTable
from pyMoist.saturation.qsat import QSat_Float, QSat_Float_Ice, QSat_Float_Liquid
from pyMoist.shared_incloud_processes import (
    cloud_effective_radius_ice,
    cloud_effective_radius_liquid,
    ice_fraction,
)


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
    TC = TE - constants.t_ice

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
            QSICE, DQSI = QSat_Float_Ice(ese, estfrz, TE, PL * 100.0, DQ_trigger=True)
            QVINC = min(QVINC, QSLIQ)  # limit to below water saturation
            # Calculate deposition onto preexisting ice

            DIFF = (
                (0.211 * 1013.25 / (PL + 0.1))
                * (((TE + 0.1) / constants.t_ice) ** 1.94)
                * 1e-4
            )  # From Seinfeld and Pandis 2006
            DENAIR = PL * 100.0 / constants.rdry / TE
            DENICE = 1000.0 * (0.9167 - 1.75e-4 * TC - 5.0e-7 * TC * TC)  # From PK 97
            LHcorr = 1.0 + DQSI * constants.latent_heat_sublimation / (
                constants.rdry / constants.kappa
            )  # must be ice deposition

            if NIX > 1.0 and QILS > 1.0e-10:
                DC = max(
                    (QILS / (NIX * DENICE * constants.PI)) ** 0.333, 20.0e-6
                )  # Assumme monodisperse size dsitribution
            else:
                DC = 20.0e-6

            TEFF = (
                NIX * DENAIR * 2.0 * constants.PI * DIFF * DC / LHcorr
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


def initial_calc(
    EIS: FloatFieldIJ,
    dw_land: Float,
    dw_ocean: Float,
    TURNRHCRIT_PARAM: Float,
    minrhcrit: FloatField,
    PLmb_at_klcl: FloatFieldIJ,
    PLmb: FloatField,
    PLEmb_top: FloatField,
    AREA: FloatFieldIJ,
    ALPHA: FloatField,
):
    with computation(PARALLEL), interval(...):
        # Send the condensates through the pdf after convection
        facEIS = max(0.0, min(1.0, EIS / 10.0)) ** 2
        # determine combined minrhcrit in stable/unstable regimes
        minrhcrit = (1.0 - dw_ocean) * (1.0 - facEIS) + (1.0 - dw_land) * facEIS
        # determine the turn pressure using the LCL
        if TURNRHCRIT_PARAM <= 0:
            turnrhcrit = (
                PLmb_at_klcl - 250
            )  # implementation of for loop needs to go here
        else:
            turnrhcrit = TURNRHCRIT_PARAM

    # lower turn from maxrhcrit=1.0 # implemented in multiple "with" statements
    # to deal with hybrid indexing
    with computation(PARALLEL), interval(0, -1):
        # Use Slingo-Ritter (1985) formulation for critical rel ative humidity
        RHCRIT = 1.0
        if PLmb <= turnrhcrit:
            RHCRIT = minrhcrit
        else:
            RHCRIT = minrhcrit + (1.0 - minrhcrit) / (19.0) * (
                (
                    atan(
                        (2.0 * (PLmb - turnrhcrit) / (PLEmb_top - turnrhcrit) - 1.0)
                        * tan(20.0 * constants.PI / 21.0 - 0.5 * constants.PI)
                    )
                    + 0.5 * constants.PI
                )
                * 21.0
                / constants.PI
                - 1.0
            )
    with computation(PARALLEL), interval(-1, None):
        # lower turn from maxrhcrit=1.0
        if PLmb <= turnrhcrit:
            RHCRIT = minrhcrit
        else:
            RHCRIT = 1.0
    with computation(PARALLEL), interval(...):
        # include grid cell area scaling and limit RHcrit to > 70%\
        ALPHA = max(0.0, min(0.30, (1.0 - RHCRIT) * sqrt(sqrt(AREA / 1.0e10))))


def hystpdf(
    DT_MOIST: Float,
    ALPHA: FloatField,
    PDFSHAPE: Float,
    cnv_frc: FloatFieldIJ,
    srf_type: FloatFieldIJ,
    PL: FloatField,
    Q: FloatField,
    QLLS: FloatField,
    QLCN: FloatField,
    QILS: FloatField,
    QICN: FloatField,
    TE: FloatField,
    CLLS: FloatField,
    CLCN: FloatField,
    NL: FloatField,
    NI: FloatField,
    RHX: FloatField,
    ese: FloatField_VaporSaturationTable,
    esw: FloatField_VaporSaturationTable,
    esx: FloatField_VaporSaturationTable,
    estfrz: Float,
    estlqu: Float,
):
    from __externals__ import use_bergeron

    # Reference Fortran: Process_Library.F90: subroutine hystpdf
    # with PDFSHAPE = 1, use_bergeron = True, and SC_ICE = False
    with computation(PARALLEL), interval(...):
        scice = 1.0  # don't ask, I don't know
        if CLCN < 1.0:
            tmpARR = 1.0 / (1.0 - CLCN)
        else:
            tmpARR = 0.0
        if (
            CLCN > 0.0
        ):  # need to make sure this is equivilant to fortran CLCN > tiny(0.0)
            QAx = (QLCN + QICN) / CLCN
        else:
            QAx = 0.0
        CFn = CLLS * tmpARR
        QCn = (QLLS + QILS) * tmpARR
        QCi = QILS * tmpARR
        TEn = TE
        QSx, _ = QSat_Float(ese, esx, TE, PL)
        QVn = (Q - QSx * CLCN) * tmpARR

        QT = QCn + QVn  # Total LS water after microphysics

        count = 1
        while count <= 20:
            QVp = QVn
            QCp = QCn
            CFp = CFn
            TEp = TEn
            QSn, DQS = QSat_Float(ese, esx, TEn, PL, DQSAT_trigger=True)

            if PDFSHAPE < 3:  # 1 = top-hat 2 = triangulat
                sigmaqt1 = ALPHA * QSn
                sigmaqt2 = ALPHA * QSn
            elif PDFSHAPE == 4:  # lognormal (sigma is dimensionless)
                sigmaqt1 = max(ALPHA / sqrt(3.0), 0.001)

            if PDFSHAPE < 5:
                CFn = pdffrac(PDFSHAPE, QT, sigmaqt1, sigmaqt2, QSn, CFn)
                QCn = pdfcondensate(PDFSHAPE, QT, sigmaqt1, sigmaqt2, QSn, QCn)

            if use_bergeron:
                DQALL = QCn - QCp
                Nfac = (
                    100.0 * PL * constants.R_AIR / TEn
                )  # density times conversion factor
                NIv = NI / Nfac
                fQi, DQALL = bergeron_partition(
                    DT_MOIST,
                    PL,
                    TEn,
                    QT,
                    QILS,
                    QICN,
                    QLLS,
                    QLCN,
                    NIv,
                    DQALL,
                    cnv_frc,
                    srf_type,
                    ese,
                    esw,
                    estfrz,
                    estlqu,
                    count,
                )
            else:
                fQi = ice_fraction(TEn, cnv_frc, srf_type)

            latent_heat_factor = (1.0 - fQi) * (
                constants.latent_heat_vaporization / constants.cp
            ) + fQi * (constants.latent_heat_sublimation / constants.cp)
            if PDFSHAPE == 1:
                QCn = QCp + (QCn - QCp) / (
                    1.0 - (CFn * (ALPHA - 1.0) - (QCn / QSn)) * DQS * latent_heat_factor
                )
            elif PDFSHAPE == 5:
                QCn = QCp + 0.5 * (QCn - QCp)

            if CLCN > 0:
                QAo = QAx
            else:
                QAo = 0

            QVn = QVp - (QCn - QCp)
            TEn = (
                TEp
                + (1.0 - fQi)
                * constants.latent_heat_vaporization
                / constants.cp
                * ((QCn - QCp) * (1.0 - CLCN) + (QAo - QAx) * CLCN)
                + fQi
                * constants.latent_heat_sublimation
                / constants.cp
                * ((QCn - QCp) * (1.0 - CLCN) + (QAo - QAx) * CLCN)
            )

            PDFITERS = count
            # TODO Differences in constants between fortran and python cause
            # calculation of TEn to be slightly off (at order 10^-5), causing
            # the following conditional to occasionsally trigger incorrectly,
            # leading to occasional errors in various fields at the end of the
            # PDF stencil.
            if abs(TEn - TEp) < 0.00001:
                count = 21  # break out of loop
            else:
                count = count + 1

        if CLCN < 1.0:
            CLLS = CFn * (1.0 - CLCN)
            QCn = QCn * (1.0 - CLCN)
            QAo = QAo * CLCN
        else:
            # Special case CLCN=1, i.e., box filled with anvil.
            # - Note: no guarantee QV_box > QS_box
            CLLS = 0.0  # Remove any LS cloud
            QAo = QLCN + QICN + QLLS + QILS  # Add all LS condensate to anvil type
            QCn = 0.0  # Remove same from new LS
            QT = QAo + Q  # Update total water
            # Now set anvil condensate to any excess of total water
            # over QSx (saturation value at top)
            QAo = max(QT - QSx, 0.0)

        # Now take {\em New} condensate and partition into ice and liquid

        # large scale
        QCx = QCn - (QLLS + QILS)
        if QCx < 0.0:  # net evaporation
            dQLLS = max(QCx, -QLLS)  # Water evaporates first
            dQILS = max(QCx - dQLLS, -QILS)  # Then sublimation
        else:
            dQLLS = (1.0 - fQi) * QCx
            dQILS = fQi * QCx

        # convective
        QAx = QAo - (QLCN + QICN)
        if QAx < 0.0:  # net evaporation
            dQLCN = max(QAx, -QLCN)  # Water evaporates first
            dQICN = max(QAx - dQLCN, -QICN)  # Then sublimation
        else:
            dQLCN = (1.0 - fQi) * QAx
            dQICN = fQi * QAx

        # Clean-up cloud if fractions are too small
        if CLCN < 1.0e-5:
            dQICN = -QICN
            dQLCN = -QLCN
        if CLLS < 1.0e-5:
            dQILS = -QILS
            dQLLS = -QLLS

        QICN = QICN + dQICN
        QLCN = QLCN + dQLCN
        QILS = QILS + dQILS
        QLLS = QLLS + dQLLS
        Q = Q - (dQICN + dQILS + dQLCN + dQLLS)
        TE = (
            TE
            + constants.latent_heat_vaporization
            / constants.cpdry
            * (dQICN + dQILS + dQLCN + dQLLS)
            + constants.latent_heat_fusion / constants.cpdry * (dQICN + dQILS)
        )

        # We need to take care of situations where QS moves past QA
        # during QSAT iteration. This should be only when QA/AF is small
        # to begin with. Effect is to make QAo negative. So, we
        # "evaporate" offending QA's
        # We get rid of anvil fraction also, although strictly
        # speaking, PDF-wise, we should not do this.

        if QAo <= 0.0:
            Q = Q + QICN + QLCN
            TE = (
                TE
                - constants.latent_heat_sublimation / constants.cpdry * QICN
                - constants.latent_heat_vaporization / constants.cpdry * QLCN
            )
            QICN = 0.0
            QLCN = 0.0
            CLCN = 0.0

        denom, _ = QSat_Float(ese, esx, TE, PL)
        RHX = Q / denom


def melt_freeze(
    dt: Float,
    cnv_frc: FloatFieldIJ,
    srf_type: FloatFieldIJ,
    T: FloatField,
    QLCN: FloatField,
    QICN: FloatField,
):
    with computation(PARALLEL), interval(...):
        if T <= constants.t_ice:
            fQi = ice_fraction(T, cnv_frc, srf_type)
            dQil = QLCN * (1.0 - exp(-dt * fQi / constants.taufrz))
            dQil = max(0.0, dQil)
            QICN = QICN + dQil
            QLCN = QLCN - dQil
            T = (
                T
                + (
                    constants.latent_heat_sublimation
                    - constants.latent_heat_vaporization
                )
                * dQil
                / constants.cp
            )
        else:
            dQil = -QICN
            dQil = min(0.0, dQil)
            QICN = QICN + dQil
            QLCN = QLCN - dQil
            T = (
                T
                + (
                    constants.latent_heat_sublimation
                    - constants.latent_heat_vaporization
                )
                * dQil
                / constants.cp
            )


def evaporate(
    DT_MOIST: Float,
    CCW_EVAP_EFF: Float,
    PLmb: FloatField,
    T: FloatField,
    Q: FloatField,
    QLCN: FloatField,
    QICN: FloatField,
    CLCN: FloatField,
    NACTL: FloatField,
    NACTI: FloatField,
    QST: FloatField,
    EVAPC: FloatField,
):
    with computation(PARALLEL), interval(...):
        EVAPC = Q
        RHCRIT = 1
        # Evaporation of cloud water. DelGenio et al formulation
        # (Eq.s 15-17, 1996, J. Clim., 9, 270-303)
        ES = (
            100.0 * PLmb * QST / (constants.epsilon + (1.0 - constants.epsilon) * QST)
        )  # (100's <-^ convert from mbar to Pa)
        RHx = min(Q / QST, 1.00)
        K1 = (
            (constants.latent_heat_vaporization ** 2)
            * constants.RHO_W
            / (constants.K_COND * constants.rvap * (T ** 2))
        )
        K2 = (
            constants.rvap
            * T
            * constants.RHO_W
            / (constants.DIFFU * (1000.0 / PLmb) * ES)
        )
        # Here, DIFFU is given for 1000 mb so 1000./PLmb accounts
        # for increased diffusivity at lower pressure
        if CLCN > 0.0 and QLCN > 0.0:
            QCm = QLCN / CLCN
        else:
            QCm = 0.0
        RADIUS = cloud_effective_radius_liquid(PLmb, T, QCm, NACTL, NACTI)
        if RHx < RHCRIT and RADIUS > 0.0:
            EVAP = (
                CCW_EVAP_EFF
                * QLCN
                * DT_MOIST
                * (RHCRIT - RHx)
                / ((K1 + K2) * RADIUS ** 2)
            )
            EVAP = min(EVAP, QLCN)
        else:
            EVAP = 0.0
        QC = QLCN + QICN
        if QC > 0.0:
            CLCN = CLCN * (QC - EVAP) / QC
        Q = Q + EVAP
        QLCN = QLCN - EVAP
        T = T - (constants.latent_heat_vaporization / constants.cpdry) * EVAP
        EVAPC = (Q - EVAPC) / DT_MOIST


def sublimate(
    DT_MOIST: Float,
    CCI_EVAP_EFF: Float,
    PLmb: FloatField,
    T: FloatField,
    Q: FloatField,
    QLCN: FloatField,
    QICN: FloatField,
    CLCN: FloatField,
    NACTL: FloatField,
    NACTI: FloatField,
    QST: FloatField,
    EVAPC: FloatField,
):
    with computation(PARALLEL), interval(...):
        SUBLC = Q
        RHCRIT = 1
        # Sublimation of cloud water. DelGenio et al formulation
        # (Eq.s 15-17, 1996, J. Clim., 9, 270-303)
        ES = (
            100.0 * PLmb * QST / (constants.epsilon + (1.0 - constants.epsilon) * QST)
        )  # (100s <-^ convert from mbar to Pa)
        RHx = min(Q / QST, 1.00)
        K1 = (
            (constants.latent_heat_vaporization ** 2)
            * constants.RHO_I
            / (constants.K_COND * constants.rvap * (T ** 2))
        )
        K2 = (
            constants.rvap
            * T
            * constants.RHO_I
            / (constants.DIFFU * (1000.0 / PLmb) * ES)
        )
        # Here, DIFFU is given for 1000 mb so 1000./PLmb accounts
        # for increased diffusivity at lower pressure
        if CLCN > 0.0 and QICN > 0.0:
            QCm = QICN / CLCN
        else:
            QCm = 0.0
        radius = cloud_effective_radius_ice(PLmb, T, QCm, NACTL, NACTI)
        if RHx < RHCRIT and radius > 0.0:
            SUBL = (
                CCI_EVAP_EFF
                * QICN
                * DT_MOIST
                * (RHCRIT - RHx)
                / ((K1 + K2) * radius ** 2)
            )
            SUBL = min(SUBL, QICN)
        else:
            SUBL = 0.0
        QC = QLCN + QICN
        if QC > 0.0:
            CLCN = CLCN * (QC - SUBL) / QC
        Q = Q + SUBL
        QICN = QICN - SUBL
        T = T - (constants.latent_heat_sublimation / constants.cpdry) * SUBL
        SUBLC = (Q - SUBLC) / DT_MOIST
