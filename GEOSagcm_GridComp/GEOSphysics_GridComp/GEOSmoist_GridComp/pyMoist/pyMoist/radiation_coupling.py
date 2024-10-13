import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import PARALLEL, computation, interval, log10

import pyMoist.constants as constants
from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float, FloatField


@gtscript.function
def air_density(PL: Float, TE: Float) -> Float:
    """
    Calculate air density [kg/m^3]

    Parameters:
    PL (Float): Pressure level.
    TE (Float): Temperature.

    Returns:
    Float: Calculated air density.
    """
    air_density = (100.0 * PL) / (constants.MAPL_RGAS * TE)
    return air_density


@gtscript.function
def cloud_effective_radius_ice(
    PL: Float,
    TE: Float,
    QC: Float,
    NNL: Float,
    NNI: Float,
) -> Float:
    """
    Calculate the effective radius of ice clouds [m]
    Implementation of LDRADIUS4 for Ice clouds

    Parameters:
    PL (Float): Pressure level.
    TE (Float): Temperature.
    QC (Float): Ice cloud mixing ratio.
    NNL (Float): Number concentration of liquid cloud droplets.
                 Not used in function body, but included in the Fortran code.
    NNI (Float): Number concentration of ice cloud crystals.
                 Not used in function body, but included in the Fortran code.

    Returns:
    Float: Effective radius of ice clouds.
    """
    # Calculate ice water content
    WC = (
        1.0e3 * air_density(PL, TE) * QC
    )  # air density [g/m3] * ice cloud mixing ratio [kg/kg]
    # Calculate radius in meters [m]
    if constants.ICE_RADII_PARAM == 1:
        # Ice cloud effective radius -- [klaus wyser, 1998]
        if TE > constants.MAPL_TICE or QC <= 0.0:
            BB = -2.0
        else:
            BB = -2.0 + log10(WC / 50.0) * (1.0e-3 * (constants.MAPL_TICE - TE) ** 1.5)
        BB = min(max(BB, -6.0), -2.0)
        RADIUS = 377.4 + 203.3 * BB + 37.91 * BB ** 2 + 2.3696 * BB ** 3
        RADIUS = min(150.0e-6, max(5.0e-6, 1.0e-6 * RADIUS))
    else:
        # Ice cloud effective radius ----- [Sun, 2001]
        TC = TE - constants.MAPL_TICE
        ZFSR = 1.2351 + 0.0105 * TC
        AA = 45.8966 * (WC ** 0.2214)
        BB = 0.79570 * (WC ** 0.2535)
        RADIUS = ZFSR * (AA + BB * (TE - 83.15))
        RADIUS = min(150.0e-6, max(5.0e-6, 1.0e-6 * RADIUS * 0.64952))
    return RADIUS


@gtscript.function
def cloud_effective_radius_liquid(
    PL: Float,
    TE: Float,
    QC: Float,
    NNL: Float,
    NNI: Float,
) -> Float:
    """
    Calculate the effective radius of liquid clouds [m]
    Implementation of LDRADIUS4 for liquid clouds

    Parameters:
    PL (Float): Pressure level.
    TE (Float): Temperature.
    QC (Float): Liquid cloud mixing ratio.
    NNL (Float): Number concentration of liquid cloud droplets.
    NNI (Float): Number concentration of ice cloud crystals. Not used in function body.

    Returns:
    Float: Effective radius of liquid clouds.
    """
    # Calculate liquid water content
    WC = (
        1.0e3 * air_density(PL, TE) * QC
    )  # air density [g/m3] * liquid cloud mixing ratio [kg/kg]
    # Calculate cloud drop number concentration from the aerosol model + ....
    NNX = max(NNL * 1.0e-6, 10.0)
    # Calculate Radius in meters [m]
    if constants.LIQ_RADII_PARAM == 1:
        # Jason Version
        RADIUS = min(
            60.0e-6,
            max(
                2.5e-6,
                1.0e-6
                * constants.BX
                * (WC / NNX) ** constants.R13BBETA
                * constants.ABETA
                * 6.92,
            ),
        )
    else:
        # [liu&daum, 2000 and 2005. liu et al 2008]
        RADIUS = min(
            60.0e-6,
            max(2.5e-6, 1.0e-6 * constants.LBX * (WC / NNX) ** constants.LBE),
        )
    return RADIUS


def _fix_up_clouds_stencil(
    QV: FloatField,
    TE: FloatField,
    QLC: FloatField,
    QIC: FloatField,
    CF: FloatField,
    QLA: FloatField,
    QIA: FloatField,
    AF: FloatField,
) -> None:
    """
    Fix cloud variables to ensure physical consistency.

    Parameters:
    QV (3D inout): Water vapor mixing ratio.
    TE (3D inout): Temperature.
    QLC (3D inout): Liquid cloud mixing ratio.
    QIC (3D inout): Ice cloud mixing ratio.
    CF (3D inout): Cloud fraction.
    QLA (3D inout): Anvil liquid cloud mixing ratio.
    QIA (3D inout): Anvil ice cloud mixing ratio.
    AF (3D inout): Anvil cloud fraction.
    """
    with computation(PARALLEL), interval(...):
        # Fix if Anvil cloud fraction too small
        if AF < 1.0e-5:
            QV = QV + QLA + QIA
            TE = TE - (constants.ALHLBCP) * QLA - (constants.ALHSBCP) * QIA
            AF = 0.0
            QLA = 0.0
            QIA = 0.0
        # Fix if LS cloud fraction too small
        if CF < 1.0e-5:
            QV = QV + QLC + QIC
            TE = TE - (constants.ALHLBCP) * QLC - (constants.ALHSBCP) * QIC
            CF = 0.0
            QLC = 0.0
            QIC = 0.0
        # LS LIQUID too small
        if QLC < 1.0e-8:
            QV = QV + QLC
            TE = TE - (constants.ALHLBCP) * QLC
            QLC = 0.0
        # LS ICE too small
        if QIC < 1.0e-8:
            QV = QV + QIC
            TE = TE - (constants.ALHSBCP) * QIC
            QIC = 0.0
        # Anvil LIQUID too small
        if QLA < 1.0e-8:
            QV = QV + QLA
            TE = TE - (constants.ALHLBCP) * QLA
            QLA = 0.0
        # Anvil ICE too small
        if QIA < 1.0e-8:
            QV = QV + QIA
            TE = TE - (constants.ALHSBCP) * QIA
            QIA = 0.0
        # Fix ALL cloud quants if Anvil cloud LIQUID+ICE too small
        if (QLA + QIA) < 1.0e-8:
            QV = QV + QLA + QIA
            TE = TE - (constants.ALHLBCP) * QLA - (constants.ALHSBCP) * QIA
            AF = 0.0
            QLA = 0.0
            QIA = 0.0
        # Fix ALL cloud quants if LS cloud LIQUID+ICE too small
        if (QLC + QIC) < 1.0e-8:
            QV = QV + QLC + QIC
            TE = TE - (constants.ALHLBCP) * QLC - (constants.ALHSBCP) * QIC
            CF = 0.0
            QLC = 0.0
            QIC = 0.0


def _radcouple_stencil(
    TE: FloatField,
    PL: FloatField,
    CF: FloatField,
    AF: FloatField,
    QV: FloatField,
    QClLS: FloatField,
    QCiLS: FloatField,
    QClAN: FloatField,
    QCiAN: FloatField,
    QRN_ALL: FloatField,
    QSN_ALL: FloatField,
    QGR_ALL: FloatField,
    NL: FloatField,
    NI: FloatField,
    RAD_QV: FloatField,
    RAD_QL: FloatField,
    RAD_QI: FloatField,
    RAD_QR: FloatField,
    RAD_QS: FloatField,
    RAD_QG: FloatField,
    RAD_CF: FloatField,
    RAD_RL: FloatField,
    RAD_RI: FloatField,
    FAC_RL: Float,
    MIN_RL: Float,
    MAX_RL: Float,
    FAC_RI: Float,
    MIN_RI: Float,
    MAX_RI: Float,
) -> None:
    """
    Couple radiation with cloud variables to ensure physical consistency.

    Parameters:
    TE (3D in): Temperature.
    PL (3D in): Pressure level.
    CF (3D in): Cloud fraction.
    AF (3D in): Anvil cloud fraction.
    QV (3D in): Water vapor mixing ratio.
    QClLS (3D in): Liquid cloud mixing ratio (large-scale).
    QCiLS (3D in): Ice cloud mixing ratio (large-scale).
    QClAN (3D in): Liquid cloud mixing ratio (anvil).
    QCiAN (3D in): Ice cloud mixing ratio (anvil).
    QRN_ALL (3D in): Rain mixing ratio.
    QSN_ALL (3D in): Snow mixing ratio.
    QGR_ALL (3D in): Graupel mixing ratio.
    NL (3D in): Number concentration of liquid cloud droplets.
    NI (3D in): Number concentration of ice cloud crystals.
    RAD_QV (3D inout): Radiation water vapor mixing ratio.
    RAD_QL (3D inout): Radiation liquid cloud mixing ratio.
    RAD_QI (3D inout): Radiation ice cloud mixing ratio.
    RAD_QR (3D inout): Radiation rain mixing ratio.
    RAD_QS (3D inout): Radiation snow mixing ratio.
    RAD_QG (3D inout): Radiation graupel mixing ratio.
    RAD_CF (3D inout): Radiation cloud fraction.
    RAD_RL (3D out): Radiation liquid effective radius.
    RAD_RI (3D out): Radiation ice effective radius.
    FAC_RL (Float): Factor for liquid effective radius.
    MIN_RL (Float): Minimum liquid effective radius.
    MAX_RL (Float): Maximum liquid effective radius.
    FAC_RI (Float): Factor for ice effective radius.
    MIN_RI (Float): Minimum ice effective radius.
    MAX_RI (Float): Maximum ice effective radius.
    """
    with computation(PARALLEL), interval(...):
        # water vapor
        RAD_QV = QV

        # total cloud fraction
        RAD_CF = max(min(CF + AF, 1.0), 0.0)
        if RAD_CF >= 1.0e-5:
            RAD_QL = (QClLS + QClAN) / RAD_CF if (QClLS + QClAN) >= 1.0e-8 else 0.0
            RAD_QI = (QCiLS + QCiAN) / RAD_CF if (QCiLS + QCiAN) >= 1.0e-8 else 0.0
            RAD_QR = QRN_ALL / RAD_CF if QRN_ALL >= 1.0e-8 else 0.0
            RAD_QS = QSN_ALL / RAD_CF if QSN_ALL >= 1.0e-8 else 0.0
            RAD_QG = QGR_ALL / RAD_CF if QGR_ALL >= 1.0e-8 else 0.0
        else:
            RAD_CF = 0.0
            RAD_QL = 0.0
            RAD_QI = 0.0
            RAD_QR = 0.0
            RAD_QS = 0.0
            RAD_QG = 0.0

        # Cap the high end of condensates
        RAD_QL = min(RAD_QL, 0.01)
        RAD_QI = min(RAD_QI, 0.01)
        RAD_QR = min(RAD_QR, 0.01)
        RAD_QS = min(RAD_QS, 0.01)
        RAD_QG = min(RAD_QG, 0.01)

        # Liquid radii - Brams formulation with limits
        RAD_RL = max(
            MIN_RL,
            min(cloud_effective_radius_liquid(PL, TE, RAD_QL, NL, NI) * FAC_RL, MAX_RL),
        )
        # Ice radii - Brams formulation with limits
        RAD_RI = max(
            MIN_RI,
            min(cloud_effective_radius_ice(PL, TE, RAD_QI, NL, NI) * FAC_RI, MAX_RI),
        )


class RadiationCoupling:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        do_qa: bool,
    ) -> None:
        """
        Initialize the RadiationCoupling class.

        Parameters:
        stencil_factory (StencilFactory): Factory to create stencils.
        quantity_factory (QuantityFactory): Factory to create quantities.
        do_qa (bool): Flag to indicate if QA should be performed.
        """
        self._fix_up_clouds = stencil_factory.from_dims_halo(
            func=_fix_up_clouds_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )
        self._radcouple = stencil_factory.from_dims_halo(
            func=_radcouple_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )
        self.do_qa = do_qa
        # GEOS_GFDL_1M_InterfaceMod.F90:866 not implemented.
        # Implement QSAT0 and QSAT3 logic if diagnostics are needed,
        # found in GEOS/src/Shared/@GMAO_Shared/GEOS_Shared/GEOS_Utilities.F90.
        if self.do_qa:
            # RHX = Q/GEOS_QSAT( T, PLmb)
            raise NotImplementedError(
                "[Radiation Coupling] Diagnostic (do_qa) not implemented."
                "(GEOS_QSAT missing)"
            )

    def __call__(
        self,
        Q: FloatField,
        T: FloatField,
        QLLS: FloatField,
        QILS: FloatField,
        CLLS: FloatField,
        QLCN: FloatField,
        QICN: FloatField,
        CLCN: FloatField,
        PLmb: FloatField,
        QRAIN: FloatField,
        QSNOW: FloatField,
        QGRAUPEL: FloatField,
        NACTL: FloatField,
        NACTI: FloatField,
        RAD_QV: FloatField,
        RAD_QL: FloatField,
        RAD_QI: FloatField,
        RAD_QR: FloatField,
        RAD_QS: FloatField,
        RAD_QG: FloatField,
        RAD_CF: FloatField,
        CLDREFFL: FloatField,
        CLDREFFI: FloatField,
        FAC_RL: Float,
        MIN_RL: Float,
        MAX_RL: Float,
        FAC_RI: Float,
        MIN_RI: Float,
        MAX_RI: Float,
    ):
        """
        Perform the radiation coupling process.

        Parameters:
        Q (FloatField): Water vapor mixing ratio.
        T (FloatField): Temperature.
        QLLS (FloatField): Liquid cloud mixing ratio (large-scale).
        QILS (FloatField): Ice cloud mixing ratio (large-scale).
        CLLS (FloatField): Cloud fraction (large-scale).
        QLCN (FloatField): Liquid cloud mixing ratio (anvil).
        QICN (FloatField): Ice cloud mixing ratio (anvil).
        CLCN (FloatField): Cloud fraction (anvil).
        PLmb (FloatField): Pressure level in millibars.
        QRAIN (FloatField): Rain mixing ratio.
        QSNOW (FloatField): Snow mixing ratio.
        QGRAUPEL (FloatField): Graupel mixing ratio.
        NACTL (FloatField): Number concentration of liquid cloud droplets.
        NACTI (FloatField): Number concentration of ice cloud crystals.
        RAD_QV (FloatField): Radiation water vapor mixing ratio.
        RAD_QL (FloatField): Radiation liquid cloud mixing ratio.
        RAD_QI (FloatField): Radiation ice cloud mixing ratio.
        RAD_QR (FloatField): Radiation rain mixing ratio.
        RAD_QS (FloatField): Radiation snow mixing ratio.
        RAD_QG (FloatField): Radiation graupel mixing ratio.
        RAD_CF (FloatField): Radiation cloud fraction.
        CLDREFFL (FloatField): Radiation liquid effective radius.
        CLDREFFI (FloatField): Radiation ice effective radius.
        FAC_RL (Float): Factor for liquid effective radius.
        MIN_RL (Float): Minimum liquid effective radius.
        MAX_RL (Float): Maximum liquid effective radius.
        FAC_RI (Float): Factor for ice effective radius.
        MIN_RI (Float): Minimum ice effective radius.
        MAX_RI (Float): Maximum ice effective radius.
        """
        self._fix_up_clouds(
            QV=Q,
            TE=T,
            QLC=QLLS,
            QIC=QILS,
            CF=CLLS,
            QLA=QLCN,
            QIA=QICN,
            AF=CLCN,
        )
        self._radcouple(
            TE=T,
            PL=PLmb,
            CF=CLLS,
            AF=CLCN,
            QV=Q,
            QClLS=QLLS,
            QCiLS=QILS,
            QClAN=QLCN,
            QCiAN=QICN,
            QRN_ALL=QRAIN,
            QSN_ALL=QSNOW,
            QGR_ALL=QGRAUPEL,
            NL=NACTL,
            NI=NACTI,
            RAD_QV=RAD_QV,
            RAD_QL=RAD_QL,
            RAD_QI=RAD_QI,
            RAD_QR=RAD_QR,
            RAD_QS=RAD_QS,
            RAD_QG=RAD_QG,
            RAD_CF=RAD_CF,
            RAD_RL=CLDREFFL,
            RAD_RI=CLDREFFI,
            FAC_RL=FAC_RL,
            MIN_RL=MIN_RL,
            MAX_RL=MAX_RL,
            FAC_RI=FAC_RI,
            MIN_RI=MIN_RI,
            MAX_RI=MAX_RI,
        )
