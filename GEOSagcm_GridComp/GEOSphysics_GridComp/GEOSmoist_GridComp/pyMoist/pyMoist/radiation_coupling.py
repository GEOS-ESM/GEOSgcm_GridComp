# pylint: ignore reportInvalidTypeForm
from gt4py.cartesian import gtscript
from gt4py.cartesian.gtscript import (
    PARALLEL,
    computation,
    interval,
    log10,
)
from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float, FloatField

import pyMoist.radiation_coupling_constants as radconstants

# ruff: noqa: PLR0913


@gtscript.function
def air_density(pl: Float, te: Float) -> Float:
    """
    Calculate air density [kg/m^3]

    Parameters:
    pl (Float): Pressure level.
    te (Float): Temperature.

    Returns:
    Float: Calculated air density.

    """
    air_density = (100.0 * pl) / (radconstants.MAPL_RGAS * te)
    return air_density


@gtscript.function
def cloud_effective_radius_ice(
    pl: Float,
    te: Float,
    qc: Float,
) -> Float:
    """
    Calculate the effective radius of ice clouds [m]
    Implementation of LDRADIUS4 for Ice clouds

    Parameters:
    PL (Float): Pressure level.
    TE (Float): Temperature.
    qc (Float): Ice cloud mixing ratio.

    Returns:
    Float: Effective radius of ice clouds.
    """
    # Calculate ice water content
    wc = 1.0e3 * air_density(pl, te) * qc  # air density [g/m3] * ice cloud mixing ratio [kg/kg]
    # Calculate radius in meters [m]
    if radconstants.ICE_RADII_PARAM == 1:
        # Ice cloud effective radius -- [klaus wyser, 1998]
        if te > radconstants.MAPL_TICE or qc <= 0.0:
            bb = -2.0
        else:
            bb = -2.0 + log10(wc / 50.0) * (1.0e-3 * (radconstants.MAPL_TICE - te) ** 1.5)
        bb = min(max(bb, -6.0), -2.0)
        radius = 377.4 + 203.3 * bb + 37.91 * bb**2 + 2.3696 * bb**3
        radius = min(150.0e-6, max(5.0e-6, 1.0e-6 * radius))
    else:
        # Ice cloud effective radius ----- [Sun, 2001]
        tc = te - radconstants.MAPL_TICE
        zfsr = 1.2351 + 0.0105 * tc
        aa = 45.8966 * (wc**0.2214)
        bb = 0.79570 * (wc**0.2535)
        radius = zfsr * (aa + bb * (te - 83.15))
        radius = min(150.0e-6, max(5.0e-6, 1.0e-6 * radius * 0.64952))
    return radius


@gtscript.function
def cloud_effective_radius_liquid(
    pl: Float,
    te: Float,
    qc: Float,
    nnl: Float,
) -> Float:
    """
    Calculate the effective radius of liquid clouds [m]
    Implementation of LDRADIUS4 for liquid clouds

    Parameters:
    pl (Float): Pressure level.
    te (Float): Temperature.
    qc (Float): Liquid cloud mixing ratio.
    nnl (Float): Number concentration of liquid cloud droplets.

    Returns:
    Float: Effective radius of liquid clouds.
    """
    # Calculate liquid water content
    liquid_water_content = (
        1.0e3 * air_density(pl, te) * qc
    )  # air density [g/m3] * liquid cloud mixing ratio [kg/kg]
    # Calculate cloud drop number concentration from the aerosol model + ....
    nnx = max(nnl * 1.0e-6, 10.0)
    # Calculate Radius in meters [m]
    if radconstants.LIQ_RADII_PARAM == 1:
        # Jason Version
        radius = min(
            60.0e-6,
            max(
                2.5e-6,
                1.0e-6
                * radconstants.BX
                * (liquid_water_content / nnx) ** radconstants.R13BBETA
                * radconstants.ABETA
                * 6.92,
            ),
        )
    else:
        # [liu&daum, 2000 and 2005. liu et al 2008]
        radius = min(
            60.0e-6,
            max(
                2.5e-6,
                1.0e-6 * radconstants.LBX * (liquid_water_content / nnx) ** radconstants.LBE,
            ),
        )
    return radius


def _fix_up_clouds_stencil(
    qv: FloatField,
    te: FloatField,
    q_lc: FloatField,
    q_ic: FloatField,
    cf: FloatField,
    q_la: FloatField,
    q_ia: FloatField,
    af: FloatField,
) -> None:
    """
    Fix cloud variables to ensure physical consistency.

    Parameters:
    qv (FloatField): Water vapor mixing ratio.
    te (FloatField): Temperature.
    q_lc (FloatField): Liquid cloud mixing ratio.
    q_ic (FloatField): Ice cloud mixing ratio.
    cf (FloatField): Cloud fraction.
    q_la (FloatField): Anvil liquid cloud mixing ratio.
    q_ia (FloatField): Anvil ice cloud mixing ratio.
    af (FloatField): Anvil cloud fraction.
    """
    with computation(PARALLEL), interval(...):
        # Fix if Anvil cloud fraction too small
        if af < 1.0e-5:
            qv = qv + q_la + q_ia
            te = te - (radconstants.ALHLBCP) * q_la - (radconstants.ALHSBCP) * q_ia
            af = 0.0
            q_la = 0.0
            q_ia = 0.0
        # Fix if LS cloud fraction too small
        if cf < 1.0e-5:
            qv = qv + q_lc + q_ic
            te = te - (radconstants.ALHLBCP) * q_lc - (radconstants.ALHSBCP) * q_ic
            cf = 0.0
            q_lc = 0.0
            q_ic = 0.0
        # LS LIQUID too small
        if q_lc < 1.0e-8:
            qv = qv + q_lc
            te = te - (radconstants.ALHLBCP) * q_lc
            q_lc = 0.0
        # LS ICE too small
        if q_ic < 1.0e-8:
            qv = qv + q_ic
            te = te - (radconstants.ALHSBCP) * q_ic
            q_ic = 0.0
        # Anvil LIQUID too small
        if q_la < 1.0e-8:
            qv = qv + q_la
            te = te - (radconstants.ALHLBCP) * q_la
            q_la = 0.0
        # Anvil ICE too small
        if q_ia < 1.0e-8:
            qv = qv + q_ia
            te = te - (radconstants.ALHSBCP) * q_ia
            q_ia = 0.0
        # Fix ALL cloud quants if Anvil cloud LIQUID+ICE too small
        if (q_la + q_ia) < 1.0e-8:
            qv = qv + q_la + q_ia
            te = te - (radconstants.ALHLBCP) * q_la - (radconstants.ALHSBCP) * q_ia
            af = 0.0
            q_la = 0.0
            q_ia = 0.0
        # Fix ALL cloud quants if LS cloud LIQUID+ICE too small
        if (q_lc + q_ic) < 1.0e-8:
            qv = qv + q_lc + q_ic
            te = te - (radconstants.ALHLBCP) * q_lc - (radconstants.ALHSBCP) * q_ic
            cf = 0.0
            q_lc = 0.0
            q_ic = 0.0


def _radcouple_stencil(
    te: FloatField,
    pl: FloatField,
    cf: FloatField,
    af: FloatField,
    qv: FloatField,
    q_cl_ls: FloatField,
    q_ci_ls: FloatField,
    q_cl_an: FloatField,
    q_ci_an: FloatField,
    q_rn_all: FloatField,
    q_sn_all: FloatField,
    q_gr_all: FloatField,
    nl: FloatField,
    rad_qv: FloatField,  # pylint: ignore W0613
    rad_ql: FloatField,
    rad_qi: FloatField,
    rad_qr: FloatField,
    rad_qs: FloatField,
    rad_qg: FloatField,
    rad_cf: FloatField,
    rad_rl: FloatField,
    rad_ri: FloatField,
    fac_rl: Float,
    min_rl: Float,
    max_rl: Float,
    fac_ri: Float,
    min_ri: Float,
    max_ri: Float,
) -> None:
    """
    Couple radiation with cloud variables to ensure physical consistency.

    Parameters:
    te (FloatField): Temperature.
    pl (FloatField): Pressure level.
    cf (FloatField): Cloud fraction.
    af (FloatField): Anvil cloud fraction.
    qv (FloatField): Water vapor mixing ratio.
    q_cl_ls (FloatField): Liquid cloud mixing ratio (large-scale).
    q_ci_ls (FloatField): Ice cloud mixing ratio (large-scale).
    q_cl_an (FloatField): Liquid cloud mixing ratio (anvil).
    q_ci_an (FloatField): Ice cloud mixing ratio (anvil).
    q_rn_all (FloatField): Rain mixing ratio.
    q_sn_all (FloatField): Snow mixing ratio.
    q_gr_all (FloatField): Graupel mixing ratio.
    nl (FloatField): Number concentration of liquid cloud droplets.
    rad_qv (FloatField): Radiation water vapor mixing ratio.
    rad_ql (FloatField): Radiation liquid cloud mixing ratio.
    rad_qi (FloatField): Radiation ice cloud mixing ratio.
    rad_qr (FloatField): Radiation rain mixing ratio.
    rad_qs (FloatField): Radiation snow mixing ratio.
    rad_qg (FloatField): Radiation graupel mixing ratio.
    rad_cf (FloatField): Radiation cloud fraction.
    rad_rl (FloatField): Radiation liquid effective radius.
    rad_ri (FloatField): Radiation ice effective radius.
    fac_rl (Float): Factor for liquid effective radius.
    min_rl (Float): Minimum liquid effective radius.
    max_rl (Float): Maximum liquid effective radius.
    fac_ri (Float): Factor for ice effective radius.
    min_ri (Float): Minimum ice effective radius.
    max_ri (Float): Maximum ice effective radius.
    """
    with computation(PARALLEL), interval(...):
        # water vapor
        rad_qv = qv

        # total cloud fraction
        rad_cf = max(min(cf + af, 1.0), 0.0)
        if rad_cf >= 1.0e-5:
            rad_ql = (q_cl_ls + q_cl_an) / rad_cf if (q_cl_ls + q_cl_an) >= 1.0e-8 else 0.0
            rad_qi = (q_ci_ls + q_ci_an) / rad_cf if (q_ci_ls + q_ci_an) >= 1.0e-8 else 0.0
            rad_qr = q_rn_all / rad_cf if q_rn_all >= 1.0e-8 else 0.0
            rad_qs = q_sn_all / rad_cf if q_sn_all >= 1.0e-8 else 0.0
            rad_qg = q_gr_all / rad_cf if q_gr_all >= 1.0e-8 else 0.0
        else:
            rad_cf = 0.0
            rad_ql = 0.0
            rad_qi = 0.0
            rad_qr = 0.0
            rad_qs = 0.0
            rad_qg = 0.0

        # Cap the high end of condensates
        rad_ql = min(rad_ql, 0.01)
        rad_qi = min(rad_qi, 0.01)
        rad_qr = min(rad_qr, 0.01)
        rad_qs = min(rad_qs, 0.01)
        rad_qg = min(rad_qg, 0.01)

        # Liquid radii - Brams formulation with limits
        rad_rl = max(
            min_rl,
            min(
                cloud_effective_radius_liquid(pl, te, rad_ql, nl) * fac_rl,
                max_rl,
            ),
        )
        # Ice radii - Brams formulation with limits
        rad_ri = max(
            min_ri,
            min(
                cloud_effective_radius_ice(pl, te, rad_qi) * fac_ri,
                max_ri,
            ),
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
                "[Radiation Coupling] Diagnostic (do_qa) not implemented." "(GEOS_QSAT missing)",
            )

    def __call__(
        self,
        q: FloatField,
        te: FloatField,
        q_lls: FloatField,
        q_ils: FloatField,
        clls: FloatField,
        q_lcn: FloatField,
        q_icn: FloatField,
        clcn: FloatField,
        pl: FloatField,
        q_rain: FloatField,
        q_snow: FloatField,
        q_graupel: FloatField,
        nactl: FloatField,
        nacti: FloatField,
        rad_q_v: FloatField,
        rad_q_l: FloatField,
        rad_q_i: FloatField,
        rad_q_r: FloatField,
        rad_q_s: FloatField,
        rad_q_g: FloatField,
        rad_c_f: FloatField,
        cldreff_l: FloatField,
        cldreff_i: FloatField,
        fac_rl: Float,
        min_rl: Float,
        max_rl: Float,
        fac_ri: Float,
        min_ri: Float,
        max_ri: Float,
    ):
        """
        Perform the radiation coupling process.

        Parameters:
        q (FloatField): Water vapor mixing ratio.
        te (FloatField): Temperature.
        q_lls (FloatField): Liquid cloud mixing ratio (large-scale).
        q_ils (FloatField): Ice cloud mixing ratio (large-scale).
        clls (FloatField): Cloud fraction (large-scale).
        q_lcn (FloatField): Liquid cloud mixing ratio (anvil).
        q_icn (FloatField): Ice cloud mixing ratio (anvil).
        clcn (FloatField): Cloud fraction (anvil).
        pl (FloatField): Pressure level in millibars.
        q_rain (FloatField): Rain mixing ratio.
        q_snow (FloatField): Snow mixing ratio.
        q_graupel (FloatField): Graupel mixing ratio.
        nactl (FloatField): Number concentration of liquid cloud droplets.
        nacti (FloatField): Number concentration of ice cloud crystals.
        rad_q_v (FloatField): Radiation water vapor mixing ratio.
        rad_q_l (FloatField): Radiation liquid cloud mixing ratio.
        rad_q_i (FloatField): Radiation ice cloud mixing ratio.
        rad_q_r (FloatField): Radiation rain mixing ratio.
        rad_q_s (FloatField): Radiation snow mixing ratio.
        rad_q_g (FloatField): Radiation graupel mixing ratio.
        rad_c_f (FloatField): Radiation cloud fraction.
        cldreff_l (FloatField): Radiation liquid effective radius.
        cldreff_i (FloatField): Radiation ice effective radius.
        fac_rl (Float): Factor for liquid effective radius.
        min_rl (Float): Minimum liquid effective radius.
        max_rl (Float): Maximum liquid effective radius.
        fac_ri (Float): Factor for ice effective radius.
        min_ri (Float): Minimum ice effective radius.
        max_ri (Float): Maximum ice effective radius.
        """
        self._fix_up_clouds(
            qv=q,
            te=te,
            q_lc=q_lls,
            q_ic=q_ils,
            cf=clls,
            q_la=q_lcn,
            q_ia=q_icn,
            af=clcn,
        )
        self._radcouple(
            te,
            pl,
            clls,
            clcn,
            q,
            q_lls,
            q_ils,
            q_lcn,
            q_icn,
            q_rain,
            q_snow,
            q_graupel,
            nactl,
            nacti,
            rad_q_v,
            rad_q_l,
            rad_q_i,
            rad_q_r,
            rad_q_s,
            rad_q_g,
            rad_c_f,
            cldreff_l,
            cldreff_i,
            fac_rl,
            min_rl,
            max_rl,
            fac_ri,
            min_ri,
            max_ri,
        )
