import copy

from gt4py.cartesian.gtscript import (
    BACKWARD,
    FORWARD,
    PARALLEL,
    THIS_K,
    computation,
    erfc,
    exp,
    f32,
    f64,
    i32,
    i64,
    interval,
    isnan,
    log,
    sqrt,
)

import pyMoist.constants as constants
from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import BoolFieldIJ, Float, FloatField, FloatFieldIJ, Int, IntField, IntFieldIJ
from pyMoist.field_types import FloatField_NTracers, FloatFieldIJ_NTracers
from pyMoist.saturation.formulation import SaturationFormulation
from pyMoist.saturation.qsat import FloatField_Extra_Dim, QSat, QSat_Float
from pyMoist.UW.config import UWConfiguration
from pyMoist.UW.temporaries import Temporaries
from pyMoist.UW.uwshcu_functions import (
    compute_alpha,
    compute_mumin2,
    compute_ppen,
    conden,
    exnerfn,
    getbuoy,
    qsinvert,
    roots,
    single_cin,
    slope_bot,
    slope_mid,
    slope_top,
    zvir,
)


def compute_uwshcu_invert_before(
    # Inputs
    k0: Int,
    pmid0_inv: FloatField,
    u0_inv: FloatField,
    v0_inv: FloatField,
    zmid0_inv: FloatField,
    exnmid0_inv: FloatField,
    dp0_inv: FloatField,
    qv0_inv: FloatField,
    ql0_inv: FloatField,
    qi0_inv: FloatField,
    t0_inv: FloatField,
    tke_inv: FloatField,
    pifc0_inv: FloatField,
    zifc0_inv: FloatField,
    exnifc0_inv: FloatField,
    kpbl_inv: IntFieldIJ,
    cnvtr: FloatFieldIJ,
    frland: FloatFieldIJ,
    tr0_inout: FloatField_NTracers,
    CNV_Tracers: FloatField_NTracers,
    # Outputs
    pmid0_in: FloatField,
    u0_in: FloatField,
    v0_in: FloatField,
    zmid0_in: FloatField,
    exnmid0_in: FloatField,
    dp0_in: FloatField,
    qv0_in: FloatField,
    ql0_in: FloatField,
    qi0_in: FloatField,
    th0_in: FloatField,
    tke_in: FloatField,
    pifc0_in: FloatField,
    zifc0_in: FloatField,
    exnifc0_in: FloatField,
    kpbl_in: IntFieldIJ,
    cnvtrmax: FloatFieldIJ,
):
    from __externals__ import k_end, ncnst

    with computation(PARALLEL), interval(...):
        # Flip mid-level variables
        k_inv = k_end - THIS_K
        pmid0_in = pmid0_inv.at(K=k_inv)
        u0_in = u0_inv.at(K=k_inv)
        v0_in = v0_inv.at(K=k_inv)
        zmid0_in = zmid0_inv.at(K=k_inv)
        exnmid0_in = exnmid0_inv.at(K=k_inv)
        dp0_in = dp0_inv.at(K=k_inv)
        qv0_in = qv0_inv.at(K=k_inv)
        ql0_in = ql0_inv.at(K=k_inv)
        qi0_in = qi0_inv.at(K=k_inv)
        th0_in = t0_inv.at(K=k_inv) / exnmid0_inv.at(K=k_inv)
        n = 0
        while n < ncnst:
            tr0_inout[0, 0, 0][n] = CNV_Tracers.at(K=k_inv, ddim=[n])
            n += 1

    with computation(FORWARD), interval(...):
        tke_in[0, 0, 1] = 0.0
        pifc0_in[0, 0, 1] = 0.0
        zifc0_in[0, 0, 1] = 0.0
        exnifc0_in[0, 0, 1] = 0.0

    with computation(FORWARD), interval(...):
        # Flip interface variables
        k_inv = k0 - THIS_K
        tke_in[0, 0, 1] = tke_inv.at(K=k_inv)
        tke_in = tke_inv.at(K=k_inv)
        pifc0_in[0, 0, 1] = pifc0_inv.at(K=k_inv)
        pifc0_in = pifc0_inv.at(K=k_inv)
        zifc0_in[0, 0, 1] = zifc0_inv.at(K=k_inv)
        zifc0_in = zifc0_inv.at(K=k_inv)
        exnifc0_in[0, 0, 1] = exnifc0_inv.at(K=k_inv)
        exnifc0_in = exnifc0_inv.at(K=k_inv)

        kpbl_in = kpbl_inv

    with computation(FORWARD), interval(...):
        cnvtrmax = min(1e-5, max(0.0, cnvtr))
        if frland > 0.5:
            cnvtrmax = 0.0
        if isnan(cnvtrmax):
            cnvtrmax = 0.0


def compute_thermodynamic_variables(
    pmid0_in: FloatField,
    zmid0_in: FloatField,
    exnmid0_in: FloatField,
    dp0_in: FloatField,
    u0_in: FloatField,
    v0_in: FloatField,
    u0: FloatField,
    v0: FloatField,
    qv0_in: FloatField,
    ql0_in: FloatField,
    qi0_in: FloatField,
    th0_in: FloatField,
    tr0_inout: FloatField_NTracers,
    cush_inout: FloatFieldIJ,
    cush: FloatFieldIJ,
    umf_out: FloatField,
    shfx: FloatFieldIJ,
    evap: FloatFieldIJ,
    qtflx_out: FloatField,
    slflx_out: FloatField,
    uflx_out: FloatField,
    vflx_out: FloatField,
    qt0: FloatField,
    t0: FloatField,
    qv0: FloatField,
    qi0: FloatField,
    pmid0: FloatField,
    tr0: FloatField_NTracers,
    tr0_temp: FloatField,
    sstr0: FloatField_NTracers,
    ssthl0: FloatField,
    ssqt0: FloatField,
    thl0: FloatField,
    dotransport: Int,
    ssu0: FloatField,
    ssv0: FloatField,
    tscaleh: FloatFieldIJ,
    fer_out: FloatField,
    fdr_out: FloatField,
    tpert_out: FloatFieldIJ,
    qpert_out: FloatFieldIJ,
):
    """
    Start of Main Calculation.

    Stencil to compute basic thermodynamic variables directly from
    input variables for each column.

    Inputs:
    pmid0_in [FloatField]: Environmental pressure at the layer mid-point [Pa]
    zmid0_in [FloatField]: Environmental height at the layer mid-point [m]
    exnmid0_in [FloatField]: Exner function at the layer mid-point
    dp0_in [FloatField]: Environmental layer pressure thickness [Pa] > 0
    u0_in [FloatField]: Environmental zonal wind [m/s]
    v0_in [FloatField]: Environmental meridional wind [m/s]
    qv0_in [FloatField]: Environmental water vapor specific humidity [kg/kg]
    ql0_in [FloatField]: Environmental liquid water specific humidity [kg/kg]
    qi0_in [FloatField]: Environmental ice specific humidity [kg/kg]
    th0_in [FloatField]: Environmental potential temperature [K]
    tr0_inout [FloatField_NTracers]: Environmental tracers [#, kg/kg]
    cush_inout [FloatFieldIJ]: Convective scale height [m]
    dotransport [Int]: Transport tracers [1 true]

    Outputs:
    u0 [FloatField]: Environmental zonal wind [m/s]
    v0 [FloatField]: Environmental meridional wind [m/s]
    cush [FloatFieldIJ]: Convective scale height [m]
    umf_out [FloatField]: Updraft mass flux at the interfaces [kg/m2/s]
    shfx [FloatFieldIJ]: Surface sensible heat [J]
    evap [FloatFieldIJ]: Surface evaporation [kg/m^2/s]
    qtflx_out [FloatField]: Mixing ratio flux [?]
    slflx_out [FloatField]: Sensible heat flux [?]
    uflx_out [FloatField]: Zonal wind flux [?]
    vflx_out [FloatField]: Meridional wind flux [?]
    qt0 [FloatField]: Mixing ratio [?]
    t0 [FloatField]: Environmental temperature [K]
    qv0 [FloatField]: Environmental specific humidity
    qi0 [FloatField]: Environmental ice specific humidity
    pmid0 [FloatField]: Environmental pressure at the layer mid-point [Pa]
    tr0 [FloatField_NTracers]: Environmental tracers [#, kg/kg]
    tr0_temp [FloatField]: Environmental tracers [#, kg/kg]
    sstr0 [FloatField_NTracers]: Convective tracer [?]
    ssthl0 [FloatField]: Temperature [?]
    ssqt0 [FloatField]: [?]
    thl0 [FloatField]: Temperature [?]
    ssu0 [FloatField]: [?]
    ssv0 [FloatField]: [?]
    tscaleh [FloatFieldIJ]: [?]
    fer_out [FloatField]: Fractional lateral entrainment rate [1/Pa]
    fdr_out [FloatField]: Fractional lateral detrainment rate [1/Pa]
    tpert_out [FloatFieldIJ]: Temperature perturbation
    qpert_out [FloatFieldIJ]: Humidity perturbation
    """
    from __externals__ import k_end, ncnst

    with computation(FORWARD), interval(...):

        # Initialize output variables defined for all grid points
        umf_out[0, 0, 1] = 0.0
        dcm_out = 0.0
        cufrc_out = 0.0
        fer_out = constants.MAPL_UNDEF
        fdr_out = constants.MAPL_UNDEF
        qldet_out = 0.0
        qidet_out = 0.0
        qlsub_out = 0.0
        qisub_out = 0.0
        ndrop_out = 0.0
        nice_out = 0.0
        qtflx_out[0, 0, 1] = 0.0
        slflx_out[0, 0, 1] = 0.0
        uflx_out[0, 0, 1] = 0.0
        vflx_out[0, 0, 1] = 0.0
        tpert_out = 0.0
        qpert_out = 0.0

        # Initialize variable that are calculated from inputs
        zmid0 = zmid0_in
        dp0 = dp0_in
        cush = cush_inout
        tscaleh = cush_inout

    with computation(PARALLEL), interval(...):
        pmid0 = pmid0_in
        u0 = u0_in
        v0 = v0_in
        qv0 = qv0_in
        ql0 = ql0_in
        qi0 = qi0_in
        qt0 = qv0 + ql0 + qi0
        exnmid0 = exnmid0_in
        t0 = th0_in * exnmid0
        thl0 = (
            t0
            - ((constants.MAPL_ALHL * ql0) / constants.MAPL_CP)
            - ((constants.MAPL_ALHS * qi0) / constants.MAPL_CP)
        ) / exnmid0

        if dotransport == 1:
            n = 0
            # Loop over tracers
            while n < ncnst:
                tr0[0, 0, 0][n] = tr0_inout[0, 0, 0][n]
                n += 1

    with computation(PARALLEL), interval(0, 1):
        # Compute slopes of environmental variables at bottom layer
        ssthl0 = slope_bot(thl0, pmid0)
        ssqt0 = slope_bot(qt0, pmid0)
        ssu0 = slope_bot(u0, pmid0)
        ssv0 = slope_bot(v0, pmid0)

        if dotransport == 1:
            n = 0
            while n < ncnst:
                tr0_temp = tr0[0, 0, 0][n]
                sstr0[0, 0, 0][n] = slope_bot(tr0_temp, pmid0)
                n += 1

    with computation(PARALLEL), interval(1, -1):
        # Compute slopes of environmental variables at mid layers
        ssthl0 = slope_mid(k_end, thl0, pmid0)
        ssqt0 = slope_mid(k_end, qt0, pmid0)
        ssu0 = slope_mid(k_end, u0, pmid0)
        ssv0 = slope_mid(k_end, v0, pmid0)

        if dotransport == 1:
            n = 0
            while n < ncnst:
                tr0_temp = tr0[0, 0, 0][n]
                sstr0[0, 0, 0][n] = slope_mid(k_end, tr0_temp, pmid0)
                n += 1

    with computation(PARALLEL), interval(-1, None):
        # Compute slopes of environmental variables at top layer
        ssthl0 = slope_top(k_end, thl0, pmid0)
        ssqt0 = slope_top(k_end, qt0, pmid0)
        ssu0 = slope_top(k_end, u0, pmid0)
        ssv0 = slope_top(k_end, v0, pmid0)

        if dotransport == 1:
            n = 0
            while n < ncnst:
                tr0_temp = tr0[0, 0, 0][n]
                sstr0[0, 0, 0][n] = slope_top(k_end, tr0_temp, pmid0)
                n += 1


def compute_thv0_thvl0(
    pmid0_in: FloatField,
    exnmid0_in: FloatField,
    qv0_in: FloatField,
    ql0_in: FloatField,
    qi0_in: FloatField,
    th0_in: FloatField,
    zmid0: FloatField,
    pifc0_in: FloatField,
    ssthl0: FloatField,
    ssqt0: FloatField,
    ese: FloatField_Extra_Dim,
    esx: FloatField_Extra_Dim,
    umf_out: FloatField,
    qtflx_out: FloatField,
    slflx_out: FloatField,
    uflx_out: FloatField,
    vflx_out: FloatField,
    u0_in: FloatField,
    v0_in: FloatField,
    k0: Int,
    condensation: BoolFieldIJ,
    dotransport: Int,
    ssu0: FloatField,
    ssv0: FloatField,
    tr0: FloatField_NTracers,
    sstr0: FloatField_NTracers,
    tr0_o: FloatField_NTracers,
    sstr0_o: FloatField_NTracers,
    trflx: FloatField_NTracers,
    trten: FloatField_NTracers,
    tru: FloatField_NTracers,
    tru_emf: FloatField_NTracers,
    umf_zint: FloatField,
    emf: FloatField,
    slflx: FloatField,
    qtflx: FloatField,
    uflx: FloatField,
    vflx: FloatField,
    thlu: FloatField,
    qtu: FloatField,
    uu: FloatField,
    vu: FloatField,
    wu: FloatField,
    thvu: FloatField,
    thlu_emf: FloatField,
    qtu_emf: FloatField,
    uu_emf: FloatField,
    vu_emf: FloatField,
    uemf: FloatField,
    thvl0bot: FloatField,
    thvl0top: FloatField,
    thvl0: FloatField,
    qt0: FloatField,
    t0: FloatField,
    qv0: FloatField,
    thl0: FloatField,
    ql0: FloatField,
    qi0: FloatField,
    thv0bot: FloatField,
    thv0top: FloatField,
    uten: FloatField,
    vten: FloatField,
    s0: FloatField,
    qcu: FloatField,
    qlu: FloatField,
    qiu: FloatField,
    cufrc: FloatField,
    ufrc: FloatField,
    qlten_det: FloatField,
    qiten_det: FloatField,
    qlten_sink: FloatField,
    qiten_sink: FloatField,
    sten: FloatField,
    slten: FloatField,
    qiten: FloatField,
    qv0_o: FloatField,
    ql0_o: FloatField,
    qi0_o: FloatField,
    t0_o: FloatField,
    s0_o: FloatField,
    u0_o: FloatField,
    v0_o: FloatField,
    qt0_o: FloatField,
    thl0_o: FloatField,
    thvl0_o: FloatField,
    ssthl0_o: FloatField,
    ssqt0_o: FloatField,
    thv0bot_o: FloatField,
    thv0top_o: FloatField,
    thvl0bot_o: FloatField,
    thvl0top_o: FloatField,
    ssu0_o: FloatField,
    ssv0_o: FloatField,
    cush_inout: FloatFieldIJ,
):
    """
    Stencil to compute internal environmental variables

    The 'condensation' mask is introduced in this stencil.
    This mask is used to indicate if the condensation process has occurred.
    It is set to False by default and becomes True if condensation occurs.
    Once condensation occurs, the UW shallow convection scheme is done
    computing at that column.

    Variables ending in '_o' indicate variables used in the second
    iteration of the implicit CIN calculation.

    Inputs:
    pmid0_in [FloatField]: Environmental pressure at midpoints [Pa]
    exnmid0_in [FloatField]: Exner function at the layer mid-point
    qv0_in [FloatField]: Environmental specific humidity
    ql0_in [FloatField]: Environmental liquid water specific humidity [kg/kg]
    qi0_in [FloatField]: Environmental ice specific humidity [kg/kg]
    th0_in [FloatField]: Environmental potential temperature [K]
    zmid0 [FloatField]: Environmental height at the layer mid-point [m]
    pifc0_in [FloatField]: Environmental pressure at interfaces [Pa]
    ese [FloatField_Extra_Dim]: Used in QSat_Float [n/a]
    esx [FloatField_Extra_Dim]: Used in QSat_Float [n/a]
    u0_in [FloatField]: Environmental zonal wind [m/s]
    v0_in [FloatField]: Environmental meridional wind [m/s]
    k0 [Int]: Number of levels
    dotransport [Int]: Transport tracers [1 true]
    ssu0 [FloatField]: [?]
    ssv0 [FloatField]: [?]
    tr0 [FloatField_NTracers]: Environmental tracers [#, kg/kg]
    sstr0 [FloatField_NTracers]: Convective tracer [?]
    t0 [FloatField]: Environmental temperature [K]

    Outputs:
    ssthl0 [FloatField]: Temperature [?]
    ssqt0 [FloatField]: Temperature [?]
    umf_out [FloatField]: Updraft mass flux at the interfaces [kg/m2/s]
    qtflx_out [FloatField]: Mixing ratio flux [?]
    slflx_out [FloatField]: Sensible heat flux [?]
    uflx_out [FloatField]: Zonal wind flux [?]
    vflx_out [FloatField]: Meridional wind flux [?]
    condensation [BoolFieldIJ]: Mask that indicates if condensation has occurred
    tr0_o [FloatField_NTracers]: Environmental tracers [#, kg/kg]
    sstr0_o [FloatField_NTracers]: Convective tracer [?]
    trflx [FloatField_NTracers]: Tracer PBL flux [?]
    trten [FloatField_NTracers]: Tendency of [?]
    tru [FloatField_NTracers]: Updraft tracers [#, kg/kg]
    tru_emf [FloatField_NTracers]: Penetrative Downdraft tracers at entraining interfaces [#, kg/kg]
    umf_zint [FloatField]: Updraft mass flux at the interfaces [kg/m2/s]
    emf [FloatField]: Penetrative [?]
    slflx [FloatField]: Sensible heat flux [?]
    qtflx [FloatField]: Mixing ratio flux [?]
    uflx [FloatField]: Zonal wind flux [?]
    vflx [FloatField]: Meridional wind flux [?]
    thlu [FloatField]: Updraft liquid potential temperature at the interface [K]
    qtu [FloatField]: Updraft total specific humidity at the interface [kg/kg]
    uu [FloatField]: Updraft zonal wind at the interface [m/s]
    vu [FloatField]: Updraft meridional wind at the interface [m/s]
    wu [FloatField]: Updraft vertical velocity at the interface [m/s]
    thvu [FloatField]: Updraft virtual potential temperature at the interface [m/s]
    thlu_emf [FloatField]: Penetrative downdraft liquid potential temperature at entraining interfaces [K]
    qtu_emf [FloatField]: Penetrative downdraft total water at entraining interfaces [kg/kg]
    uu_emf [FloatField]: Penetrative downdraft zonal wind at entraining interfaces [m/s]
    vu_emf [FloatField]: Penetrative downdraft meridional wind at entraining interfaces [m/s]
    uemf [FloatField]: Net updraft mass flux at the interface ( emf + umf ) [kg/m2/s]
    thvl0bot [FloatField]: Temperature at bottom [?]
    thvl0top [FloatField]: Temperature at top [?]
    thvl0 [FloatField]: Temperature [?]
    qt0 [FloatField]: Mixing ratio [?]
    qv0 [FloatField]: Environmental specific humidity
    thl0 [FloatField]: Temperature [?]
    ql0 [FloatField]: Environmental liquid water specific humidity
    qi0 [FloatField]: Environmental ice specific humidity
    thv0bot [FloatField]: Temperature at bottom [?]
    thv0top [FloatField]: Potential temperature at top [?]
    uten [FloatField]: Tendency of zonal wind [m/s2]
    vten [FloatField]: Tendency of meridional wind [m/s2]
    s0 [FloatField]: Environmental dry static energy [J/kg]
    qcu [FloatField]: Condensate water specific humidity within cumulus updraft at the layer mid-point [kg/kg]
    qlu [FloatField]: Liquid water specific humidity within cumulus updraft at the layer mid-point [kg/kg]
    qiu [FloatField]: Ice specific humidity within cumulus updraft at the layer mid-point [kg/kg]
    cufrc [FloatField]: Shallow cumulus cloud fraction at the layer mid-point [fraction]
    ufrc [FloatField]: Cumulus updraft fraction [fraction]
    qlten_det [FloatField]: [?]
    qiten_det [FloatField]: [?]
    qlten_sink [FloatField]: Liquid condensate tendency by compensating subsidence/upwelling [kg/kg/s]
    qiten_sink [FloatField]: Ice condensate tendency by compensating subsidence/upwelling [kg/kg/s]
    sten [FloatField]: Tendency of dry static energy [J/kg/s]
    slten [FloatField]: Tendency of [?]
    qiten [FloatField]: Tendency of ice specific humidity [kg/kg/s]
    qv0_o [FloatField]: Environmental water vapor specific humidity [kg/kg]
    ql0_o [FloatField]: Environmental liquid water specific humidity [kg/kg]
    qi0_o [FloatField]: Environmental ice specific humidity [kg/kg]
    t0_o [FloatField]: Environmental temperature [K]
    s0_o [FloatField]: Environmental dry static energy [J/kg]
    u0_o [FloatField]: Environmental zonal wind [m/s]
    v0_o [FloatField]: Environmental meridional wind [m/s]
    qt0_o [FloatField]: Mixing ratio [?]
    thl0_o [FloatField]: Temperature [?]
    thvl0_o [FloatField]: Temperature [?]
    ssthl0_o [FloatField]: Temperature [?]
    ssqt0_o [FloatField]: Temperature [?]
    thv0bot_o [FloatField]: Temperature at bottom [?]
    thv0top_o [FloatField]: Temperature at top [?]
    thvl0bot_o [FloatField]: Temperature at bottom [?]
    thvl0top_o [FloatField]: Temperature at top [?]
    ssu0_o [FloatField]: [?]
    ssv0_o [FloatField]: [?]
    cush_inout [FloatFieldIJ]: Convective scale height [m]
    """
    from __externals__ import ncnst

    with computation(PARALLEL), interval(...):
        pmid0 = pmid0_in
        exnmid0 = exnmid0_in
        qv0 = qv0_in
        ql0 = ql0_in
        qi0 = qi0_in
        qt0 = qv0 + ql0 + qi0
        t0 = th0_in * exnmid0
        s0 = constants.MAPL_GRAV * zmid0 + constants.MAPL_CP * t0
        thl0 = (
            t0
            - constants.MAPL_LATENT_HEAT_VAPORIZATION * ql0 / constants.MAPL_CP
            - constants.MAPL_LATENT_HEAT_SUBLIMATION * qi0 / constants.MAPL_CP
        ) / exnmid0
        thvl0 = (1.0 + zvir * qt0) * thl0
        pifc0 = pifc0_in

        thl0bot = thl0 + ssthl0 * (pifc0 - pmid0)
        qt0bot = qt0 + ssqt0 * (pifc0 - pmid0)

        thj, qvj, qlj, qij, qse, id_check = conden(pifc0, thl0bot, qt0bot, ese, esx)

    with computation(FORWARD), interval(...):
        if id_check == 1:
            condensation = (
                True  # Indicates if condensation has occurred (e.g., Fortran go to 333)
            )
            umf_out[0, 0, 1] = 0.0
            dcm_out = 0.0
            qvten_out = 0.0
            qlten_out = 0.0
            qiten_out = 0.0
            sten_out = 0.0
            uten_out = 0.0
            vten_out = 0.0
            qrten_out = 0.0
            qsten_out = 0.0
            cufrc_out = 0.0
            cush_inout = -1.0
            qldet_out = 0.0
            qidet_out = 0.0
            qtflx_out[0, 0, 1] = 0.0
            slflx_out[0, 0, 1] = 0.0
            uflx_out[0, 0, 1] = 0.0
            vflx_out[0, 0, 1] = 0.0
            fer_out = constants.MAPL_UNDEF
            fdr_out = constants.MAPL_UNDEF

    with computation(PARALLEL), interval(...):
        if not condensation:
            thv0bot = thj * (1.0 + zvir * qvj - qlj - qij)
            thvl0bot = thl0bot * (1.0 + zvir * qt0bot)

            thl0top = thl0 + ssthl0 * (pifc0_in[0, 0, 1] - pmid0)
            qt0top = qt0 + ssqt0 * (pifc0_in[0, 0, 1] - pmid0)

    with computation(PARALLEL), interval(0, -1):
        if not condensation:
            thj, qvj, qlj, qij, qse, id_check = conden(
                pifc0_in[0, 0, 1], thl0top, qt0top, ese, esx
            )

    with computation(FORWARD), interval(0, -1):
        if not condensation:
            if id_check == 1:
                condensation = True
                umf_out[0, 0, 1] = 0.0
                dcm_out = 0.0
                qvten_out = 0.0
                qlten_out = 0.0
                qiten_out = 0.0
                sten_out = 0.0
                uten_out = 0.0
                vten_out = 0.0
                qrten_out = 0.0
                qsten_out = 0.0
                cufrc_out = 0.0
                cush_inout = -1.0
                qldet_out = 0.0
                qidet_out = 0.0
                qtflx_out[0, 0, 1] = 0.0
                slflx_out[0, 0, 1] = 0.0
                uflx_out[0, 0, 1] = 0.0
                vflx_out[0, 0, 1] = 0.0
                fer_out = constants.MAPL_UNDEF
                fdr_out = constants.MAPL_UNDEF

    with computation(PARALLEL), interval(0, -1):
        if not condensation:
            thv0top = thj * (1.0 + zvir * qvj - qlj - qij)
            thvl0top = thl0top * (1.0 + zvir * qt0top)

    with computation(PARALLEL), interval(-1, None):
        if not condensation:
            thv0top = thv0bot
            thvl0top = thvl0bot

    with computation(PARALLEL), interval(...):
        if not condensation:
            # Save input and related environmental thermodynamic variables
            # for use at "iter_cin=2" when "del_CIN >= 0"
            qv0_o = qv0
            ql0_o = ql0
            qi0_o = qi0
            t0_o = t0
            s0_o = s0
            u0_o = u0_in
            v0_o = v0_in
            qt0_o = qt0
            thl0_o = thl0
            thvl0_o = thvl0
            ssthl0_o = ssthl0
            ssqt0_o = ssqt0
            thv0bot_o = thv0bot
            thv0top_o = thv0top
            thvl0bot_o = thvl0bot
            thvl0top_o = thvl0top
            ssu0_o = ssu0
            ssv0_o = ssv0

            if dotransport == 1:
                n = 0
                while n < ncnst:
                    tr0_o[0, 0, 0][n] = tr0[0, 0, 0][n]
                    sstr0_o[0, 0, 0][n] = sstr0[0, 0, 0][n]
                    n += 1

    with computation(FORWARD), interval(...):
        if not condensation:
            # Initialize output variables at each grid point
            umf_zint[0, 0, 1] = 0.0
            dcm = 0.0
            emf[0, 0, 1] = 0.0
            slflx[0, 0, 1] = 0.0
            qtflx[0, 0, 1] = 0.0
            uflx[0, 0, 1] = 0.0
            vflx[0, 0, 1] = 0.0
            qvten = 0.0
            qlten = 0.0
            qiten = 0.0
            sten = 0.0
            uten = 0.0
            vten = 0.0
            qrten = 0.0
            qsten = 0.0
            dwte = 0.0
            diten = 0.0
            cufrc = 0.0
            qcu = 0.0
            qlu = 0.0
            qiu = 0.0
            fer = 0.0
            fdr = 0.0
            xco = 0.0
            cin = 0.0
            cinlcl = 0.0
            cbmf = 0.0
            qc = 0.0
            qc_l = 0.0
            qc_i = 0.0
            cnt = f32(k0)
            cnb = 0.0
            qtten = 0.0
            slten = 0.0
            ufrc = 0.0

            thlu = constants.MAPL_UNDEF
            thlu[0, 0, 1] = constants.MAPL_UNDEF
            qtu = constants.MAPL_UNDEF
            qtu[0, 0, 1] = constants.MAPL_UNDEF
            uu[0, 0, 1] = constants.MAPL_UNDEF
            uu = constants.MAPL_UNDEF
            vu[0, 0, 1] = constants.MAPL_UNDEF
            vu = constants.MAPL_UNDEF
            wu[0, 0, 1] = constants.MAPL_UNDEF
            thvu[0, 0, 1] = constants.MAPL_UNDEF
            thlu_emf[0, 0, 1] = constants.MAPL_UNDEF
            thlu_emf = constants.MAPL_UNDEF
            qtu_emf[0, 0, 1] = constants.MAPL_UNDEF
            qtu_emf = constants.MAPL_UNDEF
            uu_emf[0, 0, 1] = constants.MAPL_UNDEF
            uu_emf = constants.MAPL_UNDEF
            vu_emf[0, 0, 1] = constants.MAPL_UNDEF
            vu_emf = constants.MAPL_UNDEF

            ufrclcl = 0.0
            wlcl = 0.0
            cbmflimit = 0.0

            uemf[0, 0, 1] = 0.0
            comsub = 0.0
            qlten_sink = 0.0
            qiten_sink = 0.0
            qlten_det = 0.0
            qiten_det = 0.0

            if dotransport == 1:
                n = 0
                while n < ncnst:
                    trflx[0, 0, 1][n] = 0.0
                    trten[0, 0, 0][n] = 0.0
                    tru[0, 0, 1][n] = 0.0
                    tru_emf[0, 0, 1][n] = 0.0
                    n += 1


def find_pbl_height(
    iteration: i32,
    kpbl_in: IntFieldIJ,
    k0: Int,
    condensation: BoolFieldIJ,
    umf_out: FloatField,
    qtflx_out: FloatField,
    slflx_out: FloatField,
    uflx_out: FloatField,
    vflx_out: FloatField,
    kinv: IntField,
    cush: FloatFieldIJ,
    cush_inout: FloatFieldIJ,
    tscaleh: FloatFieldIJ,
):
    """
    Stencil to find PBL top height interface index, 'kinv-1' where 'kinv' is the layer
    index with PBLH in it. When PBLH is exactly at interface, 'kinv' is the
    layer index having PBLH as a lower interface.

    Inputs:
    iteration [i32]: Iteration of implicit CIN loop (i.e., 0 or 1)
    kpbl_in [IntFieldIJ]: Height of PBL [m]
    k0 [Int]: Number of levels

    Inouts:
    condensation [BoolFieldIJ]: Mask that indicates if condensation has occurred
    cush [FloatFieldIJ]: Convective scale height [m]

    Outputs:
    umf_out [FloatField]: Updraft mass flux at the interfaces [kg/m2/s]
    qtflx_out [FloatField]: Mixing ratio flux [?]
    slflx_out [FloatField]: Sensible heat flux [?]
    uflx_out [FloatField]: Zonal wind flux [?]
    vflx_out [FloatField]: Meridional wind flux [?]
    kinv [IntField]: Inversion layer with PBL top interface as lower interface
    cush_inout [FloatFieldIJ]: Convective scale height [m]
    tscaleh [FloatFieldIJ]: [?]
    """
    with computation(FORWARD), interval(...):
        if not condensation:
            if iteration != i32(0):
                tscaleh = cush

    with computation(FORWARD), interval(...):
        if not condensation:
            # Cumulus scale height
            # In contrast to the premitive code, cumulus scale height is iteratively
            # calculated at each time step, and at each iterative cin step.
            # It is not clear whether I should locate below two lines within or out
            # of the iterative cin loop.

            cush = -1.0
            qtavg = 0.0

            # Find PBL top height interface index, 'kinv-1' where 'kinv' is the layer
            # index with PBLH in it. When PBLH is exactly at interface, 'kinv' is the
            # layer index having PBLH as a lower interface.
            # In the previous code, I set the lower limit of 'kinv' by 2  in order to
            # be consistent with the other parts of the code. However in the modified
            # code, I allowed 'kinv' to be 1 & if 'kinv = 1', I just exit the program
            # without performing cumulus convection. This new approach seems to be
            # more reasonable: if PBL height is within 'kinv=1' layer, surface is STL
            # interface (bflxs <= 0) and interface just above the surface should be
            # either non-turbulent (Ri>0.19) or stably turbulent (0<=Ri<0.19 but this
            # interface is identified as a base external interface of upperlying CL.
            # Thus, when 'kinv=1', PBL scheme guarantees 'bflxs <= 0'.  For this case
            # it is reasonable to assume that cumulus convection does not happen.
            # When these is SBCL, PBL height from the PBL scheme is likely to be very
            # close at 'kinv-1' interface, but not exactly, since 'zi' information is
            # changed between two model time steps. In order to ensure correct identi
            # fication of 'kinv' for general case including SBCL, I imposed an offset
            # of 5 [m] in the below 'kinv' finding block.

            # Invert kpbl index
            if kpbl_in > i64(k0 / 2):
                kinv = k0 - kpbl_in + 1
            else:
                kinv = 5

            if kinv <= i64(1):
                condensation = True
                umf_out[0, 0, 1] = 0.0
                dcm_out = 0.0
                qvten_out = 0.0
                qlten_out = 0.0
                qiten_out = 0.0
                sten_out = 0.0
                uten_out = 0.0
                vten_out = 0.0
                qrten_out = 0.0
                qsten_out = 0.0
                cufrc_out = 0.0
                cush_inout = -1.0
                qldet_out = 0.0
                qidet_out = 0.0
                qtflx_out[0, 0, 1] = 0.0
                slflx_out[0, 0, 1] = 0.0
                uflx_out[0, 0, 1] = 0.0
                vflx_out[0, 0, 1] = 0.0
                fer_out = constants.MAPL_UNDEF
                fdr_out = constants.MAPL_UNDEF

            if not condensation:
                if kinv >= i64(k0 / 4):
                    condensation = True
                    umf_out[0, 0, 1] = 0.0
                    dcm_out = 0.0
                    qvten_out = 0.0
                    qlten_out = 0.0
                    qiten_out = 0.0
                    sten_out = 0.0
                    uten_out = 0.0
                    vten_out = 0.0
                    qrten_out = 0.0
                    qsten_out = 0.0
                    cufrc_out = 0.0
                    cush_inout = -1.0
                    qldet_out = 0.0
                    qidet_out = 0.0
                    qtflx_out[0, 0, 1] = 0.0
                    slflx_out[0, 0, 1] = 0.0
                    uflx_out[0, 0, 1] = 0.0
                    vflx_out[0, 0, 1] = 0.0
                    fer_out = constants.MAPL_UNDEF
                    fdr_out = constants.MAPL_UNDEF


def find_pbl_averages(
    condensation: BoolFieldIJ,
    thvl0bot: FloatField,
    thvl0top: FloatField,
    kinv: IntField,
    pifc0: FloatField,
    tke_in: FloatField,
    u0: FloatField,
    v0: FloatField,
    thvl0: FloatField,
    zmid0: FloatField,
    qtsrchgt: Float,
    qt0: FloatField,
    thvlmin: FloatField,
    tkeavg: FloatField,
    uavg: FloatField,
    vavg: FloatField,
    thvlavg: FloatField,
    qtavg: FloatField,
    iteration: i32,
):
    with computation(FORWARD), interval(0, 1):
        if not condensation:
            thvlmin = 1000.0
            thvlmin = min(
                thvlmin,
                min(thvl0bot, thvl0top),
            )

    with computation(FORWARD), interval(1, None):
        if not condensation:
            if THIS_K <= (kinv - 2):
                thvlmin = min(
                    thvlmin[0, 0, -1],
                    min(thvl0bot, thvl0top),
                )

    with computation(FORWARD), interval(...):
        if not condensation:

            # Find PBL averaged tke ('tkeavg') and minimum 'thvl' ('thvlmin') in the PBL
            # In the current code, 'tkeavg' is obtained by averaging all interfacial TKE
            # within the PBL. However, in order to be conceptually consistent with   PBL
            # scheme, 'tkeavg' should be calculated by considering surface buoyancy flux
            # If surface buoyancy flux is positive ( bflxs >0 ), surface interfacial TKE
            # should be included in calculating 'tkeavg', while if bflxs <= 0,   surface
            # interfacial TKE should not be included in calculating 'tkeavg'.   I should
            # modify the code when 'bflxs' is available as an input of cumulus scheme.
            # 'thvlmin' is a minimum 'thvl' within PBL obtained by comparing top &  base
            # interface values of 'thvl' in each layers within the PBL.

            dpsum = 0.0
            thvlavg = 0.0
            tkeavg = 0.0
            uavg = 0.0
            vavg = 0.0

            lev = 0
            while lev < kinv:
                kabove = lev + 1
                dpi = pifc0.at(K=lev) - pifc0.at(K=kabove)
                dpsum = dpsum + dpi
                tkeavg = tkeavg + dpi * tke_in.at(K=kabove)
                uavg = uavg + dpi * u0.at(K=lev)
                vavg = vavg + dpi * v0.at(K=lev)
                thvlavg = thvlavg + dpi * thvl0.at(K=lev)
                lev += 1

            tkeavg = tkeavg / dpsum
            uavg = uavg / dpsum
            vavg = vavg / dpsum
            thvlavg = thvlavg / dpsum

            # Interpolate qt to specified height
            lev = 0
            while zmid0.at(K=lev) < qtsrchgt:
                lev += 1
            if lev > 0:
                kbelow = lev - 1
                qtavg = qt0.at(K=kbelow) * (zmid0.at(K=lev) - qtsrchgt) + qt0.at(
                    K=lev
                ) * (qtsrchgt - zmid0.at(K=kbelow))
                qtavg = qtavg / (zmid0.at(K=lev) - zmid0.at(K=kbelow))
            else:
                qtavg = qt0.at(K=0)


def find_cumulus_characteristics(
    condensation: BoolFieldIJ,
    windsrcavg: Int,
    pifc0: FloatField,
    t0: FloatField,
    qv0: FloatField,
    shfx: FloatFieldIJ,
    evap: FloatFieldIJ,
    thlsrc_fac: Float,
    qtsrc_fac: Float,
    qt0: FloatField,
    qtavg: FloatField,
    thvlmin: FloatField,
    uavg: FloatField,
    vavg: FloatField,
    kinv: IntField,
    u0: FloatField,
    v0: FloatField,
    ssu0: FloatField,
    ssv0: FloatField,
    pmid0: FloatField,
    dotransport: Int,
    tr0: FloatField_NTracers,
    trsrc: FloatFieldIJ_NTracers,
    qtsrc: FloatField,
    thvlsrc: FloatField,
    thlsrc: FloatField,
    usrc: FloatField,
    vsrc: FloatField,
    tpert_out: FloatFieldIJ,
    qpert_out: FloatFieldIJ,
    iteration: i32,
):
    """
    Stencil to find characteristics of cumulus source air:
    qtsrc,thlsrc,usrc,vsrc.

    Note that 'thlsrc' was concocted using 'thvlsrc' and 'qtsrc'.
    'qtsrc' is defined as the lowest layer mid-point value;   'thlsrc'
    is from 'qtsrc' and 'thvlmin=thvlsrc'; 'usrc' & 'vsrc' are defined
    as the values just below the PBL top interface.

    Inputs:
    condensation [BoolFieldIJ]: Mask that indicates if condensation has occurred
    windsrcavg [Int]: Source air uses PBL mean momentum
    pifc0 [FloatField]: Environmental pressure at the interfaces [Pa]
    t0 [FloatField]: Environmental temperature [K]
    qv0 [FloatField]: Environmental specific humidity
    shfx [FloatFieldIJ]: Surface sensible heat [J]
    evap [FloatFieldIJ]: Surface evaporation [kg/m^2/s]
    thlsrc_fac [Float]: Scaling factor for thlsrc perturbation
    qtsrc_fac [Float]: Scaling factor for qtsrc perturbation
    qt0 [FloatField]: Mixing ratio [?]
    qtavg [FloatField]: Average qt over the PBL [?]
    thvlmin [FloatField]: Minimum 'thvl' within PBL, obtained by comparing top & base
    uavg [FloatField]: Average u over the PBL [?]
    vavg [FloatField]: Average v over the PBL [?]
    kinv [IntField]: Inversion layer with PBL top interface as lower interface
    u0 [FloatField]: Environmental zonal wind [m/s]
    v0 [FloatField]: Environmental meridional wind [m/s]
    ssu0 [FloatField]: [?]
    ssv0 [FloatField]: [?]
    pmid0 [FloatField]: Environmental pressure at the layer mid-point [Pa]
    dotransport [Int]: Transport tracers [1 true]
    tr0 [FloatField_NTracers]: Environmental tracers [#, kg/kg]
    iteration [i32]: Iteration of implicit CIN loop (i.e., 0 or 1)

    Outputs:
    trsrc [FloatFieldIJ_NTracers]: Tracers of cumulus source air [?]
    qtsrc [FloatField]: Mixing ratio of cumulus source air [?]
    thvlsrc [FloatField]: Temperature of cumulus source air [K] [?]
    thlsrc [FloatField]: Temperature of cumulus source air [K] [?]
    usrc [FloatField]: Zonal wind of cumulus source air [m/s] [?]
    vsrc [FloatField]: Meridional wind of cumulus source air [m/s] [?]
    tpert_out [FloatFieldIJ]: Temperature perturbation
    qpert_out [FloatFieldIJ]: Humidity perturbation
    """
    from __externals__ import ncnst

    with computation(FORWARD), interval(...):
        if not condensation:
            if windsrcavg == 1:
                # Caution: This code has not been tested, since windsrcavg is 0
                # An error will be raised if windsrcavg == 1
                zrho = pifc0.at(K=0) / (
                    287.04 * (t0.at(K=0) * (1.0 + 0.608 * qv0.at(K=0)))
                )
                buoyflx = (
                    -shfx / constants.MAPL_CP - 0.608 * t0.at(K=0) * evap
                ) / zrho  # K m s-1
                delzg = (50.0) * constants.MAPL_GRAV  # assume 50m surface scale
                wstar = max(0.0, 0.001 - 0.41 * buoyflx * delzg / t0.at(K=0))  # m3 s-3
                qpert_out = 0.0
                tpert_out = 0.0
                if wstar > 0.001:
                    wstar = 1.0 * wstar**0.3333
                    tpert_out = (
                        thlsrc_fac * shfx / (zrho * wstar * constants.MAPL_CP)
                    )  # K
                    qpert_out = qtsrc_fac * evap / (zrho * wstar)  # kg kg-1
                    qpert_out = max(
                        min(qpert_out, 0.02 * qt0.at(K=0)), 0.0
                    )  # limit to 1% of QT
                    tpert_out = 0.1 + max(min(tpert_out, 1.0), 0.0)  # limit to 1K
                    qtsrc = qtavg + qpert_out
                    thvlsrc = thvlmin + tpert_out * (1.0 + zvir * qtsrc)
                    thlsrc = thvlsrc / (1.0 + zvir * qtsrc)
                    usrc = uavg
                    vsrc = vavg
            else:
                qtsrc = qt0.at(K=0)
                thvlsrc = thvlmin.at(K=kinv - 2)
                thlsrc = thvlsrc / (1.0 + zvir * qtsrc)
                kbelow = kinv - 2
                kbelow_zdim = kinv - 1
                usrc = u0.at(K=kbelow) + ssu0.at(K=kbelow) * (
                    pifc0.at(K=kbelow_zdim) - pmid0.at(K=kbelow)
                )
                vsrc = v0.at(K=kbelow) + ssv0.at(K=kbelow) * (
                    pifc0.at(K=kbelow_zdim) - pmid0.at(K=kbelow)
                )

            if dotransport == 1:
                n = 0
                while n < ncnst:
                    trsrc[0, 0][n] = tr0.at(K=0, ddim=[n])
                    n += 1


def find_klcl(
    condensation: BoolFieldIJ,
    pifc0: FloatField,
    umf_out: FloatField,
    qtflx_out: FloatField,
    slflx_out: FloatField,
    uflx_out: FloatField,
    vflx_out: FloatField,
    qtsrc: FloatField,
    thlsrc: FloatField,
    ese: FloatField_Extra_Dim,
    esx: FloatField_Extra_Dim,
    k0: Int,
    thl0: FloatField,
    ssthl0: FloatField,
    pmid0: FloatField,
    qt0: FloatField,
    ssqt0: FloatField,
    plcl: FloatField,
    klcl: IntField,
    thl0lcl: FloatField,
    qt0lcl: FloatField,
    thv0lcl: FloatField,
    cush_inout: FloatFieldIJ,
):
    """
    Find LCL of the source air and a layer index containing LCL (klcl)
    When the LCL is exactly at the interface, 'klcl' is a layer index
    having 'plcl' as the lower interface similar to the 'kinv' case.

    In the previous code, I assumed that if LCL is located within the
    lowest model layer ( 1 ) or the top model layer ( k0 ), then  no
    convective adjustment is performed and just exited.   However, in
    the revised code, I relaxed the first constraint and  even though
    LCL is at the lowest model layer, I allowed cumulus convection to
    be initiated. For this case, cumulus convection should be started
    from the PBL top height, as shown in the following code.
    When source air is already saturated even at the surface, klcl is
    set to 1.

    Inputs:
    pifc0 [FloatField]: Environmental pressure at the interfaces [Pa]
    qtsrc [FloatField]: Mixing ratio of cumulus source air [?]
    thlsrc [FloatField]: Temperature of cumulus source air [K] [?]
    ese [FloatField_Extra_Dim]: Used in QSat_Float [n/a]
    esx [FloatField_Extra_Dim]: Used in QSat_Float [n/a]
    k0 [Int]: Number of levels
    thl0 [FloatField]: Temperature [?]
    ssthl0 [FloatField]: Temperature [?]
    pmid0 [FloatField]: Environmental pressure at the layer mid-point [Pa]
    qt0 [FloatField]: Mixing ratio [?]
    ssqt0 [FloatField]: [?]

    Inouts:
    condensation [BoolFieldIJ]: Mask that indicates if condensation has occurred

    Outputs:
    umf_out [FloatField]: Updraft mass flux at the interfaces [kg/m2/s]
    qtflx_out [FloatField]: Mixing ratio flux [?]
    slflx_out [FloatField]: Sensible heat flux [?]
    uflx_out [FloatField]: Zonal wind flux [?]
    vflx_out [FloatField]: Meridional wind flux [?]
    plcl [FloatField]: LCL of source air [Pa]
    klcl [IntField]: Layer containing LCL of source air
    thl0lcl [FloatField]: Temperature at the LCL [?]
    qt0lcl [FloatField]: Mixing ratio at LCL [?]
    thv0lcl [FloatField]: Temperature at LCL [?]
    cush_inout [FloatFieldIJ]: Convective scale height [m]
    """
    with computation(FORWARD), interval(...):
        if not condensation:
            if pifc0.at(K=0) < 70000 or pifc0.at(K=0) > 115000.0:
                condensation = True
                umf_out[0, 0, 1] = 0.0
                dcm_out = 0.0
                qvten_out = 0.0
                qlten_out = 0.0
                qiten_out = 0.0
                sten_out = 0.0
                uten_out = 0.0
                vten_out = 0.0
                qrten_out = 0.0
                qsten_out = 0.0
                cufrc_out = 0.0
                cush_inout = -1.0
                qldet_out = 0.0
                qidet_out = 0.0
                qtflx_out[0, 0, 1] = 0.0
                slflx_out[0, 0, 1] = 0.0
                uflx_out[0, 0, 1] = 0.0
                vflx_out[0, 0, 1] = 0.0
                fer_out = constants.MAPL_UNDEF
                fdr_out = constants.MAPL_UNDEF

            if not condensation:
                if qtsrc > 0.1 or qtsrc < 1e-8:
                    condensation = True
                    umf_out[0, 0, 1] = 0.0
                    dcm_out = 0.0
                    qvten_out = 0.0
                    qlten_out = 0.0
                    qiten_out = 0.0
                    sten_out = 0.0
                    uten_out = 0.0
                    vten_out = 0.0
                    qrten_out = 0.0
                    qsten_out = 0.0
                    cufrc_out = 0.0
                    cush_inout = -1.0
                    qldet_out = 0.0
                    qidet_out = 0.0
                    qtflx_out[0, 0, 1] = 0.0
                    slflx_out[0, 0, 1] = 0.0
                    uflx_out[0, 0, 1] = 0.0
                    vflx_out[0, 0, 1] = 0.0
                    fer_out = constants.MAPL_UNDEF
                    fdr_out = constants.MAPL_UNDEF

            if not condensation:
                if thlsrc > 400.0 or thlsrc < 100.0:
                    condensation = True
                    umf_out[0, 0, 1] = 0.0
                    dcm_out = 0.0
                    qvten_out = 0.0
                    qlten_out = 0.0
                    qiten_out = 0.0
                    sten_out = 0.0
                    uten_out = 0.0
                    vten_out = 0.0
                    qrten_out = 0.0
                    qsten_out = 0.0
                    cufrc_out = 0.0
                    cush_inout = -1.0
                    qldet_out = 0.0
                    qidet_out = 0.0
                    qtflx_out[0, 0, 1] = 0.0
                    slflx_out[0, 0, 1] = 0.0
                    uflx_out[0, 0, 1] = 0.0
                    vflx_out[0, 0, 1] = 0.0
                    fer_out = constants.MAPL_UNDEF
                    fdr_out = constants.MAPL_UNDEF

            if not condensation:
                plcl = qsinvert(qtsrc, thlsrc, pifc0.at(K=0), ese, esx)

                lev = 0
                klcl_flag = 0.0
                while lev < k0 + 1 and klcl_flag == 0.0:
                    kidx = lev
                    if pifc0.at(K=kidx) < plcl:
                        klcl = lev
                        klcl_flag = 1.0
                    lev += 1

                if klcl_flag == 0.0:
                    klcl = 0

                klcl = max(0, klcl)
                klcl = klcl - 1  # Adjust klcl by 1

                if plcl < 60000.0:
                    condensation = True
                    umf_out[0, 0, 1] = 0.0
                    dcm_out = 0.0
                    qvten_out = 0.0
                    qlten_out = 0.0
                    qiten_out = 0.0
                    sten_out = 0.0
                    uten_out = 0.0
                    vten_out = 0.0
                    qrten_out = 0.0
                    qsten_out = 0.0
                    cufrc_out = 0.0
                    cush_inout = -1.0
                    qldet_out = 0.0
                    qidet_out = 0.0
                    qtflx_out[0, 0, 1] = 0.0
                    slflx_out[0, 0, 1] = 0.0
                    uflx_out[0, 0, 1] = 0.0
                    vflx_out[0, 0, 1] = 0.0
                    fer_out = constants.MAPL_UNDEF

                if not condensation:
                    # Calculate environmental virtual potential temperature at LCL,
                    # 'thv0lcl' which is solely used in the 'cin' calculation. Note
                    # that 'thv0lcl' is calculated first by calculating  'thl0lcl'
                    # and 'qt0lcl' at the LCL, and performing 'conden' afterward,
                    # in fully consistent with the other parts of the code.
                    thl0lcl = thl0.at(K=klcl) + ssthl0.at(K=klcl) * (
                        plcl - pmid0.at(K=klcl)
                    )
                    qt0lcl = qt0.at(K=klcl) + ssqt0.at(K=klcl) * (
                        plcl - pmid0.at(K=klcl)
                    )
                    thj, qvj, qlj, qij, qse, id_check = conden(
                        plcl, thl0lcl, qt0lcl, ese, esx
                    )

                    if id_check == 1:
                        condensation = True
                        umf_out[0, 0, 1] = 0.0
                        dcm_out = 0.0
                        qvten_out = 0.0
                        qlten_out = 0.0
                        qiten_out = 0.0
                        sten_out = 0.0
                        uten_out = 0.0
                        vten_out = 0.0
                        qrten_out = 0.0
                        qsten_out = 0.0
                        cufrc_out = 0.0
                        cush_inout = -1.0
                        qldet_out = 0.0
                        qidet_out = 0.0
                        qtflx_out[0, 0, 1] = 0.0
                        slflx_out[0, 0, 1] = 0.0
                        uflx_out[0, 0, 1] = 0.0
                        vflx_out[0, 0, 1] = 0.0
                        fer_out = constants.MAPL_UNDEF
                        fdr_out = constants.MAPL_UNDEF

    with computation(FORWARD), interval(...):
        if not condensation:
            thv0lcl = thj * (1.0 + zvir * qvj - qlj - qij)


def compute_cin_cinlcl(
    condensation: BoolFieldIJ,
    stop_cin: BoolFieldIJ,
    klcl: IntField,
    kinv: IntField,
    thvlsrc: FloatField,
    pifc0: FloatField,
    thv0bot: FloatField,
    thv0top: FloatField,
    plcl: FloatField,
    thv0lcl: FloatField,
    thlsrc: FloatField,
    qtsrc: FloatField,
    ese: FloatField_Extra_Dim,
    esx: FloatField_Extra_Dim,
    umf_out: FloatField,
    qtflx_out: FloatField,
    slflx_out: FloatField,
    uflx_out: FloatField,
    vflx_out: FloatField,
    cin_IJ: FloatFieldIJ,
    cinlcl_IJ: FloatFieldIJ,
    plfc_IJ: FloatFieldIJ,
    klfc_IJ: IntFieldIJ,
    plfc: FloatField,
    klfc: IntField,
    cin: FloatField,
    thvubot: FloatField,
    thvutop: FloatField,
    k0: Int,
    iteration: i32,
    rbuoy: Float,
    rkfre: FloatFieldIJ,
    tkeavg: FloatField,
    epsvarw: Float,
    thvlmin: FloatField,
    usrc: FloatField,
    vsrc: FloatField,
    dotransport: Int,
    trsrc: FloatFieldIJ_NTracers,
    trsrc_o: FloatFieldIJ_NTracers,
    cin_i: FloatFieldIJ,
    cinlcl_i: FloatFieldIJ,
    ke: FloatFieldIJ,
    kinv_o: IntField,
    klcl_o: IntField,
    klfc_o: IntField,
    plcl_o: FloatField,
    plfc_o: FloatField,
    tkeavg_o: FloatField,
    thvlmin_o: FloatField,
    qtsrc_o: FloatField,
    thvlsrc_o: FloatField,
    thlsrc_o: FloatField,
    usrc_o: FloatField,
    vsrc_o: FloatField,
    thv0lcl_o: FloatField,
    cinlcl: FloatField,
    cush_inout: FloatFieldIJ,
):
    """
    Compute Convective Inhibition, 'cin' & 'cinlcl' [J/kg]=[m2/s2] TKE unit.

    'cin' (cinlcl) is computed from the PBL top interface to LFC (LCL) using
    piecewisely reconstructed environmental profiles, assuming environmental
    buoyancy profile within each layer ( or from LCL to upper interface in
    each layer ) is simply a linear profile. For the purpose of cin (cinlcl)
    calculation, we simply assume that lateral entrainment does not occur in
    updrafting cumulus plume, i.e., cumulus source air property is conserved.

    Below explains some rules used in the calculations of cin (cinlcl).   In
    general, both 'cin' and 'cinlcl' are calculated from a PBL top interface
    to LCL and LFC, respectively :

    1. If LCL is lower than the PBL height, cinlcl = 0 and cin is calculated
    from PBL height to LFC.

    2. If LCL is higher than PBL height,   'cinlcl' is calculated by summing
    both positive and negative cloud buoyancy up to LCL using 'single_cin'
    From the LCL to LFC, however, only negative cloud buoyancy is counted
    to calculate final 'cin' upto LFC.

    3. If either 'cin' or 'cinlcl' is negative, they are set to be zero.

    In the below code, 'klfc' is the layer index containing 'LFC' similar to
    'kinv' and 'klcl'.
    """
    from __externals__ import ncnst

    with computation(FORWARD), interval(...):
        if iteration != i32(0):
            stop_cin = (
                False  # Indicates if CIN calcuation is done (e.g., Fortran go to 35)
            )
            cin_IJ = 0.0
            cinlcl_IJ = 0.0
            plfc_IJ = 0.0
            klfc_IJ = 0.0

    with computation(FORWARD), interval(1, -1):
        # Case 1. LCL height is higher than PBL interface ( 'pLCL <=ps0(kinv-1)' )
        cin = 0.0
        cinlcl = 0.0
        plfc = 0.0
        if not condensation and klcl >= kinv - 1 and not stop_cin:
            if not stop_cin and THIS_K >= kinv - 1 and THIS_K < klcl:
                thvubot = thvlsrc
                thvutop = thvlsrc
                cin = cin[0, 0, -1] + single_cin(
                    pifc0,
                    thv0bot,
                    pifc0[0, 0, 1],
                    thv0top,
                    thvubot,
                    thvutop,
                )

            elif not stop_cin and THIS_K == klcl:
                # ----- Bottom to LCL
                thvubot = thvlsrc
                thvutop = thvlsrc
                cin = cin[0, 0, -1] + single_cin(
                    pifc0,
                    thv0bot,
                    plcl,
                    thv0lcl,
                    thvubot,
                    thvutop,
                )
                cinlcl = max(cin, 0.0)
                cin = cinlcl

                # ----- LCL to Top
                (
                    thj,
                    qvj,
                    qlj,
                    qij,
                    qse,
                    id_check,
                ) = conden(
                    pifc0[0, 0, 1],
                    thlsrc,
                    qtsrc,
                    ese,
                    esx,
                )

                if id_check == 1:
                    condensation = True
                    umf_out[0, 0, 1] = 0.0
                    dcm_out = 0.0
                    qvten_out = 0.0
                    qlten_out = 0.0
                    qiten_out = 0.0
                    sten_out = 0.0
                    uten_out = 0.0
                    vten_out = 0.0
                    qrten_out = 0.0
                    qsten_out = 0.0
                    cufrc_out = 0.0
                    cush_inout = -1.0
                    qldet_out = 0.0
                    qidet_out = 0.0
                    qtflx_out[0, 0, 1] = 0.0
                    slflx_out[0, 0, 1] = 0.0
                    uflx_out[0, 0, 1] = 0.0
                    vflx_out[0, 0, 1] = 0.0
                    fer_out = constants.MAPL_UNDEF
                    fdr_out = constants.MAPL_UNDEF

                if not condensation:
                    thvutop = thj * (1.0 + zvir * qvj - qlj - qij)

                    plfc, cin = getbuoy(
                        plcl,
                        thv0lcl,
                        pifc0[0, 0, 1],
                        thv0top,
                        thvubot,
                        thvutop,
                        cin,
                        plfc,
                    )

                    if plfc > 0.0:
                        klfc = THIS_K
                        stop_cin = True

            else:
                if not condensation and THIS_K > klcl and not stop_cin:
                    thvubot = thvutop[0, 0, -1]
                    (
                        thj,
                        qvj,
                        qlj,
                        qij,
                        qse,
                        id_check,
                    ) = conden(
                        pifc0[0, 0, 1],
                        thlsrc,
                        qtsrc,
                        ese,
                        esx,
                    )

                    if id_check == 1:
                        condensation = True
                        umf_out[0, 0, 1] = 0.0
                        dcm_out = 0.0
                        qvten_out = 0.0
                        qlten_out = 0.0
                        qiten_out = 0.0
                        sten_out = 0.0
                        uten_out = 0.0
                        vten_out = 0.0
                        qrten_out = 0.0
                        qsten_out = 0.0
                        cufrc_out = 0.0
                        cush_inout = -1.0
                        qldet_out = 0.0
                        qidet_out = 0.0
                        qtflx_out[0, 0, 1] = 0.0
                        slflx_out[0, 0, 1] = 0.0
                        uflx_out[0, 0, 1] = 0.0
                        vflx_out[0, 0, 1] = 0.0
                        fer_out = constants.MAPL_UNDEF
                        fdr_out = constants.MAPL_UNDEF

                    if not condensation and not stop_cin:
                        thvutop = thj * (1.0 + zvir * qvj - qlj - qij)
                        plfc, cin = getbuoy(
                            pifc0,
                            thv0bot,
                            pifc0[0, 0, 1],
                            thv0top,
                            thvubot,
                            thvutop,
                            cin[0, 0, -1],
                            plfc[0, 0, -1],
                        )

                        if plfc > 0.0:
                            klfc = THIS_K
                            stop_cin = True

    with computation(FORWARD), interval(1, -1):
        # Case 2. LCL height is lower than PBL interface ( 'pLCL > ps0(kinv-1)')
        if not condensation and klcl < kinv - 1 and not stop_cin:
            if not stop_cin and THIS_K >= kinv - 1:
                (
                    thj,
                    qvj,
                    qlj,
                    qij,
                    qse,
                    id_check,
                ) = conden(
                    pifc0,
                    thlsrc,
                    qtsrc,
                    ese,
                    esx,
                )

                if id_check == 1 and not stop_cin:
                    condensation = True
                    umf_out[0, 0, 1] = 0.0
                    dcm_out = 0.0
                    qvten_out = 0.0
                    qlten_out = 0.0
                    qiten_out = 0.0
                    sten_out = 0.0
                    uten_out = 0.0
                    vten_out = 0.0
                    qrten_out = 0.0
                    qsten_out = 0.0
                    cufrc_out = 0.0
                    cush_inout = -1.0
                    qldet_out = 0.0
                    qidet_out = 0.0
                    qtflx_out[0, 0, 1] = 0.0
                    slflx_out[0, 0, 1] = 0.0
                    uflx_out[0, 0, 1] = 0.0
                    vflx_out[0, 0, 1] = 0.0
                    fer_out = constants.MAPL_UNDEF
                    fdr_out = constants.MAPL_UNDEF

                if not condensation and not stop_cin:
                    thvubot = thj * (1.0 + zvir * qvj - qlj - qij)
                    (
                        thj,
                        qvj,
                        qlj,
                        qij,
                        qse,
                        id_check,
                    ) = conden(
                        pifc0[0, 0, 1],
                        thlsrc,
                        qtsrc,
                        ese,
                        esx,
                    )

                if id_check == 1 and not stop_cin:
                    condensation = True
                    umf_out[0, 0, 1] = 0.0
                    dcm_out = 0.0
                    qvten_out = 0.0
                    qlten_out = 0.0
                    qiten_out = 0.0
                    sten_out = 0.0
                    uten_out = 0.0
                    vten_out = 0.0
                    qrten_out = 0.0
                    qsten_out = 0.0
                    cufrc_out = 0.0
                    cush_inout = -1.0
                    qldet_out = 0.0
                    qidet_out = 0.0
                    qtflx_out[0, 0, 1] = 0.0
                    slflx_out[0, 0, 1] = 0.0
                    uflx_out[0, 0, 1] = 0.0
                    vflx_out[0, 0, 1] = 0.0
                    fer_out = constants.MAPL_UNDEF
                    fdr_out = constants.MAPL_UNDEF

                if not condensation and not stop_cin:
                    thvutop = thj * (1.0 + zvir * qvj - qlj - qij)

                    plfc, cin = getbuoy(
                        pifc0,
                        thv0bot,
                        pifc0[0, 0, 1],
                        thv0top,
                        thvubot,
                        thvutop,
                        cin[0, 0, -1],
                        plfc[0, 0, -1],
                    )

                    if plfc > 0.0:
                        klfc = THIS_K
                        stop_cin = True

    with computation(FORWARD), interval(1, None):
        if not condensation:
            # Store cin, cinlcl, plfc, and klfc as 2D fields
            if cin == 0.0 and cin[0, 0, -1] != 0.0:
                cin_IJ = cin[0, 0, -1]
                cinlcl_IJ = cinlcl.at(K=klcl)
            if plfc == 0.0 and plfc[0, 0, -1] != 0.0:
                plfc_IJ = plfc[0, 0, -1]
            if klfc == 0.0 and klfc[0, 0, -1] != 0.0:
                klfc_IJ = klfc[0, 0, -1]

    with computation(FORWARD), interval(...):
        if not condensation:
            if klfc_IJ >= k0 - 1:
                klfc_IJ = k0 - 1
                condensation = True
                umf_out[0, 0, 1] = 0.0
                dcm_out = 0.0
                qvten_out = 0.0
                qlten_out = 0.0
                qiten_out = 0.0
                sten_out = 0.0
                uten_out = 0.0
                vten_out = 0.0
                qrten_out = 0.0
                qsten_out = 0.0
                cufrc_out = 0.0
                cush_inout = -1.0
                qldet_out = 0.0
                qidet_out = 0.0
                qtflx_out[0, 0, 1] = 0.0
                slflx_out[0, 0, 1] = 0.0
                uflx_out[0, 0, 1] = 0.0
                vflx_out[0, 0, 1] = 0.0
                fer_out = constants.MAPL_UNDEF
                fdr_out = constants.MAPL_UNDEF

    with computation(FORWARD), interval(...):
        if not condensation:
            # In order to calculate implicit 'cin' (or 'cinlcl'), save the initially
            # calculated 'cin' and 'cinlcl', and other related variables. These will
            # be restored after calculating implicit CIN.

            if iteration == i32(0):
                cin_i = cin_IJ
                cinlcl_i = cinlcl_IJ
                ke = rbuoy / (rkfre * tkeavg + epsvarw)
                kinv_o = kinv - 1
                klcl_o = klcl
                klfc_o = klfc_IJ
                plcl_o = plcl
                plfc_o = plfc_IJ
                tkeavg_o = tkeavg
                thvlmin_o = thvlmin
                qtsrc_o = qtsrc
                thvlsrc_o = thvlsrc
                thlsrc_o = thlsrc
                usrc_o = usrc
                vsrc_o = vsrc
                thv0lcl_o = thv0lcl

                if dotransport == 1:
                    n = 0
                    while n < ncnst:
                        trsrc_o[0, 0][n] = trsrc[0, 0][n]
                        n += 1


def avg_initial_and_final_cin(
    condensation: BoolFieldIJ,
    iteration: i32,
    cin_IJ: FloatFieldIJ,
    cinlcl_IJ: FloatFieldIJ,
    use_CINcin: i32,
    cin_i: FloatFieldIJ,
    cinlcl_i: FloatFieldIJ,
    ke: FloatFieldIJ,
    dotransport: Int,
    trsrc: FloatFieldIJ_NTracers,
    trsrc_o: FloatFieldIJ_NTracers,
    tr0: FloatField_NTracers,
    tr0_o: FloatField_NTracers,
    sstr0: FloatField_NTracers,
    sstr0_o: FloatField_NTracers,
    kinv_o: IntField,
    klcl_o: IntField,
    klfc_o: IntField,
    plcl_o: FloatField,
    plfc_o: FloatField,
    tkeavg_o: FloatField,
    thvlmin_o: FloatField,
    qtsrc_o: FloatField,
    thvlsrc_o: FloatField,
    thlsrc_o: FloatField,
    usrc_o: FloatField,
    vsrc_o: FloatField,
    thv0lcl_o: FloatField,
    qv0_o: FloatField,
    ql0_o: FloatField,
    qi0_o: FloatField,
    t0_o: FloatField,
    s0_o: FloatField,
    u0_o: FloatField,
    v0_o: FloatField,
    qt0_o: FloatField,
    thl0_o: FloatField,
    thvl0_o: FloatField,
    ssthl0_o: FloatField,
    ssqt0_o: FloatField,
    thv0bot_o: FloatField,
    thv0top_o: FloatField,
    thvl0bot_o: FloatField,
    thvl0top_o: FloatField,
    ssu0_o: FloatField,
    ssv0_o: FloatField,
    thvlmin_IJ: FloatFieldIJ,
    umf_zint: FloatField,
    emf: FloatField,
    slflx: FloatField,
    qtflx: FloatField,
    uflx: FloatField,
    vflx: FloatField,
    k0: Int,
    ufrc: FloatField,
    thlu: FloatField,
    qtu: FloatField,
    uu: FloatField,
    vu: FloatField,
    wu: FloatField,
    thvu: FloatField,
    thlu_emf: FloatField,
    qtu_emf: FloatField,
    uu_emf: FloatField,
    vu_emf: FloatField,
    trflx: FloatField_NTracers,
    trten: FloatField_NTracers,
    tru: FloatField_NTracers,
    tru_emf: FloatField_NTracers,
    umf_s: FloatField,
    zifc0: FloatField,
    dcm_s: FloatField,
    qvten_s: FloatField,
    qlten_s: FloatField,
    qiten_s: FloatField,
    sten_s: FloatField,
    uten_s: FloatField,
    vten_s: FloatField,
    qrten_s: FloatField,
    qsten_s: FloatField,
    qldet_s: FloatField,
    qidet_s: FloatField,
    qlsub_s: FloatField,
    qisub_s: FloatField,
    cush_s: FloatField,
    cufrc_s: FloatField,
    qtflx_out: FloatField,
    qtflx_s: FloatField,
    slflx_out: FloatField,
    slflx_s: FloatField,
    uflx_out: FloatField,
    uflx_s: FloatField,
    vflx_out: FloatField,
    vflx_s: FloatField,
    fer_s: FloatField,
    fdr_s: FloatField,
    umf_out: FloatField,
    kinv: IntField,
    klcl: IntField,
    plcl: FloatField,
    thv0bot: FloatField,
    thv0lcl: FloatField,
    thv0top: FloatField,
    thlsrc: FloatField,
    qtsrc: FloatField,
    thl0: FloatField,
    ssthl0: FloatField,
    ssqt0: FloatField,
    ssu0: FloatField,
    ssv0: FloatField,
    qt0: FloatField,
    u0: FloatField,
    v0: FloatField,
    qi0: FloatField,
    ql0: FloatField,
    qv0: FloatField,
    s0: FloatField,
    dcm: FloatField,
    qvten_out: FloatField,
    dcm_out: FloatField,
    qlten_out: FloatField,
    qiten_out: FloatField,
    del_CIN: FloatFieldIJ,
    sten_out: FloatField,
    uten_out: FloatField,
    vten_out: FloatField,
    qrten_out: FloatField,
    qsten_out: FloatField,
    cufrc_out: FloatField,
    cush_inout: FloatFieldIJ,
    qldet_out: FloatField,
    qidet_out: FloatField,
    qlsub_out: FloatField,
    qisub_out: FloatField,
    fdr_out: FloatField,
    fer_out: FloatField,
):
    """
    Calculate implicit 'cin' by averaging initial and final cins. Note that
    implicit CIN is adopted only when cumulus convection stabilized the system,
    i.e., only when 'del_CIN >0'. If 'del_CIN<=0', just use explicit CIN. Note
    also that since 'cinlcl' is set to zero whenever LCL is below the PBL top,
    (see above CIN calculation part), the use of 'implicit CIN=cinlcl'  is not
    good. Thus, when using implicit CIN, always try to only use 'implicit CIN=
    cin', not 'implicit CIN=cinlcl'. However, both 'CIN=cin' and 'CIN=cinlcl'
    are good when using explicit CIN.
    """
    from __externals__ import ncnst

    with computation(FORWARD), interval(1, None):
        if not condensation:
            if thvlmin_o == 0.0 and thvlmin_o[0, 0, -1] != 0.0:
                thvlmin_IJ = thvlmin_o[0, 0, -1]

    with computation(FORWARD), interval(...):
        if not condensation:
            cin_f = cin_IJ
            cinlcl_f = cinlcl_IJ

            if use_CINcin == 1:
                del_CIN = cin_f - cin_i
            else:
                del_CIN = cinlcl_f - cinlcl_i

    with computation(FORWARD), interval(...):
        if not condensation:
            # Calculate implicit 'cin' by averaging initial and final cins.    Note that
            # implicit CIN is adopted only when cumulus convection stabilized the system
            # i.e., only when 'del_CIN >0'. If 'del_CIN<=0', just use explicit CIN. Note
            # also that since 'cinlcl' is set to zero whenever LCL is below the PBL top,
            # (see above CIN calculation part), the use of 'implicit CIN=cinlcl'  is not
            # good. Thus, when using implicit CIN, always try to only use 'implicit CIN=
            # cin', not 'implicit CIN=cinlcl'. However, both 'CIN=cin' and 'CIN=cinlcl'
            # are good when using explicit CIN.

            if del_CIN > 0.0:
                # Calculate implicit 'cin' and 'cinlcl'. Note that when we chose
                # to use 'implicit CIN = cin', choose 'cinlcl = cinlcl_i' below:
                # because iterative CIN only aims to obtain implicit CIN,  once
                # we obtained 'implicit CIN=cin', it is good to use the original
                # profiles information for all the other variables after that.
                # Note 'cinlcl' will be explicitly used in calculating  'wlcl' &
                # 'ufrclcl' after calculating 'winv' & 'ufrcinv'  at the PBL top
                # interface later, after calculating 'cbmf'.
                alpha = compute_alpha(del_CIN, ke)
                cin_IJ = cin_i + alpha * del_CIN
                if use_CINcin == 1:
                    cinlcl_IJ = cinlcl_i

                # Restore the original values from the previous 'iter_cin' step (1)
                # to compute correct tendencies for (n+1) time step by implicit CIN

                kinv = kinv_o
                klcl = klcl_o
                klfc = klfc_o
                plcl = plcl_o
                plfc = plfc_o
                tkeavg = tkeavg_o
                thvlmin = thvlmin_IJ
                qtsrc = qtsrc_o
                thvlsrc = thvlsrc_o
                thlsrc = thlsrc_o
                usrc = usrc_o
                vsrc = vsrc_o
                thv0lcl = thv0lcl_o

                if dotransport == 1:
                    n = 0
                    while n < ncnst:
                        trsrc[0, 0][n] = trsrc_o[0, 0][n]
                        n += 1

                qv0 = qv0_o
                ql0 = ql0_o
                qi0 = qi0_o
                t0 = t0_o
                s0 = s0_o
                u0 = u0_o
                v0 = v0_o
                qt0 = qt0_o
                thl0 = thl0_o
                thvl0 = thvl0_o
                ssthl0 = ssthl0_o
                ssqt0 = ssqt0_o
                thv0bot = thv0bot_o
                thv0top = thv0top_o
                thvl0bot = thvl0bot_o
                thvl0top = thvl0top_o
                ssu0 = ssu0_o
                ssv0 = ssv0_o

                if dotransport == 1:
                    n = 0
                    while n < ncnst:
                        tr0[0, 0, 0][n] = tr0_o[0, 0, 0][n]
                        sstr0[0, 0, 0][n] = sstr0_o[0, 0, 0][n]

                        n += 1

                # Initialize all fluxes, tendencies, and other variables
                # in association with cumulus convection.
                umf_zint[0, 0, 1] = 0.0
                dcm = 0.0
                emf[0, 0, 1] = 0.0
                slflx[0, 0, 1] = 0.0
                qtflx[0, 0, 1] = 0.0
                uflx[0, 0, 1] = 0.0
                vflx[0, 0, 1] = 0.0
                qvten = 0.0
                qlten = 0.0
                qiten = 0.0
                sten = 0.0
                uten = 0.0
                vten = 0.0
                qrten = 0.0
                qsten = 0.0
                dwten = 0.0
                diten = 0.0
                cufrc = 0.0
                qcu = 0.0
                qlu = 0.0
                qiu = 0.0
                fer = 0.0
                fdr = 0.0
                xco = 0.0
                qc = 0.0
                # qldet        = 0.0
                # qidet        = 0.0
                qlten_sub = 0.0
                qiten_sub = 0.0
                qc_l = 0.0
                qc_i = 0.0
                cbmf = 0.0
                cnt = k0
                cnb = 0.0
                qtten = 0.0
                slten = 0.0
                ufrc[0, 0, 1] = 0.0

                thlu[0, 0, 1] = constants.MAPL_UNDEF
                qtu[0, 0, 1] = constants.MAPL_UNDEF
                uu[0, 0, 1] = constants.MAPL_UNDEF
                vu[0, 0, 1] = constants.MAPL_UNDEF
                wu[0, 0, 1] = constants.MAPL_UNDEF
                thvu[0, 0, 1] = constants.MAPL_UNDEF
                thlu_emf[0, 0, 1] = constants.MAPL_UNDEF
                qtu_emf[0, 0, 1] = constants.MAPL_UNDEF
                uu_emf[0, 0, 1] = constants.MAPL_UNDEF
                vu_emf[0, 0, 1] = constants.MAPL_UNDEF

                if dotransport == 1:
                    n = 0
                    while n < ncnst:
                        trflx[0, 0, 1][n] = 0.0
                        trten[0, 0, 0][n] = 0.0
                        tru[0, 0, 1][n] = 0.0
                        tru_emf[0, 0, 1][n] = 0.0
                        n += 1

                # Below are diagnostic output variables for detailed
                # analysis of cumulus scheme.
                ufrclcl = 0.0
                wlcl = 0.0
                cbmflimit = 0.0

            else:  # When 'del_CIN < 0', use explicit CIN instead of implicit CIN.

                # Identifier showing whether explicit or implicit CIN is used
                ind_delcin = 1.0

                # Restore original output values of "iter_cin = 1" and exit

                umf_out[0, 0, 1] = umf_s[0, 0, 1]
                if THIS_K >= 0 and THIS_K <= (kinv - 1):
                    umf_out = umf_s.at(K=kinv - 1) * zifc0 / zifc0.at(K=kinv - 1)

                dcm_out = dcm_s
                qvten_out = qvten_s
                qlten_out = qlten_s
                qiten_out = qiten_s
                sten_out = sten_s
                uten_out = uten_s
                vten_out = vten_s
                qrten_out = qrten_s
                qsten_out = qsten_s
                qldet_out = qldet_s
                qidet_out = qidet_s
                qlsub_out = qlsub_s
                qisub_out = qisub_s
                cush_inout = cush_s
                cufrc_out = cufrc_s
                qtflx_out[0, 0, 1] = qtflx_s[0, 0, 1]
                slflx_out[0, 0, 1] = slflx_s[0, 0, 1]
                uflx_out[0, 0, 1] = uflx_s[0, 0, 1]
                vflx_out[0, 0, 1] = vflx_s[0, 0, 1]

                # Below are diagnostic output variables for detailed analysis of cumulus
                # scheme.
                # The order of vertical index is reversed for this internal diagnostic
                # output.
                fer_out = fer_s
                fdr_out = fdr_s

    with computation(FORWARD), interval(...):
        if not condensation:
            if del_CIN <= 0.0:
                condensation = True  # Done computing at this column


def define_prel_krel(
    condensation: BoolFieldIJ,
    iteration: i32,
    klcl: IntField,
    kinv: IntField,
    pifc0: FloatField,
    thv0bot: FloatField,
    plcl: FloatField,
    thv0lcl: FloatField,
    krel: IntField,
    prel: FloatField,
    thv0rel: FloatField,
    umf_out: FloatField,
):
    with computation(FORWARD), interval(...):
        if not condensation:
            # Define a release level, 'prel' and release layer, 'krel'.
            # 'prel' is the lowest level from which buoyancy sorting occurs, and
            # 'krel' is the layer index containing 'prel' in it, similar to  the
            # previous definitions of 'kinv', 'klcl', and 'klfc'.    In order to
            # ensure that only PBL scheme works within the PBL,  if LCL is below
            # PBL top height, then 'krel = kinv', while if LCL is above  PBL top
            # height, then 'krel = klcl'.   Note however that regardless of  the
            # definition of 'krel', cumulus convection induces fluxes within PBL
            # through 'fluxbelowinv'.  We can make cumulus convection start from
            # any level, even within the PBL by appropriately defining 'krel'  &
            # 'prel' here. Then it must be accompanied by appropriate definition
            # of source air properties, CIN, and re-setting of 'fluxbelowinv', &
            # many other stuffs.
            # Note that even when 'prel' is located above the PBL top height, we
            # still have cumulus convection between PBL top height and 'prel':
            # we simply assume that no lateral mixing occurs in this range.
            if iteration != i32(0):
                kinv = kinv + 1  # Adjust kinv

            if klcl < kinv - 1:
                krel = kinv - 1
                prel = pifc0.at(K=krel)
                thv0rel = thv0bot.at(K=krel)
            else:
                krel = klcl
                prel = plcl
                thv0rel = thv0lcl


def calc_cumulus_base_mass_flux(
    condensation: BoolFieldIJ,
    iteration: i32,
    use_CINcin: i32,
    cin_IJ: FloatFieldIJ,
    rbuoy: Float,
    cinlcl_IJ: FloatFieldIJ,
    rkfre: FloatFieldIJ,
    tkeavg: FloatField,
    epsvarw: Float,
    umf_out: FloatField,
    qtflx_out: FloatField,
    slflx_out: FloatField,
    uflx_out: FloatField,
    vflx_out: FloatField,
    kinv: IntField,
    pifc0: FloatField,
    thv0top: FloatField,
    exnifc0: FloatField,
    dp0: FloatField,
    dt: Float,
    mumin1: Float,
    rmaxfrac: Float,
    winv: FloatField,
    cbmf: FloatField,
    rho0inv: FloatField,
    ufrcinv: FloatField,
    wcrit: FloatFieldIJ,
    cush_inout: FloatFieldIJ,
):
    """
    Calculate cumulus base mass flux ('cbmf'), fractional area ('ufrcinv'),
    and mean vertical velocity (winv) of cumulus updraft at PBL top interface
    Also, calculate updraft fractional area (ufrclcl) and vertical velocity at
    the LCL (wlcl). When LCL is below PBLH, cinlcl = 0 and 'ufrclcl = ufrcinv'
    and 'wlcl = winv. Only updrafts strong enough to overcome CIN can rise
    over PBL top interface.

    Thus, in order to calculate cumulus mass flux at PBL top interface, 'cbmf'
    we need to know 'CIN' ( the strength of potential energy barrier ) and
    'sigmaw' (a standard deviation of updraft vertical velocity at the PBL top
    interface, a measure of turbulentce strength in the PBL ). Naturally, the
    ratio of these two variables, 'mu' - normalized CIN by TKE- is key
    variable controlling 'cbmf'.  If 'mu' becomes large, only small fraction
    of updrafts with very strong TKE can rise over the PBL - both 'cbmf' and
    'ufrc' becomes small, but 'winv' becomes large ( this can be easily
    understood by PDF of w at PBL top ).  If 'mu' becomes small, lots of
    updraft can rise over the PBL top - both 'cbmf' and 'ufrc' becomes large,
    but 'winv' becomes small. Thus, all of the key variables associated with
    cumulus convection  at the PBL top - 'cbmf', 'ufrc', 'winv' where 'cbmf =
    rho*ufrc*winv' - are a unique functions of 'mu', normalized CIN. Although
    these are uniquely determined by 'mu', we usually impose two comstraints
    on 'cbmf' and 'ufrc': (1) because we will simply assume that subsidence
    warming and drying of 'kinv-1' layer in association with 'cbmf' at PBL top
    interface is confined only in 'kinv-1' layer, cbmf must not be larger than
    the mass within the 'kinv-1' layer. Otherwise, instability will occur due
    to the breaking of stability con. If we consider semi-Lagrangian vertical
    advection scheme and explicitly consider the extent of vertical movement
    of each layer in association with cumulus mass flux, we don't need to
    impose this constraint. However,  using a  semi-Lagrangian scheme is a
    future research subject. Note that this constraint should be applied for
    all interfaces above PBL top as well as PBL top interface.   As a result,
    this 'cbmf' constraint impose a 'lower' limit on mu - 'mumin0'. (2) in
    order for mass flux parameterization - rho*(w'a')= M*(a_c-a_e) - to be
    valid, cumulus updraft fractional area should be much smaller than 1. In
    current code, we impose 'rmaxfrac = 0.1 ~ 0.2'   through the whole
    vertical layers where cumulus convection occurs. At the PBL top interface,
    the same constraint is made by imposing another lower 'lower' limit on mu,
    'mumin1'. After that, also limit 'ufrclcl' to be smaller than 'rmaxfrac'
    by 'mumin2'.

    Calculate normalized CIN, 'mu' satisfying all the three constraints
    imposed on 'cbmf'('mumin0'), 'ufrc' at the PBL top - 'ufrcinv' - ( by
    'mumin1' from a parameter sentence), and 'ufrc' at the LCL - 'ufrclcl'
    ( by 'mumin2'). Note that 'cbmf' does not change between PBL top and LCL
    because we assume that buoyancy sorting does not occur when cumulus
    updraft is unsaturated.

    Inputs:
    iteration [i32]: Iteration of implicit CIN loop (i.e., 0 or 1)
    use_CINcin [i32]: if true, calc CIN thru L..
    cin_IJ [FloatFieldIJ]: Convective inhibition [J/kg]
    rbuoy [Float]: Non-hydro pressure effect on updraft
    cinlcl_IJ [FloatFieldIJ]: Convective inhibition at LCL [J/kg]
    rkfre [FloatFieldIJ]: Resolution dependent vertical velocity variance as fraction of tke
    tkeavg [FloatField]: Average tke over the PBL [m2/s2]
    epsvarw [Float]: Variance of PBL w by mesoscale
    kinv [IntField]: Inversion layer with PBL top interface as lower interface
    pifc0 [FloatField]: Environmental pressure at the interfaces [Pa]
    thv0top [FloatField]: Potential temperature at top [?]
    exnifc0 [FloatField]: Exner function at interfaces
    dp0 [FloatField]: Environmental layer pressure thickness [Pa] > 0
    dt [Float]: Timestep [s]
    mumin1 [Float]: Normalized CIN ('mu') corresponding to 'rmaxfrac' at the PBL top
    rmaxfrac [Float]: Maximum core updraft fraction

    Inouts:
    condensation [BoolFieldIJ]: Mask that indicates if condensation has occurred

    Outputs:
    umf_out [FloatField]: Updraft mass flux at the interfaces [kg/m2/s]
    qtflx_out [FloatField]: Mixing ratio flux [?]
    slflx_out [FloatField]: Sensible heat flux [?]
    uflx_out [FloatField]: Zonal wind flux [?]
    vflx_out [FloatField]: Meridional wind flux [?]
    winv [FloatField]: Mean vertical velocity of cumulus updraft at PBL top interface [m/s]
    cbmf [FloatField]: Cloud base mass flux [kg/m2/s]
    rho0inv [FloatField]: Density of water [?]
    ufrcinv [FloatField]: Fractional area [unitless]
    wcrit [FloatFieldIJ]: Critical mixing ratio [?]
    cush_inout [FloatFieldIJ]: Convective scale height [m]
    """
    with computation(FORWARD), interval(...):
        if not condensation:
            if use_CINcin == i32(1):
                wcrit = sqrt(2.0 * cin_IJ * rbuoy)
            else:
                wcrit = sqrt(2.0 * cinlcl_IJ * rbuoy)

            sigmaw = sqrt(rkfre * tkeavg + epsvarw)
            mu = wcrit / sigmaw / 1.4142

            if mu >= 3.0:
                condensation = True
                umf_out[0, 0, 1] = 0.0
                dcm_out = 0.0
                qvten_out = 0.0
                qlten_out = 0.0
                qiten_out = 0.0
                sten_out = 0.0
                uten_out = 0.0
                vten_out = 0.0
                qrten_out = 0.0
                qsten_out = 0.0
                cufrc_out = 0.0
                cush_inout = -1.0
                qldet_out = 0.0
                qidet_out = 0.0
                qtflx_out[0, 0, 1] = 0.0
                slflx_out[0, 0, 1] = 0.0
                uflx_out[0, 0, 1] = 0.0
                vflx_out[0, 0, 1] = 0.0
                fer_out = constants.MAPL_UNDEF
                fdr_out = constants.MAPL_UNDEF

            if not condensation:
                kbelow = kinv - 1
                rho0inv = pifc0.at(K=kbelow) / (
                    constants.MAPL_RDRY
                    * thv0top.at(K=kbelow - 1)
                    * exnifc0.at(K=kbelow)
                )
                cbmf = (rho0inv * sigmaw / 2.5066) * exp(-(mu**2))

                # 1. 'cbmf' constraint
                cbmflimit = 0.9 * dp0.at(K=kbelow - 1) / constants.MAPL_GRAV / dt
                mumin0 = 0.0
                if cbmf > cbmflimit:
                    mumin0 = sqrt(-log(2.5066 * cbmflimit / rho0inv / sigmaw))

                # # 2. 'ufrcinv' constraint
                mu = max(max(mu, mumin0), mumin1)

                # # 3. 'ufrclcl' constraint
                mulcl = sqrt(2.0 * cinlcl_IJ * rbuoy) / 1.4142 / sigmaw
                mulclstar = sqrt(
                    max(
                        0.0,
                        2.0
                        * (exp(-(mu**2)) / 2.5066) ** 2
                        * (1.0 / erfc(mu) ** 2 - 0.25 / rmaxfrac**2),
                    )
                )

                if mulcl > 1.0e-8 and mulcl > mulclstar:
                    mumin2 = compute_mumin2(mulcl, rmaxfrac, mu)
                    mu = max(mu, mumin2)

    with computation(FORWARD), interval(...):
        if not condensation:
            # Calculate final ['cbmf','ufrcinv','winv'] at the PBL top interface.
            # Note that final 'cbmf' here is obtained in such that 'ufrcinv' and
            # 'ufrclcl' are smaller than ufrcmax with no instability.

            cbmf = rkfre * (rho0inv * sigmaw / 2.5066) * exp((-(mu**2)))
            winv = sigmaw * (2.0 / 2.5066) * exp(-(mu**2)) / erfc(mu)
            ufrcinv = cbmf / winv / rho0inv


def define_updraft_properties(
    condensation: BoolFieldIJ,
    iteration: i32,
    winv: FloatField,
    cinlcl_IJ: FloatFieldIJ,
    rbuoy: Float,
    umf_out: FloatField,
    qtflx_out: FloatField,
    slflx_out: FloatField,
    uflx_out: FloatField,
    vflx_out: FloatField,
    cbmf: FloatField,
    rho0inv: FloatField,
    krel: IntField,
    ufrc: FloatField,
    ufrcinv: FloatField,
    kinv: IntField,
    umf_zint: FloatField,
    wu: FloatField,
    emf: FloatField,
    thlu: FloatField,
    qtu: FloatField,
    thlsrc: FloatField,
    qtsrc: FloatField,
    prel: FloatField,
    ese: FloatField_Extra_Dim,
    esx: FloatField_Extra_Dim,
    thvu: FloatField,
    wlcl: FloatField,
    ufrclcl: FloatField,
    cush_inout: FloatFieldIJ,
):
    """
    Calculate ['ufrclcl','wlcl'] at the LCL. When LCL is below PBL top,
    it automatically becomes 'ufrclcl = ufrcinv' & 'wlcl = winv', since
    it was already set to 'cinlcl=0' if LCL is below PBL top interface.
    Note 'cbmf' at the PBL top is the same as 'cbmf' at the LCL.  Note
    also that final 'cbmf' here is obtained in such that 'ufrcinv' and
    'ufrclcl' are smaller than ufrcmax and there is no instability.
    By construction, it must be 'wlcl > 0' but for assurance, I checked
    this again in the below block. If 'ufrclcl < 0.1%', just exit.

    Inputs:
    iteration [i32]: Iteration of implicit CIN loop (i.e., 0 or 1)
    winv [FloatField]: Mean vertical velocity of cumulus updraft at PBL top interface [m/s]
    cinlcl_IJ [FloatFieldIJ]: Convective inhibition at LCL [J/kg]
    rbuoy [Float]: Non-hydro pressure effect on updraft
    cbmf [FloatField]: Cloud base mass flux [kg/m2/s]
    rho0inv [FloatField]: Density of water [?]
    krel [IntField]: Release layer where buoyancy sorting first occurs
    kinv [IntField]: Inversion layer with PBL top interface as lower interface
    thlsrc [FloatField]: Temperature of cumulus source air [K] [?]
    qtsrc [FloatField]: Mixing ratio of cumulus source air [?]
    prel [FloatField]: Lowest level from which buoyancy sorting occurs [Pa]
    ese [FloatField_Extra_Dim]: Used in QSat_Float [n/a]
    esx [FloatField_Extra_Dim]: Used in QSat_Float [n/a]

    Inouts:
    condensation [BoolFieldIJ]: Mask that indicates if condensation has occurred

    Outputs:
    umf_out [FloatField]: Updraft mass flux at the interfaces [kg/m2/s]
    qtflx_out [FloatField]: Mixing ratio flux [?]
    slflx_out [FloatField]: Sensible heat flux [?]
    uflx_out [FloatField]: Zonal wind flux [?]
    vflx_out [FloatField]: Meridional wind flux [?]
    ufrc [FloatField]: Cumulus updraft fraction [fraction]
    umf_zint [FloatField]: Updraft mass flux at the interfaces [kg/m2/s]
    wu [FloatField]: Updraft vertical velocity at the interface [m/s]
    emf [FloatField]: Penetrative [?]
    thlu [FloatField]: Updraft liquid potential temperature at the interface [K]
    qtu [FloatField]: Updraft total specific humidity at the interface [kg/kg]
    thvu [FloatField]: Updraft virtual potential temperature at the interface [m/s]
    wlcl [FloatField]: Vertical velocity at LCL [m/s]
    ufrclcl [FloatField]: Updraft fractional area [no unit]
    cush_inout [FloatFieldIJ]: Convective scale height [m]
    """
    with computation(FORWARD), interval(...):
        if not condensation:
            wtw = (winv * winv) - (2.0 * cinlcl_IJ * rbuoy)

            if wtw <= 0.0:
                condensation = True
                umf_out[0, 0, 1] = 0.0
                dcm_out = 0.0
                qvten_out = 0.0
                qlten_out = 0.0
                qiten_out = 0.0
                sten_out = 0.0
                uten_out = 0.0
                vten_out = 0.0
                qrten_out = 0.0
                qsten_out = 0.0
                cufrc_out = 0.0
                cush_inout = -1.0
                qldet_out = 0.0
                qidet_out = 0.0
                qtflx_out[0, 0, 1] = 0.0
                slflx_out[0, 0, 1] = 0.0
                uflx_out[0, 0, 1] = 0.0
                vflx_out[0, 0, 1] = 0.0
                fer_out = constants.MAPL_UNDEF
                fdr_out = constants.MAPL_UNDEF

            if not condensation:

                wlcl = sqrt(wtw)
                ufrclcl = cbmf / wlcl / rho0inv
                wrel = wlcl

                if ufrclcl <= 0.0001:
                    condensation = True
                    umf_out[0, 0, 1] = 0.0
                    dcm_out = 0.0
                    qvten_out = 0.0
                    qlten_out = 0.0
                    qiten_out = 0.0
                    sten_out = 0.0
                    uten_out = 0.0
                    vten_out = 0.0
                    qrten_out = 0.0
                    qsten_out = 0.0
                    cufrc_out = 0.0
                    cush_inout = -1.0
                    qldet_out = 0.0
                    qidet_out = 0.0
                    qtflx_out[0, 0, 1] = 0.0
                    slflx_out[0, 0, 1] = 0.0
                    uflx_out[0, 0, 1] = 0.0
                    vflx_out[0, 0, 1] = 0.0
                    fer_out = constants.MAPL_UNDEF
                    fdr_out = constants.MAPL_UNDEF

                if not condensation:
                    if THIS_K == (krel - 1):
                        ufrc[0, 0, 1] = ufrclcl

    with computation(FORWARD), interval(...):
        if not condensation:
            # Define updraft properties at the level where buoyancy sorting starts to be
            # happening, i.e., by definition, at 'prel' level within the release layer.
            # Because no lateral entrainment occurs upto 'prel', conservative scalars of
            # cumulus updraft at release level is same as those of source air.  However,
            # horizontal momentums of source air are modified by horizontal PGF forcings
            # from PBL top interface to 'prel'.  For this case, we should add additional
            # horizontal momentum from PBL top interface to 'prel' as will be done below
            # to 'usrc' and 'vsrc'. Note that below cumulus updraft properties - umf, wu
            # thlu, qtu, thvu, uu, vu - are defined all interfaces not at the layer mid-
            # point. From the index notation of cumulus scheme, wu(k) is the cumulus up-
            # draft vertical velocity at the top interface of k layer.
            # Diabatic horizontal momentum forcing should be treated as a kind of 'body'
            # forcing without actual mass exchange between convective updraft and
            # environment, but still taking horizontal momentum from the environment to
            # the convective updrafts. Thus, diabatic convective momentum transport
            # vertically redistributes environmental horizontal momentum.

            if THIS_K >= (kinv - 2) and THIS_K <= (krel - 1):
                umf_zint = cbmf
                wu = winv

            if THIS_K == (krel - 1):
                emf = 0.0
                umf_zint = cbmf
                wu = wrel
                thlu = thlsrc
                qtu = qtsrc

            thj, qvj, qlj, qij, qse, id_check = conden(prel, thlsrc, qtsrc, ese, esx)

            if id_check == 1:
                condensation = True
                umf_out[0, 0, 1] = 0.0
                dcm_out = 0.0
                qvten_out = 0.0
                qlten_out = 0.0
                qiten_out = 0.0
                sten_out = 0.0
                uten_out = 0.0
                vten_out = 0.0
                qrten_out = 0.0
                qsten_out = 0.0
                cufrc_out = 0.0
                cush_inout = -1.0
                qldet_out = 0.0
                qidet_out = 0.0
                qtflx_out[0, 0, 1] = 0.0
                slflx_out[0, 0, 1] = 0.0
                uflx_out[0, 0, 1] = 0.0
                vflx_out[0, 0, 1] = 0.0
                fer_out = constants.MAPL_UNDEF
                fdr_out = constants.MAPL_UNDEF

            if not condensation:
                if THIS_K == (krel - 1):
                    thvu = thj * (1.0 + zvir * qvj - qlj - qij)


def define_env_properties(
    condensation: BoolFieldIJ,
    iteration: i32,
    krel: IntField,
    kinv: IntField,
    PGFc: Float,
    ssu0: FloatField,
    ssv0: FloatField,
    prel: FloatField,
    pifc0: FloatField,
    uu: FloatField,
    vu: FloatField,
    usrc: FloatField,
    vsrc: FloatField,
    dotransport: Int,
    tru: FloatField_NTracers,
    trsrc: FloatFieldIJ_NTracers,
    thv0rel: FloatField,
    thl0: FloatField,
    ssthl0: FloatField,
    pmid0: FloatField,
    qt0: FloatField,
    ssqt0: FloatField,
    u0: FloatField,
    v0: FloatField,
    tre: FloatFieldIJ_NTracers,
    tr0: FloatField_NTracers,
    sstr0: FloatField_NTracers,
    uplus: FloatFieldIJ,
    vplus: FloatFieldIJ,
    uplus_3D: FloatField,
    vplus_3D: FloatField,
    qsat_pe: FloatField,
    pe: FloatFieldIJ,
    thle: FloatFieldIJ,
    qte: FloatFieldIJ,
    dpe: FloatFieldIJ,
    exne: FloatFieldIJ,
    thvebot: FloatFieldIJ,
    ue: FloatFieldIJ,
    ve: FloatFieldIJ,
):
    """
    Define environmental properties at the level where buoyancy sorting occurs
    ('pe', normally, layer midpoint except in the 'krel' layer). In the 'krel'
    layer where buoyancy sorting starts to occur, however, 'pe' is defined
    differently because LCL is regarded as lower interface for mixing purpose.

    Inputs:
    condensation [BoolFieldIJ]: Mask that indicates if condensation has occurred
    iteration [i32]: Iteration of implicit CIN loop (i.e., 0 or 1)
    krel [IntField]: Release layer where buoyancy sorting first occurs
    kinv [IntField]: Inversion layer with PBL top interface as lower interface
    PGFc [Float]: Pressure gradient force
    ssu0 [FloatField]: [?]
    ssv0 [FloatField]: [?]
    prel [FloatField]: Lowest level from which buoyancy sorting occurs [Pa]
    pifc0 [FloatField]: Environmental pressure at the interfaces [Pa]
    usrc [FloatField]: Zonal wind of cumulus source air [m/s] [?]
    vsrc [FloatField]: Meridional wind of cumulus source air [m/s] [?]
    dotransport [Int]: Transport tracers [1 true]
    trsrc [FloatFieldIJ_NTracers]: Tracers of cumulus source air [?]
    thv0rel [FloatField]: [?]
    thl0 [FloatField]: Temperature [?]
    ssthl0 [FloatField]: [?]
    pmid0 [FloatField]: Environmental pressure at the layer mid-point [Pa]
    qt0 [FloatField]: Mixing ratio [?]
    ssqt0 [FloatField]: [?]
    u0 [FloatField]: Environmental zonal wind [m/s]
    v0 [FloatField]: Environmental meridional wind [m/s]
    tr0 [FloatField_NTracers]: Environmental tracers [#, kg/kg]
    sstr0 [FloatField_NTracers]: Convective tracer [?]

    Outputs:
    uu [FloatField]: Updraft zonal wind at the interface [m/s]
    vu [FloatField]: Updraft meridional wind at the interface [m/s]
    tru [FloatField_NTracers]: Updraft tracers [#, kg/kg]
    tre [FloatFieldIJ_NTracers]: [?]
    uplus [FloatFieldIJ]: [?]
    vplus [FloatFieldIJ]: [?]
    uplus_3D [FloatField]: [?]
    vplus_3D [FloatField]: [?]
    qsat_pe [FloatField]: Saturation vapor pressure at level where buoyancy occurs [?]
    pe [FloatFieldIJ]: Level where buoyancy occurs [Pa]
    thle [FloatFieldIJ]: Temperature at level where buoyancy occurs [?]
    qte [FloatFieldIJ]: [?]
    dpe [FloatFieldIJ]: [?]
    exne [FloatFieldIJ]: Exner function at level where buoyancy occurs [?]
    thvebot [FloatFieldIJ]: [?]
    ue [FloatFieldIJ]: Zonal wind at level where buoyancy occurs [?]
    ve [FloatFieldIJ]: Meridional wind at level where buoyancy occurs [?]
    """
    from __externals__ import ncnst

    with computation(FORWARD), interval(...):
        if not condensation:
            uplus = 0.0
            vplus = 0.0

            if krel == kinv - 1:
                uplus = PGFc * ssu0.at(K=kinv - 1) * (prel - pifc0.at(K=kinv - 1))
                vplus = PGFc * ssv0.at(K=kinv - 1) * (prel - pifc0.at(K=kinv - 1))

            else:
                if THIS_K >= kinv - 1 and THIS_K <= max(krel - 1, kinv - 1):
                    uplus_3D = uplus + PGFc * ssu0 * (pifc0[0, 0, 1] - pifc0)
                    vplus_3D = vplus + PGFc * ssv0 * (pifc0[0, 0, 1] - pifc0)

                uplus = uplus_3D.at(K=max(krel - 1, kinv - 1))
                vplus = vplus_3D.at(K=max(krel - 1, kinv - 1))

                uplus = uplus + PGFc * ssu0.at(K=krel) * (prel - pifc0.at(K=krel))
                vplus = vplus + PGFc * ssv0.at(K=krel) * (prel - pifc0.at(K=krel))

            if THIS_K == (krel - 1):
                uu = usrc + uplus
                vu = vsrc + vplus

            if dotransport == 1:
                n = 0
                while n < ncnst:
                    if THIS_K == krel - 1:
                        tru[0, 0, 0][n] = trsrc[0, 0][n]
                    n += 1

            pe = 0.5 * (prel + pifc0.at(K=krel + 1))
            qsat_pe = 0.5 * (prel + pifc0.at(K=krel + 1))
            dpe = prel - pifc0.at(K=krel + 1)
            exne = exnerfn(pe)
            thvebot = thv0rel
            thle = thl0.at(K=krel) + ssthl0.at(K=krel) * (pe - pmid0.at(K=krel))
            qte = qt0.at(K=krel) + ssqt0.at(K=krel) * (pe - pmid0.at(K=krel))
            ue = u0.at(K=krel) + ssu0.at(K=krel) * (pe - pmid0.at(K=krel))
            ve = v0.at(K=krel) + ssv0.at(K=krel) * (pe - pmid0.at(K=krel))
            if dotransport == 1:
                n = 0
                while n < ncnst:
                    tre[0, 0][n] = tr0.at(K=krel, ddim=[n]) + sstr0.at(
                        K=krel, ddim=[n]
                    ) * (pe - pmid0.at(K=krel))
                    n += 1


def buoyancy_sorting(
    condensation: BoolFieldIJ,
    tscaleh: FloatFieldIJ,
    krel: IntField,
    wlcl: FloatField,
    prel: FloatField,
    pifc0: FloatField,
    thv0rel: FloatField,
    thl0: FloatField,
    ssthl0: FloatField,
    pmid0: FloatField,
    qt0: FloatField,
    ssqt0: FloatField,
    u0: FloatField,
    v0: FloatField,
    ssu0: FloatField,
    ssv0: FloatField,
    dotransport: Int,
    tre: FloatFieldIJ_NTracers,
    tr0: FloatField_NTracers,
    sstr0: FloatField_NTracers,
    k0: Int,
    thlu: FloatField,
    qtu: FloatField,
    wu: FloatField,
    niter_xc: Int,
    ese: FloatField_Extra_Dim,
    esx: FloatField_Extra_Dim,
    umf_out: FloatField,
    qtflx_out: FloatField,
    slflx_out: FloatField,
    uflx_out: FloatField,
    vflx_out: FloatField,
    qsat_pe: FloatField,
    criqc: Float,
    cridist_opt: Int,
    rle: Float,
    zifc0: FloatField,
    rbuoy: Float,
    mixscale: Float,
    rkm: Float,
    zmid0: FloatField,
    detrhgt: Float,
    thlue: FloatField,
    qtue: FloatField,
    wue: FloatField,
    wtwb: FloatFieldIJ,
    dp0: FloatField,
    dt: Float,
    thv0bot: FloatField,
    exnmid0: FloatField,
    rmaxfrac: Float,
    thv0top: FloatField,
    exnifc0: FloatField,
    use_self_detrain: Int,
    rdrag: Float,
    PGFc: Float,
    tru: FloatField_NTracers,
    emf: FloatField,
    thvu: FloatField,
    umf_zint: FloatField,
    rei: FloatField,
    uu: FloatField,
    vu: FloatField,
    ufrc: FloatField,
    pe: FloatFieldIJ,
    thle: FloatFieldIJ,
    qte: FloatFieldIJ,
    dpe: FloatFieldIJ,
    exne: FloatFieldIJ,
    thvebot: FloatFieldIJ,
    ue: FloatFieldIJ,
    ve: FloatFieldIJ,
    drage: FloatFieldIJ,
    bogbot: FloatFieldIJ,
    bogtop: FloatFieldIJ,
    kpen_IJ: IntFieldIJ,
    kbup_IJ: IntFieldIJ,
    rhomid0j: FloatFieldIJ,
    fer: FloatField,
    dwten: FloatField,
    diten: FloatField,
    fdr: FloatField,
    dcm: FloatField,
    xco: FloatField,
    stop_buoyancy_sort: BoolFieldIJ,
    iteration: i32,
    cush_inout: FloatFieldIJ,
):
    """
    Buoyancy-Sorting Mixing

    In order to complete buoyancy-sorting mixing at layer mid-point, and so
    calculate 'updraft mass flux, updraft w velocity, conservative scalars'
    at the upper interface of each layer, we need following 3 information.

    1. Pressure where mixing occurs ('pe'), and temperature at 'pe' which is
         necessary to calculate various thermodynamic coefficients at pe. This
         temperature is obtained by undiluted cumulus properties lifted to pe.
    2. Undiluted updraft properties at pe - conservative scalar and vertical
         velocity -which are assumed to be the same as the properties at lower
         interface only for calculation of fractional lateral entrainment  and
         detrainment rate ( fer(k) and fdr(k) [Pa-1] ), respectively.    Final
         values of cumulus conservative scalars and w at the top interface are
         calculated afterward after obtaining fer(k) & fdr(k).
    3. Environmental properties at pe.
    """
    from __externals__ import ncnst

    with computation(FORWARD), interval(...):
        # Define cumulus scale height.
        # Cumulus scale height is defined as the maximum height cumulus can reach.
        # In case of premitive code, cumulus scale height ('cush')  at the current
        # time step was assumed to be the same as 'cush' of previous time step.
        # However, I directly calculated cush at each time step using an iterative
        # method. Note that within the cumulus scheme, 'cush' information is  used
        # only at two places during buoyancy-sorting process:
        # (1) Even negatively buoyancy mixtures with strong vertical velocity
        #       enough to rise up to 'rle*scaleh' (rle = 0.1) from pe are entrained
        #       into cumulus updraft,
        # (2) The amount of mass that is involved in buoyancy-sorting mixing
        #        process at pe is rei(k) = rkm/scaleh/rho*g [Pa-1]
        # In terms of (1), I think critical stopping distance might be replaced by
        # layer thickness. In future, we will use rei(k) = (0.5*rkm/z0(k)/rho/g).
        # In the premitive code,  'scaleh' was largely responsible for the jumping
        # variation of precipitation amount.

        if not condensation:
            scaleh = tscaleh

            if tscaleh <= 0.0:
                scaleh = 1000

            # Save time : Set iter_scaleh = 1. This will automatically use 'cush' from
            # the previous time step at the first implicit iteration. At the second
            # implicit iteration, it will use the updated 'cush' by the first implicit
            # cin. So, this updating has an effect of doing one iteration for cush
            # calculation, which is good. So, only this setting of 'iter_scaleh = 1'
            # is sufficient-enough to save computation time.

            iter_scaleh = 1

            # Initialization of 'kbup' and 'kpen'
            # 'kbup' is the top-most layer in which cloud buoyancy is positive
            # both at the top and bottom interface of the layer. 'kpen' is the
            # layer upto which cumulus panetrates ,i.e., cumulus w at the base
            # interface is positive, but becomes negative at the top interface.
            # Here, we initialize 'kbup' and 'kpen'. These initializations are
            # not trivial but important, expecially   in calculating turbulent
            # fluxes without confliction among several physics as explained in
            # detail in the part of turbulent fluxes calculation later.   Note
            # that regardless of whether 'kbup' and 'kpen' are updated or  not
            # during updraft motion,  penetrative entrainments are dumped down
            # across the top interface of 'kbup' later.      More specifically,
            # penetrative entrainment heat and moisture fluxes are  calculated
            # from the top interface of 'kbup' layer  to the base interface of
            # 'kpen' layer. Because of this, initialization of 'kbup' & 'kpen'
            # influence the convection system when there are not updated.  The
            # below initialization of 'kbup = krel' assures  that  penetrative
            # entrainment fluxes always occur at interfaces above the PBL  top
            # interfaces (i.e., only at interfaces k >=kinv ), which seems  to
            # be attractable considering that the most correct fluxes  at  the
            # PBL top interface can be ontained from the 'fluxbelowinv'  using
            # reconstructed PBL height.

            kbup_IJ = krel
            kpen_IJ = krel

            # Since 'wtw' is continuously updated during vertical motion,
            # I need below initialization command within this 'iter_scaleh'
            # do loop. Similarily, I need initializations of environmental
            # properties at 'krel' layer as below.

            wtw = wlcl * wlcl
            pe = 0.5 * (prel + pifc0.at(K=krel + 1))
            dpe = prel - pifc0.at(K=krel + 1)
            exne = exnerfn(pe)
            thvebot = thv0rel
            thle = thl0.at(K=krel) + ssthl0.at(K=krel) * (pe - pmid0.at(K=krel))
            qte = qt0.at(K=krel) + ssqt0.at(K=krel) * (pe - pmid0.at(K=krel))
            ue = u0.at(K=krel) + ssu0.at(K=krel) * (pe - pmid0.at(K=krel))
            ve = v0.at(K=krel) + ssv0.at(K=krel) * (pe - pmid0.at(K=krel))

            if dotransport == 1:
                n = 0
                while n < ncnst:
                    tre[0, 0][n] = tr0.at(K=krel, ddim=[n]) + sstr0.at(
                        K=krel, ddim=[n]
                    ) * (pe - pmid0.at(K=krel, ddim=[n]))
                    n += 1

            # Cumulus rises upward from 'prel' ( or base interface of  'krel' layer )
            # until updraft vertical velocity becomes zero.
            # Buoyancy sorting is performed via two stages. (1) Using cumulus updraft
            # properties at the base interface of each layer,perform buoyancy sorting
            # at the layer mid-point, 'pe',  and update cumulus properties at the top
            # interface, and then  (2) by averaging updated cumulus properties at the
            # top interface and cumulus properties at the base interface,   calculate
            # cumulus updraft properties at pe that will be used  in buoyancy sorting
            # mixing - thlue, qtue and, wue.  Using this averaged properties, perform
            # buoyancy sorting again at pe, and re-calculate fer(k) and fdr(k). Using
            # this recalculated fer(k) and fdr(k),  finally calculate cumulus updraft
            # properties at the top interface - thlu, qtu, thvu, uu, vu. In the below,
            # 'iter_xc = 1' performs the first stage, while 'iter_xc= 2' performs the
            # second stage. We can increase the number of iterations, 'nter_xc'.as we
            # want, but a sample test indicated that about 3 - 5 iterations  produced
            # satisfactory converent solution. Finally, identify 'kbup' and 'kpen'.

            if iteration != i32(0):
                stop_buoyancy_sort = (
                    False  # Indicates that buoyancy sorting processes are done
                )
                # (e.g., Fortran go to 45)

    with computation(FORWARD), interval(1, -2):
        if (
            THIS_K >= krel
            and THIS_K < k0 - 1
            and not stop_buoyancy_sort
            and not condensation
        ):

            thlue = thlu[0, 0, -1]
            qtue = qtu[0, 0, -1]
            wue = wu[0, 0, -1]
            wtwb = wtw

            iter_xc = 1
            while iter_xc <= niter_xc and not condensation:
                wtw = wu[0, 0, -1] * wu[0, 0, -1]

                # Calculate environmental and cumulus saturation 'excess' at 'pe'.
                # Note that in order to calculate saturation excess, we should use
                # liquid water temperature instead of temperature  as the argument
                # of "qsat". But note normal argument of "qsat" is temperature.

                thj, qvj, qlj, qij, qse, id_check = conden(pe, thle, qte, ese, esx)

                if id_check == 1:
                    condensation = True
                    umf_out[0, 0, 1] = 0.0
                    dcm_out = 0.0
                    qvten_out = 0.0
                    qlten_out = 0.0
                    qiten_out = 0.0
                    sten_out = 0.0
                    uten_out = 0.0
                    vten_out = 0.0
                    qrten_out = 0.0
                    qsten_out = 0.0
                    cufrc_out = 0.0
                    cush_inout = -1.0
                    qldet_out = 0.0
                    qidet_out = 0.0
                    qtflx_out[0, 0, 1] = 0.0
                    slflx_out[0, 0, 1] = 0.0
                    uflx_out[0, 0, 1] = 0.0
                    vflx_out[0, 0, 1] = 0.0
                    fer_out = constants.MAPL_UNDEF
                    fdr_out = constants.MAPL_UNDEF

                if not condensation:
                    thv0j = thj * (1.0 + zvir * qvj - qlj - qij)
                    rhomid0j = pe / (constants.MAPL_RDRY * thv0j * exne)
                    qsat_arg = thle * exne
                    qs, _ = QSat_Float(ese, esx, qsat_arg, qsat_pe / 100.0)
                    excess0 = qte - qs

                    thj, qvj, qlj, qij, qse, id_check = conden(
                        pe, thlue, qtue, ese, esx
                    )

                    if id_check == 1:
                        condensation = True
                        umf_out[0, 0, 1] = 0.0
                        dcm_out = 0.0
                        qvten_out = 0.0
                        qlten_out = 0.0
                        qiten_out = 0.0
                        sten_out = 0.0
                        uten_out = 0.0
                        vten_out = 0.0
                        qrten_out = 0.0
                        qsten_out = 0.0
                        cufrc_out = 0.0
                        cush_inout = -1.0
                        qldet_out = 0.0
                        qidet_out = 0.0
                        qtflx_out[0, 0, 1] = 0.0
                        slflx_out[0, 0, 1] = 0.0
                        uflx_out[0, 0, 1] = 0.0
                        vflx_out[0, 0, 1] = 0.0
                        fer_out = constants.MAPL_UNDEF
                        fdr_out = constants.MAPL_UNDEF

                    if not condensation:
                        # Detrain excessive condensate larger than 'criqc' from the
                        # cumulus updraft before performing buoyancy sorting. All I
                        # should to do is to update 'thlue' &  'que' here. Below
                        # modification is completely compatible with the other part of
                        # the code since 'thule' & 'qtue' are used only for buoyancy
                        # sorting. I found that as long as I use 'niter_xc >= 2',
                        # detraining excessive condensate before buoyancy sorting has
                        # negligible influence on the buoyancy sorting results.

                        if (qlj + qij) > criqc:
                            exql = ((qlj + qij) - criqc) * qlj / (qlj + qij)
                            exqi = ((qlj + qij) - criqc) * qij / (qlj + qij)
                            qtue = qtue - exql - exqi
                            thlue = (
                                thlue
                                + (
                                    (
                                        constants.MAPL_LATENT_HEAT_VAPORIZATION
                                        / constants.MAPL_CP
                                        / exne
                                    )
                                    * exql
                                )
                                + (
                                    (
                                        constants.MAPL_LATENT_HEAT_SUBLIMATION
                                        / constants.MAPL_CP
                                        / exne
                                    )
                                    * exqi
                                )
                            )

                        thj, qvj, qlj, qij, qse, id_check = conden(
                            pe, thlue, qtue, ese, esx
                        )

                        if id_check == 1:
                            condensation = True
                            umf_out[0, 0, 1] = 0.0
                            dcm_out = 0.0
                            qvten_out = 0.0
                            qlten_out = 0.0
                            qiten_out = 0.0
                            sten_out = 0.0
                            uten_out = 0.0
                            vten_out = 0.0
                            qrten_out = 0.0
                            qsten_out = 0.0
                            cufrc_out = 0.0
                            cush_inout = -1.0
                            qldet_out = 0.0
                            qidet_out = 0.0
                            qtflx_out[0, 0, 1] = 0.0
                            slflx_out[0, 0, 1] = 0.0
                            uflx_out[0, 0, 1] = 0.0
                            vflx_out[0, 0, 1] = 0.0
                            fer_out = constants.MAPL_UNDEF
                            fdr_out = constants.MAPL_UNDEF

                        if not condensation:
                            thvj = thj * (1.0 + zvir * qvj - qlj - qij)
                            tj = (
                                thj * exne
                            )  # This 'tj' is used for computing thermo. coeffs. below
                            qsat_arg = thlue * exne
                            qs, _ = QSat_Float(ese, esx, qsat_arg, qsat_pe / 100.0)
                            excessu = qtue - qs

                            # Calculate critical mixing fraction, 'xc'. Mixture with
                            # mixing ratio smaller than 'xc' will be entrained into
                            # cumulus updraft. Both the saturated updrafts with
                            # 'positive buoyancy' or 'negative buoyancy + strong
                            # vertical velocity enough to rise certain threshold
                            # distance' are kept into the updraft in the below program.
                            # If the core updraft is unsaturated, we can set 'xc = 0'
                            # and let the cumulus  convection still works or we may
                            # exit. Current below code does not entrain unsaturated
                            # mixture. However it should be modified such that it also
                            # entrain unsaturated mixture.

                            # cridis : Critical stopping distance for buoyancy sorting
                            # purpose. scaleh is only used here.

                            if cridist_opt == 0:
                                cridis = rle * scaleh  # Original code
                            else:
                                cridis = rle * (zifc0[0, 0, 1] - zifc0)
                                # New code

                            # Buoyancy Sorting
                            # Case 1 : When both cumulus and env. are unsaturated or
                            # saturated.
                            xsat = 0.0

                            if (excessu <= 0.0 and excess0 <= 0.0) or (
                                excessu >= 0.0 and excess0 >= 0.0
                            ):
                                xc = min(
                                    1.0,
                                    max(
                                        0.0,
                                        1.0
                                        - 2.0
                                        * rbuoy
                                        * constants.MAPL_GRAV
                                        * cridis
                                        / wue**2.0
                                        * (1.0 - thvj / thv0j),
                                    ),
                                )
                                aquad = 0.0
                                bquad = 0.0
                                cquad = 0.0
                                if excessu > 0.0:
                                    xsat = 1.0
                                else:
                                    xsat = 0.0

                            else:
                                # Case 2 : When either cumulus or env. is saturated. !
                                xsat = excessu / (excessu - excess0)
                                thlxsat = thlue + xsat * (thle - thlue)
                                qtxsat = qtue + xsat * (qte - qtue)

                                thj, qvj, qlj, qij, qse, id_check = conden(
                                    pe, thlxsat, qtxsat, ese, esx
                                )

                                if id_check == 1:
                                    condensation = True
                                    umf_out[0, 0, 1] = 0.0
                                    dcm_out = 0.0
                                    qvten_out = 0.0
                                    qlten_out = 0.0
                                    qiten_out = 0.0
                                    sten_out = 0.0
                                    uten_out = 0.0
                                    vten_out = 0.0
                                    qrten_out = 0.0
                                    qsten_out = 0.0
                                    cufrc_out = 0.0
                                    cush_inout = -1.0
                                    qldet_out = 0.0
                                    qidet_out = 0.0
                                    qtflx_out[0, 0, 1] = 0.0
                                    slflx_out[0, 0, 1] = 0.0
                                    uflx_out[0, 0, 1] = 0.0
                                    vflx_out[0, 0, 1] = 0.0
                                    fer_out = constants.MAPL_UNDEF
                                    fdr_out = constants.MAPL_UNDEF

                                if not condensation:
                                    thvxsat = thj * (1.0 + zvir * qvj - qlj - qij)

                                    # kk=1 : Cumulus segment, kk=2 : Environment segment
                                    kk = 1
                                    while kk <= 2:
                                        if xsat == 1.0:
                                            xsat = 1.0 + 1e-6
                                        if kk == 1:
                                            thv_x0 = thvj
                                            thv_x1 = (1.0 - 1.0 / xsat) * thvj + (
                                                1.0 / xsat
                                            ) * thvxsat
                                        else:
                                            thv_x1 = thv0j
                                            thv_x0 = ((xsat / (xsat - 1.0)) * thv0j) + (
                                                (1.0 / (1.0 - xsat)) * thvxsat
                                            )

                                        aquad = wue**2
                                        bquad = (
                                            2.0
                                            * rbuoy
                                            * constants.MAPL_GRAV
                                            * cridis
                                            * (thv_x1 - thv_x0)
                                            / thv0j
                                            - 2.0 * wue**2
                                        )
                                        cquad = (
                                            2.0
                                            * rbuoy
                                            * constants.MAPL_GRAV
                                            * cridis
                                            * (thv_x0 - thv0j)
                                            / thv0j
                                            + wue**2
                                        )
                                        if kk == 1:
                                            if (bquad**2 - 4.0 * aquad * cquad) >= 0.0:
                                                xs1, xs2, status = roots(
                                                    aquad, bquad, cquad
                                                )
                                                x_cu = min(
                                                    1.0,
                                                    max(
                                                        0.0,
                                                        min(xsat, min(xs1, xs2)),
                                                    ),
                                                )
                                            else:
                                                x_cu = xsat

                                        else:
                                            if (bquad**2 - 4.0 * aquad * cquad) >= 0.0:
                                                xs1, xs2, status = roots(
                                                    aquad, bquad, cquad
                                                )
                                                x_en = min(
                                                    1.0,
                                                    max(
                                                        0.0,
                                                        max(xsat, min(xs1, xs2)),
                                                    ),
                                                )
                                            else:
                                                x_en = 1.0

                                        kk += 1

                                    if x_cu == xsat:
                                        xc = max(x_cu, x_en)
                                    else:
                                        xc = x_cu

                            # Compute fractional lateral entrainment & detrainment rate
                            # in each layers. The unit of rei(k), fer(k), and fdr(k) is
                            # [Pa-1]. Alternative choice of 'rei(k)' is also shown
                            # below, where coefficient 0.5 was from approximate tuning
                            # against the BOMEX case. In order to prevent the onset of
                            # instability in association with cumulus induced
                            # subsidence advection, cumulus mass flux at the top
                            # interface  in any layer should be smaller than (90% of)
                            # total mass within that layer. I imposed limits on
                            # 'rei(k)' as below,  in such that stability condition is
                            # always satisfied. Below limiter of 'rei(k)' becomes
                            # negative for some cases, causing error. So, for the time
                            # being, I came back to the original limiter.

                            if not condensation:
                                ee2 = xc**2
                                ud2 = 1.0 - 2.0 * xc + xc**2
                                if min(scaleh, mixscale) != 0.0:
                                    rei = (
                                        (
                                            rkm
                                            + max(
                                                0.0,
                                                (zmid0 - detrhgt) / 200.0,
                                            )
                                        )
                                        / min(scaleh, mixscale)
                                        / constants.MAPL_GRAV
                                        / rhomid0j
                                    )

                                else:
                                    rei = (
                                        (0.5 * rkm)
                                        / zmid0
                                        / constants.MAPL_GRAV
                                        / rhomid0j
                                    )

                                if xc > 0.5:
                                    rei = min(
                                        rei,
                                        0.9
                                        * log(
                                            dp0
                                            / constants.MAPL_GRAV
                                            / dt
                                            / umf_zint[0, 0, -1]
                                            + 1.0
                                        )
                                        / dpe
                                        / (2.0 * xc - 1.0),
                                    )

                                fer = rei * ee2
                                fdr = rei * ud2
                                xco = xc

                                # Iteration Start due to 'maxufrc' constraint
                                # Calculate cumulus updraft mass flux and penetrative
                                # entrainment mass flux. Note that  non-zero
                                # penetrative entrainment mass flux will be asigned
                                # only to interfaces from the top interface of 'kbup'
                                # layer to the base interface of 'kpen' layer as will
                                # be shown later.

                                umf_zint = umf_zint[0, 0, -1] * exp(dpe * (fer - fdr))

                                emf = 0.0

                                dcm = (
                                    0.5
                                    * (umf_zint + umf_zint[0, 0, -1])
                                    * rei
                                    * dpe
                                    * min(1.0, max(0.0, xsat - xc))
                                )

                                # Compute cumulus updraft properties at the top
                                # interface. Also use Tayler expansion in order to
                                # treat limiting case

                                if fer * dpe < 1.0e-4:
                                    thlu = (
                                        thlu[0, 0, -1]
                                        + (thle + ssthl0 * dpe / 2.0 - thlu[0, 0, -1])
                                        * fer
                                        * dpe
                                    )

                                    qtu = (
                                        qtu[0, 0, -1]
                                        + (qte + ssqt0 * dpe / 2.0 - qtu[0, 0, -1])
                                        * fer
                                        * dpe
                                    )

                                    uu = (
                                        uu[0, 0, -1]
                                        + (ue + ssu0 * dpe / 2.0 - uu[0, 0, -1])
                                        * fer
                                        * dpe
                                        - PGFc * ssu0 * dpe
                                    )

                                    vu = (
                                        vu[0, 0, -1]
                                        + (ve + ssv0 * dpe / 2.0 - vu[0, 0, -1])
                                        * fer
                                        * dpe
                                        - PGFc * ssv0 * dpe
                                    )

                                    if dotransport == 1:
                                        n = 0
                                        while n < ncnst:
                                            tru[0, 0, 0][n] = (
                                                tru[0, 0, -1][n]
                                                + (
                                                    tre[0, 0][n]
                                                    + sstr0.at(K=THIS_K, ddim=[n])
                                                    * dpe
                                                    / 2.0
                                                    - tru[0, 0, -1][n]
                                                )
                                                * fer
                                                * dpe
                                            )
                                            n += 1

                                else:
                                    thlu = (
                                        thle + ssthl0 / fer - ssthl0 * dpe / 2.0
                                    ) - (
                                        thle
                                        + ssthl0 * dpe / 2.0
                                        - thlu[0, 0, -1]
                                        + ssthl0 / fer
                                    ) * exp(
                                        -fer * dpe
                                    )

                                    qtu = (qte + ssqt0 / fer - ssqt0 * dpe / 2.0) - (
                                        qte
                                        + ssqt0 * dpe / 2.0
                                        - qtu[0, 0, -1]
                                        + ssqt0 / fer
                                    ) * exp(-fer * dpe)

                                    uu = (
                                        ue
                                        + (1.0 - PGFc) * ssu0 / fer
                                        - ssu0 * dpe / 2.0
                                    ) - (
                                        ue
                                        + ssu0 * dpe / 2.0
                                        - uu[0, 0, -1]
                                        + (1.0 - PGFc) * ssu0 / fer
                                    ) * exp(
                                        -fer * dpe
                                    )
                                    vu = (
                                        ve
                                        + (1.0 - PGFc) * ssv0 / fer
                                        - ssv0 * dpe / 2.0
                                    ) - (
                                        ve
                                        + ssv0 * dpe / 2.0
                                        - vu[0, 0, -1]
                                        + (1.0 - PGFc) * ssv0 / fer
                                    ) * exp(
                                        -fer * dpe
                                    )

                                    if dotransport == 1:
                                        n = 0
                                        while n < ncnst:
                                            tru[0, 0, 0][n] = (
                                                tre[0, 0][n]
                                                + sstr0[0, 0, 0][n] / fer
                                                - sstr0[0, 0, 0][n] * dpe / 2.0
                                            ) - (
                                                tre[0, 0][n]
                                                + sstr0[0, 0, 0][n] * dpe / 2.0
                                                - tru[0, 0, -1][n]
                                                + sstr0[0, 0, 0][n] / fer
                                            ) * exp(
                                                -fer * dpe
                                            )
                                            n += 1

                                # Expel some of cloud water and ice from cumulus
                                # updraft at the top interface.  Note that this
                                # is not 'detrainment' term  but a 'sink' term of
                                # cumulus updraft qt ( or one part of 'source' term
                                # of  mean environmental qt ). At this stage, as
                                # the most simplest choice, if condensate amount
                                # within cumulus updraft is larger than a critical
                                # value, 'criqc', expels the surplus condensate from
                                # cumulus updraft to the environment. A certain
                                # fraction ( e.g., 'frc_sus' ) of this expelled
                                # condesnate will be in a form that can be suspended
                                # in the layer k where it was formed, while the
                                # other fraction, '1-frc_sus' will be in a form of
                                # precipitatble (e.g.,can potentially fall down
                                # across the base interface of layer k ). In turn we
                                # should describe subsequent falling of
                                # precipitable condensate ('1-frc_sus') across the
                                # base interface of the layer k, &  evaporation of
                                # precipitating water in the below layer k-1 and
                                # associated evaporative cooling of the later, k-1,
                                # and falling of 'non-evaporated precipitating water
                                # (which was initially formed in layer k ) and a
                                # newly-formed precipitable water in the layer, k-1',
                                # across the base interface of the lower layer k-1.
                                # Cloud microphysics should correctly describe all
                                # of these process.  In a near future, I should
                                # significantly modify this cloud microphysics,
                                # including precipitation-induced downdraft also.

                                thj, qvj, qlj, qij, qse, id_check = conden(
                                    pifc0[0, 0, 1],
                                    thlu,
                                    qtu,
                                    ese,
                                    esx,
                                )

                                if id_check == 1:
                                    condensation = True
                                    umf_out[0, 0, 1] = 0.0
                                    dcm_out = 0.0
                                    qvten_out = 0.0
                                    qlten_out = 0.0
                                    qiten_out = 0.0
                                    sten_out = 0.0
                                    uten_out = 0.0
                                    vten_out = 0.0
                                    qrten_out = 0.0
                                    qsten_out = 0.0
                                    cufrc_out = 0.0
                                    cush_inout = -1.0
                                    qldet_out = 0.0
                                    qidet_out = 0.0
                                    qtflx_out[0, 0, 1] = 0.0
                                    slflx_out[0, 0, 1] = 0.0
                                    uflx_out[0, 0, 1] = 0.0
                                    vflx_out[0, 0, 1] = 0.0
                                    fer_out = constants.MAPL_UNDEF
                                    fdr_out = constants.MAPL_UNDEF

                                if not condensation:
                                    if (qlj + qij) > criqc:
                                        exql = ((qlj + qij) - criqc) * qlj / (qlj + qij)
                                        exqi = ((qlj + qij) - criqc) * qij / (qlj + qij)

                                        # It is very important to re-update 'qtu'
                                        # and 'thlu' at the upper interface after
                                        # expelling condensate from cumulus updraft
                                        # at the top interface of the layer. As
                                        # mentioned above, this is a 'sink' of
                                        # cumulus qt (or equivalently, a 'source'
                                        # of environmentasl qt), not a regular
                                        # convective'detrainment'.

                                        qtu = qtu - exql - exqi
                                        thlu = (
                                            thlu
                                            + (
                                                constants.MAPL_LATENT_HEAT_VAPORIZATION
                                                / exnifc0[0, 0, 1]
                                                / constants.MAPL_CP
                                            )
                                            * exql
                                            + (
                                                constants.MAPL_LATENT_HEAT_SUBLIMATION
                                                / exnifc0[0, 0, 1]
                                                / constants.MAPL_CP
                                            )
                                            * exqi
                                        )

                                        # Expelled cloud condensate into the
                                        # environment from the updraft. After
                                        # all the calculation later, 'dwten' and
                                        # 'diten' will have a unit of [ kg/kg/s ],
                                        # because it is a tendency of qt.
                                        # Restoration of 'dwten' and 'diten' to
                                        # this correct unit through  multiplying
                                        # 'umf(k)*g/dp0(k)' will be performed
                                        # later after finally updating 'umf' using
                                        # a 'rmaxfrac' constraint near the end of
                                        # this updraft buoyancy sorting loop.

                                        dwten = exql
                                        diten = exqi

                                    else:
                                        dwten = 0.0
                                        diten = 0.0

                                    # Update 'thvu(k)' after detraining condensate
                                    # from cumulus updraft.
                                    thj, qvj, qlj, qij, qse, id_check = conden(
                                        pifc0[0, 0, 1],
                                        thlu,
                                        qtu,
                                        ese,
                                        esx,
                                    )

                                    if id_check == 1:
                                        condensation = True
                                        umf_out[0, 0, 1] = 0.0
                                        dcm_out = 0.0
                                        qvten_out = 0.0
                                        qlten_out = 0.0
                                        qiten_out = 0.0
                                        sten_out = 0.0
                                        uten_out = 0.0
                                        vten_out = 0.0
                                        qrten_out = 0.0
                                        qsten_out = 0.0
                                        cufrc_out = 0.0
                                        cush_inout = -1.0
                                        qldet_out = 0.0
                                        qidet_out = 0.0
                                        qtflx_out[0, 0, 1] = 0.0
                                        slflx_out[0, 0, 1] = 0.0
                                        uflx_out[0, 0, 1] = 0.0
                                        vflx_out[0, 0, 1] = 0.0
                                        fer_out = constants.MAPL_UNDEF
                                        fdr_out = constants.MAPL_UNDEF

                                    if not condensation:

                                        thvu = thj * (1.0 + zvir * qvj - qlj - qij)

                                        # Calculate updraft vertical velocity at the
                                        # upper interface. In order to calculate
                                        # 'wtw' at the upper interface, we use 'wtw'
                                        # at the lower interface. Note  'wtw'  is
                                        # continuously updated as cumulus updraft
                                        # rises.

                                        bogbot = rbuoy * (
                                            thvu[0, 0, -1] / thvebot - 1.0
                                        )  # Cloud buoyancy at base interface
                                        bogtop = rbuoy * (
                                            thvu / thv0top - 1.0
                                        )  # Cloud buoyancy at top  interface

                                        delbog = bogtop - bogbot
                                        drage = fer * (1.0 + rdrag)
                                        expfac = exp(-2.0 * drage * dpe)

                                        wtwb = wtw

                                        if drage * dpe > 1.0e-3:
                                            wtw = wtw * expfac + (
                                                delbog
                                                + (1.0 - expfac)
                                                * (
                                                    bogbot
                                                    + delbog / (-2.0 * drage * dpe)
                                                )
                                            ) / (rhomid0j * drage)
                                        else:
                                            wtw = (
                                                wtw + dpe * (bogbot + bogtop) / rhomid0j
                                            )

                                        # Force the plume rise at least to klfc of the
                                        # undiluted plume. Because even the below is
                                        # not complete, I decided not to include this.

                                        # Repeat 'iter_xc' iteration loop until
                                        # 'iter_xc = niter_xc'. Also treat the case
                                        # even when wtw < 0 at the 'kpen' interface.

                                        if wtw > 0.0:
                                            thlue = 0.5 * (thlu[0, 0, -1] + thlu)
                                            qtue = 0.5 * (qtu[0, 0, -1] + qtu)
                                            wue = 0.5 * sqrt(max(wtwb + wtw, 0.0))

                                        else:
                                            iter_xc = (
                                                niter_xc + 1
                                            )  # Break out of iter_xc loop

                iter_xc += 1  # end iter_xc loop

            # Add the contribution of self-detrainment  to vertical variations of
            # cumulus updraft mass flux. The reason why we are trying to include
            # self-detrainment is as follows.  In current scheme,  vertical
            # variation of updraft mass flux is not fully consistent with the
            # vertical variation of updraft vertical w. For example, within a
            # given layer, let's assume that  cumulus w is positive at the base
            # interface, while negative at the top interface. This means that
            # cumulus updraft cannot reach to the top interface of the layer.
            # However, cumulus updraft mass flux at the top interface is not zero
            # according to the vertical tendency equation of cumulus mass flux.
            # Ideally, cumulus updraft mass flux at the top interface should be
            # zero for this case. In order to assures that cumulus updraft mass
            # flux goes to zero when cumulus updraft vertical velocity goes to
            # zero, we are imposing self-detrainment term as below by considering
            # layer-mean cloud buoyancy and cumulus updraft vertical velocity
            # square at the top interface. Use of auto-detrainment term will be
            # determined by setting 'use_self_detrain=.true.' in the parameter
            # sentence.

            if not condensation:
                if use_self_detrain == 1:
                    autodet = min(
                        0.5
                        * constants.MAPL_GRAV
                        * (bogbot + bogtop)
                        / (max(wtw, 0.0) + 1.0e-4),
                        0.0,
                    )

                    umf_zint = umf_zint * exp(
                        0.637 * (dpe / rhomid0j / constants.MAPL_GRAV) * autodet
                    )

                if umf_zint == 0.0:
                    wtw = -1.0

                # 'kbup' is the upper most layer in which cloud buoyancy  is positive
                # both at the base and top interface.  'kpen' is the upper most layer
                # up to cumulus can reach. Usually, 'kpen' is located higher than the
                # 'kbup'. Note we initialized these by 'kbup = krel' & 'kpen = krel'.
                # As explained before, it is possible that only 'kpen' is updated,
                # while 'kbup' keeps its initialization value. For this case, current
                # scheme will simply turns-off penetrative entrainment fluxes and use
                # normal buoyancy-sorting fluxes for 'kbup <= k <= kpen-1' interfaces,
                # in order to describe shallow continental cumulus convection.

                if bogtop > 0.0 and wtw > 0.0:
                    kbup_IJ = THIS_K

                if wtw <= 0.0:
                    kpen_IJ = THIS_K
                    stop_buoyancy_sort = (
                        True  # Done calculating updraft properties, break out of loop
                    )

                if not stop_buoyancy_sort:
                    wu = sqrt(wtw)

                    if wu > 100.0:
                        condensation = True
                        umf_out[0, 0, 1] = 0.0
                        dcm_out = 0.0
                        qvten_out = 0.0
                        qlten_out = 0.0
                        qiten_out = 0.0
                        sten_out = 0.0
                        uten_out = 0.0
                        vten_out = 0.0
                        qrten_out = 0.0
                        qsten_out = 0.0
                        cufrc_out = 0.0
                        cush_inout = -1.0
                        qldet_out = 0.0
                        qidet_out = 0.0
                        qtflx_out[0, 0, 1] = 0.0
                        slflx_out[0, 0, 1] = 0.0
                        uflx_out[0, 0, 1] = 0.0
                        vflx_out[0, 0, 1] = 0.0
                        fer_out = constants.MAPL_UNDEF
                        fdr_out = constants.MAPL_UNDEF

                    if not condensation:

                        # Iteration end due to 'rmaxfrac' constraint

                        # Calculate updraft fractional area at the upper interface and
                        # set upper limit to 'ufrc' by 'rmaxfrac'. In order to keep the
                        # consistency  among ['ufrc','umf','wu (or wtw)'], if ufrc is
                        # limited by 'rmaxfrac', either 'umf' or 'wu' should be changed.
                        # Although both 'umf' and 'wu (wtw)' at the current upper
                        # interface are used for updating 'umf' & 'wu'  at the next
                        # upper interface, 'umf' is a passive variable not influencing
                        # the buoyancy sorting process in contrast to 'wtw'. This is a
                        # reason why we adjusted 'umf' instead of 'wtw'. In turn we
                        # updated 'fdr' here instead of 'fer',  which guarantees  that
                        # all previously updated thermodynamic variables at the upper
                        # interface before applying 'rmaxfrac' constraint are already
                        # internally consistent,  even though 'ufrc'  is  limited by
                        # 'rmaxfrac'. Thus, we don't need to go through interation loop
                        # again. If we update 'fer' however, we should go through above
                        # iteration loop.

                        rhoifc0j = pifc0[0, 0, 1] / (
                            constants.MAPL_RDRY
                            * 0.5
                            * (thv0bot[0, 0, 1] + thv0top)
                            * exnifc0[0, 0, 1]
                        )

                        ufrc[0, 0, 1] = umf_zint / (rhoifc0j * wu)

                        if ufrc[0, 0, 1] > rmaxfrac:
                            ufrc[0, 0, 1] = rmaxfrac
                            umf_zint = rmaxfrac * rhoifc0j * wu
                            fdr = fer - log(umf_zint / umf_zint[0, 0, -1]) / dpe

                        # Update environmental properties for at the mid-point of next
                        # upper layer for use in buoyancy sorting.

                        pe = pmid0[0, 0, 1]
                        dpe = dp0[0, 0, 1]
                        exne = exnmid0[0, 0, 1]
                        thvebot = thv0bot[0, 0, 1]
                        thle = thl0[0, 0, 1]
                        qte = qt0[0, 0, 1]
                        ue = u0[0, 0, 1]
                        ve = v0[0, 0, 1]
                        if dotransport == 1:
                            n = 0
                            while n < ncnst:
                                tre[0, 0][n] = tr0[0, 0, 1][n]
                                n += 1


def calc_ppen(
    condensation: BoolFieldIJ,
    drage: FloatFieldIJ,
    bogbot: FloatFieldIJ,
    bogtop: FloatFieldIJ,
    pifc0: FloatField,
    kpen_IJ: IntFieldIJ,
    kpen: IntField,
    wu: FloatField,
    rhomid0j: FloatFieldIJ,
    dp0: FloatField,
    wtwb: FloatFieldIJ,
    ppen: FloatFieldIJ,
):
    """
    Calculate 'ppen( < 0 )', updraft penetrative distance from the lower
    interface of 'kpen' layer. Note that bogbot & bogtop at the 'kpen' layer
    either when fer is zero or non-zero was already calculated above.
    It seems that below qudarature solving formula is valid only when
    bogbot < 0.

    Note: The suffix '_IJ' is used to indicate 2D variables (e.g., kpen_IJ)

    Inputs:
    condensation [BoolFieldIJ]: Mask that indicates if condensation has occurred
    drage [FloatFieldIJ]: Drag coefficient [no unit] [?]
    bogbot [FloatFieldIJ]: [?]
    bogtop [FloatFieldIJ]: [?]
    pifc0 [FloatField]: Environmental pressure at the interfaces [Pa]
    kpen_IJ [IntFieldIJ]: Highest layer with positive updraft velocity
    wu [FloatField]: Updraft vertical velocity at the interface [m/s]
    rhomid0j [FloatFieldIJ]: Density of water at mid layers [?]
    dp0 [FloatField]: Environmental layer pressure thickness [Pa] > 0
    wtwb [FloatFieldIJ]: [?]

    Outputs:
    kpen [IntField]: Highest layer with positive updraft velocity
    ppen [FloatFieldIJ]: Highest interface evel where Cu w = 0 [Pa]
    """
    with computation(FORWARD), interval(...):
        if not condensation:
            # Up to this point, we finished all of buoyancy sorting processes from the
            # 'krel' layer to 'kpen' layer: at the top interface of individual layers,
            # we calculated updraft and penetrative mass fluxes [umf(k) & emf(k) = 0],
            # updraft fractional area [ ufrc(k) ],  updraft vertical velocity [wu(k)],
            # updraft  thermodynamic variables [thlu(k),qtu(k),uu(k),vu(k),thvu(k)].
            # In the layer,we also calculated fractional entrainment-detrainment rate
            # [ fer(k), fdr(k) ], and detrainment tendency of water and ice from
            # cumulus updraft [ dwten(k), diten(k) ]. In addition, we updated and
            # identified 'krel' and 'kpen' layer index, if any.  In the 'kpen' layer,
            # we calculated everything mentioned above except the 'wu(k)' and 'ufrc(k)'
            # since a real value of updraft vertical velocity is not defined at the kpen
            # top interface (note 'ufrc' at the top interface of layer is calculated
            # from 'umf(k)' and 'wu(k)'). As mentioned before, special treatment is
            # required when 'kbup' is not updated and so 'kbup = krel'.

            # During the 'iter_scaleh' iteration loop, non-physical (with non-zero
            # values) values can remain in the variable arrays above (also 'including'
            # in case of wu and ufrc at the top interface) the 'kpen' layer. This can
            # happen when the kpen layer index identified from the 'iter_scaleh = 1'
            # iteration loop is located at above the kpen layer index identified from
            # 'iter_scaleh = 3' iteration loop. Thus, in the following calculations,
            # we should only use the values in each variables only up to finally
            # identified 'kpen' layer & 'kpen' interface except 'wu' and 'ufrc' at the
            # top interface of 'kpen' layer.    Note that in order to prevent any
            # problems due to these non-physical values, I re-initialized the values
            # of [ umf(kpen:k0), emf(kpen:k0), dwten(kpen+1:k0), diten(kpen+1:k0),
            # fer(kpen:k0), fdr(kpen+1:k0), ufrc(kpen:k0) ] to be zero after
            # 'iter_scaleh' do loop.

            # Calculate 'ppen( < 0 )', updraft penetrative distance from the lower
            # interface of 'kpen' layer. Note that bogbot & bogtop at the 'kpen' layer
            # either when fer is zero or non-zero was already calculated above.
            # It seems that below qudarature solving formula is valid only when
            # bogbot < 0.

            # Below solving equation is clearly wrong ! I should revise this !

            kpen = kpen_IJ  # Convert kpen_IJ to a 3D field

            if drage == 0.0:
                aquad = (bogtop - bogbot) / (pifc0.at(K=kpen + 1) - pifc0.at(K=kpen))
                bquad = 2.0 * bogbot
                cquad = -1 * wu.at(K=kpen - 1) ** 2 * rhomid0j
                xc1, xc2, status = roots(aquad, bquad, cquad)
                if status == 0:
                    if xc1 <= 0.0 and xc2 <= 0.0:
                        ppen = max(xc1, xc2)
                        ppen = min(0.0, max(-dp0.at(K=kpen), ppen))

                    elif xc1 > 0.0 and xc2 > 0.0:
                        ppen = -1 * dp0.at(K=kpen)
                    else:
                        ppen = min(xc1, xc2)
                        ppen = min(0.0, max(-dp0.at(K=kpen), ppen))

                else:
                    ppen = -1 * dp0.at(K=kpen)

            else:
                ppen = compute_ppen(
                    wtwb, drage, bogbot, bogtop, rhomid0j, dp0.at(K=kpen)
                )


def recalc_condensate(
    condensation: BoolFieldIJ,
    fer: FloatField,
    kpen: IntField,
    ppen: FloatFieldIJ,
    thlu: FloatField,
    thl0: FloatField,
    ssthl0: FloatField,
    qtu: FloatField,
    qt0: FloatField,
    ssqt0: FloatField,
    pifc0: FloatField,
    ese: FloatField_Extra_Dim,
    esx: FloatField_Extra_Dim,
    umf_out: FloatField,
    qtflx_out: FloatField,
    slflx_out: FloatField,
    uflx_out: FloatField,
    vflx_out: FloatField,
    criqc: Float,
    thv0bot: FloatField,
    thv0top: FloatField,
    exnifc0: FloatField,
    zifc0: FloatField,
    kbup_IJ: IntFieldIJ,
    kbup: IntField,
    krel: IntField,
    k0: Int,
    umf_zint: FloatField,
    emf: FloatField,
    ufrc: FloatField,
    dwten: FloatField,
    diten: FloatField,
    dwten_temp: FloatField,
    diten_temp: FloatField,
    thlu_top: FloatFieldIJ,
    qtu_top: FloatFieldIJ,
    cldhgt: FloatFieldIJ,
    umf_temp: FloatField,
    fdr: FloatField,
    xco: FloatField,
    cush: FloatFieldIJ,
    cush_inout: FloatFieldIJ,
    iteration: i32,
):
    """
    Re-calculate the amount of expelled condensate from cloud updraft
    at the cumulus top. This is necessary for refined calculations of
    bulk cloud microphysics at the cumulus top. Note that ppen < 0.
    In the below, I explicitly calculate 'thlu_top' & 'qtu_top' by
    using non-zero 'fer(kpen)'.

    Inputs:
    fer [FloatField]: Fractional lateral entrainment rate [1/Pa]
    kpen [IntField]: Highest layer with positive updraft velocity
    ppen [FloatFieldIJ]: Highest interface evel where Cu w = 0 [Pa]
    thlu [FloatField]: Updraft liquid potential temperature at the interface [K]
    thl0 [FloatField]: Temperature [?]
    ssthl0 [FloatField]: [?]
    qtu [FloatField]: Updraft total specific humidity at the interface [kg/kg]
    qt0 [FloatField]: Mixing ratio [?]
    ssqt0 [FloatField]: Mixing ratio [?]
    pifc0 [FloatField]: Environmental pressure at the interfaces [Pa]
    ese [FloatField_Extra_Dim]: Used in QSat_Float [n/a]
    esx [FloatField_Extra_Dim]: Used in QSat_Float [n/a]
    criqc [Float]: Maximum condensate that can be hold by cumulus updraft [kg/kg]
    thv0bot [FloatField]: Temperature at bottom [?]
    thv0top [FloatField]: Temperature at top [?]
    exnifc0 [FloatField]: Exner function at interfaces
    zifc0 [FloatField]: Environmental height at interfaces [m]
    kbup_IJ [IntFieldIJ]: Top layer in which buoyancy is positive at top interface
    krel [IntField]: Release layer where buoyancy sorting first occurs
    k0 [Int]: Number of levels
    umf_zint [FloatField]: Updraft mass flux at the interfaces [kg/m2/s]
    emf [FloatField]: [?]
    ufrc [FloatField]: Cumulus updraft fraction [fraction]
    iteration [i32]: Iteration of implicit CIN loop (i.e., 0 or 1)

    Inouts:
    condensation [BoolFieldIJ]: Mask that indicates if condensation has occurred

    Outputs:
    umf_out [FloatField]: Updraft mass flux at the interfaces [kg/m2/s]
    qtflx_out [FloatField]: Mixing ratio flux [?]
    slflx_out [FloatField]: Sensible heat flux [?]
    uflx_out [FloatField]: Zonal wind flux [?]
    vflx_out [FloatField]: Meridional wind flux [?]
    kbup [IntField]: Top layer in which buoyancy is positive at top interface
    dwten [FloatField]: Detrained [?]
    diten [FloatField]: Detrained [?]
    dwten_temp [FloatField]: Detrained [?]
    diten_temp [FloatField]: Detrained [?]
    thlu_top [FloatFieldIJ]: Updraft liquid potential temperature at top [K]
    qtu_top [FloatFieldIJ]: Updraft total specific humidity at top [kg/kg]
    cldhgt [FloatFieldIJ]: Cloud height [?]
    umf_temp [FloatField]: Updraft mass flux at the interfaces [kg/m2/s]
    fdr [FloatField]: Fractional lateral detrainment rate [1/Pa]
    xco [FloatField]: [?]
    cush [FloatFieldIJ]: Convective scale height [m]
    cush_inout [FloatFieldIJ]: Convective scale height [m]
    """
    with computation(FORWARD), interval(...):
        if not condensation:
            if fer.at(K=kpen) * (-ppen) < 1.0e-4:
                thlu_top = thlu.at(K=kpen - 1) + (
                    thl0.at(K=kpen)
                    + ssthl0.at(K=kpen) * (-ppen) / 2.0
                    - thlu.at(K=kpen - 1)
                ) * fer.at(K=kpen) * (-ppen)
                qtu_top = qtu.at(K=kpen - 1) + (
                    qt0.at(K=kpen)
                    + ssqt0.at(K=kpen) * (-ppen) / 2.0
                    - qtu.at(K=kpen - 1)
                ) * fer.at(K=kpen) * (-ppen)
            else:
                thlu_top = (
                    thl0.at(K=kpen)
                    + ssthl0.at(K=kpen) / fer.at(K=kpen)
                    - ssthl0.at(K=kpen) * (-ppen) / 2.0
                ) - (
                    thl0.at(K=kpen)
                    + ssthl0.at(K=kpen) * (-ppen) / 2.0
                    - thlu.at(K=kpen - 1)
                    + ssthl0.at(K=kpen) / fer.at(K=kpen)
                ) * exp(
                    -fer.at(K=kpen) * (-ppen)
                )
                qtu_top = (
                    qt0.at(K=kpen)
                    + ssqt0.at(K=kpen) / fer.at(K=kpen)
                    - ssqt0.at(K=kpen) * (-ppen) / 2.0
                ) - (
                    qt0.at(K=kpen)
                    + ssqt0.at(K=kpen) * (-ppen) / 2.0
                    - qtu.at(K=kpen - 1)
                    + ssqt0.at(K=kpen) / fer.at(K=kpen)
                ) * exp(
                    -fer.at(K=kpen) * (-ppen)
                )

            thj, qvj, qlj, qij, qse, id_check = conden(
                pifc0.at(K=kpen) + ppen, thlu_top, qtu_top, ese, esx
            )

            if id_check == 1:
                condensation = True
                umf_out[0, 0, 1] = 0.0
                dcm_out = 0.0
                qvten_out = 0.0
                qlten_out = 0.0
                qiten_out = 0.0
                sten_out = 0.0
                uten_out = 0.0
                vten_out = 0.0
                qrten_out = 0.0
                qsten_out = 0.0
                cufrc_out = 0.0
                cush_inout = -1.0
                qldet_out = 0.0
                qidet_out = 0.0
                qtflx_out[0, 0, 1] = 0.0
                slflx_out[0, 0, 1] = 0.0
                uflx_out[0, 0, 1] = 0.0
                vflx_out[0, 0, 1] = 0.0
                fer_out = constants.MAPL_UNDEF
                fdr_out = constants.MAPL_UNDEF

            if not condensation:
                p00 = 1e5
                rovcp = constants.MAPL_RGAS / constants.MAPL_CP
                exntop = ((pifc0.at(K=kpen) + ppen) / p00) ** rovcp
                if (qlj + qij) > criqc:
                    if THIS_K == kpen:
                        dwten = ((qlj + qij) - criqc) * qlj / (qlj + qij)
                        diten = ((qlj + qij) - criqc) * qij / (qlj + qij)

                    qtu_top = qtu_top - dwten.at(K=kpen) - diten.at(K=kpen)
                    thlu_top = (
                        thlu_top
                        + (
                            (
                                constants.MAPL_LATENT_HEAT_VAPORIZATION
                                / constants.MAPL_CP
                                / exntop
                            )
                            * dwten.at(K=kpen)
                        )
                        + (
                            (
                                constants.MAPL_LATENT_HEAT_SUBLIMATION
                                / constants.MAPL_CP
                                / exntop
                            )
                            * diten.at(K=kpen)
                        )
                    )
                else:
                    if THIS_K == kpen:
                        dwten = 0.0
                        diten = 0.0

                # Calculate cumulus scale height as the top height that cumulus can
                # reach.
                rhoifc0j = pifc0.at(K=kpen) / (
                    constants.MAPL_RDRY
                    * 0.5
                    * (thv0bot.at(K=kpen) + thv0top.at(K=kpen - 1))
                    * exnifc0.at(K=kpen)
                )
                cush = zifc0.at(K=kpen) - ppen / rhoifc0j / constants.MAPL_GRAV
                scaleh = cush

                # The 'forcedCu' is logical identifier saying whether cumulus updraft
                # overcome the buoyancy barrier just above the PBL top. If it is true,
                # cumulus did not overcome the barrier -  this is a shallow convection
                # with negative cloud buoyancy, mimicking  shallow continental cumulus
                # convection. Depending on 'forcedCu' parameter, treatment of heat  &
                # moisture fluxes at the entraining interfaces, 'kbup <= k < kpen - 1'
                # will be set up in a different ways, as will be shown later.

                kbup = kbup_IJ  # Convert kbup into 3D field

                if kbup == krel:
                    forcedCu = True
                else:
                    forcedCu = False

                # Filtering of unerasonable cumulus adjustment here.  This is a very
                # important process which should be done cautiously. Various ways of
                # filtering are possible depending on cases mainly using the indices
                # of key layers - 'klcl','kinv','krel','klfc','kbup','kpen'. At this
                # stage, the followings are all possible : 'kinv >= 2', 'klcl >= 1',
                # 'krel >= kinv', 'kbup >= krel', 'kpen >= krel'. I must design this
                # filtering very cautiously, in such that none of  realistic cumulus
                # convection is arbitrarily turned-off. Potentially, I might turn-off
                # cumulus convection if layer-mean 'ql > 0' in the 'kinv-1' layer,in
                # order to suppress cumulus convection growing, based at the Sc top.
                # This is one of potential future modifications. Note that ppen < 0.

                cldhgt = pifc0.at(K=kpen) + ppen

                if forcedCu:
                    condensation = True
                    umf_out[0, 0, 1] = 0.0
                    dcm_out = 0.0
                    qvten_out = 0.0
                    qlten_out = 0.0
                    qiten_out = 0.0
                    sten_out = 0.0
                    uten_out = 0.0
                    vten_out = 0.0
                    qrten_out = 0.0
                    qsten_out = 0.0
                    cufrc_out = 0.0
                    cush_inout = -1.0
                    qldet_out = 0.0
                    qidet_out = 0.0
                    qtflx_out[0, 0, 1] = 0.0
                    slflx_out[0, 0, 1] = 0.0
                    uflx_out[0, 0, 1] = 0.0
                    vflx_out[0, 0, 1] = 0.0
                    fer_out = constants.MAPL_UNDEF
                    fdr_out = constants.MAPL_UNDEF

                if not condensation:
                    # Re-initializing some key variables above the 'kpen' layer in
                    # order to suppress the influence of non-physical values above
                    # 'kpen', in association with the use of 'iter_scaleh' loop.
                    # Note that umf, emf,  ufrc are defined at the interfaces (0:k0),
                    # while 'dwten','diten', 'fer', 'fdr' are defined at layer
                    # mid-points. Initialization of 'fer' and 'fdr' is for correct
                    # writing purpose of diagnostic output. Note that we set
                    # umf(kpen)=emf(kpen)=ufrc(kpen)=0, in consistent  with wtw < 0
                    # at the top interface of 'kpen' layer. However, we still have
                    # non-zero expelled cloud condensate in the 'kpen' layer.
                    if THIS_K >= kpen - 1 and THIS_K <= k0:
                        umf_zint[0, 0, 1] = 0.0
                        emf[0, 0, 1] = 0.0
                    if THIS_K >= kpen and THIS_K <= k0:
                        ufrc[0, 0, 1] = 0.0
                    if THIS_K >= kpen + 1 and THIS_K < k0:
                        dwten = 0.0
                        diten = 0.0
                        fer = 0.0
                        fdr = 0.0
                        xco = 0.0

                    # Update output variables as needed
                    dwten_temp = dwten
                    diten_temp = diten
                    umf_temp[0, 0, 1] = umf_zint


def calc_entrainment_mass_flux(
    condensation: BoolFieldIJ,
    k0: Int,
    thlu: FloatField,
    qtu: FloatField,
    uu: FloatField,
    vu: FloatField,
    tru: FloatField_NTracers,
    dotransport: Int,
    tru_emf: FloatField_NTracers,
    kpen: IntField,
    kbup: IntField,
    pifc0: FloatField,
    thv0bot: FloatField,
    thv0top: FloatField,
    exnifc0: FloatField,
    umf_zint: FloatField,
    ppen: FloatFieldIJ,
    rei: FloatField,
    rpen: Float,
    dp0: FloatField,
    dt: Float,
    thl0: FloatField,
    ssthl0: FloatField,
    pmid0: FloatField,
    qt0: FloatField,
    ssqt0: FloatField,
    u0: FloatField,
    ssu0: FloatField,
    v0: FloatField,
    ssv0: FloatField,
    tr0: FloatField_NTracers,
    sstr0: FloatField_NTracers,
    use_cumpenent: Int,
    thlu_emf: FloatField,
    qtu_emf: FloatField,
    uu_emf: FloatField,
    vu_emf: FloatField,
    emf: FloatField,
    iteration: i32,
):
    """
    Calculate downward penetrative entrainment mass flux, 'emf(k) < 0',  and
    thermodynamic properties of penetratively entrained airs at   entraining
    interfaces. emf(k) is defined from the top interface of the  layer  kbup
    to the bottom interface of the layer 'kpen'. Note even when  kbup = krel,
    i.e.,even when 'kbup' was not updated in the above buoyancy  sorting  do
    loop (i.e., 'kbup' remains as the initialization value),   below do loop
    of penetrative entrainment flux can be performed without  any conceptual
    or logical problems, because we have already computed all  the variables
    necessary for performing below penetrative entrainment block.
    In the below 'do' loop, 'k' is an interface index at which non-zero 'emf'
    (penetrative entrainment mass flux) is calculated. Since cumulus updraft
    is negatively buoyant in the layers between the top interface of 'kbup'
    layer (interface index, kbup) and the top interface of 'kpen' layer, the
    fractional lateral entrainment, fer(k) within these layers will be close
    to zero - so it is likely that only strong lateral detrainment occurs in
    thses layers. Under this situation,we can easily calculate the amount of
    detrainment cumulus air into these negatively buoyanct layers by  simply
    comparing cumulus updraft mass fluxes between the base and top interface
    of each layer: emf(k) = emf(k-1)*exp(-fdr(k)*dp0(k))
                          ~ emf(k-1)*(1-rei(k)*dp0(k))
                   emf(k-1)-emf(k) ~ emf(k-1)*rei(k)*dp0(k)
    Current code assumes that about 'rpen~10' times of these detrained  mass
    are penetratively re-entrained down into the 'k-1' interface. And all of
    these detrained masses are finally dumped down into the top interface of
    'kbup' layer. Thus, the amount of penetratively entrained air across the
    top interface of 'kbup' layer with 'rpen~10' becomes too large.
    Note that this penetrative entrainment part can be completely turned-off
    and we can simply use normal buoyancy-sorting involved turbulent  fluxes
    by modifying 'penetrative entrainment fluxes' part below.

    Calculate entrainment mass flux and conservative scalars of entraining
    free air at interfaces of 'kbup <= k < kpen - 1'

    Inputs:
    condensation [BoolFieldIJ]: Mask that indicates if condensation has occurred
    k0 [Int]: Number of levels
    thlu [FloatField]: Updraft liquid potential temperature at the interface [K]
    qtu [FloatField]: Updraft total specific humidity at the interface [kg/kg]
    uu [FloatField]: Updraft zonal wind at the interface [m/s]
    vu [FloatField]: Updraft meridional wind at the interface [m/s]
    tru [FloatField_NTracers]: Updraft tracers [#, kg/kg]
    dotransport [Int]: Transport tracers [1 true]
    kpen [IntField]: Highest layer with positive updraft velocity
    kbup [IntField]: Top layer in which buoyancy is positive at top interface
    pifc0 [FloatField]: Environmental pressure at the interfaces [Pa]
    thv0bot [FloatField]: Temperature at bottom [?]
    thv0top [FloatField]: Temperature at top [?]
    exnifc0 [FloatField]: Exner function at interfaces
    umf_zint [FloatField]: Updraft mass flux at the interfaces [kg/m2/s]
    ppen [FloatFieldIJ]: Highest interface evel where Cu w = 0 [Pa]
    rei [FloatField]: Updraft fractional mixing rate with the environment [1/Pa]
    rpen [Float]: Penentrative entrainment factor
    dp0 [FloatField]: Environmental layer pressure thickness [Pa] > 0
    dt [Float]: Timestep [s]
    thl0 [FloatField]: Temperature [?]
    ssthl0 [FloatField]: [?]
    pmid0 [FloatField]: Environmental pressure at the layer mid-point [Pa]
    qt0 [FloatField]: Mixing ratio [?]
    ssqt0 [FloatField]: [?]
    u0 [FloatField]: Environmental zonal wind [m/s]
    ssu0 [FloatField]: [?]
    v0 [FloatField]: Environmental meridional wind [m/s]
    ssv0 [FloatField]: [?]
    tr0 [FloatField_NTracers]: Environmental tracers [#, kg/kg]
    sstr0 [FloatField_NTracers]: Convective tracer [?]
    use_cumpenent [Int]: Cumulative penetrative entrainment
    iteration [i32]: Iteration of implicit CIN loop (i.e., 0 or 1)

    Outputs:
    tru_emf [FloatField_NTracers]: Penetrative Downdraft tracers at entraining interfaces [#, kg/kg]
    thlu_emf [FloatField]: Penetrative downdraft liquid potential temperature at entraining interfaces [K]
    qtu_emf [FloatField]: Penetrative downdraft total water at entraining interfaces [kg/kg]
    uu_emf [FloatField]: Penetrative downdraft zonal wind at entraining interfaces [m/s]
    vu_emf [FloatField]: Penetrative downdraft meridional wind at entraining interfaces [m/s]
    emf [FloatField]: [?]
    """
    from __externals__ import ncnst

    with computation(FORWARD), interval(...):
        if not condensation:
            thlu_emf[0, 0, 1] = thlu
            qtu_emf[0, 0, 1] = qtu
            uu_emf[0, 0, 1] = uu
            vu_emf[0, 0, 1] = vu
            if dotransport == 1:
                n = 0
                while n < ncnst:
                    tru_emf[0, 0, 1][n] = tru[0, 0, 0][n]
                    n += 1

    with computation(BACKWARD), interval(...):
        if not condensation:
            if (
                THIS_K <= kpen - 1 and THIS_K >= kbup
            ):  # Here, 'k' is an interface index at which
                rhoifc0j = pifc0[0, 0, 1] / (
                    constants.MAPL_RDRY
                    * 0.5
                    * (thv0bot[0, 0, 1] + thv0top)
                    * exnifc0[0, 0, 1]
                )

                if THIS_K == (kpen - 1):
                    # Note that 'ppen' has already been calculated in the above
                    # 'iter_scaleh' loop assuming zero lateral entrainmentin the
                    # layer 'kpen'.

                    # Calculate returning mass flux, emf ( < 0 )
                    # Current penetrative entrainment rate with 'rpen~10' is too
                    # large and future refinement is necessary including the
                    # definition of 'thl','qt' of penetratively entrained air.
                    # Penetratively entrained airs across the 'kpen-1' interface
                    # is assumed to have the properties of the base interface of
                    # 'kpen' layer. Note that 'emf ~ - umf/ufrc = - w * rho'. Thus,
                    # below limit sets an upper limit of |emf| to be ~ 10cm/s,
                    # which is very loose constraint. Here, I used more restricted
                    # constraint on the limit of emf, assuming 'emf' cannot exceed
                    # a net mass within the layer above the interface. Similar to
                    # the case of warming and drying due to cumulus updraft induced
                    # compensating subsidence, penetrative entrainment induces
                    # compensating upwelling -     in order to prevent numerical
                    # instability in association with compensating upwelling, we
                    # should similarily limit the amount of penetrative entrainment
                    # at the interface by the amount of masses within the layer just
                    # above the penetratively entraining interface.

                    emf = max(
                        max(
                            umf_zint * ppen * rei.at(K=kpen) * rpen,
                            -0.1 * rhoifc0j,
                        ),
                        -0.9 * dp0.at(K=kpen) / constants.MAPL_GRAV / dt,
                    )
                    thlu_emf = thl0.at(K=kpen) + ssthl0.at(K=kpen) * (
                        pifc0[0, 0, 1] - pmid0.at(K=kpen)
                    )
                    qtu_emf = qt0.at(K=kpen) + ssqt0.at(K=kpen) * (
                        pifc0[0, 0, 1] - pmid0.at(K=kpen)
                    )
                    uu_emf = u0.at(K=kpen) + ssu0.at(K=kpen) * (
                        pifc0[0, 0, 1] - pmid0.at(K=kpen)
                    )
                    vu_emf = v0.at(K=kpen) + ssv0.at(K=kpen) * (
                        pifc0[0, 0, 1] - pmid0.at(K=kpen)
                    )

                    if dotransport == 1:
                        n = 0
                        while n < ncnst:
                            tru_emf[0, 0, 0][n] = tr0.at(K=kpen, ddim=[n]) + sstr0.at(
                                K=kpen, ddim=[n]
                            ) * (pifc0[0, 0, 1] - pmid0.at(K=kpen))
                            n += 1

                else:

                    # Note we are coming down from the higher interfaces to the lower
                    # interfaces. Also note that 'emf < 0'. So, below operation is a
                    # summing not subtracting. In order to ensure numerical stability,
                    # I imposed a modified correct limit of '-0.9*dp0(k+1)/g/dt' on
                    # emf(k).

                    if (
                        use_cumpenent == 1
                    ):  # Original Cumulative Penetrative Entrainment

                        emf = max(
                            max(
                                emf[0, 0, 1]
                                - umf_zint * dp0[0, 0, 1] * rei[0, 0, 1] * rpen,
                                -0.1 * rhoifc0j,
                            ),
                            -0.9 * dp0[0, 0, 1] / constants.MAPL_GRAV / dt,
                        )
                        if abs(emf) > abs(emf[0, 0, 1]):
                            thlu_emf = (
                                thlu_emf[0, 0, 1] * emf[0, 0, 1]
                                + thl0[0, 0, 1] * (emf - emf[0, 0, 1])
                            ) / emf
                            qtu_emf = (
                                qtu_emf[0, 0, 1] * emf[0, 0, 1]
                                + qt0[0, 0, 1] * (emf - emf[0, 0, 1])
                            ) / emf
                            uu_emf = (
                                uu_emf[0, 0, 1] * emf[0, 0, 1]
                                + u0[0, 0, 1] * (emf - emf[0, 0, 1])
                            ) / emf
                            vu_emf = (
                                vu_emf[0, 0, 1] * emf[0, 0, 1]
                                + v0[0, 0, 1] * (emf - emf[0, 0, 1])
                            ) / emf

                            if dotransport == 1:
                                n = 0
                                while n < ncnst:
                                    tru_emf[0, 0, 0][n] = (
                                        tru_emf[0, 0, 1][n] * emf[0, 0, 1]
                                        + tr0[0, 0, 1][n] * (emf - emf[0, 0, 1])
                                    ) / emf
                                    n += 1
                        else:
                            thlu_emf = thl0[0, 0, 1]
                            qtu_emf = qt0[0, 0, 1]
                            uu_emf = u0[0, 0, 1]
                            vu_emf = v0[0, 0, 1]
                            if dotransport == 1:
                                n = 0
                                while n < ncnst:
                                    tru_emf[0, 0, 0][n] = tr0[0, 0, 1][n]
                                    n += 1

                    else:  # Alternative Non-Cumulative Penetrative Entrainment
                        emf = max(
                            max(
                                -umf_zint * dp0[0, 0, 1] * rei[0, 0, 1] * rpen,
                                -0.1 * rhoifc0j,
                            ),
                            -0.9 * dp0[0, 0, 1] / constants.MAPL_GRAV / dt,
                        )
                        thlu_emf = thl0[0, 0, 1]
                        qtu_emf = qt0[0, 0, 1]
                        uu_emf = u0[0, 0, 1]
                        vu_emf = v0[0, 0, 1]
                        if dotransport == 1:
                            n = 0
                            while n < ncnst:
                                tru_emf[0, 0, 0][n] = tr0[0, 0, 1][n]
                                n += 1


def calc_pbl_fluxes(
    condensation: BoolFieldIJ,
    qtsrc: FloatField,
    qt0: FloatField,
    ssqt0: FloatField,
    pifc0: FloatField,
    pmid0: FloatField,
    kinv: IntField,
    cbmf: FloatField,
    dt: Float,
    xflx: FloatField,
    qtflx: FloatField,
    uflx: FloatField,
    vflx: FloatField,
    slflx: FloatField,
    thlsrc: FloatField,
    thl0: FloatField,
    ssthl0: FloatField,
    exnifc0: FloatField,
    usrc: FloatField,
    u0: FloatField,
    ssu0: FloatField,
    vsrc: FloatField,
    v0: FloatField,
    ssv0: FloatField,
    dotransport: Int,
    trsrc: FloatFieldIJ_NTracers,
    tr0: FloatField_NTracers,
    sstr0: FloatField_NTracers,
    trflx: FloatField_NTracers,
    xflx_ndim: FloatField_NTracers,
    iteration: i32,
):
    from __externals__ import ncnst

    with computation(FORWARD), interval(...):
        # Compute turbulent heat, moisture, momentum flux at all interfaces

        # 1. PBL fluxes :  0 <= k <= kinv - 1
        # All the information necessary to reconstruct PBL
        # height are passed to 'fluxbelowinv'.
        if not condensation:
            kinv = kinv - 1  # Adjust kinv by 1
            xsrc = qtsrc
            xmean = qt0.at(K=kinv)
            xtop = qt0.at(K=kinv + 1) + ssqt0.at(K=kinv + 1) * (
                pifc0.at(K=kinv + 1) - pmid0.at(K=kinv + 1)
            )
            xbot = qt0.at(K=kinv - 1) + ssqt0.at(K=kinv - 1) * (
                pifc0.at(K=kinv) - pmid0.at(K=kinv - 1)
            )

    with computation(FORWARD), interval(...):
        if not condensation:
            xflx[0, 0, 1] = 0.0
            xflx = 0.0

    with computation(FORWARD), interval(...):
        if not condensation:

            k_below = kinv - 1
            dp = pifc0.at(K=k_below + 1) - pifc0.at(K=kinv + 1)

            # Compute reconstructed inversion height
            xtop_ori = xtop
            xbot_ori = xbot
            rcbmf = (
                cbmf * constants.MAPL_GRAV * dt
            ) / dp  # Can be larger than 1 : 'OK'

            if xbot >= xtop:
                rpeff = (xmean - xtop) / max(1.0e-20, xbot - xtop)
            else:
                rpeff = (xmean - xtop) / min(-1.0e-20, xbot - xtop)

            rpeff = min(max(0.0, rpeff), 1.0)  # As of this, 0<= rpeff <= 1
            if rpeff == 0.0 or rpeff == 1.0:
                xbot = xmean
                xtop = xmean

            rr = rpeff / rcbmf
            pinv = pifc0.at(K=k_below + 1) - rpeff * dp  # "pinv" before detraining mass
            pinv_eff = (
                pifc0.at(K=k_below + 1) + (rcbmf - rpeff) * dp
            )  # Effective "pinv" after detraining mass

            # Compute turbulent fluxes.
            # Below two cases exactly converges at 'kinv-1' interface when rr = 1.
            if THIS_K < k_below + 1:
                xflx[0, 0, 1] = (
                    cbmf
                    * (xsrc - xbot)
                    * (pifc0.at(K=0) - pifc0[0, 0, 1])
                    / (pifc0.at(K=0) - pinv)
                )

            if THIS_K == k_below + 1:
                if rr <= 1:
                    xflx = xflx - (1.0 - rr) * cbmf * (xtop_ori - xbot_ori)

    with computation(FORWARD), interval(...):
        if not condensation:
            if THIS_K <= kinv:
                qtflx = xflx

            xsrc = thlsrc
            xmean = thl0.at(K=kinv)
            xtop = thl0.at(K=kinv + 1) + ssthl0.at(K=kinv + 1) * (
                pifc0.at(K=kinv + 1) - pmid0.at(K=kinv + 1)
            )
            xbot = thl0.at(K=kinv - 1) + ssthl0.at(K=kinv - 1) * (
                pifc0.at(K=kinv) - pmid0.at(K=kinv - 1)
            )

    with computation(FORWARD), interval(...):
        if not condensation:
            xflx[0, 0, 1] = 0.0
            xflx = 0.0

    with computation(FORWARD), interval(...):
        if not condensation:

            k_below = kinv - 1
            dp = pifc0.at(K=k_below + 1) - pifc0.at(K=kinv + 1)

            # Compute reconstructed inversion height
            xtop_ori = xtop
            xbot_ori = xbot
            rcbmf = (
                cbmf * constants.MAPL_GRAV * dt
            ) / dp  # Can be larger than 1 : 'OK'

            if xbot >= xtop:
                rpeff = (xmean - xtop) / max(1.0e-20, xbot - xtop)
            else:
                rpeff = (xmean - xtop) / min(-1.0e-20, xbot - xtop)

            rpeff = min(max(0.0, rpeff), 1.0)  # As of this, 0<= rpeff <= 1
            if rpeff == 0.0 or rpeff == 1.0:
                xbot = xmean
                xtop = xmean

            rr = rpeff / rcbmf
            pinv = pifc0.at(K=k_below + 1) - rpeff * dp  # "pinv" before detraining mass
            pinv_eff = (
                pifc0.at(K=k_below + 1) + (rcbmf - rpeff) * dp
            )  # Effective "pinv" after detraining mass

            # Compute turbulent fluxes.
            # Below two cases exactly converges at 'kinv-1' interface when rr = 1.
            if THIS_K < k_below + 1:
                xflx[0, 0, 1] = (
                    cbmf
                    * (xsrc - xbot)
                    * (pifc0.at(K=0) - pifc0[0, 0, 1])
                    / (pifc0.at(K=0) - pinv)
                )

            if THIS_K == k_below + 1:
                if rr <= 1:
                    xflx = xflx - (1.0 - rr) * cbmf * (xtop_ori - xbot_ori)

    with computation(FORWARD), interval(...):
        if not condensation:
            if THIS_K <= kinv:
                slflx = constants.MAPL_CP * exnifc0 * xflx

            xsrc = usrc
            xmean = u0.at(K=kinv)
            xtop = u0.at(K=kinv + 1) + ssu0.at(K=kinv + 1) * (
                pifc0.at(K=kinv + 1) - pmid0.at(K=kinv + 1)
            )
            xbot = u0.at(K=kinv - 1) + ssu0.at(K=kinv - 1) * (
                pifc0.at(K=kinv) - pmid0.at(K=kinv - 1)
            )

    with computation(FORWARD), interval(...):
        if not condensation:
            xflx[0, 0, 1] = 0.0
            xflx = 0.0

    with computation(FORWARD), interval(...):
        if not condensation:

            k_below = kinv - 1
            dp = pifc0.at(K=k_below + 1) - pifc0.at(K=kinv + 1)

            # Compute reconstructed inversion height
            xtop_ori = xtop
            xbot_ori = xbot
            rcbmf = (
                cbmf * constants.MAPL_GRAV * dt
            ) / dp  # Can be larger than 1 : 'OK'

            if xbot >= xtop:
                rpeff = (xmean - xtop) / max(1.0e-20, xbot - xtop)
            else:
                rpeff = (xmean - xtop) / min(-1.0e-20, xbot - xtop)

            rpeff = min(max(0.0, rpeff), 1.0)  # As of this, 0<= rpeff <= 1
            if rpeff == 0.0 or rpeff == 1.0:
                xbot = xmean
                xtop = xmean

            rr = rpeff / rcbmf
            pinv = pifc0.at(K=k_below + 1) - rpeff * dp  # "pinv" before detraining mass
            pinv_eff = (
                pifc0.at(K=k_below + 1) + (rcbmf - rpeff) * dp
            )  # Effective "pinv" after detraining mass

            # Compute turbulent fluxes.
            # Below two cases exactly converges at 'kinv-1' interface when rr = 1.
            if THIS_K < k_below + 1:
                xflx[0, 0, 1] = (
                    cbmf
                    * (xsrc - xbot)
                    * (pifc0.at(K=0) - pifc0[0, 0, 1])
                    / (pifc0.at(K=0) - pinv)
                )

            if THIS_K == k_below + 1:
                if rr <= 1:
                    xflx = xflx - (1.0 - rr) * cbmf * (xtop_ori - xbot_ori)
    with computation(FORWARD), interval(...):
        if not condensation:
            if THIS_K <= kinv:
                uflx = xflx

            xsrc = vsrc
            xmean = v0.at(K=kinv)
            xtop = v0.at(K=kinv + 1) + ssv0.at(K=kinv + 1) * (
                pifc0.at(K=kinv + 1) - pmid0.at(K=kinv + 1)
            )
            xbot = v0.at(K=kinv - 1) + ssv0.at(K=kinv - 1) * (
                pifc0.at(K=kinv) - pmid0.at(K=kinv - 1)
            )

    with computation(FORWARD), interval(...):
        if not condensation:
            xflx[0, 0, 1] = 0.0
            xflx = 0.0

    with computation(FORWARD), interval(...):
        if not condensation:

            k_below = kinv - 1
            dp = pifc0.at(K=k_below + 1) - pifc0.at(K=kinv + 1)

            # Compute reconstructed inversion height
            xtop_ori = xtop
            xbot_ori = xbot
            rcbmf = (
                cbmf * constants.MAPL_GRAV * dt
            ) / dp  # Can be larger than 1 : 'OK'

            if xbot >= xtop:
                rpeff = (xmean - xtop) / max(1.0e-20, xbot - xtop)
            else:
                rpeff = (xmean - xtop) / min(-1.0e-20, xbot - xtop)

            rpeff = min(max(0.0, rpeff), 1.0)  # As of this, 0<= rpeff <= 1
            if rpeff == 0.0 or rpeff == 1.0:
                xbot = xmean
                xtop = xmean

            rr = rpeff / rcbmf
            pinv = pifc0.at(K=k_below + 1) - rpeff * dp  # "pinv" before detraining mass
            pinv_eff = (
                pifc0.at(K=k_below + 1) + (rcbmf - rpeff) * dp
            )  # Effective "pinv" after detraining mass

            # Compute turbulent fluxes.
            # Below two cases exactly converges at 'kinv-1' interface when rr = 1.
            if THIS_K < k_below + 1:
                xflx[0, 0, 1] = (
                    cbmf
                    * (xsrc - xbot)
                    * (pifc0.at(K=0) - pifc0[0, 0, 1])
                    / (pifc0.at(K=0) - pinv)
                )

            if THIS_K == k_below + 1:
                if rr <= 1:
                    xflx = xflx - (1.0 - rr) * cbmf * (xtop_ori - xbot_ori)

    with computation(FORWARD), interval(...):
        if not condensation:
            if THIS_K <= kinv:
                vflx = xflx

    with computation(FORWARD), interval(...):
        if not condensation:
            if dotransport == 1:
                n = 0
                while n < ncnst:
                    xflx_ndim[0, 0, 1][n] = 0.0
                    xflx_ndim[0, 0, 0][n] = 0.0
                    n += 1

    with computation(FORWARD), interval(...):
        if not condensation:
            if dotransport == 1:
                n = 0
                while n < ncnst:
                    xsrc = trsrc[0, 0][n]
                    xmean = tr0.at(K=kinv, ddim=[n])
                    xtop = tr0.at(K=kinv + 1, ddim=[n]) + sstr0.at(
                        K=kinv + 1, ddim=[n]
                    ) * (pifc0.at(K=kinv + 1) - pmid0.at(K=kinv + 1))
                    xbot = tr0.at(K=kinv - 1, ddim=[n]) + sstr0.at(
                        K=kinv - 1, ddim=[n]
                    ) * (pifc0.at(K=kinv) - pmid0.at(K=kinv - 1))

                    k_below = kinv - 1
                    dp = pifc0.at(K=k_below + 1) - pifc0.at(K=kinv + 1)

                    # Compute reconstructed inversion height
                    xtop_ori = xtop
                    xbot_ori = xbot
                    rcbmf = (
                        cbmf * constants.MAPL_GRAV * dt
                    ) / dp  # Can be larger than 1 : 'OK'

                    if xbot >= xtop:
                        rpeff = (xmean - xtop) / max(1.0e-20, xbot - xtop)
                    else:
                        rpeff = (xmean - xtop) / min(-1.0e-20, xbot - xtop)

                    rpeff = min(max(0.0, rpeff), 1.0)  # As of this, 0<= rpeff <= 1
                    if rpeff == 0.0 or rpeff == 1.0:
                        xbot = xmean
                        xtop = xmean

                    rr = rpeff / rcbmf
                    pinv = (
                        pifc0.at(K=k_below + 1) - rpeff * dp
                    )  # "pinv" before detraining mass
                    pinv_eff = (
                        pifc0.at(K=k_below + 1) + (rcbmf - rpeff) * dp
                    )  # Effective "pinv" after detraining mass

                    # Compute turbulent fluxes.
                    # Below two cases exactly converges at 'kinv-1' interface when
                    # rr = 1.
                    if THIS_K < k_below + 1:
                        xflx_ndim[0, 0, 1][n] = (
                            cbmf
                            * (xsrc - xbot)
                            * (pifc0.at(K=0) - pifc0[0, 0, 1])
                            / (pifc0.at(K=0) - pinv)
                        )

                    if THIS_K == k_below + 1:
                        if rr <= 1:
                            xflx_ndim[0, 0, 0][n] = xflx_ndim[0, 0, 0][n] - (
                                1.0 - rr
                            ) * cbmf * (xtop_ori - xbot_ori)

                    if THIS_K <= kinv:
                        trflx[0, 0, 0][n] = xflx_ndim[0, 0, 0][n]

                    n += 1


def non_buoyancy_sorting_fluxes(
    condensation: BoolFieldIJ,
    kinv: IntField,
    krel: IntField,
    cbmf: FloatField,
    qtsrc: FloatField,
    qt0: FloatField,
    ssqt0: FloatField,
    pifc0: FloatField,
    pmid0: FloatField,
    thlsrc: FloatField,
    thl0: FloatField,
    ssthl0: FloatField,
    PGFc: Float,
    exnifc0: FloatField,
    ssu0: FloatField,
    ssv0: FloatField,
    u0: FloatField,
    v0: FloatField,
    usrc: FloatField,
    vsrc: FloatField,
    dotransport: Int,
    trflx: FloatField_NTracers,
    trsrc: FloatFieldIJ_NTracers,
    tr0: FloatField_NTracers,
    sstr0: FloatField_NTracers,
    uflx: FloatField,
    vflx: FloatField,
    slflx: FloatField,
    qtflx: FloatField,
):
    """
    2. Calculate non-buoyancy sorting fluxes : kinv <= k <= krel - 1
    Note that when 'krel = kinv', below block is never executed
    as in a desirable, expected way ( but I must check  if this
    is the case ). The non-buoyancy sorting fluxes are computed
    only when 'krel > kinv'.

    Inputs:
    condensation [BoolFieldIJ]: Mask that indicates if condensation has occurred
    kinv [IntField]: Inversion layer with PBL top interface as lower interface
    krel [IntField]: Release layer where buoyancy sorting first occurs
    cbmf [FloatField]: Cloud base mass flux [kg/m2/s]
    qtsrc [FloatField]: Mixing ratio of cumulus source air [?]
    qt0 [FloatField]: Mixing ratio [?]
    ssqt0 [FloatField]: [?]
    pifc0 [FloatField]: Environmental pressure at the interfaces [Pa]
    pmid0 [FloatField]: Environmental pressure at the layer mid-point [Pa]
    thlsrc [FloatField]: Temperature of cumulus source air [K] [?]
    thl0 [FloatField]: Temperature [?]
    ssthl0 [FloatField]: [?]
    PGFc [Float]: Pressure gradient force
    exnifc0 [FloatField]: Exner function at interfaces
    ssu0 [FloatField]: [?]
    ssv0 [FloatField]: [?]
    u0 [FloatField]: Environmental zonal wind [m/s]
    v0 [FloatField]: Environmental meridional wind [m/s]
    usrc [FloatField]: Zonal wind of cumulus source air [m/s] [?]
    vsrc [FloatField]: Meridional wind of cumulus source air [m/s] [?]
    dotransport [Int]: Transport tracers [1 true]
    trsrc [FloatFieldIJ_NTracers]: Tracers of cumulus source air [?]
    tr0 [FloatField_NTracers]: Environmental tracers [#, kg/kg]
    sstr0 [FloatField_NTracers]: Convective tracer [?]

    Outputs:
    trflx [FloatField_NTracers]: Updraft/pen.entrainment tracer flux [#/m2/s, kg/kg/m2/s]
    uflx [FloatField]: Zonal wind flux [?]
    vflx [FloatField]: Meridional wind flux [?]
    slflx [FloatField]: Sensible heat flux [?]
    qtflx [FloatField]: Mixing ratio flux [?]
    """
    from __externals__ import ncnst

    with computation(FORWARD), interval(...):
        if not condensation:
            uplus = 0.0
            vplus = 0.0
            if THIS_K >= kinv and THIS_K < (krel - 1):
                qtflx[0, 0, 1] = cbmf * (
                    qtsrc
                    - (
                        qt0[0, 0, 1]
                        + ssqt0[0, 0, 1] * (pifc0[0, 0, 1] - pmid0[0, 0, 1])
                    )
                )
                slflx[0, 0, 1] = (
                    cbmf
                    * (
                        thlsrc
                        - (
                            thl0[0, 0, 1]
                            + ssthl0[0, 0, 1] * (pifc0[0, 0, 1] - pmid0[0, 0, 1])
                        )
                    )
                    * constants.MAPL_CP
                    * exnifc0[0, 0, 1]
                )
                uplus = uplus + PGFc * ssu0 * (pifc0[0, 0, 1] - pifc0)
                vplus = vplus + PGFc * ssv0 * (pifc0[0, 0, 1] - pifc0)
                uflx[0, 0, 1] = cbmf * (
                    usrc
                    + uplus
                    - (u0[0, 0, 1] + ssu0[0, 0, 1] * (pifc0[0, 0, 1] - pmid0[0, 0, 1]))
                )
                vflx[0, 0, 1] = cbmf * (
                    vsrc
                    + vplus
                    - (v0[0, 0, 1] + ssv0[0, 0, 1] * (pifc0[0, 0, 1] - pmid0[0, 0, 1]))
                )
                if dotransport == 1:
                    n = 0
                    while n < ncnst:
                        trflx[0, 0, 1][n] = cbmf * (
                            trsrc[0, 0][n]
                            - (
                                tr0[0, 0, 1][n]
                                + sstr0[0, 0, 1][n] * (pifc0[0, 0, 1] - pmid0[0, 0, 1])
                            )
                        )
                        n += 1


def buoyancy_sorting_fluxes(
    condensation: BoolFieldIJ,
    kbup: IntField,
    krel: IntField,
    exnifc0: FloatField,
    umf_zint: FloatField,
    thlu: FloatField,
    thl0: FloatField,
    ssthl0: FloatField,
    pifc0: FloatField,
    pmid0: FloatField,
    qtu: FloatField,
    qt0: FloatField,
    ssqt0: FloatField,
    uu: FloatField,
    u0: FloatField,
    v0: FloatField,
    vu: FloatField,
    ssu0: FloatField,
    ssv0: FloatField,
    dotransport: Int,
    trflx: FloatField_NTracers,
    tru: FloatField_NTracers,
    tr0: FloatField_NTracers,
    sstr0: FloatField_NTracers,
    qtflx: FloatField,
    uflx: FloatField,
    vflx: FloatField,
    slflx: FloatField,
    iteration: i32,
):
    """
    3. Calculate buoyancy sorting fluxes : krel <= k <= kbup - 1
    In case that 'kbup = krel - 1 ' ( or even in case 'kbup = krel' ),
    buoyancy sorting fluxes are not calculated, which is consistent,
    desirable feature.

    Inputs:
    condensation [BoolFieldIJ]: Mask that indicates if condensation has occurred
    kbup [IntField]: Top layer in which buoyancy is positive at top interface
    krel [IntField]: Release layer where buoyancy sorting first occurs
    exnifc0 [FloatField]: Exner function at interfaces
    umf_zint [FloatField]: Updraft mass flux at the interfaces [kg/m2/s]
    thlu [FloatField]: Updraft liquid potential temperature at the interface [K]
    thl0 [FloatField]: Temperature [?]
    ssthl0 [FloatField]: [?]
    pifc0 [FloatField]: Environmental pressure at the interfaces [Pa]
    pmid0 [FloatField]: Environmental pressure at the layer mid-point [Pa]
    qtu [FloatField]: Updraft total specific humidity at the interface [kg/kg]
    qt0 [FloatField]: Mixing ratio [?]
    ssqt0 [FloatField]: [?]
    uu [FloatField]: Updraft zonal wind at the interface [m/s]
    u0 [FloatField]: Environmental zonal wind [m/s]
    v0 [FloatField]: Environmental meridional wind [m/s]
    vu [FloatField]: Updraft meridional wind at the interface [m/s]
    ssu0 [FloatField]: [?]
    ssv0 [FloatField]: [?]
    dotransport [Int]: Transport tracers [1 true]
    iteration [i32]: Iteration of implicit CIN loop (i.e., 0 or 1)
    tru [FloatField_NTracers]: Updraft tracers [#, kg/kg]
    tr0 [FloatField_NTracers]: Environmental tracers [#, kg/kg]
    sstr0 [FloatField_NTracers]: Convective tracer [?]

    Outputs:
    trflx [FloatField_NTracers]: Updraft/pen.entrainment tracer flux [#/m2/s, kg/kg/m2/s]
    qtflx [FloatField]: Mixing ratio flux [?]
    uflx [FloatField]: Zonal wind flux [?]
    vflx [FloatField]: Meridional wind flux [?]
    slflx [FloatField]: Sensible heat flux [?]
    """
    from __externals__ import ncnst

    with computation(FORWARD), interval(...):
        if not condensation:
            if THIS_K >= krel and THIS_K <= (kbup - 1):
                slflx[0, 0, 1] = (
                    constants.MAPL_CP
                    * exnifc0[0, 0, 1]
                    * umf_zint
                    * (
                        thlu
                        - (
                            thl0[0, 0, 1]
                            + ssthl0[0, 0, 1] * (pifc0[0, 0, 1] - pmid0[0, 0, 1])
                        )
                    )
                )
                qtflx[0, 0, 1] = umf_zint * (
                    qtu
                    - (
                        qt0[0, 0, 1]
                        + ssqt0[0, 0, 1] * (pifc0[0, 0, 1] - pmid0[0, 0, 1])
                    )
                )
                uflx[0, 0, 1] = umf_zint * (
                    uu
                    - (u0[0, 0, 1] + ssu0[0, 0, 1] * (pifc0[0, 0, 1] - pmid0[0, 0, 1]))
                )
                vflx[0, 0, 1] = umf_zint * (
                    vu
                    - (v0[0, 0, 1] + ssv0[0, 0, 1] * (pifc0[0, 0, 1] - pmid0[0, 0, 1]))
                )

                if dotransport == 1:
                    n = 0
                    while n < ncnst:
                        trflx[0, 0, 1][n] = umf_zint * (
                            tru[0, 0, 0][n]
                            - (
                                tr0[0, 0, 1][n]
                                + sstr0[0, 0, 1][n] * (pifc0[0, 0, 1] - pmid0[0, 0, 1])
                            )
                        )
                        n += 1


def penetrative_entrainment_fluxes(
    condensation: BoolFieldIJ,
    kbup: IntField,
    kpen: IntField,
    exnifc0: FloatField,
    emf: FloatField,
    thlu_emf: FloatField,
    thl0: FloatField,
    ssthl0: FloatField,
    pifc0: FloatField,
    pmid0: FloatField,
    qtu_emf: FloatField,
    qt0: FloatField,
    ssqt0: FloatField,
    uu_emf: FloatField,
    vu_emf: FloatField,
    u0: FloatField,
    v0: FloatField,
    ssu0: FloatField,
    ssv0: FloatField,
    dotransport: Int,
    trflx: FloatField_NTracers,
    tru_emf: FloatField_NTracers,
    tr0: FloatField_NTracers,
    sstr0: FloatField_NTracers,
    use_momenflx: Int,
    kinv: IntField,
    cbmf: FloatField,
    uflx: FloatField,
    vflx: FloatField,
    slflx: FloatField,
    qtflx: FloatField,
    uemf: FloatField,
    krel: IntField,
    umf_zint: FloatField,
    k0: Int,
    ql0: FloatField,
    qi0: FloatField,
    dt: Float,
    ese: FloatField_Extra_Dim,
    esx: FloatField_Extra_Dim,
    umf_out: FloatField,
    qtflx_out: FloatField,
    slflx_out: FloatField,
    uflx_out: FloatField,
    vflx_out: FloatField,
    qlten_sink: FloatField,
    qiten_sink: FloatField,
    cush_inout: FloatFieldIJ,
    iteration: i32,
):
    """
    4. Penetrative entrainment fluxes : kbup <= k <= kpen - 1
    The only confliction that can happen is when 'kbup = kinv-1'. For this
    case, turbulent flux at kinv-1 is calculated  both from 'fluxbelowinv'
    and here as penetrative entrainment fluxes.  Since penetrative flux is
    calculated later, flux at 'kinv - 1 ' will be that of penetrative flux.
    However, turbulent flux calculated at 'kinv - 1' from penetrative entr.
    is less attractable,  since more reasonable turbulent flux at 'kinv-1'
    should be obtained from 'fluxbelowinv', by considering  re-constructed
    inversion base height. This conflicting problem can be solved if we can
    initialize 'kbup = krel', instead of kbup = krel - 1. This choice seems
    to be more reasonable since it is not conflicted with 'fluxbelowinv' in
    calculating fluxes at 'kinv - 1' ( for this case, flux at 'kinv-1' is
    always from 'fluxbelowinv' ), and flux at 'krel-1' is calculated from
    the non-buoyancy sorting flux without being competed with penetrative
    entrainment fluxes. Even when we use normal cumulus flux instead of
    penetrative entrainment fluxes at 'kbup <= k <= kpen-1' interfaces,
    the initialization of kbup=krel perfectly works without any conceptual
    confliction. Thus it seems to be much better to choose 'kbup = krel'
    initialization of 'kbup', which is current choice.
    Note that below formula uses conventional updraft cumulus fluxes for
    shallow cumulus which did not overcome the first buoyancy barrier above
    PBL top while uses penetrative entrainment fluxes for the other cases
    'kbup <= k <= kpen-1' interfaces. Depending on cases, however, I can
    selelct different choice.

    Inputs:
    kbup [IntField]: Top layer in which buoyancy is positive at top interface
    kpen [IntField]: Highest layer with positive updraft velocity
    exnifc0 [FloatField]: Exner function at interfaces
    emf [FloatField]: Penetrative [?]
    thlu_emf [FloatField]: Penetrative downdraft liquid potential temperature at entraining interfaces [K]
    thl0 [FloatField]: Temperature [?]
    ssthl0 [FloatField]: [?]
    pifc0 [FloatField]: Environmental pressure at the interfaces [Pa]
    pmid0 [FloatField]: Environmental pressure at the layer mid-point [Pa]
    qtu_emf [FloatField]: Penetrative downdraft total water at entraining interfaces [kg/kg]
    qt0 [FloatField]: Mixing ratio [?]
    ssqt0 [FloatField]: [?]
    uu_emf [FloatField]: Penetrative downdraft zonal wind at entraining interfaces [m/s]
    vu_emf [FloatField]: Penetrative downdraft meridional wind at entraining interfaces [m/s]
    u0 [FloatField]: Environmental zonal wind [m/s]
    v0 [FloatField]: Environmental meridional wind [m/s]
    ssu0 [FloatField]: [?]
    ssv0 [FloatField]: [?]
    dotransport [Int]: Transport tracers [1 true]
    tru_emf [FloatField_NTracers]: Penetrative Downdraft tracers at entraining interfaces [#, kg/kg]
    tr0 [FloatField_NTracers]: Environmental tracers [#, kg/kg]
    sstr0 [FloatField_NTracers]: Convective tracer [?]
    use_momenflx [Int]: Perform momentum transport
    kinv [IntField]: Inversion layer with PBL top interface as lower interface
    cbmf [FloatField]: Cloud base mass flux [kg/m2/s]
    krel [IntField]: Release layer where buoyancy sorting first occurs
    umf_zint [FloatField]: Updraft mass flux at the interfaces [kg/m2/s]
    k0 [Int]: Number of levels
    ql0 [FloatField]: Environmental liquid water specific humidity
    qi0 [FloatField]: Environmental ice specific humidity
    dt [Float]: Timestep [s]
    ese [FloatField_Extra_Dim]: Used in QSat_Float [n/a]
    esx [FloatField_Extra_Dim]: Used in QSat_Float [n/a]
    iteration [i32]: Iteration of implicit CIN loop (i.e., 0 or 1)

    Inouts:
    condensation [BoolFieldIJ]: Mask that indicates if condensation has occurred

    Outputs:
    trflx [FloatField_NTracers]: Updraft/pen.entrainment tracer flux [#/m2/s, kg/kg/m2/s]
    uflx [FloatField]: Zonal wind flux [?]
    vflx [FloatField]: Meridional wind flux [?]
    slflx [FloatField]: Sensible heat flux [?]
    qtflx [FloatField]: Mixing ratio flux [?]
    uemf [FloatField]: Net updraft mass flux at the interface ( emf + umf ) [kg/m2/s]
    umf_out [FloatField]: Updraft mass flux at the interfaces [kg/m2/s]
    qtflx_out [FloatField]: Mixing ratio flux [?]
    slflx_out [FloatField]: Sensible heat flux [?]
    uflx_out [FloatField]: Zonal wind flux [?]
    vflx_out [FloatField]: Meridional wind flux [?]
    qlten_sink [FloatField]: Liquid condensate tendency by compensating subsidence/upwelling [kg/kg/s]
    qiten_sink [FloatField]: Ice condensate tendency by compensating subsidence/upwelling [kg/kg/s]
    cush_inout [FloatFieldIJ]: Convective scale height [m]
    """
    from __externals__ import ncnst

    with computation(FORWARD), interval(...):
        if not condensation:
            if THIS_K >= kbup and THIS_K <= kpen - 1:
                slflx[0, 0, 1] = (
                    constants.MAPL_CP
                    * exnifc0[0, 0, 1]
                    * emf
                    * (thlu_emf - (thl0 + ssthl0 * (pifc0[0, 0, 1] - pmid0)))
                )
                qtflx[0, 0, 1] = emf * (
                    qtu_emf - (qt0 + ssqt0 * (pifc0[0, 0, 1] - pmid0))
                )
                uflx[0, 0, 1] = emf * (uu_emf - (u0 + ssu0 * (pifc0[0, 0, 1] - pmid0)))
                vflx[0, 0, 1] = emf * (vu_emf - (v0 + ssv0 * (pifc0[0, 0, 1] - pmid0)))
                if dotransport == 1:
                    n = 0
                    while n < ncnst:
                        trflx[0, 0, 1][n] = emf * (
                            tru_emf[0, 0, 0][n]
                            - (
                                tr0[0, 0, 0][n]
                                + sstr0[0, 0, 0][n] * (pifc0[0, 0, 1] - pmid0)
                            )
                        )
                        n += 1

            # Turn-off cumulus momentum flux as an option
            if use_momenflx == 0:
                uflx = 0.0
                vflx = 0.0
                uflx[0, 0, 1] = 0.0
                vflx[0, 0, 1] = 0.0

            # Condensate tendency by compensating subsidence/upwelling
            uemf[0, 0, 1] = 0.0

            if THIS_K >= 0 and THIS_K <= (
                kinv - 2
            ):  # Assume linear updraft mass flux within the PBL.
                uemf[0, 0, 1] = (
                    cbmf
                    * (pifc0.at(K=0) - pifc0[0, 0, 1])
                    / (pifc0.at(K=0) - pifc0.at(K=kinv))
                )

            if THIS_K >= (kinv - 1) and THIS_K <= (krel - 1):
                uemf[0, 0, 1] = cbmf
            if THIS_K >= krel and THIS_K <= (kbup - 1):
                uemf[0, 0, 1] = umf_zint
            if THIS_K >= kbup and THIS_K <= (kpen - 1):
                uemf[0, 0, 1] = (
                    emf  # Only use penetrative entrainment flux consistently.
                )

            comsub = 0.0

            if THIS_K <= kpen:
                comsub = 0.5 * (uemf[0, 0, 1] + uemf)  # comsub defined on interfaces

    with computation(FORWARD), interval(0, 1):
        if not condensation:
            if THIS_K <= kpen:
                if comsub < 0.0:
                    thlten_sub = 0.0
                    qtten_sub = 0.0
                    qlten_sub = 0.0
                    qiten_sub = 0.0
                    nlten_sub = 0.0
                    niten_sub = 0.0

    with computation(FORWARD), interval(0, -1):
        if not condensation:
            if THIS_K <= kpen:
                if comsub >= 0.0:
                    thlten_sub = (
                        constants.MAPL_GRAV
                        * comsub
                        * (thl0[0, 0, 1] - thl0)
                        / (pmid0 - pmid0[0, 0, 1])
                    )
                    qtten_sub = (
                        constants.MAPL_GRAV
                        * comsub
                        * (qt0[0, 0, 1] - qt0)
                        / (pmid0 - pmid0[0, 0, 1])
                    )
                    qlten_sub = (
                        constants.MAPL_GRAV
                        * comsub
                        * (ql0[0, 0, 1] - ql0)
                        / (pmid0 - pmid0[0, 0, 1])
                    )
                    qiten_sub = (
                        constants.MAPL_GRAV
                        * comsub
                        * (qi0[0, 0, 1] - qi0)
                        / (pmid0 - pmid0[0, 0, 1])
                    )

    with computation(FORWARD), interval(1, None):
        if not condensation:
            if THIS_K <= kpen:
                if comsub < 0.0:
                    thlten_sub = (
                        constants.MAPL_GRAV
                        * comsub
                        * (thl0 - thl0[0, 0, -1])
                        / (pmid0[0, 0, -1] - pmid0)
                    )
                    qtten_sub = (
                        constants.MAPL_GRAV
                        * comsub
                        * (qt0 - qt0[0, 0, -1])
                        / (pmid0[0, 0, -1] - pmid0)
                    )
                    qlten_sub = (
                        constants.MAPL_GRAV
                        * comsub
                        * (ql0 - ql0[0, 0, -1])
                        / (pmid0[0, 0, -1] - pmid0)
                    )
                    qiten_sub = (
                        constants.MAPL_GRAV
                        * comsub
                        * (qi0 - qi0[0, 0, -1])
                        / (pmid0[0, 0, -1] - pmid0)
                    )

    with computation(FORWARD), interval(...):
        if not condensation:
            if THIS_K <= kpen:
                if comsub >= 0.0:
                    if THIS_K == k0 - 1:
                        thlten_sub = 0.0
                        qtten_sub = 0.0
                        qlten_sub = 0.0
                        qiten_sub = 0.0
                        nlten_sub = 0.0
                        niten_sub = 0.0

    with computation(FORWARD), interval(...):
        if not condensation:
            if THIS_K <= kpen:
                thl_prog = thl0 + thlten_sub * dt
                qt_prog = max(qt0 + qtten_sub * dt, 1.0e-12)
                thj, qvj, qlj, qij, qse, id_check = conden(
                    pmid0, thl_prog, qt_prog, ese, esx
                )

                if id_check == 1:
                    condensation = True
                    umf_out[0, 0, 1] = 0.0
                    dcm_out = 0.0
                    qvten_out = 0.0
                    qlten_out = 0.0
                    qiten_out = 0.0
                    sten_out = 0.0
                    uten_out = 0.0
                    vten_out = 0.0
                    qrten_out = 0.0
                    qsten_out = 0.0
                    cufrc_out = 0.0
                    cush_inout = -1.0
                    qldet_out = 0.0
                    qidet_out = 0.0
                    qtflx_out[0, 0, 1] = 0.0
                    slflx_out[0, 0, 1] = 0.0
                    uflx_out[0, 0, 1] = 0.0
                    vflx_out[0, 0, 1] = 0.0
                    fer_out = constants.MAPL_UNDEF
                    fdr_out = constants.MAPL_UNDEF

                if not condensation:
                    qlten_sink = max(
                        qlten_sub, -ql0 / dt
                    )  # For consistency with prognostic macrophysics scheme
                    qiten_sink = max(
                        qiten_sub, -qi0 / dt
                    )  # For consistency with prognostic macrophysics scheme


def calc_momentum_tendency(
    condensation: BoolFieldIJ,
    kpen: IntField,
    uflx: FloatField,
    vflx: FloatField,
    dp0: FloatField,
    u0: FloatField,
    v0: FloatField,
    dt: Float,
    uf: FloatField,
    vf: FloatField,
    uten: FloatField,
    vten: FloatField,
):
    """
    Stencil to calculate momentum tendency at each layer.

    Inputs:
    kpen [IntField]: Highest layer with positive updraft velocity
    uflx [FloatField]: Zonal wind flux [?]
    vflx [FloatField]: Meridional wind flux [?]
    dp0 [FloatField]: Environmental layer pressure thickness [Pa] > 0
    u0 [FloatField]: Environmental zonal wind [m/s]
    v0 [FloatField]: Environmental meridional wind [m/s]
    dt [Float]: Timestep [s]

    Inouts:
    condensation [BoolFieldIJ]: Mask that indicates if condensation has occurred

    Outputs:
    uf [FloatField]: [?]
    vf [FloatField]: [?]
    uten [FloatField]: Tendency of zonal wind [m/s2]
    vten [FloatField]: Tendency of meridional wind [m/s2]
    """
    with computation(FORWARD), interval(...):
        if not condensation:
            # Calculate convective tendencies at each layer
            # Momentum tendency
            if THIS_K <= kpen:
                uten = (uflx - uflx[0, 0, 1]) * constants.MAPL_GRAV / dp0
                vten = (vflx - vflx[0, 0, 1]) * constants.MAPL_GRAV / dp0
                uf = u0 + uten * dt
                vf = v0 + vten * dt


def calc_thermodynamic_tendencies(
    condensation: BoolFieldIJ,
    kpen: IntField,
    umf_zint: FloatField,
    dp0: FloatField,
    frc_rasn: Float,
    slflx: FloatField,
    uflx: FloatField,
    vflx: FloatField,
    qtflx: FloatField,
    u0: FloatField,
    v0: FloatField,
    uf: FloatField,
    vf: FloatField,
    dwten: FloatField,
    diten: FloatField,
    dwten_temp: FloatField,
    diten_temp: FloatField,
    umf_temp: FloatField,
    krel: IntField,
    prel: FloatField,
    thlu: FloatField,
    qtu: FloatField,
    ese: FloatField_Extra_Dim,
    esx: FloatField_Extra_Dim,
    umf_out: FloatField,
    qtflx_out: FloatField,
    slflx_out: FloatField,
    uflx_out: FloatField,
    vflx_out: FloatField,
    pifc0: FloatField,
    ppen: FloatFieldIJ,
    thlu_top: FloatFieldIJ,
    qtu_top: FloatFieldIJ,
    qlubelow: FloatFieldIJ,
    qiubelow: FloatFieldIJ,
    qlj_2D: FloatFieldIJ,
    qij_2D: FloatFieldIJ,
    kbup: IntField,
    fdr: FloatField,
    ql0: FloatField,
    qi0: FloatField,
    pmid0: FloatField,
    thlu_emf: FloatField,
    qtu_emf: FloatField,
    emf: FloatField,
    dt: Float,
    qlten_sink: FloatField,
    qiten_sink: FloatField,
    qrten: FloatField,
    qsten: FloatField,
    qvten: FloatField,
    qlten: FloatField,
    sten: FloatField,
    qiten: FloatField,
    qc: FloatField,
    qlten_det: FloatField,
    qiten_det: FloatField,
    slten: FloatField,
    cush_inout: FloatFieldIJ,
):
    """
    Tendencies of thermodynamic variables.
    This part requires a careful treatment of bulk cloud microphysics.
    Relocations of 'precipitable condensates' either into the surface
    or into the tendency of 'krel' layer will be performed just after
    finishing the below 'do-loop'.

    Inputs:
    kpen [IntField]: Highest layer with positive updraft velocity
    dp0 [FloatField]: Environmental layer pressure thickness [Pa] > 0
    frc_rasn [Float]: Precipitation fraction of expelled condensate
    slflx [FloatField]: Sensible heat flux [?]
    uflx [FloatField]: Zonal wind flux [?]
    vflx [FloatField]: Meridional wind flux [?]
    qtflx [FloatField]: Mixing ratio flux [?]
    u0 [FloatField]: Environmental zonal wind [m/s]
    v0 [FloatField]: Environmental meridional wind [m/s]
    uf [FloatField]: [?]
    vf [FloatField]: [?]
    umf_temp [FloatField]: Updraft mass flux at the interfaces [kg/m2/s]
    krel [IntField]: Release layer where buoyancy sorting first occurs
    prel [FloatField]: Lowest level from which buoyancy sorting occurs [Pa]
    thlu [FloatField]: Updraft liquid potential temperature at the interface [K]
    qtu [FloatField]: Updraft total specific humidity at the interface [kg/kg]
    ese [FloatField_Extra_Dim]: Used in QSat_Float [n/a]
    esx [FloatField_Extra_Dim]: Used in QSat_Float [n/a]
    pifc0 [FloatField]: Environmental pressure at the interfaces [Pa]
    ppen [FloatFieldIJ]: Highest interface evel where Cu w = 0 [Pa]
    thlu_top [FloatFieldIJ]: Updraft liquid potential temperature at top [K]
    qtu_top [FloatFieldIJ]: Updraft total specific humidity at top [kg/kg]
    kbup [IntField]: Top layer in which buoyancy is positive at top interface
    fdr [FloatField]: Fractional lateral detrainment rate [1/Pa]
    ql0 [FloatField]: Environmental liquid water specific humidity
    qi0 [FloatField]: Environmental ice specific humidity
    pmid0 [FloatField]: Environmental pressure at the layer mid-point [Pa]
    thlu_emf [FloatField]: Penetrative downdraft liquid potential temperature at entraining interfaces [K]
    qtu_emf [FloatField]: Penetrative downdraft total water at entraining interfaces [kg/kg]
    emf [FloatField]: Penetrative [?]
    dt [Float]: Timestep [s]
    qlten_sink [FloatField]: Liquid condensate tendency by compensating subsidence/upwelling [kg/kg/s]
    qiten_sink [FloatField]: Ice condensate tendency by compensating subsidence/upwelling [kg/kg/s]

    Inouts:
    condensation [BoolFieldIJ]: Mask that indicates if condensation has occurred

    Outputs:
    umf_zint [FloatField]: Updraft mass flux at the interfaces [kg/m2/s]
    dwten [FloatField]: Detrained [?]
    diten [FloatField]: Detrained [?]
    dwten_temp [FloatField]: Detrained [?]
    diten_temp [FloatField]: Detrained [?]
    umf_out [FloatField]: Updraft mass flux at the interfaces [kg/m2/s]
    qtflx_out [FloatField]: Mixing ratio flux [?]
    slflx_out [FloatField]: Sensible heat flux [?]
    uflx_out [FloatField]: Zonal wind flux [?]
    vflx_out [FloatField]: Meridional wind flux [?]
    qlubelow [FloatFieldIJ]: [?]
    qiubelow [FloatFieldIJ]: [?]
    qlj_2D [FloatFieldIJ]: [?]
    qij_2D [FloatFieldIJ]: [?]
    qrten [FloatField]: Tendency of rain water specific humidity [kg/kg/s]
    qsten [FloatField]: Tendency of snow specific humidity [kg/kg/s]
    qvten [FloatField]: Tendency of water vapor specific humidity [kg/kg/s]
    qlten [FloatField]: Tendency of liquid water specific humidity [kg/kg/s]
    sten [FloatField]: Tendency of dry static energy [J/kg/s]
    qiten [FloatField]: Tendency of ice specific humidity [kg/kg/s]
    qc [FloatField]: Tendency due [?]
    qlten_det [FloatField]: [?]
    qiten_det [FloatField]: [?]
    slten [FloatField]: Tendency of [?]
    cush_inout [FloatFieldIJ]: Convective scale height [m]
    """
    with computation(FORWARD), interval(...):
        if not condensation:

            umf_zint[0, 0, 1] = umf_temp[0, 0, 1]  # Update umf

            if THIS_K <= kpen:
                # Compute 'slten', 'qtten', 'qvten', 'qlten', 'qiten', and 'sten'

                # Key assumptions made in this 'cumulus scheme' are :
                # 1. Cumulus updraft expels condensate into the environment at the top
                #   interface of each layer. Note that in addition to this expel process
                #   ('source' term), cumulus updraft can modify layer mean condensate
                #   through normal detrainment forcing or compensating subsidence.
                # 2. Expelled water can be either 'sustaining' or 'precipitating'
                #   condensate. By definition, 'suataining condensate' will remain in
                #   the layer where it was formed, while 'precipitating condensate'
                #   will fall across the base of the layer where it was formed.
                # 3. All precipitating condensates are assumed to fall into the release
                #   layer or ground as soon as it was formed without being evaporated
                #   during the falling process down to the desinated layer ( either
                #   release layer of surface ).

                # 'dwten(k)','diten(k)' : Production rate of condensate  within the
                # layer k [ kg/kg/s ] by the expels of condensate from cumulus updraft.
                # It is important to note that in terms of moisture tendency equation,
                # this is a 'source' term of enviromental 'qt'.  More importantly,
                # these source are already counted in the turbulent heat and moisture
                # fluxes we computed until now, assuming all the expelled condensate
                # remain in the layer where it was formed. Thus, in calculation of
                # 'qtten' and 'slten' below, we MUST NOT add or subtract these terms
                # explicitly in order not to double or miss count, unless some expelled
                # condensates fall down out of the layer.  Note this falling-down
                # process ( i.e., precipitation process ) and  associated 'qtten' and
                # 'slten' and production of surface precipitation flux  will be treated
                # later in 'zm_conv_evap' in 'convect_shallow_tend' subroutine. In below
                # we are converting expelled cloud condensate into correct unit. I found
                # that below use of '0.5 * (umf(k-1) + umf(k))' causes conservation
                # errors at some columns in global simulation. So, I returned to
                # originals. This will cause no precipitation flux at 'kpen' layer since
                # umf(kpen)=0.

                dwten = (
                    dwten_temp
                    * 0.5
                    * (umf_zint + umf_zint[0, 0, 1])
                    * constants.MAPL_GRAV
                    / dp0
                )  # [ kg/kg/s ]

                diten = (
                    diten_temp
                    * 0.5
                    * (umf_zint + umf_zint[0, 0, 1])
                    * constants.MAPL_GRAV
                    / dp0
                )  # [ kg/kg/s ]

                # 'qrten(k)','qsten(k)' : Production rate of rain and snow within the
                # layer k [ kg/kg/s ] by cumulus expels of condensates to the
                # environment.
                qrten = frc_rasn * dwten
                qsten = frc_rasn * diten

                # 'slten(k)','qtten(k)'
                # Note that 'slflx(k)' and 'qtflx(k)' we have calculated already
                # included all the contributions of (1) expels of condensate
                # (dwten(k), diten(k)), (2) mass detrainment
                # ( delta * umf * ( qtu - qt ) ), & (3) compensating subsidence
                # ( M * dqt / dz ). Thus 'slflx(k)' and 'qtflx(k)' we computed is a
                # hybrid turbulent flux containing one part of 'source' term -
                # expel of condensate. In order to calculate 'slten' and 'qtten',
                # we should add additional 'source' term, if any. If the expelled
                # condensate falls down across the base of the layer, it will be
                # another sink (negative source) term.  Note also that we included
                # frictional heating terms in the below calculation of 'slten'.

                slten = (slflx - slflx[0, 0, 1]) * constants.MAPL_GRAV / dp0

                if THIS_K == 0:
                    slten = slten - constants.MAPL_GRAV / 4.0 / dp0 * (
                        uflx[0, 0, 1] * (uf[0, 0, 1] - uf + u0[0, 0, 1] - u0)
                        + vflx[0, 0, 1] * (vf[0, 0, 1] - vf + v0[0, 0, 1] - v0)
                    )

                elif THIS_K >= 1 and THIS_K <= (kpen - 1):
                    slten = slten - constants.MAPL_GRAV / 4.0 / dp0 * (
                        (
                            uflx[0, 0, 1] * (uf[0, 0, 1] - uf + u0[0, 0, 1] - u0)
                            + uflx
                            * (uf - uf.at(K=THIS_K - 1) + u0 - u0.at(K=THIS_K - 1))
                        )
                        + (
                            vflx[0, 0, 1] * (vf[0, 0, 1] - vf + v0[0, 0, 1] - v0)
                            + vflx
                            * (vf - vf.at(K=THIS_K - 1) + v0 - v0.at(K=THIS_K - 1))
                        )
                    )

                elif THIS_K == kpen:
                    slten = slten - constants.MAPL_GRAV / 4.0 / dp0 * (
                        uflx * (uf - uf.at(K=THIS_K - 1) + u0 - u0.at(K=THIS_K - 1))
                        + vflx * (vf - vf.at(K=THIS_K - 1) + v0 - v0.at(K=THIS_K - 1))
                    )

                qtten = (qtflx - qtflx[0, 0, 1]) * constants.MAPL_GRAV / dp0

                # Compute condensate tendency, including reserved condensate
                # We assume that eventual detachment and detrainment occurs in kbup
                # layer  due to downdraft buoyancy sorting. In the layer above the
                # kbup, only penetrative entrainment exists. Penetrative entrained
                # air is assumed not to contain any condensate.

                # Compute in-cumulus condensate at the layer mid-point.
                if THIS_K < krel or THIS_K > kpen:
                    qlu_mid = 0.0
                    qiu_mid = 0.0
                    qlj_2D = 0.0
                    qij_2D = 0.0
                elif THIS_K == krel:
                    thj, qvj, qlj_2D, qij_2D, qse, id_check = conden(
                        prel, thlu.at(K=THIS_K - 1), qtu.at(K=THIS_K - 1), ese, esx
                    )

                    if id_check == 1:
                        condensation = True
                        umf_out[0, 0, 1] = 0.0
                        dcm_out = 0.0
                        qvten_out = 0.0
                        qlten_out = 0.0
                        qiten_out = 0.0
                        sten_out = 0.0
                        uten_out = 0.0
                        vten_out = 0.0
                        qrten_out = 0.0
                        qsten_out = 0.0
                        cufrc_out = 0.0
                        cush_inout = -1.0
                        qldet_out = 0.0
                        qidet_out = 0.0
                        qtflx_out[0, 0, 1] = 0.0
                        slflx_out[0, 0, 1] = 0.0
                        uflx_out[0, 0, 1] = 0.0
                        vflx_out[0, 0, 1] = 0.0
                        fer_out = constants.MAPL_UNDEF
                        fdr_out = constants.MAPL_UNDEF
                    if not condensation:
                        qlubelow = qlj_2D
                        qiubelow = qij_2D
                        thj, qvj, qlj_2D, qij_2D, qse, id_check = conden(
                            pifc0[0, 0, 1], thlu, qtu, ese, esx
                        )

                        if id_check == 1:
                            condensation = True
                            umf_out[0, 0, 1] = 0.0
                            dcm_out = 0.0
                            qvten_out = 0.0
                            qlten_out = 0.0
                            qiten_out = 0.0
                            sten_out = 0.0
                            uten_out = 0.0
                            vten_out = 0.0
                            qrten_out = 0.0
                            qsten_out = 0.0
                            cufrc_out = 0.0
                            cush_inout = -1.0
                            qldet_out = 0.0
                            qidet_out = 0.0
                            qtflx_out[0, 0, 1] = 0.0
                            slflx_out[0, 0, 1] = 0.0
                            uflx_out[0, 0, 1] = 0.0
                            vflx_out[0, 0, 1] = 0.0
                            fer_out = constants.MAPL_UNDEF
                            fdr_out = constants.MAPL_UNDEF

                        if not condensation:
                            qlu_mid = (
                                0.5
                                * (qlubelow + qlj_2D)
                                * (prel - pifc0[0, 0, 1])
                                / (pifc0 - pifc0[0, 0, 1])
                            )
                            qiu_mid = (
                                0.5
                                * (qiubelow + qij_2D)
                                * (prel - pifc0[0, 0, 1])
                                / (pifc0 - pifc0[0, 0, 1])
                            )

                elif THIS_K == kpen:
                    if not condensation:
                        thj, qvj, qlj_2D, qij_2D, qse, id_check = conden(
                            pifc0 + ppen, thlu_top, qtu_top, ese, esx
                        )

                        if id_check == 1:
                            condensation = True
                            umf_out[0, 0, 1] = 0.0
                            dcm_out = 0.0
                            qvten_out = 0.0
                            qlten_out = 0.0
                            qiten_out = 0.0
                            sten_out = 0.0
                            uten_out = 0.0
                            vten_out = 0.0
                            qrten_out = 0.0
                            qsten_out = 0.0
                            cufrc_out = 0.0
                            cush_inout = -1.0
                            qldet_out = 0.0
                            qidet_out = 0.0
                            qtflx_out[0, 0, 1] = 0.0
                            slflx_out[0, 0, 1] = 0.0
                            uflx_out[0, 0, 1] = 0.0
                            vflx_out[0, 0, 1] = 0.0
                            fer_out = constants.MAPL_UNDEF
                            fdr_out = constants.MAPL_UNDEF
                        if not condensation:
                            qlu_mid = (
                                0.5
                                * (qlubelow + qlj_2D)
                                * (-ppen)
                                / (pifc0 - pifc0[0, 0, 1])
                            )

                            qiu_mid = (
                                0.5
                                * (qiubelow + qij_2D)
                                * (-ppen)
                                / (pifc0 - pifc0[0, 0, 1])
                            )
                            qlu_top = qlj_2D
                            qiu_top = qij_2D

                else:
                    if not condensation:
                        thj, qvj, qlj_2D, qij_2D, qse, id_check = conden(
                            pifc0[0, 0, 1], thlu, qtu, ese, esx
                        )

                        if id_check == 1:
                            condensation = True
                            umf_out[0, 0, 1] = 0.0
                            dcm_out = 0.0
                            qvten_out = 0.0
                            qlten_out = 0.0
                            qiten_out = 0.0
                            sten_out = 0.0
                            uten_out = 0.0
                            vten_out = 0.0
                            qrten_out = 0.0
                            qsten_out = 0.0
                            cufrc_out = 0.0
                            cush_inout = -1.0
                            qldet_out = 0.0
                            qidet_out = 0.0
                            qtflx_out[0, 0, 1] = 0.0
                            slflx_out[0, 0, 1] = 0.0
                            uflx_out[0, 0, 1] = 0.0
                            vflx_out[0, 0, 1] = 0.0
                            fer_out = constants.MAPL_UNDEF
                            fdr_out = constants.MAPL_UNDEF
                        if not condensation:
                            qlu_mid = 0.5 * (qlubelow + qlj_2D)
                            qiu_mid = 0.5 * (qiubelow + qij_2D)

                if not condensation:
                    qlubelow = qlj_2D
                    qiubelow = qij_2D

                    # 1. Non-precipitating portion of expelled condensate
                    qc_l = (1.0 - frc_rasn) * dwten  # [ kg/kg/s ]
                    qc_i = (1.0 - frc_rasn) * diten  # [ kg/kg/s ]

                    # 2. Detrained Condensate
                    if THIS_K <= kbup:
                        qc_l = (
                            qc_l
                            + constants.MAPL_GRAV
                            * 0.5
                            * (umf_zint + umf_zint[0, 0, 1])
                            * fdr
                            * qlu_mid
                        )  # [ kg/kg/s ]
                        qc_i = (
                            qc_i
                            + constants.MAPL_GRAV
                            * 0.5
                            * (umf_zint + umf_zint[0, 0, 1])
                            * fdr
                            * qiu_mid
                        )  # [ kg/kg/s ]
                        qc_lm = (
                            -constants.MAPL_GRAV
                            * 0.5
                            * (umf_zint + umf_zint[0, 0, 1])
                            * fdr
                            * ql0
                        )
                        qc_im = (
                            -constants.MAPL_GRAV
                            * 0.5
                            * (umf_zint + umf_zint[0, 0, 1])
                            * fdr
                            * qi0
                        )

                    else:
                        qc_lm = 0.0
                        qc_im = 0.0
                        nc_lm = 0.0
                        nc_im = 0.0

                    # 3. Detached Updraft
                    if THIS_K == kbup:
                        qc_l = qc_l + constants.MAPL_GRAV * umf_zint[
                            0, 0, 1
                        ] * qlj_2D / (
                            pifc0 - pifc0[0, 0, 1]
                        )  # [ kg/kg/s ]
                        qc_i = qc_i + constants.MAPL_GRAV * umf_zint[
                            0, 0, 1
                        ] * qij_2D / (
                            pifc0 - pifc0[0, 0, 1]
                        )  # [ kg/kg/s ]
                        qc_lm = qc_lm - constants.MAPL_GRAV * umf_zint[
                            0, 0, 1
                        ] * ql0 / (
                            pifc0 - pifc0[0, 0, 1]
                        )  # [ kg/kg/s ]
                        qc_im = qc_im - constants.MAPL_GRAV * umf_zint[
                            0, 0, 1
                        ] * qi0 / (
                            pifc0 - pifc0[0, 0, 1]
                        )  # [ kg/kg/s ]

                    # 4. Cumulative Penetrative entrainment detrained in the 'kbup'
                    # layer. Explicitly compute the properties detrained penetrative
                    # entrained airs in k = kbup layer.

                if THIS_K == kbup:
                    thj, qvj, ql_emf_kbup, qi_emf_kbup, qse, id_check = conden(
                        pmid0,
                        thlu_emf,
                        qtu_emf,
                        ese,
                        esx,
                    )

                    if id_check == 1:
                        condensation = True
                        umf_out[0, 0, 1] = 0.0
                        dcm_out = 0.0
                        qvten_out = 0.0
                        qlten_out = 0.0
                        qiten_out = 0.0
                        sten_out = 0.0
                        uten_out = 0.0
                        vten_out = 0.0
                        qrten_out = 0.0
                        qsten_out = 0.0
                        cufrc_out = 0.0
                        cush_inout = -1.0
                        qldet_out = 0.0
                        qidet_out = 0.0
                        qtflx_out[0, 0, 1] = 0.0
                        slflx_out[0, 0, 1] = 0.0
                        uflx_out[0, 0, 1] = 0.0
                        vflx_out[0, 0, 1] = 0.0
                        fer_out = constants.MAPL_UNDEF
                        fdr_out = constants.MAPL_UNDEF

                    if not condensation:
                        if ql_emf_kbup < 0.0:
                            nl_emf_kbup = 0.0

                        if qi_emf_kbup < 0.0:
                            ni_emf_kbup = 0.0

                        qc_lm = qc_lm - constants.MAPL_GRAV * emf * (
                            ql_emf_kbup - ql0
                        ) / (
                            pifc0 - pifc0[0, 0, 1]
                        )  # [ kg/kg/s ]
                        qc_im = qc_im - constants.MAPL_GRAV * emf * (
                            qi_emf_kbup - qi0
                        ) / (
                            pifc0 - pifc0[0, 0, 1]
                        )  # [ kg/kg/s ]

                if not condensation:
                    qlten_det = qc_l + qc_lm
                    qiten_det = qc_i + qc_im

                    if ((qc_lm + qlten_sink) * dt + ql0) < 0.0:
                        totsink = qc_lm + qlten_sink
                        if totsink != 0.0:
                            qc_lm = -(ql0 / dt) * qc_lm / totsink
                            qlten_sink = -(ql0 / dt) * qlten_sink / totsink
                            qlten_det = qc_l + qc_lm

                    if ((qc_im + qiten_sink) * dt + qi0) < 0.0:
                        totsink = qc_im + qiten_sink
                        if totsink != 0.0:
                            qc_im = -(qi0 / dt) * qc_im / totsink
                            qiten_sink = -(qi0 / dt) * qiten_sink / totsink
                            qiten_det = qc_i + qc_im

                    qlten = qrten + qlten_sink + qlten_det
                    qiten = qsten + qiten_sink + qiten_det

                    qvten = qtten - qlten - qiten

                    sten = (
                        slten
                        + constants.MAPL_ALHL * qlten
                        + constants.MAPL_ALHS * qiten
                    )

                    qc = qc_l + qc_i

                    qlten = qlten - qrten
                    qiten = qiten - qsten
                    qtten = qlten + qiten + qvten

                    slten = (
                        sten - constants.MAPL_ALHL * qlten - constants.MAPL_ALHS * qiten
                    )
                    slten = (
                        slten
                        + constants.MAPL_ALHL * qrten
                        + constants.MAPL_ALHS * qsten
                    )
                    sten = (
                        slten
                        + constants.MAPL_ALHL * qlten
                        + constants.MAPL_ALHS * qiten
                    )


def prevent_negative_condensate(
    condensation: BoolFieldIJ,
    qv0: FloatField,
    dt: Float,
    qvten: FloatField,
    ql0: FloatField,
    qlten: FloatField,
    qi0: FloatField,
    s0: FloatField,
    sten: FloatField,
    dp0: FloatField,
    qiten: FloatField,
    k0: Int,
    qmin: FloatField,
):

    with computation(FORWARD), interval(...):
        if not condensation:
            # Prevent the onset-of negative condensate at the next time step
            # Potentially, this block can be moved just in front of the above
            # block.

            qv0_star = qv0 + qvten * dt
            ql0_star = ql0 + qlten * dt
            qi0_star = qi0 + qiten * dt
            s0_star = s0 + sten * dt

            qmin = 0.0

    with computation(BACKWARD), interval(...):
        # Calculate positive_moisture_single
        if not condensation:
            ixcldice = 1
            ixcldliq = 2
            dql: f64 = max(f64(0.0), f64(1.0) * qmin.at(K=ixcldliq) - ql0_star)
            dqi: f64 = max(f64(0.0), f64(1.0) * qmin.at(K=ixcldice) - qi0_star)
            qlten = qlten + dql / dt
            qiten = qiten + dqi / dt
            qvten = qvten - (dql + dqi) / dt
            sten = (
                sten
                + constants.MAPL_LATENT_HEAT_VAPORIZATION * (dql / dt)
                + constants.MAPL_LATENT_HEAT_SUBLIMATION * (dqi / dt)
            )
            ql0_star = ql0_star + dql
            qi0_star = qi0_star + dqi
            qv0_star = qv0_star - dql - dqi

            s0_star = (
                s0_star
                + constants.MAPL_LATENT_HEAT_VAPORIZATION * dql
                + constants.MAPL_LATENT_HEAT_SUBLIMATION * dqi
            )

            dqv = max(0.0, 1.0 * qmin.at(K=0) - qv0_star)
            qvten = qvten + dqv / dt
            qv0_star = qv0_star + dqv

    with computation(BACKWARD), interval(1, None):
        if not condensation:
            qv0_star[0, 0, -1] = qv0_star[0, 0, -1] - dqv * dp0 / dp0[0, 0, -1]
            qvten[0, 0, -1] = qvten[0, 0, -1] - dqv * dp0 / dp0[0, 0, -1] / dt

    with computation(BACKWARD), interval(...):
        if not condensation:
            qv0_star = max(qv0_star, qmin)
            ql0_star = max(ql0_star, qmin)
            qi0_star = max(qi0_star, qmin)

    with computation(PARALLEL), interval(...):
        if not condensation:
            # Extra moisture used to satisfy 'qv(i,1)=qvmin' is proportionally
            # extracted from all the layers that has 'qv > 2*qvmin'. This fully
            # preserves column moisture.
            if dqv > f64(1.0e-20):
                sum: f64 = 0.0
                if THIS_K <= k0:
                    if qv0_star > f64(2.0) * qmin.at(K=0):
                        sum = sum + qv0_star * dp0
                aa: f64 = dqv * dp0.at(K=1) / max(f64(1.0e-20), sum)
                if aa < f64(0.5):
                    if THIS_K <= k0 - 1:
                        if qv0_star > f64(2.0) * qmin.at(K=0):
                            dum: f64 = aa * qv0_star
                            qv0_star = qv0_star - dum
                            qvten = qvten - dum / dt

                # else:
                # print('Full positive_moisture is impossible in uwshcu')

    with computation(FORWARD), interval(...):
        if not condensation:
            qtten = qvten + qlten + qiten
            slten = (
                sten
                - constants.MAPL_LATENT_HEAT_VAPORIZATION * qlten
                - constants.MAPL_LATENT_HEAT_SUBLIMATION * qiten
            )


def calc_tracer_tendencies(
    condensation: BoolFieldIJ,
    dotransport: Int,
    k0: Int,
    dt: Float,
    dp0: FloatField,
    trflx_d: FloatField_NTracers,
    trflx_u: FloatField_NTracers,
    trmin: FloatFieldIJ_NTracers,
    tr0: FloatField_NTracers,
    trflx: FloatField_NTracers,
    trten: FloatField_NTracers,
    iteration: i32,
):
    from __externals__ import ncnst

    with computation(FORWARD), interval(...):
        if not condensation:
            # Tendencies of tracers
            if dotransport == 1:
                n = 0
                while n < ncnst:
                    trmin[0, 0][n] = 0.0
                    trflx_d[0, 0, 1][n] = 0.0
                    trflx_d[0, 0, 0][n] = 0.0
                    trflx_u[0, 0, 1][n] = 0.0
                    trflx_u[0, 0, 0][n] = 0.0
                    n += 1

    with computation(FORWARD), interval(1, -1):
        if not condensation:
            if dotransport == 1:
                n = 0
                while n < ncnst:
                    pdelx = dp0
                    dum = (
                        (tr0[0, 0, 0][n] - trmin[0, 0][n])
                        * pdelx
                        / constants.MAPL_GRAV
                        / dt
                        + trflx[0, 0, 0][n]
                        - trflx[0, 0, 1][n]
                        + trflx_d[0, 0, -1][n]
                    )
                    trflx_d[0, 0, 0][n] = min(0.0, dum)
                    n += 1

    with computation(BACKWARD), interval(2, None):
        if not condensation:
            if dotransport == 1:
                n = 0
                while n < ncnst:
                    pdelx = dp0
                    dum = (
                        (tr0[0, 0, 0][n] - trmin[0, 0][n])
                        * pdelx
                        / constants.MAPL_GRAV
                        / dt
                        + trflx[0, 0, 0][n]
                        - trflx[0, 0, 1][n]
                        + trflx_d[0, 0, -1][n]
                        - trflx_d[0, 0, 0][n]
                        - trflx_u[0, 0, 0][n]
                    )
                    trflx_u[0, 0, -1][n] = max(0.0, -dum)
                    n += 1

    with computation(FORWARD), interval(0, None):
        if not condensation:
            if dotransport == 1:
                n = 0
                while n < ncnst:
                    pdelx = dp0
                    trten[0, 0, 0][n] = (
                        (
                            trflx[0, 0, 0][n]
                            - trflx[0, 0, 1][n]
                            + trflx_d[0, 0, 0][n]
                            - trflx_d[0, 0, 0][n]
                            + trflx_u[0, 0, 0][n]
                            - trflx_u[0, 0, 0][n]
                        )
                        * constants.MAPL_GRAV
                        / pdelx
                    )
                    n += 1


def compute_diagnostic_outputs(
    condensation: BoolFieldIJ,
    prel: FloatField,
    thlu: FloatField,
    qtu: FloatField,
    krel: IntField,
    ese: FloatField_Extra_Dim,
    esx: FloatField_Extra_Dim,
    umf_out: FloatField,
    qtflx_out: FloatField,
    slflx_out: FloatField,
    uflx_out: FloatField,
    vflx_out: FloatField,
    qcubelow: FloatFieldIJ,
    qlubelow: FloatFieldIJ,
    qiubelow: FloatFieldIJ,
    rcwp: FloatFieldIJ,
    rlwp: FloatFieldIJ,
    riwp: FloatFieldIJ,
    cush_inout: FloatFieldIJ,
):
    """
    Stencil to compute default diagnostic outputs.

    Note that since 'qtu(krel-1:kpen-1)' & 'thlu(krel-1:kpen-1)' has
    been adjusted after detraining cloud condensate into environment
    during cumulus updraft motion,  below calculations will  exactly
    reproduce in-cloud properties as shown in the output analysis.

    Inputs:
    prel [FloatField]: Lowest level from which buoyancy sorting occurs [Pa]
    thlu [FloatField]: Updraft liquid potential temperature at the interface [K]
    qtu [FloatField]: Updraft total specific humidity at the interface [kg/kg]
    krel [IntField]: Release layer where buoyancy sorting first occurs
    ese [FloatField_Extra_Dim]: Used in QSat_Float [n/a]
    esx [FloatField_Extra_Dim]: Used in QSat_Float [n/a]

    Inouts:
    condensation [BoolFieldIJ]: Mask that indicates if condensation has occurred

    Outputs:
    umf_out [FloatField]: Updraft mass flux at the interfaces [kg/m2/s]
    qtflx_out [FloatField]: Mixing ratio flux [?]
    slflx_out [FloatField]: Sensible heat flux [?]
    uflx_out [FloatField]: Zonal wind flux [?]
    vflx_out [FloatField]:  Meridional wind flux [?]
    qcubelow [FloatFieldIJ]: [?]
    qlubelow [FloatFieldIJ]: [?]
    qiubelow [FloatFieldIJ]: [?]
    rcwp [FloatFieldIJ]: Layer mean Cumulus LWP+IWP [kg/m2]
    rlwp [FloatFieldIJ]: Layer mean Cumulus LWP [kg/m2]
    riwp [FloatFieldIJ]: Layer mean Cumulus IWP [kg/m2]
    cush_inout [FloatFieldIJ]: Convective scale height [m]
    """
    with computation(FORWARD), interval(...):
        if not condensation:
            # Compute default diagnostic outputs
            # Note that since 'qtu(krel-1:kpen-1)' & 'thlu(krel-1:kpen-1)' has
            # been adjusted after detraining cloud condensate into environment
            # during cumulus updraft motion,  below calculations will  exactly
            # reproduce in-cloud properties as shown in the output analysis.

            thj, qvj, qlj, qij, qse, id_check = conden(
                prel, thlu.at(K=krel - 1), qtu.at(K=krel - 1), ese, esx
            )

            if id_check == 1:
                condensation = True
                umf_out[0, 0, 1] = 0.0
                dcm_out = 0.0
                qvten_out = 0.0
                qlten_out = 0.0
                qiten_out = 0.0
                sten_out = 0.0
                uten_out = 0.0
                vten_out = 0.0
                qrten_out = 0.0
                qsten_out = 0.0
                cufrc_out = 0.0
                cush_inout = -1.0
                qldet_out = 0.0
                qidet_out = 0.0
                qtflx_out[0, 0, 1] = 0.0
                slflx_out[0, 0, 1] = 0.0
                uflx_out[0, 0, 1] = 0.0
                vflx_out[0, 0, 1] = 0.0
                fer_out = constants.MAPL_UNDEF
                fdr_out = constants.MAPL_UNDEF

            if not condensation:
                qcubelow = qlj + qij
                qlubelow = qlj
                qiubelow = qij
                rcwp = 0.0
                rlwp = 0.0
                riwp = 0.0


def calc_cumulus_condensate_at_interface(
    condensation: BoolFieldIJ,
    krel: IntField,
    kpen: IntField,
    pifc0: FloatField,
    ppen: FloatFieldIJ,
    thlu_top: FloatFieldIJ,
    qtu_top: FloatFieldIJ,
    thlu: FloatField,
    qtu: FloatField,
    ese: FloatField_Extra_Dim,
    esx: FloatField_Extra_Dim,
    umf_out: FloatField,
    qtflx_out: FloatField,
    slflx_out: FloatField,
    uflx_out: FloatField,
    vflx_out: FloatField,
    ufrc: FloatField,
    ufrclcl: FloatField,
    prel: FloatField,
    criqc: Float,
    qcu: FloatField,
    qlu: FloatField,
    qiu: FloatField,
    qcubelow: FloatFieldIJ,
    qlubelow: FloatFieldIJ,
    qiubelow: FloatFieldIJ,
    rcwp: FloatFieldIJ,
    rlwp: FloatFieldIJ,
    riwp: FloatFieldIJ,
    cufrc: FloatField,
    cush_inout: FloatFieldIJ,
):
    """
    Stencil to calculate cumulus condensate at the upper interface of each layer.

    Inputs:
    krel [IntField]: Release layer where buoyancy sorting first occurs
    kpen [IntField]: Highest layer with positive updraft velocity
    pifc0 [FloatField]: Environmental pressure at the interfaces [Pa]
    ppen [FloatFieldIJ]: Highest interface evel where Cu w = 0 [Pa]
    thlu_top [FloatFieldIJ]: Updraft liquid potential temperature at top [K]
    qtu_top [FloatFieldIJ]: Updraft total specific humidity at top [kg/kg]
    thlu [FloatField]: Updraft liquid potential temperature at the interface [K]
    qtu [FloatField]: Updraft total specific humidity at the interface [kg/kg]
    ese [FloatField_Extra_Dim]: Used in QSat_Float [n/a]
    esx [FloatField_Extra_Dim]:  Used in QSat_Float [n/a]
    ufrc [FloatField]: Cumulus updraft fraction [fraction]
    ufrclcl [FloatField]: Cumulus updraft fraction at the LCL [fraction]
    prel [FloatField]: Lowest level from which buoyancy sorting occurs [Pa]
    criqc [Float]: Maximum condensate that can be hold by cumulus updraft [kg/kg]
    qcubelow [FloatFieldIJ]: [?]
    qlubelow [FloatFieldIJ]: [?]
    qiubelow [FloatFieldIJ]: [?]

    Inouts:
    condensation [BoolFieldIJ]: Mask that indicates if condensation has occurred

    Outputs:
    umf_out [FloatField]: Updraft mass flux at the interfaces [kg/m2/s]
    qtflx_out [FloatField]: Mixing ratio flux [?]
    slflx_out [FloatField]: Sensible heat flux [?]
    uflx_out [FloatField]: Zonal wind flux [?]
    vflx_out [FloatField]: Meridional wind flux [?]
    qcu [FloatField]: Condensate water specific humidity within cumulus updraft at the layer mid-point [kg/kg]
    qlu [FloatField]: Liquid water specific humidity within cumulus updraft at the layer mid-point [kg/kg]
    qiu [FloatField]: Ice specific humidity within cumulus updraft at the layer mid-point [kg/kg]
    rcwp [FloatFieldIJ]: Layer mean Cumulus LWP+IWP [kg/m2]
    rlwp [FloatFieldIJ]: Layer mean Cumulus LWP [kg/m2]
    riwp [FloatFieldIJ]: Layer mean Cumulus IWP [kg/m2]
    cufrc [FloatField]: Shallow cumulus cloud fraction at the layer mid-point [fraction]
    cush_inout [FloatFieldIJ]: Convective scale height [m]
    """
    with computation(FORWARD), interval(...):
        if not condensation:
            if THIS_K >= krel and THIS_K <= kpen:
                if not condensation:
                    # Calculate cumulus condensate at the upper interface of each layer.
                    # Note 'ppen < 0' and at 'k=kpen' layer, I used 'thlu_top'&'qtu_top'
                    # which explicitly considered zero or non-zero 'fer(kpen)'.

                    if THIS_K == kpen:
                        thj, qvj, qlj, qij, qse, id_check = conden(
                            pifc0 + ppen, thlu_top, qtu_top, ese, esx
                        )
                    else:
                        thj, qvj, qlj, qij, qse, id_check = conden(
                            pifc0[0, 0, 1], thlu, qtu, ese, esx
                        )

                    if id_check == 1:
                        condensation = True
                        umf_out[0, 0, 1] = 0.0
                        dcm_out = 0.0
                        qvten_out = 0.0
                        qlten_out = 0.0
                        qiten_out = 0.0
                        sten_out = 0.0
                        uten_out = 0.0
                        vten_out = 0.0
                        qrten_out = 0.0
                        qsten_out = 0.0
                        cufrc_out = 0.0
                        cush_inout = -1.0
                        qldet_out = 0.0
                        qidet_out = 0.0
                        qtflx_out[0, 0, 1] = 0.0
                        slflx_out[0, 0, 1] = 0.0
                        uflx_out[0, 0, 1] = 0.0
                        vflx_out[0, 0, 1] = 0.0
                        fer_out = constants.MAPL_UNDEF
                        fdr_out = constants.MAPL_UNDEF

                    if not condensation:
                        # Calculate in-cloud mean LWC ( qlu(k) ), IWC ( qiu(k) ),
                        # & layer mean cumulus fraction ( cufrc(k) ),
                        # vertically-integrated layer mean LWP and IWP. Expel some
                        # of in-cloud condensate at the upper interface if it is
                        # largr than criqc. Note cumulus cloud fraction is assumed
                        # to be twice of core updraft fractional area. Thus LWP
                        # and IWP will be twice of actual value coming from our
                        # scheme.

                        qcu = 0.5 * (qcubelow + qlj + qij)
                        qlu = 0.5 * (qlubelow + qlj)
                        qiu = 0.5 * (qiubelow + qij)
                        cufrc = ufrc + ufrc[0, 0, 1]

                        if THIS_K == krel:
                            cufrc = (
                                (ufrclcl + ufrc[0, 0, 1])
                                * (prel - pifc0[0, 0, 1])
                                / (pifc0 - pifc0[0, 0, 1])
                            )
                        elif THIS_K == kpen:
                            cufrc = (ufrc + 0.0) * (-ppen) / (pifc0 - pifc0[0, 0, 1])
                            if (qlj + qij) > criqc:
                                qcu = 0.5 * (qcubelow + criqc)
                                qlu = 0.5 * (qlubelow + criqc * qlj / (qlj + qij))
                                qiu = 0.5 * (qiubelow + criqc * qij / (qlj + qij))

                        rcwp = (
                            rcwp
                            + (qlu + qiu)
                            * (pifc0 - pifc0[0, 0, 1])
                            / constants.MAPL_GRAV
                            * cufrc
                        )

                        rlwp = (
                            rlwp
                            + qlu
                            * (pifc0 - pifc0[0, 0, 1])
                            / constants.MAPL_GRAV
                            * cufrc
                        )

                        riwp = (
                            riwp
                            + qiu
                            * (pifc0 - pifc0[0, 0, 1])
                            / constants.MAPL_GRAV
                            * cufrc
                        )

                        qcubelow = qlj + qij
                        qlubelow = qlj
                        qiubelow = qij

            if not condensation:
                # Cloud top and base interface indices
                cnt = f32(kpen)
                cnb = f32(krel - 1)

                # End of formal calculation. Below blocks are for implicit CIN
                # calculations with re-initialization and save variables at
                # iter_cin = 1.


def adjust_implicit_CIN_inputs(
    condensation: BoolFieldIJ,
    qv0: FloatField,
    qvten: FloatField,
    dt: Float,
    ql0: FloatField,
    qlten: FloatField,
    qi0: FloatField,
    qiten: FloatField,
    s0: FloatField,
    sten: FloatField,
    u0: FloatField,
    uten: FloatField,
    v0: FloatField,
    vten: FloatField,
    t0: FloatField,
    dotransport: Int,
    tr0_s: FloatField_NTracers,
    tr0: FloatField_NTracers,
    trten: FloatField_NTracers,
    umf_s: FloatField,
    umf_zint: FloatField,
    dcm: FloatField,
    qrten: FloatField,
    qsten: FloatField,
    cush: FloatFieldIJ,
    cufrc: FloatField,
    slflx_s: FloatField,
    slflx: FloatField,
    qtflx_s: FloatField,
    qtflx: FloatField,
    uflx_s: FloatField,
    uflx: FloatField,
    vflx_s: FloatField,
    vflx: FloatField,
    qcu: FloatField,
    qlu: FloatField,
    qiu: FloatField,
    fer: FloatField,
    fdr: FloatField,
    xco: FloatField,
    cin_IJ: FloatFieldIJ,
    cinlcl_IJ: FloatFieldIJ,
    cbmf: FloatField,
    qc: FloatField,
    qlten_det: FloatField,
    qiten_det: FloatField,
    qlten_sink: FloatField,
    qiten_sink: FloatField,
    ufrc_s: FloatField,
    ufrc: FloatField,
    qv0_s: FloatField,
    ql0_s: FloatField,
    qi0_s: FloatField,
    s0_s: FloatField,
    t0_s: FloatField,
    dcm_s: FloatField,
    qvten_s: FloatField,
    qlten_s: FloatField,
    qiten_s: FloatField,
    sten_s: FloatField,
    uten_s: FloatField,
    vten_s: FloatField,
    qrten_s: FloatField,
    qsten_s: FloatField,
    qldet_s: FloatField,
    qidet_s: FloatField,
    qlsub_s: FloatField,
    qisub_s: FloatField,
    cush_s: FloatField,
    cufrc_s: FloatField,
    fer_s: FloatField,
    fdr_s: FloatField,
):
    """
    Stencil to adjust the original input profiles for implicit CIN
    calculation. Save the output from "iter_cin = 0".

    These output will be writed-out if "iter_cin = 0" was not performed
    for some reason.

    Variables ending in '_s' indicate output from first iteration
    of CIN calculation.

    Inputs:
    condensation [BoolFieldIJ]: Mask that indicates if condensation has occurred
    qv0 [FloatField]: Environmental specific humidity
    qvten [FloatField]: Tendency of water vapor specific humidity [kg/kg/s]
    dt [Float]: Timestep [s]
    ql0 [FloatField]: Environmental liquid water specific humidity
    qlten [FloatField]: Tendency of liquid water specific humidity [kg/kg/s]
    qi0 [FloatField]: Environmental ice specific humidity
    qiten [FloatField]: Tendency of ice specific humidity [kg/kg/s]
    s0 [FloatField]: Environmental dry static energy [J/kg]
    sten [FloatField]: Tendency of dry static energy [J/kg/s]
    u0 [FloatField]: Environmental zonal wind [m/s]
    uten [FloatField]: Tendency of zonal wind [m/s2]
    v0 [FloatField]: Environmental meridional wind [m/s]
    vten [FloatField]: Tendency of meridional wind [m/s2]
    t0 [FloatField]: Environmental temperature [K]
    dotransport [Int]: Transport tracers [1 true]
    tr0 [FloatField_NTracers]: Environmental tracers [#, kg/kg]
    trten [FloatField_NTracers]: Tendency of tracers [#/s, kg/kg/s]
    umf_zint [FloatField]: Updraft mass flux at the interfaces [kg/m2/s]
    dcm [FloatField]: Detrained cloudy air mass
    qrten [FloatField]: Tendency of rain water specific humidity [kg/kg/s]
    qsten [FloatField]: Tendency of snow specific humidity [kg/kg/s]
    cush [FloatFieldIJ]: Convective scale height [m]
    cufrc [FloatField]: Shallow cumulus cloud fraction at the layer mid-point [fraction]
    slflx [FloatField]: Sensible heat flux [?]
    qtflx [FloatField]: Mixing ratio flux [?]
    uflx [FloatField]: Zonal wind flux [?]
    vflx [FloatField]: Meridional wind flux [?]
    qcu [FloatField]: Condensate water specific humidity within cumulus updraft at the layer mid-point [kg/kg]
    qlu [FloatField]: Liquid water specific humidity within cumulus updraft at the layer mid-point [kg/kg]
    qiu [FloatField]: Ice specific humidity within cumulus updraft at the layer mid-point [kg/kg]
    fer [FloatField]: Fractional lateral entrainment rate [1/Pa]
    fdr [FloatField]: Fractional lateral detrainment rate [1/Pa]
    xco [FloatField]: [?]
    cin_IJ [FloatFieldIJ]: Convective inhibition [J/kg]
    cinlcl_IJ [FloatFieldIJ]: Convective inhibition at LCL [J/kg]
    cbmf [FloatField]: Cloud base mass flux [kg/m2/s]
    qc [FloatField]: Tendency due [?]
    qlten_det [FloatField]: [?]
    qiten_det [FloatField]: [?]
    qlten_sink [FloatField]: Liquid condensate tendency by compensating subsidence/upwelling [kg/kg/s]
    qiten_sink [FloatField]: Ice condensate tendency by compensating subsidence/upwelling [kg/kg/s]
    ufrc [FloatField]: Cumulus updraft fraction [fraction]

    Outputs:
    tr0_s [FloatField_NTracers]: Environmental tracers [#, kg/kg]
    umf_s [FloatField]: Updraft mass flux at the interfaces [kg/m2/s]
    slflx_s [FloatField]: [?]
    qtflx_s [FloatField]: [?]
    uflx_s [FloatField]: [?]
    vflx_s [FloatField]: [?]
    ufrc_s [FloatField]: Cumulus updraft fraction [fraction]
    qv0_s [FloatField]: Environmental specific humidity
    ql0_s [FloatField]: Environmental liquid water specific humidity
    qi0_s [FloatField]: Environmental ice specific humidity
    s0_s [FloatField]: Environmental dry static energy [J/kg]
    t0_s [FloatField]: Environmental temperature [K]
    dcm_s [FloatField]: Detrained cloudy air mass
    qvten_s [FloatField]: Tendency of water vapor specific humidity [kg/kg/s]
    qlten_s [FloatField]: Tendency of liquid water specific humidity [kg/kg/s]
    qiten_s [FloatField]: Tendency of ice specific humidity [kg/kg/s]
    sten_s [FloatField]: Tendency of dry static energy [J/kg/s]
    uten_s [FloatField]: Tendency of zonal wind [m/s2]
    vten_s [FloatField]: Tendency of meridional wind [m/s2]
    qrten_s [FloatField]: Tendency of rain water specific humidity [kg/kg/s]
    qsten_s [FloatField]: Tendency of snow specific humidity [kg/kg/s]
    qldet_s [FloatField]: [?]
    qidet_s [FloatField]: [?]
    qlsub_s [FloatField]: [?]
    qisub_s [FloatField]: [?]
    cush_s [FloatField]: Convective scale height [m]
    cufrc_s [FloatField]: Shallow cumulus cloud fraction at the layer mid-point [fraction]
    fer_s [FloatField]: Fractional lateral entrainment rate [1/Pa]
    fdr_s [FloatField]: Fractional lateral detrainment rate [1/Pa]
    """
    from __externals__ import ncnst

    with computation(FORWARD), interval(...):
        if not condensation:
            # Adjust the original input profiles for implicit CIN calculation
            # Save the output from "iter_cin = 1"
            # These output will be writed-out if "iter_cin = 1" was not performed
            # for some reasons.
            qv0_s = qv0 + qvten * dt
            ql0_s = ql0 + qlten * dt
            qi0_s = qi0 + qiten * dt
            s0_s = s0 + sten * dt
            u0_s = u0 + uten * dt
            v0_s = v0 + vten * dt

            t0_s = t0 + sten * dt / constants.MAPL_CP

            if dotransport == 1:
                n = 0
                while n < ncnst:
                    tr0_s[0, 0, 0][n] = tr0[0, 0, 0][n] + trten[0, 0, 0][n] * dt
                    n += 1

            umf_s[0, 0, 1] = umf_zint[0, 0, 1]
            dcm_s = dcm
            qvten_s = qvten
            qlten_s = qlten
            qiten_s = qiten
            sten_s = sten
            uten_s = uten
            vten_s = vten
            qrten_s = qrten
            qsten_s = qsten
            cush_s = cush
            cufrc_s = cufrc
            slflx_s[0, 0, 1] = slflx[0, 0, 1]
            qtflx_s[0, 0, 1] = qtflx[0, 0, 1]
            uflx_s[0, 0, 1] = uflx[0, 0, 1]
            vflx_s[0, 0, 1] = vflx[0, 0, 1]
            qcu_s = qcu
            qlu_s = qlu
            qiu_s = qiu
            fer_s = fer
            fdr_s = fdr
            xc_s = xco
            cin_s = cin_IJ
            cinlcl_s = cinlcl_IJ
            cbmf_s = cbmf
            qc_s = qc
            qldet_s = qlten_det
            qidet_s = qiten_det
            qlsub_s = qlten_sink
            qisub_s = qiten_sink

            ufrc_s[0, 0, 1] = ufrc[0, 0, 1]


def recalc_environmental_variables(
    condensation: BoolFieldIJ,
    qv0_s: FloatField,
    ql0_s: FloatField,
    qi0_s: FloatField,
    s0_s: FloatField,
    t0_s: FloatField,
    exnmid0: FloatField,
    pmid0: FloatField,
    dotransport: Int,
    sstr0: FloatField_NTracers,
    tr0: FloatField_NTracers,
    u0: FloatField,
    v0: FloatField,
    pifc0: FloatField,
    ese: FloatField_Extra_Dim,
    esx: FloatField_Extra_Dim,
    umf_out: FloatField,
    qtflx_out: FloatField,
    slflx_out: FloatField,
    uflx_out: FloatField,
    vflx_out: FloatField,
    thvl0bot: FloatField,
    thv0bot: FloatField,
    thvl0top: FloatField,
    thv0top: FloatField,
    thl0: FloatField,
    qt0: FloatField,
    thvl0: FloatField,
    ssthl0: FloatField,
    ssu0: FloatField,
    ssv0: FloatField,
    ssqt0: FloatField,
    qv0: FloatField,
    ql0: FloatField,
    qi0: FloatField,
    s0: FloatField,
    t0: FloatField,
    tr0_temp: FloatField,
    cush_inout: FloatFieldIJ,
):
    """
    Stencil to update environmental variables after first iteration
    of implicit CIN calculation.

    The suffix '_s' (e.g., qv0_s) indicates variables used in the first
    iteration of the implicit CIN calculation.

    Inputs:
    qv0_s [FloatField]: Environmental water vapor specific humidity [kg/kg]
    ql0_s [FloatField]: Environmental liquid water specific humidity [kg/kg]
    qi0_s [FloatField]: Environmental ice specific humidity [kg/kg]
    s0_s [FloatField]: Environmental dry static energy [J/kg]
    t0_s [FloatField]: Environmental temperature [K]
    exnmid0 [FloatField]: Exner function at the layer mid-point
    pmid0 [FloatField]: Environmental pressure at the layer mid-point [Pa]
    dotransport [Int]: Transport tracers [1 true]
    tr0 [FloatField_NTracers]: Environmental tracers [#, kg/kg]
    u0 [FloatField]: Environmental zonal wind [m/s]
    v0 [FloatField]: Environmental meridional wind [m/s]
    pifc0 [FloatField]: Environmental pressure at the interfaces [Pa]
    ese [FloatField_Extra_Dim]: Used in QSat_Float [n/a]
    esx [FloatField_Extra_Dim]: Used in QSat_Float [n/a]
    tr0_temp [FloatField]: Environmental tracers [#, kg/kg]

    Inouts:
    condensation [BoolFieldIJ]: Mask that indicates if condensation has occurred
    cush_inout [FloatFieldIJ]: Convective scale height [m]

    Outputs:
    sstr0 [FloatField_NTracers]: Convective tracer [?]
    umf_out [FloatField]: Updraft mass flux at the interfaces [kg/m2/s]
    qtflx_out [FloatField]: Mixing ratio flux [?]
    slflx_out [FloatField]: Sensible heat flux [?]
    uflx_out [FloatField]: Zonal wind flux [?]
    vflx_out [FloatField]: Meridional wind flux [?]
    thvl0bot [FloatField]: Temperature at bottom [?]
    thv0bot [FloatField]: Temperature at bottom [?]
    thvl0top [FloatField]: Temperature at top [?]
    thv0top [FloatField]: Potential temperature at top [?]
    thl0 [FloatField]: Temperature [?]
    qt0 [FloatField]: Mixing ratio [?]
    thvl0 [FloatField]: Temperature [?]
    ssthl0 [FloatField]: [?]
    ssu0 [FloatField]: [?]
    ssv0 [FloatField]: [?]
    ssqt0 [FloatField]: [?]
    qv0 [FloatField]: Environmental specific humidity
    ql0 [FloatField]: Environmental liquid water specific humidity
    qi0 [FloatField]: Environmental ice specific humidity
    s0 [FloatField]: Environmental dry static energy [J/kg]
    t0 [FloatField]: Environmental temperature [K]
    """
    from __externals__ import k_end, ncnst

    with computation(FORWARD), interval(...):
        if not condensation:
            # Recalculate environmental variables for new cin calculation at
            # "iter_cin = 2" using the updated state variables. Perform only
            # for variables necessary for the new cin calculation.

            qv0 = qv0_s
            ql0 = ql0_s
            qi0 = qi0_s
            s0 = s0_s
            t0 = t0_s

            qt0 = qv0 + ql0 + qi0
            thl0 = (
                t0
                - constants.MAPL_LATENT_HEAT_VAPORIZATION * ql0 / constants.MAPL_CP
                - constants.MAPL_LATENT_HEAT_SUBLIMATION * qi0 / constants.MAPL_CP
            ) / exnmid0
            thvl0 = (1.0 + zvir * qt0) * thl0

    with computation(FORWARD), interval(0, 1):
        if not condensation:
            # Compute slopes of environmental variables at bottom layer
            ssthl0 = slope_bot(thl0, pmid0)
            ssqt0 = slope_bot(qt0, pmid0)
            ssu0 = slope_bot(u0, pmid0)
            ssv0 = slope_bot(v0, pmid0)

            if dotransport == 1:
                n = 0
                while n < ncnst:
                    tr0_temp = tr0[0, 0, 0][n]
                    sstr0[0, 0, 0][n] = slope_bot(tr0_temp, pmid0)
                    n += 1

    with computation(FORWARD), interval(1, -1):
        if not condensation:
            # Compute slopes of environmental variables at mid layers
            ssthl0 = slope_mid(k_end, thl0, pmid0)
            ssqt0 = slope_mid(k_end, qt0, pmid0)
            ssu0 = slope_mid(k_end, u0, pmid0)
            ssv0 = slope_mid(k_end, v0, pmid0)

            if dotransport == 1:
                n = 0
                while n < ncnst:
                    tr0_temp = tr0[0, 0, 0][n]
                    sstr0[0, 0, 0][n] = slope_mid(
                        k_end,
                        tr0_temp,
                        pmid0,
                    )
                    n += 1

    with computation(FORWARD), interval(-1, None):
        if not condensation:
            # Compute slopes of environmental variables at top layers
            ssthl0 = slope_top(k_end, thl0, pmid0)
            ssqt0 = slope_top(k_end, qt0, pmid0)
            ssu0 = slope_top(k_end, u0, pmid0)
            ssv0 = slope_top(k_end, v0, pmid0)

            if dotransport == 1:
                n = 0
                while n < ncnst:
                    tr0_temp = tr0[0, 0, 0][n]
                    sstr0[0, 0, 0][n] = slope_top(k_end, tr0_temp, pmid0)
                    n += 1

    with computation(FORWARD), interval(...):
        if not condensation:
            thl0bot = thl0 + ssthl0 * (pifc0 - pmid0)
            qt0bot = qt0 + ssqt0 * (pifc0 - pmid0)
            thj, qvj, qlj, qij, qse, id_check = conden(pifc0, thl0bot, qt0bot, ese, esx)
            if id_check == 1:
                condensation = True
                umf_out[0, 0, 1] = 0.0
                dcm_out = 0.0
                qvten_out = 0.0
                qlten_out = 0.0
                qiten_out = 0.0
                sten_out = 0.0
                uten_out = 0.0
                vten_out = 0.0
                qrten_out = 0.0
                qsten_out = 0.0
                cufrc_out = 0.0
                cush_inout = -1.0
                qldet_out = 0.0
                qidet_out = 0.0
                qtflx_out[0, 0, 1] = 0.0
                slflx_out[0, 0, 1] = 0.0
                uflx_out[0, 0, 1] = 0.0
                vflx_out[0, 0, 1] = 0.0
                fer_out = constants.MAPL_UNDEF
                fdr_out = constants.MAPL_UNDEF

            if not condensation:
                thv0bot = thj * (1.0 + zvir * qvj - qlj - qij)
                thvl0bot = thl0bot * (1.0 + zvir * qt0bot)

                thl0top = thl0 + ssthl0 * (pifc0[0, 0, 1] - pmid0)
                qt0top = qt0 + ssqt0 * (pifc0[0, 0, 1] - pmid0)
                thj, qvj, qlj, qij, qse, id_check = conden(
                    pifc0[0, 0, 1], thl0top, qt0top, ese, esx
                )

                if id_check == 1:
                    condensation = True
                    umf_out[0, 0, 1] = 0.0
                    dcm_out = 0.0
                    qvten_out = 0.0
                    qlten_out = 0.0
                    qiten_out = 0.0
                    sten_out = 0.0
                    uten_out = 0.0
                    vten_out = 0.0
                    qrten_out = 0.0
                    qsten_out = 0.0
                    cufrc_out = 0.0
                    cush_inout = -1.0
                    qldet_out = 0.0
                    qidet_out = 0.0
                    qtflx_out[0, 0, 1] = 0.0
                    slflx_out[0, 0, 1] = 0.0
                    uflx_out[0, 0, 1] = 0.0
                    vflx_out[0, 0, 1] = 0.0
                    fer_out = constants.MAPL_UNDEF
                    fdr_out = constants.MAPL_UNDEF

                if not condensation:
                    thv0top = thj * (1.0 + zvir * qvj - qlj - qij)
                    thvl0top = thl0top * (1.0 + zvir * qt0top)

                # End of iter loop


def update_output_variables(
    umf_zint: FloatField,
    zifc0: FloatField,
    kinv: IntField,
    dcm: FloatField,
    qvten: FloatField,
    qlten: FloatField,
    qiten: FloatField,
    sten: FloatField,
    uten: FloatField,
    vten: FloatField,
    qrten: FloatField,
    qsten: FloatField,
    cufrc: FloatField,
    cush: FloatFieldIJ,
    qlten_det: FloatField,
    qiten_det: FloatField,
    qlten_sink: FloatField,
    qiten_sink: FloatField,
    rdrop: Float,
    qtflx: FloatField,
    slflx: FloatField,
    uflx: FloatField,
    vflx: FloatField,
    dotransport: Int,
    trten: FloatField_NTracers,
    dt: Float,
    fer: FloatField,
    fdr: FloatField,
    kpen: IntField,
    # # Outputs
    umf_out: FloatField,
    dcm_out: FloatField,
    qvten_out: FloatField,
    qlten_out: FloatField,
    qiten_out: FloatField,
    sten_out: FloatField,
    uten_out: FloatField,
    vten_out: FloatField,
    qrten_out: FloatField,
    qsten_out: FloatField,
    cufrc_out: FloatField,
    cush_inout: FloatFieldIJ,
    qldet_out: FloatField,
    qidet_out: FloatField,
    qlsub_out: FloatField,
    qisub_out: FloatField,
    ndrop_out: FloatField,
    nice_out: FloatField,
    qtflx_out: FloatField,
    slflx_out: FloatField,
    uflx_out: FloatField,
    vflx_out: FloatField,
    tr0_inout: FloatField_NTracers,
    fer_out: FloatField,
    fdr_out: FloatField,
    umf_outvar: FloatField,
    dcm_outvar: FloatField,
    qvten_outvar: FloatField,
    qlten_outvar: FloatField,
    qiten_outvar: FloatField,
    sten_outvar: FloatField,
    uten_outvar: FloatField,
    vten_outvar: FloatField,
    qrten_outvar: FloatField,
    qsten_outvar: FloatField,
    cufrc_outvar: FloatField,
    cush_inoutvar: FloatFieldIJ,
    qldet_outvar: FloatField,
    qidet_outvar: FloatField,
    qlsub_outvar: FloatField,
    qisub_outvar: FloatField,
    qtflx_outvar: FloatField,
    slflx_outvar: FloatField,
    uflx_outvar: FloatField,
    vflx_outvar: FloatField,
    fer_outvar: FloatField,
    fdr_outvar: FloatField,
    tr0_inoutvar: FloatField_NTracers,
    del_CIN: FloatFieldIJ,
):
    """
    Stencil to update ComputeUwshcu output variables.

    If del_CIN <= 0.0, output variables remain unchanged from the first
    implicit CIN iteration.

    If del_CIN > 0.0, output variables are updated.

    Inputs:
    del_CIN [FloatFieldIJ]: Difference between initial and final CIN calculations [J/kg]
    umf_zint [FloatField]: Updraft mass flux at the interfaces [kg/m2/s]
    zifc0 [FloatField]: Environmental height at interfaces [m]
    kinv [IntField]: Inversion layer with PBL top interface as lower interface
    dcm [FloatField]: Detrained cloudy air mass
    qvten [FloatField]: Tendency of water vapor specific humidity [kg/kg/s]
    qlten [FloatField]: Tendency of liquid water specific humidity [kg/kg/s]
    qiten [FloatField]: Tendency of ice specific humidity [kg/kg/s]
    sten [FloatField]: Tendency of dry static energy [J/kg/s]
    uten [FloatField]: Tendency of zonal wind [m/s2]
    vten [FloatField]: Tendency of meridional wind [m/s2]
    qrten [FloatField]: Tendency of rain water specific humidity [kg/kg/s]
    qsten [FloatField]: Tendency of snow specific humidity [kg/kg/s]
    cufrc [FloatField]: Shallow cumulus cloud fraction at the layer mid-point [fraction]
    cush [FloatFieldIJ]: Convective scale height [m]
    qlten_det [FloatField]: [?]
    qiten_det [FloatField]: [?]
    qlten_sink [FloatField]: Liquid condensate tendency by compensating subsidence/upwelling [kg/kg/s]
    qiten_sink [FloatField]: Ice condensate tendency by compensating subsidence/upwelling [kg/kg/s]
    rdrop [Float]: Liquid drop radius
    qtflx [FloatField]: Mixing ratio flux [?]
    slflx [FloatField]: Sensible heat flux [?]
    uflx [FloatField]: Zonal wind flux [?]
    vflx [FloatField]: Meridional wind flux [?]
    dotransport [Int]: Transport tracers [1 true]
    trten [FloatField_NTracers]: Tendency of tracers [#/s, kg/kg/s]
    dt [Float]: Timestep [s]
    fer [FloatField]: Fractional lateral entrainment rate [1/Pa]
    fdr [FloatField]: Fractional lateral detrainment rate [1/Pa]
    kpen [IntField]: Highest layer with positive updraft velocity
    umf_out [FloatField]: Updraft mass flux at the interfaces [kg/m2/s]
    dcm_out [FloatField]: Detrained cloudy air mass
    qvten_out [FloatField]: Tendency of water vapor specific humidity [kg/kg/s]
    qlten_out [FloatField]: Tendency of liquid water specific humidity [kg/kg/s]
    qiten_out [FloatField]: Tendency of ice specific humidity [kg/kg/s]
    sten_out [FloatField]: Tendency of dry static energy [J/kg/s]
    uten_out [FloatField]: Tendency of zonal wind [m/s2]
    vten_out [FloatField]: Tendency of meridional wind [m/s2]
    qrten_out [FloatField]: Tendency of rain water specific humidity [kg/kg/s]
    qsten_out [FloatField]: Tendency of snow specific humidity [kg/kg/s]
    cufrc_out [FloatField]: Shallow cumulus cloud fraction at the layer mid-point [fraction]
    cush_inout [FloatFieldIJ]: Convective scale height [m]
    qldet_out [FloatField]: [?]
    qidet_out [FloatField]: [?]
    qlsub_out [FloatField]: [?]
    qisub_out [FloatField]: [?]
    qtflx_out [FloatField]: Mixing ratio flux [?]
    slflx_out [FloatField]: Sensible heat flux [?]
    uflx_out [FloatField]: Zonal wind flux [?]
    vflx_out [FloatField]: Meridional wind flux [?]
    tr0_inout [FloatField_NTracers]: Environmental tracers [#, kg/kg]
    fer_out [FloatField]: Fractional lateral entrainment rate [1/Pa]
    fdr_out [FloatField]: Fractional lateral detrainment rate [1/Pa]

    Outputs:
    umf_outvar [FloatField]: Updraft mass flux at the interfaces [kg/m2/s]
    dcm_outvar [FloatField]: Detrained cloudy air mass
    qvten_outvar [FloatField]: Tendency of water vapor specific humidity [kg/kg/s]
    qlten_outvar [FloatField]: Tendency of liquid water specific humidity [kg/kg/s]
    qiten_outvar [FloatField]: Tendency of ice specific humidity [kg/kg/s]
    sten_outvar [FloatField]: Tendency of dry static energy [J/kg/s]
    uten_outvar [FloatField]: Tendency of zonal wind [m/s2]
    vten_outvar [FloatField]: Tendency of meridional wind [m/s2]
    qrten_outvar [FloatField]: Tendency of rain water specific humidity [kg/kg/s]
    qsten_outvar [FloatField]: Tendency of snow specific humidity [kg/kg/s]
    cufrc_outvar [FloatField]: Shallow cumulus cloud fraction at the layer mid-point [fraction]
    cush_inoutvar [FloatFieldIJ]: Convective scale height [m]
    qldet_outvar [FloatField]: [?]
    qidet_outvar [FloatField]: [?]
    qlsub_outvar [FloatField]: [?]
    qisub_outvar [FloatField]: [?]
    qtflx_outvar [FloatField]: [?]
    slflx_outvar [FloatField]: [?]
    uflx_outvar [FloatField]: [?]
    vflx_outvar [FloatField]: [?]
    fer_outvar [FloatField]: Fractional lateral entrainment rate [1/Pa]
    fdr_outvar [FloatField]: Fractional lateral detrainment rate [1/Pa]
    ndrop_out [FloatField]: [?]
    nice_out [FloatField]: [?]
    tr0_inoutvar [FloatField_NTracers]: Environmental tracers [#, kg/kg]
    """
    from __externals__ import ncnst

    with computation(FORWARD), interval(...):

        if del_CIN <= 0.0:
            umf_outvar = umf_out
            dcm_outvar = dcm_out
            qvten_outvar = qvten_out
            qlten_outvar = qlten_out
            qiten_outvar = qiten_out
            sten_outvar = sten_out
            uten_outvar = uten_out
            vten_outvar = vten_out
            qrten_outvar = qrten_out
            qsten_outvar = qsten_out
            cufrc_outvar = cufrc_out
            cush_inoutvar = cush_inout
            qldet_outvar = qldet_out
            qidet_outvar = qidet_out
            qlsub_outvar = qlsub_out
            qisub_outvar = qisub_out
            qtflx_outvar = qtflx_out
            slflx_outvar = slflx_out
            uflx_outvar = uflx_out
            vflx_outvar = vflx_out
            fer_outvar = fer_out
            fdr_outvar = fdr_out
            if dotransport == 1:
                n = 0
                while n < ncnst:
                    tr0_inoutvar[0, 0, 0][n] = tr0_inout[0, 0, 0][n]
                    n += 1

        if del_CIN > 0.0:
            umf_outvar = umf_zint

            if THIS_K <= kinv - 1:
                umf_outvar = umf_zint.at(K=kinv) * zifc0 / zifc0.at(K=kinv)

            dcm_outvar = dcm
            qvten_outvar = qvten
            qlten_outvar = qlten
            qiten_outvar = qiten
            sten_outvar = sten
            uten_outvar = uten
            vten_outvar = vten
            qrten_outvar = qrten
            qsten_outvar = qsten
            cufrc_outvar = cufrc
            cush_inoutvar = cush
            qldet_outvar = qlten_det
            qidet_outvar = qiten_det
            qlsub_outvar = qlten_sink
            qisub_outvar = qiten_sink
            ndrop_out = qlten_det / (4188.787 * rdrop**3)
            nice_out = qiten_det / (3.0e-10)
            qtflx_outvar = qtflx
            slflx_outvar = slflx
            uflx_outvar = uflx
            vflx_outvar = vflx

            if dotransport == 1:
                n = 0
                while n < ncnst:
                    tr0_inoutvar[0, 0, 0][n] = (
                        tr0_inout[0, 0, 0][n] + trten[0, 0, 0][n] * dt
                    )
                    n += 1

            # # Below are specific diagnostic output for detailed
            # # analysis of cumulus scheme
            fer_outvar = constants.MAPL_UNDEF
            fdr_outvar = constants.MAPL_UNDEF
            if THIS_K <= kpen:
                fer_outvar = fer
                fdr_outvar = fdr


def compute_uwshcu_invert_after(
    # Inputs
    k0: Int,
    umf_outvar: FloatField,
    qtflx_outvar: FloatField,
    slflx_outvar: FloatField,
    uflx_outvar: FloatField,
    vflx_outvar: FloatField,
    dcm_outvar: FloatField,
    qvten_outvar: FloatField,
    qlten_outvar: FloatField,
    qiten_outvar: FloatField,
    sten_outvar: FloatField,
    uten_outvar: FloatField,
    vten_outvar: FloatField,
    qrten_outvar: FloatField,
    qsten_outvar: FloatField,
    cufrc_outvar: FloatField,
    qldet_outvar: FloatField,
    qidet_outvar: FloatField,
    qlsub_outvar: FloatField,
    qisub_outvar: FloatField,
    fer_outvar: FloatField,
    fdr_outvar: FloatField,
    ndrop_out: FloatField,
    nice_out: FloatField,
    tr0: FloatField_NTracers,
    tr0_inoutvar: FloatField_NTracers,
    # Outputs
    umf_inv: FloatField,
    dcm_inv: FloatField,
    qtflx_inv: FloatField,
    slflx_inv: FloatField,
    uflx_inv: FloatField,
    vflx_inv: FloatField,
    qvten_inv: FloatField,
    qlten_inv: FloatField,
    qiten_inv: FloatField,
    tten_inv: FloatField,
    uten_inv: FloatField,
    vten_inv: FloatField,
    qrten_inv: FloatField,
    qsten_inv: FloatField,
    cufrc_inv: FloatField,
    fer_inv: FloatField,
    fdr_inv: FloatField,
    ndrop_inv: FloatField,
    nice_inv: FloatField,
    qldet_inv: FloatField,
    qlsub_inv: FloatField,
    qidet_inv: FloatField,
    qisub_inv: FloatField,
    CNV_Tracers: FloatField_NTracers,
    dotransport: Int,
):
    from __externals__ import k_end, ncnst

    with computation(FORWARD), interval(...):
        # Revert interface variables
        k_inv = k0 - THIS_K
        umf_inv = umf_outvar.at(K=k_inv)
        qtflx_inv = qtflx_outvar.at(K=k_inv)
        slflx_inv = slflx_outvar.at(K=k_inv)
        uflx_inv = uflx_outvar.at(K=k_inv)
        vflx_inv = vflx_outvar.at(K=k_inv)

    with computation(FORWARD), interval(...):
        # Revert mid-level variables
        k_inv = k0 - THIS_K - 1
        dcm_inv = dcm_outvar.at(K=k_inv)
        qvten_inv = qvten_outvar.at(K=k_inv)
        qlten_inv = qlten_outvar.at(K=k_inv)
        qiten_inv = qiten_outvar.at(K=k_inv)
        tten_inv = sten_outvar.at(K=k_inv) / constants.MAPL_CP
        uten_inv = uten_outvar.at(K=k_inv)
        vten_inv = vten_outvar.at(K=k_inv)
        qrten_inv = qrten_outvar.at(K=k_inv)
        qsten_inv = qsten_outvar.at(K=k_inv)
        cufrc_inv = cufrc_outvar.at(K=k_inv)
        fer_inv = fer_outvar.at(K=k_inv)
        fdr_inv = fdr_outvar.at(K=k_inv)
        qldet_inv = qldet_outvar.at(K=k_inv)
        qidet_inv = qidet_outvar.at(K=k_inv)
        qlsub_inv = qlsub_outvar.at(K=k_inv)
        qisub_inv = qisub_outvar.at(K=k_inv)
        ndrop_inv = ndrop_out.at(K=k_inv)
        nice_inv = nice_out.at(K=k_inv)

        mintracer = 1.1754944e-38  # mintracer = tiny(1.0)
        if dotransport == 1:
            n = 0
            while n < ncnst:
                tr0[0, 0, 0][n] = max(mintracer, tr0_inoutvar[0, 0, 0][n])
                CNV_Tracers[0, 0, 0][n] = tr0.at(K=k_inv, ddim=[n])
                n += 1

        if THIS_K == k_end:
            dcm_inv = 0.0

        # Re-scale liquid/ice water sub-tendencies to enforce conservation
        if abs(qldet_inv + qlsub_inv) > 1e-12:
            tmp2d = qlten_inv / (qldet_inv + qlsub_inv)
            qldet_inv = tmp2d * qldet_inv
            qlsub_inv = tmp2d * qlsub_inv

        if abs(qidet_inv + qisub_inv) > 1e-12:
            tmp2d = qiten_inv / (qidet_inv + qisub_inv)
            qidet_inv = tmp2d * qidet_inv
            qisub_inv = tmp2d * qisub_inv


class ComputeUwshcuInv:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        UW_config: UWConfiguration,
        formulation: SaturationFormulation = SaturationFormulation.Staars,
    ) -> None:

        # Initialize the ComputeUwshcu class

        # Parameters:
        # stencil_factory (StencilFactory): Factory for creating stencil computations.
        # quantity_factory (QuantityFactory): Factory for creating quantities.
        # UW_config (dataclass): Data class containing configuration dependent
        # constants.
        # formulation: Saturation Formulation used for QSat.

        self.UW_config = UW_config
        self.temporaries = Temporaries.make(quantity_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory
        grid_indexing = stencil_factory.grid_indexing

        if constants.NCNST != self.UW_config.NCNST:
            raise NotImplementedError(
                f"Coding limitation: {constants.NCNST} tracers are expected, "
                f"getting {self.UW_config.NCNST}"
            )

        if self.UW_config.k0 < 5:
            raise NotImplementedError(
                f"Coding limitation: Only {self.UW_config.k0} k-levels are "
                f"available, atleast 5 are expected"
            )

        if self.UW_config.windsrcavg != 0:
            raise NotImplementedError(
                f"Coding limitation: windsrcavg is {self.UW_config.windsrcavg}"
                f"expected 0"
            )

        self._compute_uwshcu_invert_before = self.stencil_factory.from_dims_halo(
            func=compute_uwshcu_invert_before,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ncnst": UW_config.NCNST},
        )

        self._compute_thermodynamic_variables = self.stencil_factory.from_dims_halo(
            func=compute_thermodynamic_variables,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ncnst": UW_config.NCNST},
        )

        self._compute_thv0_thvl0 = self.stencil_factory.from_dims_halo(
            func=compute_thv0_thvl0,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ncnst": UW_config.NCNST},
        )

        self._find_pbl_height = self.stencil_factory.from_dims_halo(
            func=find_pbl_height,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._find_pbl_averages = self.stencil_factory.from_dims_halo(
            func=find_pbl_averages,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._find_cumulus_characteristics = self.stencil_factory.from_dims_halo(
            func=find_cumulus_characteristics,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ncnst": UW_config.NCNST},
        )

        self._find_klcl = self.stencil_factory.from_dims_halo(
            func=find_klcl,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._compute_cin_cinlcl = self.stencil_factory.from_dims_halo(
            func=compute_cin_cinlcl,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ncnst": UW_config.NCNST},
        )

        self._avg_initial_and_final_cin = self.stencil_factory.from_dims_halo(
            func=avg_initial_and_final_cin,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ncnst": UW_config.NCNST},
        )

        self._define_prel_krel = self.stencil_factory.from_dims_halo(
            func=define_prel_krel,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._calc_cumulus_base_mass_flux = self.stencil_factory.from_dims_halo(
            func=calc_cumulus_base_mass_flux,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._define_updraft_properties = self.stencil_factory.from_dims_halo(
            func=define_updraft_properties,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._define_env_properties = self.stencil_factory.from_dims_halo(
            func=define_env_properties,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ncnst": UW_config.NCNST},
        )

        self._buoyancy_sorting = self.stencil_factory.from_dims_halo(
            func=buoyancy_sorting,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ncnst": UW_config.NCNST},
        )

        self._calc_ppen = self.stencil_factory.from_dims_halo(
            func=calc_ppen,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._recalc_condensate = self.stencil_factory.from_dims_halo(
            func=recalc_condensate,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._calc_entrainment_mass_flux = self.stencil_factory.from_dims_halo(
            func=calc_entrainment_mass_flux,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ncnst": UW_config.NCNST},
        )

        self._calc_pbl_fluxes = self.stencil_factory.from_dims_halo(
            func=calc_pbl_fluxes,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ncnst": UW_config.NCNST},
        )

        self._non_buoyancy_sorting_fluxes = self.stencil_factory.from_dims_halo(
            func=non_buoyancy_sorting_fluxes,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ncnst": UW_config.NCNST},
        )

        self._buoyancy_sorting_fluxes = self.stencil_factory.from_dims_halo(
            func=buoyancy_sorting_fluxes,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ncnst": UW_config.NCNST},
        )

        self._penetrative_entrainment_fluxes = self.stencil_factory.from_dims_halo(
            func=penetrative_entrainment_fluxes,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ncnst": UW_config.NCNST},
        )

        self._calc_momentum_tendency = self.stencil_factory.from_dims_halo(
            func=calc_momentum_tendency,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._calc_thermodynamic_tendencies = self.stencil_factory.from_dims_halo(
            func=calc_thermodynamic_tendencies,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._prevent_negative_condensate = self.stencil_factory.from_dims_halo(
            func=prevent_negative_condensate,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._calc_tracer_tendencies = self.stencil_factory.from_dims_halo(
            func=calc_tracer_tendencies,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ncnst": UW_config.NCNST},
        )

        self._compute_diagnostic_outputs = self.stencil_factory.from_dims_halo(
            func=compute_diagnostic_outputs,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._calc_cumulus_condensate_at_interfaces = (
            self.stencil_factory.from_dims_halo(
                func=calc_cumulus_condensate_at_interface,
                compute_dims=[X_DIM, Y_DIM, Z_DIM],
            )
        )

        self._adjust_implicit_CIN_inputs = self.stencil_factory.from_dims_halo(
            func=adjust_implicit_CIN_inputs,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ncnst": UW_config.NCNST},
        )

        self._recalc_environmental_variables = self.stencil_factory.from_dims_halo(
            func=recalc_environmental_variables,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ncnst": UW_config.NCNST},
        )

        self._update_output_variables = self.stencil_factory.from_dims_halo(
            func=update_output_variables,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ncnst": UW_config.NCNST},
        )

        self._compute_uwshcu_invert_after = self.stencil_factory.from_dims_halo(
            func=compute_uwshcu_invert_after,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ncnst": UW_config.NCNST},
        )

        # Create masks
        # Mask that indicates if condensation has occurred (e.g., Fortran goto 333)
        self.condensation = self.quantity_factory.zeros(
            [X_DIM, Y_DIM], "n/a", dtype=bool
        )
        self.stop_cin = self.quantity_factory.zeros([X_DIM, Y_DIM], "n/a", dtype=bool)
        self.stop_buoyancy_sort = self.quantity_factory.zeros(
            [X_DIM, Y_DIM], "n/a", dtype=bool
        )

        self.ntracers_quantity_factory = self.make_ntracers_quantity_factory(
            self.quantity_factory,
        )

        # Create tracer fields
        self.trsrc = self.ntracers_quantity_factory.zeros(
            [X_DIM, Y_DIM, "ntracers"], "n/a"
        )
        self.trsrc_o = self.ntracers_quantity_factory.zeros(
            [X_DIM, Y_DIM, "ntracers"], "n/a"
        )
        self.tre = self.ntracers_quantity_factory.zeros(
            [X_DIM, Y_DIM, "ntracers"], "n/a"
        )
        self.trmin = self.ntracers_quantity_factory.zeros(
            [X_DIM, Y_DIM, "ntracers"], "n/a"
        )
        self.sstr0 = self.ntracers_quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_DIM, "ntracers"], "n/a"
        )
        self.tr0 = self.ntracers_quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_DIM, "ntracers"], "n/a"
        )
        self.tr0_o = self.ntracers_quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_DIM, "ntracers"], "n/a"
        )
        self.sstr0_o = self.ntracers_quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_DIM, "ntracers"], "n/a"
        )
        self.trten = self.ntracers_quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_DIM, "ntracers"], "n/a"
        )
        self.tr0_s = self.ntracers_quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_DIM, "ntracers"], "n/a"
        )
        self.tr0_inoutvar = self.ntracers_quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_DIM, "ntracers"], "n/a"
        )
        self.tr0_inout = self.ntracers_quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_DIM, "ntracers"], "n/a"
        )
        self.trflx = self.ntracers_quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_INTERFACE_DIM, "ntracers"], "n/a"
        )
        self.tru = self.ntracers_quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_INTERFACE_DIM, "ntracers"], "n/a"
        )
        self.tru_emf = self.ntracers_quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_INTERFACE_DIM, "ntracers"], "n/a"
        )
        self.xflx_ndim = self.ntracers_quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_INTERFACE_DIM, "ntracers"], "n/a"
        )
        self.trflx_d = self.ntracers_quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_INTERFACE_DIM, "ntracers"], "n/a"
        )
        self.trflx_u = self.ntracers_quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_INTERFACE_DIM, "ntracers"], "n/a"
        )

    @staticmethod
    def make_ntracers_quantity_factory(
        ijk_quantity_factory: QuantityFactory,
    ):
        ntracers_quantity_factory = copy.deepcopy(ijk_quantity_factory)
        ntracers_quantity_factory.set_extra_dim_lengths(
            **{
                "ntracers": constants.NCNST,
            }
        )
        return ntracers_quantity_factory

    def __call__(
        self,
        # Field inputs
        pifc0_inv: FloatField,
        zifc0_inv: FloatField,
        pmid0_inv: FloatField,
        zmid0_inv: FloatField,
        kpbl_inv: IntFieldIJ,
        exnmid0_inv: FloatField,
        exnifc0_inv: FloatField,
        dp0_inv: FloatField,
        u0_inv: FloatField,
        v0_inv: FloatField,
        qv0_inv: FloatField,
        ql0_inv: FloatField,
        qi0_inv: FloatField,
        t0_inv: FloatField,
        frland: FloatFieldIJ,
        tke_inv: FloatField,
        rkfre: FloatFieldIJ,
        cush: FloatFieldIJ,
        shfx: FloatFieldIJ,
        evap: FloatFieldIJ,
        cnvtr: FloatFieldIJ,
        CNV_Tracers: FloatField_NTracers,
        # Float/Int inputs
        dotransport: Int,
        k0: Int,
        windsrcavg: Int,
        qtsrchgt: Float,
        qtsrc_fac: Float,
        thlsrc_fac: Float,
        frc_rasn: Float,
        rbuoy: Float,
        epsvarw: Float,
        use_CINcin: i32,
        mumin1: Float,
        rmaxfrac: Float,
        PGFc: Float,
        dt: Float,
        niter_xc: Int,
        criqc: Float,
        rle: Float,
        cridist_opt: Int,
        mixscale: Float,
        rdrag: Float,
        rkm: Float,
        use_self_detrain: Int,
        detrhgt: Float,
        use_cumpenent: Int,
        rpen: Float,
        use_momenflx: Int,
        rdrop: Float,
        iter_cin: i32,
        # Outputs
        umf_inv: FloatField,
        dcm_inv: FloatField,
        qtflx_inv: FloatField,
        slflx_inv: FloatField,
        uflx_inv: FloatField,
        vflx_inv: FloatField,
        qvten_inv: FloatField,
        qlten_inv: FloatField,
        qiten_inv: FloatField,
        tten_inv: FloatField,
        uten_inv: FloatField,
        vten_inv: FloatField,
        qrten_inv: FloatField,
        qsten_inv: FloatField,
        cufrc_inv: FloatField,
        fer_inv: FloatField,
        fdr_inv: FloatField,
        ndrop_inv: FloatField,
        nice_inv: FloatField,
        qldet_inv: FloatField,
        qlsub_inv: FloatField,
        qidet_inv: FloatField,
        qisub_inv: FloatField,
        tpert_out: FloatFieldIJ,
        qpert_out: FloatFieldIJ,
        formulation: SaturationFormulation = SaturationFormulation.Staars,
    ):

        # University of Washington Shallow Convection Scheme

        # Described in Park and Bretherton. 2008. J. Climate :

        # 'The University of Washington shallow convection and
        # moist turbulent schemes and their impact on climate
        # simulations with the Community Atmosphere Model'
        # Coded in CESM by Sungsu Park. Oct.2005. May.2008.
        # Coded in GEOS by Nathan Arnold. July 2016.
        # NDSL Port by Katrina Fandrich. May 2025.

        # For general questions, email sungsup@ucar.edu or sungsu@atmos.washington.edu
        # For GEOS-specific questions, email nathan.arnold@nasa.gov
        # For NDSL-specific questions, email katrina.fandrich@nasa.gov

        ##############################################################################

        # General Notes on ComputeUwshcu:

        # Masks:

        # - There are several masks used throughout ComputeUwshcu to replicate
        #  some of the 'go tos' in the original Fortran code (e.g., go to 333).
        # - These masks are set to False by default, but can become True if the
        # calculation converges.

        # - The most frequently used mask is the 'condensation' mask, which
        # indicates if the condensation process has occurred.
        # - If condensation occurs, the mask becomes True and UW Shallow
        # Convection scheme is done for that column.

        # Special naming of variables:

        # - The suffix '_s' is used to indicate variables computed during the
        # first iteration of the iterative CIN calculation (e.g, s0_s, u0_s,
        # v0_s)

        # - The suffix '_o' is used to indicate variables computed during the
        # second iteration of the iterative CIN calculation (e.g, s0_o, u0_o,
        # v0_o)

        # - The suffix '_IJ' is used to indicate variables that were 3D at an
        # earlier stage in the
        #   computation, but have been converted to 2D (e.g., cin_IJ,
        # cinlcl_IJ).

        ######################################################################

        self.qsat = QSat(
            self.stencil_factory,
            self.quantity_factory,
            formulation=formulation,
        )

        # Initialize masks, default for all masks is False.
        self.condensation.view[:, :] = False
        self.stop_cin.view[:, :] = False
        self.stop_buoyancy_sort.view[:, :] = False

        self._compute_uwshcu_invert_before(
            # Inputs
            k0=k0,
            pmid0_inv=pmid0_inv,
            u0_inv=u0_inv,
            v0_inv=v0_inv,
            zmid0_inv=zmid0_inv,
            exnmid0_inv=exnmid0_inv,
            dp0_inv=dp0_inv,
            qv0_inv=qv0_inv,
            ql0_inv=ql0_inv,
            qi0_inv=qi0_inv,
            t0_inv=t0_inv,
            tke_inv=tke_inv,
            pifc0_inv=pifc0_inv,
            zifc0_inv=zifc0_inv,
            exnifc0_inv=exnifc0_inv,
            kpbl_inv=kpbl_inv,
            cnvtr=cnvtr,
            frland=frland,
            CNV_Tracers=CNV_Tracers,
            tr0_inout=self.tr0_inout,
            # Outputs
            pmid0_in=self.temporaries.pmid0_in,
            u0_in=self.temporaries.u0_in,
            v0_in=self.temporaries.v0_in,
            zmid0_in=self.temporaries.zmid0_in,
            exnmid0_in=self.temporaries.exnmid0_in,
            dp0_in=self.temporaries.dp0_in,
            qv0_in=self.temporaries.qv0_in,
            ql0_in=self.temporaries.ql0_in,
            qi0_in=self.temporaries.qi0_in,
            th0_in=self.temporaries.th0_in,
            tke_in=self.temporaries.tke_in,
            pifc0_in=self.temporaries.pifc0_in,
            zifc0_in=self.temporaries.zifc0_in,
            exnifc0_in=self.temporaries.exnifc0_in,
            kpbl_in=self.temporaries.kpbl_in,
            cnvtrmax=self.temporaries.cnvtrmax,
        )

        self._compute_thermodynamic_variables(
            pmid0_in=self.temporaries.pmid0_in,
            zmid0_in=self.temporaries.zmid0_in,
            exnmid0_in=self.temporaries.exnmid0_in,
            dp0_in=self.temporaries.dp0_in,
            u0_in=self.temporaries.u0_in,
            v0_in=self.temporaries.v0_in,
            u0=self.temporaries.u0,
            v0=self.temporaries.v0,
            qv0_in=self.temporaries.qv0_in,
            ql0_in=self.temporaries.ql0_in,
            qi0_in=self.temporaries.qi0_in,
            th0_in=self.temporaries.th0_in,
            tr0_inout=self.tr0_inout,
            cush_inout=cush,
            cush=cush,
            umf_out=self.temporaries.umf_out,
            shfx=shfx,
            evap=evap,
            qtflx_out=self.temporaries.qtflx_out,
            slflx_out=self.temporaries.slflx_out,
            uflx_out=self.temporaries.uflx_out,
            vflx_out=self.temporaries.vflx_out,
            qt0=self.temporaries.qt0,
            t0=self.temporaries.t0,
            qv0=self.temporaries.qv0,
            qi0=self.temporaries.qi0,
            pmid0=self.temporaries.pmid0,
            tr0=self.tr0,
            tr0_temp=self.temporaries.tr0_temp,
            sstr0=self.sstr0,
            thl0=self.temporaries.thl0,
            ssthl0=self.temporaries.ssthl0,
            ssqt0=self.temporaries.ssqt0,
            ssu0=self.temporaries.ssu0,
            ssv0=self.temporaries.ssv0,
            tscaleh=self.temporaries.tscaleh,
            dotransport=dotransport,
            fer_out=self.temporaries.fer_out,
            fdr_out=self.temporaries.fdr_out,
            tpert_out=tpert_out,
            qpert_out=qpert_out,
        )

        self._compute_thv0_thvl0(
            pmid0_in=self.temporaries.pmid0_in,
            exnmid0_in=self.temporaries.exnmid0_in,
            qv0_in=self.temporaries.qv0_in,
            ql0_in=self.temporaries.ql0_in,
            qi0_in=self.temporaries.qi0_in,
            th0_in=self.temporaries.th0_in,
            zmid0=self.temporaries.zmid0_in,
            pifc0_in=self.temporaries.pifc0_in,
            ssthl0=self.temporaries.ssthl0,
            ssqt0=self.temporaries.ssqt0,
            ese=self.qsat.ese,
            esx=self.qsat.esx,
            umf_out=self.temporaries.umf_out,
            qtflx_out=self.temporaries.qtflx_out,
            slflx_out=self.temporaries.slflx_out,
            uflx_out=self.temporaries.uflx_out,
            vflx_out=self.temporaries.vflx_out,
            u0_in=self.temporaries.u0_in,
            v0_in=self.temporaries.v0_in,
            k0=k0,
            condensation=self.condensation,
            dotransport=dotransport,
            ssu0=self.temporaries.ssu0,
            ssv0=self.temporaries.ssv0,
            tr0=self.tr0,
            sstr0=self.sstr0,
            tr0_o=self.tr0_o,
            sstr0_o=self.sstr0_o,
            trflx=self.trflx,
            trten=self.trten,
            tru=self.tru,
            tru_emf=self.tru_emf,
            umf_zint=self.temporaries.umf_zint,
            emf=self.temporaries.emf,
            slflx=self.temporaries.slflx,
            qtflx=self.temporaries.qtflx,
            uflx=self.temporaries.uflx,
            vflx=self.temporaries.vflx,
            thlu=self.temporaries.thlu,
            qtu=self.temporaries.qtu,
            uu=self.temporaries.uu,
            vu=self.temporaries.vu,
            wu=self.temporaries.wu,
            thvu=self.temporaries.thvu,
            thlu_emf=self.temporaries.thlu_emf,
            qtu_emf=self.temporaries.qtu_emf,
            uu_emf=self.temporaries.uu_emf,
            vu_emf=self.temporaries.vu_emf,
            uemf=self.temporaries.uemf,
            thvl0bot=self.temporaries.thvl0bot,
            thvl0top=self.temporaries.thvl0top,
            thvl0=self.temporaries.thvl0,
            qt0=self.temporaries.qt0,
            t0=self.temporaries.t0,
            qv0=self.temporaries.qv0,
            ql0=self.temporaries.ql0,
            qi0=self.temporaries.qi0,
            thl0=self.temporaries.thl0,
            thv0bot=self.temporaries.thv0bot,
            thv0top=self.temporaries.thv0top,
            uten=self.temporaries.uten,
            vten=self.temporaries.vten,
            s0=self.temporaries.s0,
            qcu=self.temporaries.qcu,
            qlu=self.temporaries.qlu,
            qiu=self.temporaries.qiu,
            cufrc=self.temporaries.cufrc,
            ufrc=self.temporaries.ufrc,
            qlten_det=self.temporaries.qlten_det,
            qiten_det=self.temporaries.qiten_det,
            qlten_sink=self.temporaries.qlten_sink,
            qiten_sink=self.temporaries.qiten_sink,
            sten=self.temporaries.sten,
            slten=self.temporaries.slten,
            qiten=self.temporaries.qiten,
            qv0_o=self.temporaries.qv0_o,
            ql0_o=self.temporaries.ql0_o,
            qi0_o=self.temporaries.qi0_o,
            t0_o=self.temporaries.t0_o,
            s0_o=self.temporaries.s0_o,
            u0_o=self.temporaries.u0_o,
            v0_o=self.temporaries.v0_o,
            qt0_o=self.temporaries.qt0_o,
            thl0_o=self.temporaries.thl0_o,
            thvl0_o=self.temporaries.thvl0_o,
            ssthl0_o=self.temporaries.ssthl0_o,
            ssqt0_o=self.temporaries.ssqt0_o,
            thv0bot_o=self.temporaries.thv0bot_o,
            thv0top_o=self.temporaries.thv0top_o,
            thvl0bot_o=self.temporaries.thvl0bot_o,
            thvl0top_o=self.temporaries.thvl0top_o,
            ssu0_o=self.temporaries.ssu0_o,
            ssv0_o=self.temporaries.ssv0_o,
            cush_inout=self.temporaries.cush_inout,
        )

        # This 'iteration' loop is for implicit CIN closure

        # It is important to note that this iterative cin loop is located at the outer
        # shell of the code. Thus, source air properties can also be changed during the
        # iterative cin calculation, because cumulus convection induces non-zero fluxes
        # even at interfaces below PBL top height through 'fluxbelowinv' calculation.

        for iteration in range(iter_cin):
            iteration = i32(iteration)

            self._find_pbl_height(
                iteration=iteration,
                kpbl_in=self.temporaries.kpbl_in,
                k0=k0,
                condensation=self.condensation,
                umf_out=self.temporaries.umf_out,
                qtflx_out=self.temporaries.qtflx_out,
                slflx_out=self.temporaries.slflx_out,
                uflx_out=self.temporaries.uflx,
                vflx_out=self.temporaries.vflx,
                kinv=self.temporaries.kinv,
                cush=cush,
                tscaleh=self.temporaries.tscaleh,
                cush_inout=self.temporaries.cush_inout,
            )

            self._find_pbl_averages(
                condensation=self.condensation,
                thvl0bot=self.temporaries.thvl0bot,
                thvl0top=self.temporaries.thvl0top,
                kinv=self.temporaries.kinv,
                pifc0=self.temporaries.pifc0_in,
                tke_in=self.temporaries.tke_in,
                u0=self.temporaries.u0_in,
                v0=self.temporaries.v0_in,
                thvl0=self.temporaries.thvl0,
                zmid0=self.temporaries.zmid0_in,
                qtsrchgt=qtsrchgt,
                qt0=self.temporaries.qt0,
                thvlmin=self.temporaries.thvlmin,
                tkeavg=self.temporaries.tkeavg,
                uavg=self.temporaries.uavg,
                vavg=self.temporaries.vavg,
                thvlavg=self.temporaries.thvlavg,
                qtavg=self.temporaries.qtavg,
                iteration=iteration,
            )

            self._find_cumulus_characteristics(
                condensation=self.condensation,
                windsrcavg=windsrcavg,
                pifc0=self.temporaries.pifc0_in,
                t0=self.temporaries.t0,
                qv0=self.temporaries.qv0,
                shfx=shfx,
                evap=evap,
                thlsrc_fac=thlsrc_fac,
                qtsrc_fac=qtsrc_fac,
                qt0=self.temporaries.qt0,
                qtavg=self.temporaries.qtavg,
                thvlmin=self.temporaries.thvlmin,
                uavg=self.temporaries.uavg,
                vavg=self.temporaries.vavg,
                kinv=self.temporaries.kinv,
                u0=self.temporaries.u0_in,
                v0=self.temporaries.v0_in,
                ssu0=self.temporaries.ssu0,
                ssv0=self.temporaries.ssv0,
                pmid0=self.temporaries.pmid0_in,
                dotransport=dotransport,
                tr0=self.tr0,
                trsrc=self.trsrc,
                qtsrc=self.temporaries.qtsrc,
                thvlsrc=self.temporaries.thvlsrc,
                thlsrc=self.temporaries.thlsrc,
                usrc=self.temporaries.usrc,
                vsrc=self.temporaries.vsrc,
                tpert_out=tpert_out,
                qpert_out=qpert_out,
                iteration=iteration,
            )

            self._find_klcl(
                condensation=self.condensation,
                pifc0=self.temporaries.pifc0_in,
                umf_out=self.temporaries.umf_out,
                qtflx_out=self.temporaries.qtflx_out,
                slflx_out=self.temporaries.slflx_out,
                uflx_out=self.temporaries.uflx_out,
                vflx_out=self.temporaries.vflx_out,
                qtsrc=self.temporaries.qtsrc,
                thlsrc=self.temporaries.thlsrc,
                ese=self.qsat.ese,
                esx=self.qsat.esx,
                k0=k0,
                thl0=self.temporaries.thl0,
                ssthl0=self.temporaries.ssthl0,
                pmid0=self.temporaries.pmid0_in,
                qt0=self.temporaries.qt0,
                ssqt0=self.temporaries.ssqt0,
                plcl=self.temporaries.plcl,
                klcl=self.temporaries.klcl,
                thl0lcl=self.temporaries.thl0lcl,
                qt0lcl=self.temporaries.qt0lcl,
                thv0lcl=self.temporaries.thv0lcl,
                cush_inout=self.temporaries.cush_inout,
            )

            self._compute_cin_cinlcl(
                condensation=self.condensation,
                stop_cin=self.stop_cin,
                klcl=self.temporaries.klcl,
                kinv=self.temporaries.kinv,
                thvlsrc=self.temporaries.thvlsrc,
                pifc0=self.temporaries.pifc0_in,
                thv0bot=self.temporaries.thv0bot,
                thv0top=self.temporaries.thv0top,
                plcl=self.temporaries.plcl,
                thv0lcl=self.temporaries.thv0lcl,
                thlsrc=self.temporaries.thlsrc,
                qtsrc=self.temporaries.qtsrc,
                ese=self.qsat.ese,
                esx=self.qsat.esx,
                umf_out=self.temporaries.umf_out,
                qtflx_out=self.temporaries.qtflx_out,
                slflx_out=self.temporaries.slflx_out,
                uflx_out=self.temporaries.uflx_out,
                vflx_out=self.temporaries.vflx_out,
                cin_IJ=self.temporaries.cin_IJ,
                cinlcl_IJ=self.temporaries.cinlcl_IJ,
                plfc_IJ=self.temporaries.plfc_IJ,
                klfc_IJ=self.temporaries.klfc_IJ,
                plfc=self.temporaries.plfc,
                klfc=self.temporaries.klfc,
                cin=self.temporaries.cin,
                thvubot=self.temporaries.thvubot,
                thvutop=self.temporaries.thvutop,
                k0=k0,
                iteration=iteration,
                rbuoy=rbuoy,
                rkfre=rkfre,
                tkeavg=self.temporaries.tkeavg,
                epsvarw=epsvarw,
                thvlmin=self.temporaries.thvlmin,
                usrc=self.temporaries.usrc,
                vsrc=self.temporaries.vsrc,
                dotransport=dotransport,
                trsrc=self.trsrc,
                trsrc_o=self.trsrc_o,
                cin_i=self.temporaries.cin_i,
                cinlcl_i=self.temporaries.cinlcl_i,
                ke=self.temporaries.ke,
                kinv_o=self.temporaries.kinv_o,
                klcl_o=self.temporaries.klcl_o,
                klfc_o=self.temporaries.klfc_o,
                plcl_o=self.temporaries.plcl_o,
                plfc_o=self.temporaries.plfc_o,
                tkeavg_o=self.temporaries.tkeavg_o,
                thvlmin_o=self.temporaries.thvlmin_o,
                qtsrc_o=self.temporaries.qtsrc_o,
                thvlsrc_o=self.temporaries.thvlsrc_o,
                thlsrc_o=self.temporaries.thlsrc_o,
                usrc_o=self.temporaries.usrc_o,
                vsrc_o=self.temporaries.vsrc_o,
                thv0lcl_o=self.temporaries.thv0lcl_o,
                cinlcl=self.temporaries.cinlcl,
                cush_inout=self.temporaries.cush_inout,
            )
            if iteration != i32(0):
                self._avg_initial_and_final_cin(
                    condensation=self.condensation,
                    iteration=iteration,
                    cin_IJ=self.temporaries.cin_IJ,
                    cinlcl_IJ=self.temporaries.cinlcl_IJ,
                    use_CINcin=use_CINcin,
                    cin_i=self.temporaries.cin_i,
                    cinlcl_i=self.temporaries.cinlcl_i,
                    ke=self.temporaries.ke,
                    dotransport=dotransport,
                    trsrc=self.trsrc,
                    trsrc_o=self.trsrc_o,
                    tr0=self.tr0,
                    tr0_o=self.tr0_o,
                    sstr0=self.sstr0,
                    sstr0_o=self.sstr0_o,
                    kinv_o=self.temporaries.kinv_o,
                    klcl_o=self.temporaries.klcl_o,
                    klfc_o=self.temporaries.klfc_o,
                    plcl_o=self.temporaries.plcl_o,
                    plfc_o=self.temporaries.plfc_o,
                    tkeavg_o=self.temporaries.tkeavg_o,
                    thvlmin_o=self.temporaries.thvlmin_o,
                    qtsrc_o=self.temporaries.qtsrc_o,
                    thvlsrc_o=self.temporaries.thvlsrc_o,
                    thlsrc_o=self.temporaries.thlsrc_o,
                    usrc_o=self.temporaries.usrc_o,
                    vsrc_o=self.temporaries.vsrc_o,
                    thv0lcl_o=self.temporaries.thv0lcl_o,
                    qv0_o=self.temporaries.qv0_o,
                    ql0_o=self.temporaries.ql0_o,
                    qi0_o=self.temporaries.qi0_o,
                    t0_o=self.temporaries.t0_o,
                    s0_o=self.temporaries.s0_o,
                    u0_o=self.temporaries.u0_o,
                    v0_o=self.temporaries.v0_o,
                    qt0_o=self.temporaries.qt0_o,
                    thl0_o=self.temporaries.thl0_o,
                    thvl0_o=self.temporaries.thvl0_o,
                    ssthl0_o=self.temporaries.ssthl0_o,
                    ssqt0_o=self.temporaries.ssqt0_o,
                    thv0bot_o=self.temporaries.thv0bot_o,
                    thv0top_o=self.temporaries.thv0top_o,
                    thvl0bot_o=self.temporaries.thvl0bot_o,
                    thvl0top_o=self.temporaries.thvl0top_o,
                    ssu0_o=self.temporaries.ssu0_o,
                    ssv0_o=self.temporaries.ssv0_o,
                    thvlmin_IJ=self.temporaries.thvlmin_IJ,
                    umf_zint=self.temporaries.umf_zint,
                    emf=self.temporaries.emf,
                    slflx=self.temporaries.slflx,
                    qtflx=self.temporaries.qtflx,
                    uflx=self.temporaries.uflx,
                    vflx=self.temporaries.vflx,
                    k0=k0,
                    ufrc=self.temporaries.ufrc,
                    thlu=self.temporaries.thlu,
                    qtu=self.temporaries.qtu,
                    uu=self.temporaries.uu,
                    vu=self.temporaries.vu,
                    wu=self.temporaries.wu,
                    thvu=self.temporaries.thvu,
                    thlu_emf=self.temporaries.thlu_emf,
                    qtu_emf=self.temporaries.qtu_emf,
                    uu_emf=self.temporaries.uu_emf,
                    vu_emf=self.temporaries.vu_emf,
                    trflx=self.trflx,
                    trten=self.trten,
                    tru=self.tru,
                    tru_emf=self.tru_emf,
                    umf_s=self.temporaries.umf_s,
                    zifc0=self.temporaries.zifc0_in,
                    dcm_s=self.temporaries.dcm_s,
                    qvten_s=self.temporaries.qvten_s,
                    qlten_s=self.temporaries.qlten_s,
                    qiten_s=self.temporaries.qiten_s,
                    sten_s=self.temporaries.sten_s,
                    uten_s=self.temporaries.uten_s,
                    vten_s=self.temporaries.vten_s,
                    qrten_s=self.temporaries.qrten_s,
                    qsten_s=self.temporaries.qsten_s,
                    qldet_s=self.temporaries.qldet_s,
                    qidet_s=self.temporaries.qidet_s,
                    qlsub_s=self.temporaries.qlsub_s,
                    qisub_s=self.temporaries.qisub_s,
                    cush_s=self.temporaries.cush_s,
                    cufrc_s=self.temporaries.cufrc_s,
                    qtflx_out=self.temporaries.qtflx_out,
                    qtflx_s=self.temporaries.qtflx_s,
                    slflx_out=self.temporaries.slflx_out,
                    slflx_s=self.temporaries.slflx_s,
                    uflx_out=self.temporaries.uflx_out,
                    uflx_s=self.temporaries.uflx_s,
                    vflx_out=self.temporaries.vflx_out,
                    vflx_s=self.temporaries.vflx_s,
                    fer_s=self.temporaries.fer_s,
                    fdr_s=self.temporaries.fdr_s,
                    umf_out=self.temporaries.umf_out,
                    kinv=self.temporaries.kinv,
                    klcl=self.temporaries.klcl,
                    plcl=self.temporaries.plcl,
                    thv0bot=self.temporaries.thv0bot,
                    thv0lcl=self.temporaries.thv0lcl,
                    thv0top=self.temporaries.thv0top,
                    thlsrc=self.temporaries.thlsrc,
                    qtsrc=self.temporaries.qtsrc,
                    thl0=self.temporaries.thl0,
                    ssthl0=self.temporaries.ssthl0,
                    ssqt0=self.temporaries.ssqt0,
                    ssu0=self.temporaries.ssu0,
                    ssv0=self.temporaries.ssv0,
                    qt0=self.temporaries.qt0,
                    u0=self.temporaries.u0_in,
                    v0=self.temporaries.v0_in,
                    qi0=self.temporaries.qi0,
                    ql0=self.temporaries.ql0,
                    qv0=self.temporaries.qv0,
                    s0=self.temporaries.s0,
                    qvten_out=self.temporaries.qvten_out,
                    dcm_out=self.temporaries.dcm_out,
                    qlten_out=self.temporaries.qlten_out,
                    qiten_out=self.temporaries.qiten_out,
                    del_CIN=self.temporaries.del_CIN,
                    dcm=self.temporaries.dcm,
                    sten_out=self.temporaries.sten_out,
                    uten_out=self.temporaries.uten_out,
                    vten_out=self.temporaries.vten_out,
                    qrten_out=self.temporaries.qrten_out,
                    qsten_out=self.temporaries.qsten_out,
                    cufrc_out=self.temporaries.cufrc_out,
                    cush_inout=self.temporaries.cush_inout,
                    qldet_out=self.temporaries.qldet_out,
                    qidet_out=self.temporaries.qidet_out,
                    qlsub_out=self.temporaries.qlsub_out,
                    qisub_out=self.temporaries.qisub_out,
                    fdr_out=self.temporaries.fdr_out,
                    fer_out=self.temporaries.fer_out,
                )

            self._define_prel_krel(
                condensation=self.condensation,
                iteration=iteration,
                klcl=self.temporaries.klcl,
                kinv=self.temporaries.kinv,
                pifc0=self.temporaries.pifc0_in,
                thv0bot=self.temporaries.thv0bot,
                plcl=self.temporaries.plcl,
                thv0lcl=self.temporaries.thv0lcl,
                krel=self.temporaries.krel,
                prel=self.temporaries.prel,
                thv0rel=self.temporaries.thv0rel,
                umf_out=self.temporaries.umf_out,
            )

            self._calc_cumulus_base_mass_flux(
                condensation=self.condensation,
                iteration=iteration,
                use_CINcin=use_CINcin,
                cin_IJ=self.temporaries.cin_IJ,
                rbuoy=rbuoy,
                cinlcl_IJ=self.temporaries.cinlcl_IJ,
                rkfre=rkfre,
                tkeavg=self.temporaries.tkeavg,
                epsvarw=epsvarw,
                umf_out=self.temporaries.umf_out,
                qtflx_out=self.temporaries.qtflx_out,
                slflx_out=self.temporaries.slflx_out,
                uflx_out=self.temporaries.uflx_out,
                vflx_out=self.temporaries.vflx_out,
                kinv=self.temporaries.kinv,
                pifc0=self.temporaries.pifc0_in,
                thv0top=self.temporaries.thv0top,
                exnifc0=self.temporaries.exnifc0_in,
                dp0=self.temporaries.dp0_in,
                dt=dt,
                mumin1=mumin1,
                rmaxfrac=rmaxfrac,
                winv=self.temporaries.winv,
                cbmf=self.temporaries.cbmf,
                rho0inv=self.temporaries.rho0inv,
                ufrcinv=self.temporaries.ufrcinv,
                wcrit=self.temporaries.wcrit,
                cush_inout=self.temporaries.cush_inout,
            )

            self._define_updraft_properties(
                condensation=self.condensation,
                iteration=iteration,
                winv=self.temporaries.winv,
                cinlcl_IJ=self.temporaries.cinlcl_IJ,
                rbuoy=rbuoy,
                umf_out=self.temporaries.umf_out,
                qtflx_out=self.temporaries.qtflx_out,
                slflx_out=self.temporaries.slflx_out,
                uflx_out=self.temporaries.uflx_out,
                vflx_out=self.temporaries.vflx_out,
                cbmf=self.temporaries.cbmf,
                rho0inv=self.temporaries.rho0inv,
                krel=self.temporaries.krel,
                ufrc=self.temporaries.ufrc,
                ufrcinv=self.temporaries.ufrcinv,
                kinv=self.temporaries.kinv,
                umf_zint=self.temporaries.umf_zint,
                wu=self.temporaries.wu,
                emf=self.temporaries.emf,
                thlu=self.temporaries.thlu,
                qtu=self.temporaries.qtu,
                thlsrc=self.temporaries.thlsrc,
                qtsrc=self.temporaries.qtsrc,
                prel=self.temporaries.prel,
                ese=self.qsat.ese,
                esx=self.qsat.esx,
                thvu=self.temporaries.thvu,
                wlcl=self.temporaries.wlcl,
                ufrclcl=self.temporaries.ufrclcl,
                cush_inout=self.temporaries.cush_inout,
            )

            self._define_env_properties(
                condensation=self.condensation,
                iteration=iteration,
                krel=self.temporaries.krel,
                kinv=self.temporaries.kinv,
                PGFc=PGFc,
                ssu0=self.temporaries.ssu0,
                ssv0=self.temporaries.ssv0,
                prel=self.temporaries.prel,
                pifc0=self.temporaries.pifc0_in,
                uu=self.temporaries.uu,
                vu=self.temporaries.vu,
                usrc=self.temporaries.usrc,
                vsrc=self.temporaries.vsrc,
                dotransport=dotransport,
                tru=self.tru,
                trsrc=self.trsrc,
                thv0rel=self.temporaries.thv0rel,
                thl0=self.temporaries.thl0,
                ssthl0=self.temporaries.ssthl0,
                pmid0=self.temporaries.pmid0_in,
                qt0=self.temporaries.qt0,
                ssqt0=self.temporaries.ssqt0,
                u0=self.temporaries.u0_in,
                v0=self.temporaries.v0_in,
                tre=self.tre,
                tr0=self.tr0,
                sstr0=self.sstr0,
                uplus=self.temporaries.uplus,
                vplus=self.temporaries.vplus,
                uplus_3D=self.temporaries.uplus_3D,
                vplus_3D=self.temporaries.vplus_3D,
                qsat_pe=self.temporaries.qsat_pe,
                pe=self.temporaries.pe,
                thle=self.temporaries.thle,
                qte=self.temporaries.qte,
                dpe=self.temporaries.dpe,
                exne=self.temporaries.exne,
                thvebot=self.temporaries.thvebot,
                ue=self.temporaries.ue,
                ve=self.temporaries.ve,
            )

            self._buoyancy_sorting(
                condensation=self.condensation,
                tscaleh=self.temporaries.tscaleh,
                krel=self.temporaries.krel,
                wlcl=self.temporaries.wlcl,
                prel=self.temporaries.prel,
                pifc0=self.temporaries.pifc0_in,
                thv0rel=self.temporaries.thv0rel,
                thl0=self.temporaries.thl0,
                ssthl0=self.temporaries.ssthl0,
                pmid0=self.temporaries.pmid0_in,
                qt0=self.temporaries.qt0,
                ssqt0=self.temporaries.ssqt0,
                u0=self.temporaries.u0_in,
                v0=self.temporaries.v0_in,
                ssu0=self.temporaries.ssu0,
                ssv0=self.temporaries.ssv0,
                dotransport=dotransport,
                tre=self.tre,
                tr0=self.tr0,
                sstr0=self.sstr0,
                k0=k0,
                thlu=self.temporaries.thlu,
                qtu=self.temporaries.qtu,
                wu=self.temporaries.wu,
                niter_xc=niter_xc,
                ese=self.qsat.ese,
                esx=self.qsat.esx,
                umf_out=self.temporaries.umf_out,
                qtflx_out=self.temporaries.qtflx_out,
                slflx_out=self.temporaries.slflx_out,
                uflx_out=self.temporaries.uflx_out,
                vflx_out=self.temporaries.vflx_out,
                qsat_pe=self.temporaries.qsat_pe,
                criqc=criqc,
                cridist_opt=cridist_opt,
                rle=rle,
                zifc0=self.temporaries.zifc0_in,
                rbuoy=rbuoy,
                mixscale=mixscale,
                rkm=rkm,
                zmid0=self.temporaries.zmid0_in,
                detrhgt=detrhgt,
                thlue=self.temporaries.thlue,
                qtue=self.temporaries.qtue,
                wue=self.temporaries.wue,
                wtwb=self.temporaries.wtwb,
                dp0=self.temporaries.dp0_in,
                dt=dt,
                thv0bot=self.temporaries.thv0bot,
                exnmid0=self.temporaries.exnmid0_in,
                rmaxfrac=rmaxfrac,
                thv0top=self.temporaries.thv0top,
                exnifc0=self.temporaries.exnifc0_in,
                use_self_detrain=use_self_detrain,
                rdrag=rdrag,
                PGFc=PGFc,
                tru=self.tru,
                umf_zint=self.temporaries.umf_zint,
                emf=self.temporaries.emf,
                thvu=self.temporaries.thvu,
                rei=self.temporaries.rei,
                uu=self.temporaries.uu,
                vu=self.temporaries.vu,
                ufrc=self.temporaries.ufrc,
                pe=self.temporaries.pe,
                thle=self.temporaries.thle,
                qte=self.temporaries.qte,
                dpe=self.temporaries.dpe,
                exne=self.temporaries.exne,
                thvebot=self.temporaries.thvebot,
                ue=self.temporaries.ue,
                ve=self.temporaries.ve,
                drage=self.temporaries.drage,
                bogbot=self.temporaries.bogbot,
                bogtop=self.temporaries.bogtop,
                kpen_IJ=self.temporaries.kpen_IJ,
                rhomid0j=self.temporaries.rhomid0j,
                kbup_IJ=self.temporaries.kbup_IJ,
                fer=self.temporaries.fer,
                fdr=self.temporaries.fdr,
                dwten=self.temporaries.dwten,
                diten=self.temporaries.diten,
                dcm=self.temporaries.dcm,
                xco=self.temporaries.xco,
                cush_inout=self.temporaries.cush_inout,
                stop_buoyancy_sort=self.stop_buoyancy_sort,
                iteration=iteration,
            )

            self._calc_ppen(
                condensation=self.condensation,
                drage=self.temporaries.drage,
                bogbot=self.temporaries.bogbot,
                bogtop=self.temporaries.bogtop,
                pifc0=self.temporaries.pifc0_in,
                kpen_IJ=self.temporaries.kpen_IJ,
                kpen=self.temporaries.kpen,
                wu=self.temporaries.wu,
                rhomid0j=self.temporaries.rhomid0j,
                dp0=self.temporaries.dp0_in,
                wtwb=self.temporaries.wtwb,
                ppen=self.temporaries.ppen,
            )

            self._recalc_condensate(
                condensation=self.condensation,
                fer=self.temporaries.fer,
                kpen=self.temporaries.kpen,
                ppen=self.temporaries.ppen,
                thlu=self.temporaries.thlu,
                thl0=self.temporaries.thl0,
                ssthl0=self.temporaries.ssthl0,
                qtu=self.temporaries.qtu,
                qt0=self.temporaries.qt0,
                ssqt0=self.temporaries.ssqt0,
                pifc0=self.temporaries.pifc0_in,
                ese=self.qsat.ese,
                esx=self.qsat.esx,
                umf_out=self.temporaries.umf_out,
                qtflx_out=self.temporaries.qtflx_out,
                slflx_out=self.temporaries.slflx_out,
                uflx_out=self.temporaries.uflx_out,
                vflx_out=self.temporaries.vflx_out,
                criqc=criqc,
                thv0bot=self.temporaries.thv0bot,
                thv0top=self.temporaries.thv0top,
                exnifc0=self.temporaries.exnifc0_in,
                zifc0=self.temporaries.zifc0_in,
                kbup_IJ=self.temporaries.kbup_IJ,
                kbup=self.temporaries.kbup,
                krel=self.temporaries.krel,
                k0=k0,
                umf_zint=self.temporaries.umf_zint,
                emf=self.temporaries.emf,
                ufrc=self.temporaries.ufrc,
                dwten=self.temporaries.dwten,
                diten=self.temporaries.diten,
                dwten_temp=self.temporaries.dwten_temp,
                diten_temp=self.temporaries.diten_temp,
                thlu_top=self.temporaries.thlu_top,
                qtu_top=self.temporaries.qtu_top,
                cldhgt=self.temporaries.cldhgt,
                fdr=self.temporaries.fdr,
                umf_temp=self.temporaries.umf_temp,
                xco=self.temporaries.xco,
                cush=cush,
                cush_inout=self.temporaries.cush_inout,
                iteration=iteration,
            )

            self._calc_entrainment_mass_flux(
                condensation=self.condensation,
                k0=k0,
                thlu=self.temporaries.thlu,
                qtu=self.temporaries.qtu,
                uu=self.temporaries.uu,
                vu=self.temporaries.vu,
                tru=self.tru,
                dotransport=dotransport,
                tru_emf=self.tru_emf,
                kpen=self.temporaries.kpen,
                kbup=self.temporaries.kbup,
                pifc0=self.temporaries.pifc0_in,
                thv0bot=self.temporaries.thv0bot,
                thv0top=self.temporaries.thv0top,
                exnifc0=self.temporaries.exnifc0_in,
                umf_zint=self.temporaries.umf_zint,
                ppen=self.temporaries.ppen,
                rei=self.temporaries.rei,
                rpen=rpen,
                dp0=self.temporaries.dp0_in,
                dt=dt,
                thl0=self.temporaries.thl0,
                ssthl0=self.temporaries.ssthl0,
                pmid0=self.temporaries.pmid0_in,
                qt0=self.temporaries.qt0,
                ssqt0=self.temporaries.ssqt0,
                u0=self.temporaries.u0_in,
                ssu0=self.temporaries.ssu0,
                v0=self.temporaries.v0_in,
                ssv0=self.temporaries.ssv0,
                tr0=self.tr0,
                sstr0=self.sstr0,
                use_cumpenent=use_cumpenent,
                thlu_emf=self.temporaries.thlu_emf,
                qtu_emf=self.temporaries.qtu_emf,
                uu_emf=self.temporaries.uu_emf,
                vu_emf=self.temporaries.vu_emf,
                emf=self.temporaries.emf,
                iteration=iteration,
            )

            self._calc_pbl_fluxes(
                condensation=self.condensation,
                qtsrc=self.temporaries.qtsrc,
                qt0=self.temporaries.qt0,
                ssqt0=self.temporaries.ssqt0,
                pifc0=self.temporaries.pifc0_in,
                pmid0=self.temporaries.pmid0_in,
                kinv=self.temporaries.kinv,
                cbmf=self.temporaries.cbmf,
                dt=dt,
                xflx=self.temporaries.xflx,
                qtflx=self.temporaries.qtflx,
                thlsrc=self.temporaries.thlsrc,
                thl0=self.temporaries.thl0,
                ssthl0=self.temporaries.ssthl0,
                exnifc0=self.temporaries.exnifc0_in,
                usrc=self.temporaries.usrc,
                u0=self.temporaries.u0_in,
                ssu0=self.temporaries.ssu0,
                vsrc=self.temporaries.vsrc,
                v0=self.temporaries.v0_in,
                ssv0=self.temporaries.ssv0,
                dotransport=dotransport,
                trsrc=self.trsrc,
                tr0=self.tr0,
                sstr0=self.sstr0,
                trflx=self.trflx,
                uflx=self.temporaries.uflx,
                vflx=self.temporaries.vflx,
                slflx=self.temporaries.slflx,
                iteration=iteration,
                xflx_ndim=self.xflx_ndim,
            )

            self._non_buoyancy_sorting_fluxes(
                condensation=self.condensation,
                kinv=self.temporaries.kinv,
                krel=self.temporaries.krel,
                cbmf=self.temporaries.cbmf,
                qtsrc=self.temporaries.qtsrc,
                qt0=self.temporaries.qt0,
                ssqt0=self.temporaries.ssqt0,
                pifc0=self.temporaries.pifc0_in,
                pmid0=self.temporaries.pmid0_in,
                thlsrc=self.temporaries.thlsrc,
                thl0=self.temporaries.thl0,
                ssthl0=self.temporaries.ssthl0,
                PGFc=PGFc,
                exnifc0=self.temporaries.exnifc0_in,
                ssu0=self.temporaries.ssu0,
                ssv0=self.temporaries.ssv0,
                u0=self.temporaries.u0_in,
                v0=self.temporaries.v0_in,
                usrc=self.temporaries.usrc,
                vsrc=self.temporaries.vsrc,
                dotransport=dotransport,
                trflx=self.trflx,
                trsrc=self.trsrc,
                tr0=self.tr0,
                sstr0=self.sstr0,
                uflx=self.temporaries.uflx,
                vflx=self.temporaries.vflx,
                slflx=self.temporaries.slflx,
                qtflx=self.temporaries.qtflx,
            )

            self._buoyancy_sorting_fluxes(
                condensation=self.condensation,
                kbup=self.temporaries.kbup,
                krel=self.temporaries.krel,
                exnifc0=self.temporaries.exnifc0_in,
                umf_zint=self.temporaries.umf_zint,
                thlu=self.temporaries.thlu,
                thl0=self.temporaries.thl0,
                ssthl0=self.temporaries.ssthl0,
                pifc0=self.temporaries.pifc0_in,
                pmid0=self.temporaries.pmid0_in,
                qtu=self.temporaries.qtu,
                qt0=self.temporaries.qt0,
                ssqt0=self.temporaries.ssqt0,
                uu=self.temporaries.uu,
                u0=self.temporaries.u0_in,
                v0=self.temporaries.v0_in,
                vu=self.temporaries.vu,
                ssu0=self.temporaries.ssu0,
                ssv0=self.temporaries.ssv0,
                dotransport=dotransport,
                trflx=self.trflx,
                tru=self.tru,
                tr0=self.tr0,
                sstr0=self.sstr0,
                qtflx=self.temporaries.qtflx,
                vflx=self.temporaries.vflx,
                uflx=self.temporaries.uflx,
                slflx=self.temporaries.slflx,
                iteration=iteration,
            )

            self._penetrative_entrainment_fluxes(
                condensation=self.condensation,
                kbup=self.temporaries.kbup,
                kpen=self.temporaries.kpen,
                exnifc0=self.temporaries.exnifc0_in,
                emf=self.temporaries.emf,
                thlu_emf=self.temporaries.thlu_emf,
                thl0=self.temporaries.thl0,
                ssthl0=self.temporaries.ssthl0,
                pifc0=self.temporaries.pifc0_in,
                pmid0=self.temporaries.pmid0_in,
                qtu_emf=self.temporaries.qtu_emf,
                qt0=self.temporaries.qt0,
                ssqt0=self.temporaries.ssqt0,
                uu_emf=self.temporaries.uu_emf,
                vu_emf=self.temporaries.vu_emf,
                u0=self.temporaries.u0_in,
                v0=self.temporaries.v0_in,
                ssu0=self.temporaries.ssu0,
                ssv0=self.temporaries.ssv0,
                dotransport=dotransport,
                trflx=self.trflx,
                tru_emf=self.tru_emf,
                tr0=self.tr0,
                sstr0=self.sstr0,
                use_momenflx=use_momenflx,
                cbmf=self.temporaries.cbmf,
                uflx=self.temporaries.uflx,
                vflx=self.temporaries.vflx,
                qtflx=self.temporaries.qtflx,
                slflx=self.temporaries.slflx,
                uemf=self.temporaries.uemf,
                kinv=self.temporaries.kinv,
                krel=self.temporaries.krel,
                umf_zint=self.temporaries.umf_zint,
                k0=k0,
                ql0=self.temporaries.ql0,
                qi0=self.temporaries.qi0,
                dt=dt,
                ese=self.qsat.ese,
                esx=self.qsat.esx,
                umf_out=self.temporaries.umf_out,
                qtflx_out=self.temporaries.qtflx_out,
                slflx_out=self.temporaries.slflx_out,
                uflx_out=self.temporaries.uflx_out,
                vflx_out=self.temporaries.vflx_out,
                qlten_sink=self.temporaries.qlten_sink,
                qiten_sink=self.temporaries.qiten_sink,
                cush_inout=self.temporaries.cush_inout,
                iteration=iteration,
            )

            self._calc_momentum_tendency(
                condensation=self.condensation,
                kpen=self.temporaries.kpen,
                uflx=self.temporaries.uflx,
                vflx=self.temporaries.vflx,
                dp0=self.temporaries.dp0_in,
                u0=self.temporaries.u0_in,
                v0=self.temporaries.v0_in,
                dt=dt,
                uten=self.temporaries.uten,
                vten=self.temporaries.vten,
                uf=self.temporaries.uf,
                vf=self.temporaries.vf,
            )

            self._calc_thermodynamic_tendencies(
                condensation=self.condensation,
                kpen=self.temporaries.kpen,
                umf_zint=self.temporaries.umf_zint,
                dp0=self.temporaries.dp0_in,
                frc_rasn=frc_rasn,
                slflx=self.temporaries.slflx,
                uflx=self.temporaries.uflx,
                vflx=self.temporaries.vflx,
                u0=self.temporaries.u0_in,
                v0=self.temporaries.v0_in,
                uf=self.temporaries.uf,
                vf=self.temporaries.vf,
                dwten=self.temporaries.dwten,
                diten=self.temporaries.diten,
                dwten_temp=self.temporaries.dwten_temp,
                diten_temp=self.temporaries.diten_temp,
                umf_temp=self.temporaries.umf_temp,
                qtflx=self.temporaries.qtflx,
                krel=self.temporaries.krel,
                prel=self.temporaries.prel,
                thlu=self.temporaries.thlu,
                qtu=self.temporaries.qtu,
                ese=self.qsat.ese,
                esx=self.qsat.esx,
                umf_out=self.temporaries.umf_out,
                qtflx_out=self.temporaries.qtflx_out,
                slflx_out=self.temporaries.slflx_out,
                uflx_out=self.temporaries.uflx_out,
                vflx_out=self.temporaries.vflx_out,
                pifc0=self.temporaries.pifc0_in,
                ppen=self.temporaries.ppen,
                thlu_top=self.temporaries.thlu_top,
                qtu_top=self.temporaries.qtu_top,
                qlubelow=self.temporaries.qlubelow,
                qiubelow=self.temporaries.qiubelow,
                qlj_2D=self.temporaries.qlj_2D,
                qij_2D=self.temporaries.qij_2D,
                fdr=self.temporaries.fdr,
                ql0=self.temporaries.ql0,
                qi0=self.temporaries.qi0,
                kbup=self.temporaries.kbup,
                pmid0=self.temporaries.pmid0_in,
                thlu_emf=self.temporaries.thlu_emf,
                qtu_emf=self.temporaries.qtu_emf,
                emf=self.temporaries.emf,
                qlten_sink=self.temporaries.qlten_sink,
                qiten_sink=self.temporaries.qiten_sink,
                dt=dt,
                qrten=self.temporaries.qrten,
                qsten=self.temporaries.qsten,
                qvten=self.temporaries.qvten,
                qlten=self.temporaries.qlten,
                sten=self.temporaries.sten,
                qiten=self.temporaries.qiten,
                qc=self.temporaries.qc,
                slten=self.temporaries.slten,
                qlten_det=self.temporaries.qlten_det,
                qiten_det=self.temporaries.qiten_det,
                cush_inout=self.temporaries.cush_inout,
            )

            self._prevent_negative_condensate(
                condensation=self.condensation,
                qv0=self.temporaries.qv0,
                dt=dt,
                qvten=self.temporaries.qvten,
                ql0=self.temporaries.ql0,
                qlten=self.temporaries.qlten,
                qi0=self.temporaries.qi0,
                s0=self.temporaries.s0,
                sten=self.temporaries.sten,
                dp0=self.temporaries.dp0_in,
                qiten=self.temporaries.qiten,
                k0=k0,
                qmin=self.temporaries.qmin,
            )

            self._calc_tracer_tendencies(
                condensation=self.condensation,
                dotransport=dotransport,
                k0=k0,
                dt=dt,
                dp0=self.temporaries.dp0_in,
                trflx_d=self.trflx_d,
                trflx_u=self.trflx_u,
                trmin=self.trmin,
                tr0=self.tr0,
                trflx=self.trflx,
                trten=self.trten,
                iteration=iteration,
            )

            self._compute_diagnostic_outputs(
                condensation=self.condensation,
                prel=self.temporaries.prel,
                thlu=self.temporaries.thlu,
                qtu=self.temporaries.qtu,
                krel=self.temporaries.krel,
                ese=self.qsat.ese,
                esx=self.qsat.esx,
                umf_out=self.temporaries.umf_out,
                qtflx_out=self.temporaries.qtflx_out,
                slflx_out=self.temporaries.slflx_out,
                uflx_out=self.temporaries.uflx_out,
                vflx_out=self.temporaries.vflx_out,
                qcubelow=self.temporaries.qcubelow,
                qlubelow=self.temporaries.qlubelow,
                qiubelow=self.temporaries.qiubelow,
                rcwp=self.temporaries.rcwp,
                rlwp=self.temporaries.rlwp,
                riwp=self.temporaries.riwp,
                cush_inout=self.temporaries.cush_inout,
            )

            self._calc_cumulus_condensate_at_interfaces(
                condensation=self.condensation,
                krel=self.temporaries.krel,
                kpen=self.temporaries.kpen,
                pifc0=self.temporaries.pifc0_in,
                ppen=self.temporaries.ppen,
                thlu_top=self.temporaries.thlu_top,
                qtu_top=self.temporaries.qtu_top,
                thlu=self.temporaries.thlu,
                qtu=self.temporaries.qtu,
                ese=self.qsat.ese,
                esx=self.qsat.esx,
                umf_out=self.temporaries.umf_out,
                qtflx_out=self.temporaries.qtflx_out,
                slflx_out=self.temporaries.slflx_out,
                uflx_out=self.temporaries.uflx_out,
                vflx_out=self.temporaries.vflx_out,
                ufrc=self.temporaries.ufrc,
                ufrclcl=self.temporaries.ufrclcl,
                prel=self.temporaries.prel,
                criqc=criqc,
                qcubelow=self.temporaries.qcubelow,
                qlubelow=self.temporaries.qlubelow,
                qiubelow=self.temporaries.qiubelow,
                qcu=self.temporaries.qcu,
                qlu=self.temporaries.qlu,
                qiu=self.temporaries.qiu,
                rcwp=self.temporaries.rcwp,
                rlwp=self.temporaries.rlwp,
                riwp=self.temporaries.riwp,
                cufrc=self.temporaries.cufrc,
                cush_inout=self.temporaries.cush_inout,
            )

            if iteration != iter_cin:
                self._adjust_implicit_CIN_inputs(
                    condensation=self.condensation,
                    qv0=self.temporaries.qv0,
                    qvten=self.temporaries.qvten,
                    dt=dt,
                    ql0=self.temporaries.ql0,
                    qlten=self.temporaries.qlten,
                    qi0=self.temporaries.qi0,
                    qiten=self.temporaries.qiten,
                    s0=self.temporaries.s0,
                    sten=self.temporaries.sten,
                    u0=self.temporaries.u0_in,
                    uten=self.temporaries.uten,
                    v0=self.temporaries.v0_in,
                    vten=self.temporaries.vten,
                    t0=self.temporaries.t0,
                    dotransport=dotransport,
                    tr0_s=self.tr0_s,
                    tr0=self.tr0,
                    trten=self.trten,
                    umf_s=self.temporaries.umf_s,
                    umf_zint=self.temporaries.umf_zint,
                    dcm=self.temporaries.dcm,
                    qrten=self.temporaries.qrten,
                    qsten=self.temporaries.qsten,
                    cush=cush,
                    cufrc=self.temporaries.cufrc,
                    slflx_s=self.temporaries.slflx_s,
                    slflx=self.temporaries.slflx,
                    qtflx_s=self.temporaries.qtflx_s,
                    qtflx=self.temporaries.qtflx,
                    uflx_s=self.temporaries.uflx_s,
                    uflx=self.temporaries.uflx,
                    vflx_s=self.temporaries.vflx_s,
                    vflx=self.temporaries.vflx,
                    qcu=self.temporaries.qcu,
                    qlu=self.temporaries.qlu,
                    qiu=self.temporaries.qiu,
                    fer=self.temporaries.fer,
                    fdr=self.temporaries.fdr,
                    xco=self.temporaries.xco,
                    cin_IJ=self.temporaries.cin_IJ,
                    cinlcl_IJ=self.temporaries.cinlcl_IJ,
                    cbmf=self.temporaries.cbmf,
                    qc=self.temporaries.qc,
                    qlten_det=self.temporaries.qlten_det,
                    qiten_det=self.temporaries.qiten_det,
                    qlten_sink=self.temporaries.qlten_sink,
                    qiten_sink=self.temporaries.qiten_sink,
                    ufrc_s=self.temporaries.ufrc_s,
                    ufrc=self.temporaries.ufrc,
                    qv0_s=self.temporaries.qv0_s,
                    ql0_s=self.temporaries.ql0_s,
                    qi0_s=self.temporaries.qi0_s,
                    s0_s=self.temporaries.s0_s,
                    t0_s=self.temporaries.t0_s,
                    dcm_s=self.temporaries.dcm_s,
                    qvten_s=self.temporaries.qvten_s,
                    qlten_s=self.temporaries.qlten_s,
                    qiten_s=self.temporaries.qiten_s,
                    sten_s=self.temporaries.sten_s,
                    uten_s=self.temporaries.uten_s,
                    vten_s=self.temporaries.vten_s,
                    qrten_s=self.temporaries.qrten_s,
                    qsten_s=self.temporaries.qsten_s,
                    qldet_s=self.temporaries.qldet_s,
                    qidet_s=self.temporaries.qidet_s,
                    qlsub_s=self.temporaries.qlsub_s,
                    qisub_s=self.temporaries.qisub_s,
                    cush_s=self.temporaries.cush_s,
                    cufrc_s=self.temporaries.cufrc_s,
                    fer_s=self.temporaries.fer_s,
                    fdr_s=self.temporaries.fdr_s,
                )

                self._recalc_environmental_variables(
                    condensation=self.condensation,
                    qv0_s=self.temporaries.qv0_s,
                    ql0_s=self.temporaries.ql0_s,
                    qi0_s=self.temporaries.qi0_s,
                    s0_s=self.temporaries.s0_s,
                    t0_s=self.temporaries.t0_s,
                    exnmid0=self.temporaries.exnmid0_in,
                    pmid0=self.temporaries.pmid0_in,
                    dotransport=dotransport,
                    sstr0=self.sstr0,
                    tr0=self.tr0,
                    u0=self.temporaries.u0_in,
                    v0=self.temporaries.v0_in,
                    pifc0=self.temporaries.pifc0_in,
                    ese=self.qsat.ese,
                    esx=self.qsat.esx,
                    umf_out=self.temporaries.umf_out,
                    slflx_out=self.temporaries.slflx_out,
                    qtflx_out=self.temporaries.qtflx_out,
                    uflx_out=self.temporaries.uflx_out,
                    vflx_out=self.temporaries.vflx_out,
                    thvl0bot=self.temporaries.thvl0bot,
                    thv0bot=self.temporaries.thv0bot,
                    thvl0top=self.temporaries.thvl0top,
                    thv0top=self.temporaries.thv0top,
                    thl0=self.temporaries.thl0,
                    qt0=self.temporaries.qt0,
                    thvl0=self.temporaries.thvl0,
                    ssthl0=self.temporaries.ssthl0,
                    ssu0=self.temporaries.ssu0,
                    ssv0=self.temporaries.ssv0,
                    ssqt0=self.temporaries.ssqt0,
                    qv0=self.temporaries.qv0,
                    ql0=self.temporaries.ql0,
                    qi0=self.temporaries.qi0,
                    s0=self.temporaries.s0,
                    t0=self.temporaries.t0,
                    tr0_temp=self.temporaries.tr0_temp,
                    cush_inout=self.temporaries.cush_inout,
                )

        self._update_output_variables(
            del_CIN=self.temporaries.del_CIN,
            umf_zint=self.temporaries.umf_zint,
            zifc0=self.temporaries.zifc0_in,
            kinv=self.temporaries.kinv,
            dcm=self.temporaries.dcm,
            qvten=self.temporaries.qvten,
            qlten=self.temporaries.qlten,
            qiten=self.temporaries.qiten,
            sten=self.temporaries.sten,
            uten=self.temporaries.uten,
            vten=self.temporaries.vten,
            qrten=self.temporaries.qrten,
            qsten=self.temporaries.qsten,
            cufrc=self.temporaries.cufrc,
            cush=cush,
            qlten_det=self.temporaries.qlten_det,
            qiten_det=self.temporaries.qiten_det,
            qlten_sink=self.temporaries.qlten_sink,
            qiten_sink=self.temporaries.qiten_sink,
            rdrop=rdrop,
            qtflx=self.temporaries.qtflx,
            slflx=self.temporaries.slflx,
            uflx=self.temporaries.uflx,
            vflx=self.temporaries.vflx,
            dotransport=dotransport,
            trten=self.trten,
            dt=dt,
            fer=self.temporaries.fer,
            fdr=self.temporaries.fdr,
            kpen=self.temporaries.kpen,
            umf_out=self.temporaries.umf_out,
            dcm_out=self.temporaries.dcm_out,
            qvten_out=self.temporaries.qvten_out,
            qlten_out=self.temporaries.qlten_out,
            qiten_out=self.temporaries.qiten_out,
            sten_out=self.temporaries.sten_out,
            uten_out=self.temporaries.uten_out,
            vten_out=self.temporaries.vten_out,
            qrten_out=self.temporaries.qrten_out,
            qsten_out=self.temporaries.qsten_out,
            cufrc_out=self.temporaries.cufrc_out,
            cush_inout=self.temporaries.cush_inout,
            qldet_out=self.temporaries.qldet_out,
            qidet_out=self.temporaries.qidet_out,
            qlsub_out=self.temporaries.qlsub_out,
            qisub_out=self.temporaries.qisub_out,
            ndrop_out=self.temporaries.ndrop_out,
            nice_out=self.temporaries.nice_out,
            qtflx_out=self.temporaries.qtflx_out,
            slflx_out=self.temporaries.slflx_out,
            uflx_out=self.temporaries.uflx_out,
            vflx_out=self.temporaries.vflx_out,
            tr0_inout=self.tr0_inout,
            fer_out=self.temporaries.fer_out,
            umf_outvar=self.temporaries.umf_outvar,
            dcm_outvar=self.temporaries.dcm_outvar,
            qvten_outvar=self.temporaries.qvten_outvar,
            qlten_outvar=self.temporaries.qlten_outvar,
            qiten_outvar=self.temporaries.qiten_outvar,
            sten_outvar=self.temporaries.sten_outvar,
            uten_outvar=self.temporaries.uten_outvar,
            vten_outvar=self.temporaries.vten_outvar,
            qrten_outvar=self.temporaries.qrten_outvar,
            qsten_outvar=self.temporaries.qsten_outvar,
            cufrc_outvar=self.temporaries.cufrc_outvar,
            cush_inoutvar=self.temporaries.cush_inoutvar,
            qldet_outvar=self.temporaries.qldet_outvar,
            qidet_outvar=self.temporaries.qidet_outvar,
            qlsub_outvar=self.temporaries.qlsub_outvar,
            qisub_outvar=self.temporaries.qisub_outvar,
            qtflx_outvar=self.temporaries.qtflx_outvar,
            slflx_outvar=self.temporaries.slflx_outvar,
            uflx_outvar=self.temporaries.uflx_outvar,
            vflx_outvar=self.temporaries.vflx_outvar,
            fdr_out=self.temporaries.fdr_out,
            fer_outvar=self.temporaries.fer_outvar,
            fdr_outvar=self.temporaries.fdr_outvar,
            tr0_inoutvar=self.tr0_inoutvar,
        )

        self._compute_uwshcu_invert_after(
            # Inputs
            k0=k0,
            umf_outvar=self.temporaries.umf_outvar,
            qtflx_outvar=self.temporaries.qtflx_outvar,
            slflx_outvar=self.temporaries.slflx_outvar,
            uflx_outvar=self.temporaries.uflx_outvar,
            vflx_outvar=self.temporaries.vflx_outvar,
            dcm_outvar=self.temporaries.dcm_outvar,
            qvten_outvar=self.temporaries.qvten_outvar,
            qlten_outvar=self.temporaries.qlten_outvar,
            qiten_outvar=self.temporaries.qiten_outvar,
            sten_outvar=self.temporaries.sten_outvar,
            uten_outvar=self.temporaries.uten_outvar,
            vten_outvar=self.temporaries.vten_outvar,
            qrten_outvar=self.temporaries.qrten_outvar,
            qsten_outvar=self.temporaries.qsten_outvar,
            cufrc_outvar=self.temporaries.cufrc_outvar,
            qldet_outvar=self.temporaries.qldet_outvar,
            qidet_outvar=self.temporaries.qidet_outvar,
            qlsub_outvar=self.temporaries.qlsub_outvar,
            qisub_outvar=self.temporaries.qisub_outvar,
            fer_outvar=self.temporaries.fer_outvar,
            fdr_outvar=self.temporaries.fdr_outvar,
            ndrop_out=self.temporaries.ndrop_out,
            nice_out=self.temporaries.nice_out,
            tr0=self.tr0,
            tr0_inoutvar=self.tr0_inoutvar,
            CNV_Tracers=CNV_Tracers,
            # Outputs
            umf_inv=umf_inv,
            dcm_inv=dcm_inv,
            qtflx_inv=qtflx_inv,
            slflx_inv=slflx_inv,
            uflx_inv=uflx_inv,
            vflx_inv=vflx_inv,
            qvten_inv=qvten_inv,
            qlten_inv=qlten_inv,
            qiten_inv=qiten_inv,
            tten_inv=tten_inv,
            uten_inv=uten_inv,
            vten_inv=vten_inv,
            qrten_inv=qrten_inv,
            qsten_inv=qsten_inv,
            cufrc_inv=cufrc_inv,
            fer_inv=fer_inv,
            fdr_inv=fdr_inv,
            ndrop_inv=ndrop_inv,
            nice_inv=nice_inv,
            qldet_inv=qldet_inv,
            qlsub_inv=qlsub_inv,
            qidet_inv=qidet_inv,
            qisub_inv=qisub_inv,
            dotransport=dotransport,
        )
