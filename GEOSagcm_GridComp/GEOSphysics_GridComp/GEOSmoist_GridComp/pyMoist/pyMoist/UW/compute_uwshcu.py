import copy
import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import (
    computation,
    interval,
    PARALLEL,
    FORWARD,
    BACKWARD,
    exp,
    sqrt,
    log,
    erfc,
    THIS_K,
    f32,
    i64,
    i32,
    f64,
)
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import (
    FloatField,
    Float,
    Int,
    IntFieldIJ,
    FloatFieldIJ,
    BoolField,
    IntField,
)
from ndsl import StencilFactory, QuantityFactory
from pyMoist.field_types import FloatField_NTracers
from pyMoist.saturation.qsat import QSat, QSat_Float, FloatField_Extra_Dim
from pyMoist.saturation.formulation import SaturationFormulation
import pyMoist.constants as constants
from pyMoist.UW.uwshcu_functions import (
    slope,
    conden,
    qsinvert,
    single_cin,
    getbuoy,
    compute_alpha,
    compute_mumin2,
    compute_ppen,
    exnerfn,
    roots,
    zvir,
)


def compute_thermodynamic_variables(
    k0: Int,
    dt: Float,
    ncnst: Int,
    pifc0_in: FloatField,
    zifc0_in: FloatField,
    exnifc0_in: FloatField,
    pmid0_in: FloatField,
    zmid0_in: FloatField,
    exnmid0_in: FloatField,
    dp0_in: FloatField,
    u0_in: FloatField,
    v0_in: FloatField,
    qv0_in: FloatField,
    ql0_in: FloatField,
    qi0_in: FloatField,
    th0_in: FloatField,
    tr0_inout: FloatField_NTracers,
    kpbl_in: IntFieldIJ,
    frland_in: FloatFieldIJ,
    # tke_in: FloatField,
    rkfre: FloatFieldIJ,
    cush_inout: FloatFieldIJ,
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
    fer_out: FloatField,
    fdr_out: FloatField,
    # qldet_out: FloatField,
    # qidet_out: FloatField,
    # qlsub_out: FloatField,
    # qisub_out: FloatField,
    # ndrop_out: FloatField,
    # nice_out: FloatField,
    shfx: FloatFieldIJ,
    evap: FloatFieldIJ,
    # cnvtr: FloatFieldIJ,
    # tpert_out: FloatFieldIJ,
    # qpert_out: FloatFieldIJ,
    qtflx_out: FloatField,
    slflx_out: FloatField,
    uflx_out: FloatField,
    vflx_out: FloatField,
    # cush: FloatFieldIJ,
    # pifc0: FloatField,
    # zifc0: FloatField,
    # exnifc0: FloatField,
    # tke: FloatField,
    # u0: FloatField,
    # v0: FloatField,
    # dp0: FloatField,
    thvl0: FloatField,
    thvl0bot: FloatField,
    thv0bot: FloatField,
    # thl0bot: FloatFieldIJ,
    # thvl0top: FloatField,
    zmid0: FloatField,
    qt0: FloatField,
    t0: FloatField,
    qv0: FloatField,
    pmid0: FloatField,
    # zvir: Float,
    tr0: FloatField_NTracers,
    ssu0: FloatField,
    ssv0: FloatField,
    ssthl0: FloatField,
    ssqt0: FloatField,
    sstr0: FloatField_NTracers,
    thl0: FloatField,
    thl0top: FloatField,
    qt0top: FloatField,
    # qv0_o: FloatField,
    # ql0_o: FloatField,
    # qi0_o: FloatField,
    # t0_o: FloatField,
    # s0_o: FloatField,
    # u0_o: FloatField,
    # v0_o: FloatField,
    # qt0_o: FloatField,
    # thl0_o: FloatField,
    # thvl0_o: FloatField,
    # ssthl0_o: FloatField,
    # ssqt0_o: FloatField,
    # thv0bot_o: FloatField,
    # thv0top_o: FloatField,
    # thvl0bot_o: FloatField,
    # thvl0top_o: FloatField,
    # ssu0_o: FloatField,
    # ssv0_o: FloatField,
    # tr0_o: FloatField,
    # sstr0_o: FloatField,
    thj: FloatField,
    qij: FloatField,
    qlj: FloatField,
    qvj: FloatField,
    qse: FloatField,
    id_check: IntField,
    thv0top: FloatField,
    thvl0top: FloatField,
    tr0_o: FloatField_NTracers,
    sstr0_o: FloatField_NTracers,
    trflx: FloatField_NTracers,
    trten: FloatField_NTracers,
    tru: FloatField_NTracers,
    tru_emf: FloatField_NTracers,
    dotransport: Float,
    k_idx: FloatField,
    ese: FloatField_Extra_Dim,
    esx: FloatField_Extra_Dim,
    id_exit: BoolField,
):
    """
    University of Washington Shallow Convection Scheme

    Described in Park and Bretherton. 2008. J. Climate :

    'The University of Washington shallow convection and
    moist turbulent schemes and their impact on climate
    simulations with the Community Atmosphere Model'
    Coded in CESM by Sungsu Park. Oct.2005. May.2008.
    Coded in GEOS by Nathan Arnold. July 2016.
    For general questions, email sungsup@ucar.edu or
    sungsu@atmos.washington.edu
    For GEOS-specific questions, email nathan.arnold@nasa.gov
    """

    with computation(FORWARD), interval(...):
        ### TEMPORARY TEST SETTINGS
        ixcldice = 1
        ixcldliq = 2
        ixnumliq = 3
        ixnumice = 4

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

        # Initialize more variables
        exit_UWCu = 0.0
        exit_conden = 0.0
        exit_klclk0 = 0.0
        exit_klfck0 = 0.0
        exit_ufrc = 0.0
        exit_wtw = 0.0
        exit_drycore = 0.0
        exit_wu = 0.0
        exit_cufilter = 0.0
        exit_kinv1 = 0.0
        exit_rei = 0.0

        limit_shcu = 0.0
        limit_negcon = 0.0
        limit_ufrc = 0.0
        limit_ppen = 0.0
        limit_emf = 0.0
        limit_cinlcl = 0.0
        limit_cin = 0.0
        limit_cbmf = 0.0
        limit_rei = 0.0

        ind_delcin = 0.0

        # Initialize variable that are calculated from inputs
        # pifc0[0, 0, 1] = pifc0_in  # uncomment later
        # zifc0[0, 0, 1] = zifc0_in  # uncomment later
        pmid0 = pmid0_in
        zmid0 = zmid0_in
        dp0 = dp0_in
        u0 = u0_in
        v0 = v0_in
        qv0 = qv0_in
        ql0 = ql0_in
        qi0 = qi0_in
        # tke = tke_in[0, 0, 1]
        # pblh           = pblh_in
        cush = cush_inout

    # Start Main Calculation
    # Compute basic thermodynamic variables directly from
    # input variables for each column
    with computation(PARALLEL), interval(0, 1):
        # Compute interval environmental variables
        pmid0 = pmid0_in
        pmid0_above = pmid0_in[0, 0, 1]
        u0 = u0_in
        u0_above = u0_in[0, 0, 1]
        v0 = v0_in
        v0_above = v0_in[0, 0, 1]
        qv0 = qv0_in
        qv0_above = qv0_in[0, 0, 1]
        ql0 = ql0_in
        ql0_above = ql0_in[0, 0, 1]
        qi0 = qi0_in
        qi0_above = qi0_in[0, 0, 1]
        qt0 = qv0 + ql0 + qi0
        qt0_above = qv0_above + ql0_above + qi0_above
        exnmid0 = exnmid0_in
        exnmid0_above = exnmid0_in[0, 0, 1]
        t0 = th0_in * exnmid0
        t0_above = th0_in[0, 0, 1] * exnmid0_above
        thl0 = (
            t0
            - ((constants.MAPL_ALHL * ql0) / constants.MAPL_CP)
            - ((constants.MAPL_ALHS * qi0) / constants.MAPL_CP)
        ) / exnmid0
        thl0_above = (
            t0_above
            - ((constants.MAPL_ALHL * ql0_above) / constants.MAPL_CP)
            - ((constants.MAPL_ALHS * qi0_above) / constants.MAPL_CP)
        ) / exnmid0_above

        if dotransport == 1.0:
            n = 0
            # Loop over tracers
            while n < ncnst:
                tr0[0, 0, 0][n] = tr0_inout[0, 0, 0][n]
                n += 1

        # Compute slopes of environmental variables in each layer
        ssthl0 = slope(
            k_idx, thl0, thl0_above, thl0_above, pmid0, pmid0_above, pmid0_above
        )
        ssqt0 = slope(k_idx, qt0, qt0_above, qt0_above, pmid0, pmid0_above, pmid0_above)
        ssu0 = slope(k_idx, u0, u0_above, u0_above, pmid0, pmid0_above, pmid0_above)
        ssv0 = slope(k_idx, v0, v0_above, v0_above, pmid0, pmid0_above, pmid0_above)

        if dotransport == 1.0:
            n = 0
            while n < ncnst:
                sstr0[0, 0, 0][n] = slope(
                    k_idx,
                    tr0[0, 0, 0][n],
                    tr0[0, 0, 1][n],
                    tr0[0, 0, 1][n],
                    pmid0,
                    pmid0_above,
                    pmid0_above,
                )
                n += 1

    with computation(PARALLEL), interval(1, -1):
        # Compute interval environmental variables
        pmid0 = pmid0_in
        pmid0_above = pmid0_in[0, 0, 1]
        pmid0_below = pmid0_in[0, 0, -1]
        u0 = u0_in
        u0_above = u0_in[0, 0, 1]
        u0_below = u0_in[0, 0, -1]
        v0 = v0_in
        v0_above = v0_in[0, 0, 1]
        v0_below = v0_in[0, 0, -1]
        qv0 = qv0_in
        qv0_above = qv0_in[0, 0, 1]
        qv0_below = qv0_in[0, 0, -1]
        ql0 = ql0_in
        ql0_above = ql0_in[0, 0, 1]
        ql0_below = ql0_in[0, 0, -1]
        qi0 = qi0_in
        qi0_above = qi0_in[0, 0, 1]
        qi0_below = qi0_in[0, 0, -1]

        if dotransport == 1.0:
            n = 0
            # Loop over tracers
            while n < ncnst:
                tr0[0, 0, 0][n] = tr0_inout[0, 0, 0][n]
                n += 1

        exnmid0 = exnmid0_in
        exnmid0_above = exnmid0_in[0, 0, 1]
        exnmid0_below = exnmid0_in[0, 0, -1]
        t0 = th0_in * exnmid0
        t0_above = th0_in[0, 0, 1] * exnmid0_above
        t0_below = th0_in[0, 0, -1] * exnmid0_below
        qt0: f32 = f32(qv0) + f32(ql0) + f32(qi0)
        qt0_above = qv0_above + ql0_above + qi0_above
        qt0_below = qv0_below + ql0_below + qi0_below
        thl0 = (
            t0
            - ((constants.MAPL_ALHL * ql0) / constants.MAPL_CP)
            - ((constants.MAPL_ALHS * qi0) / constants.MAPL_CP)
        ) / exnmid0
        thl0_above = (
            t0_above
            - ((constants.MAPL_ALHL * ql0_above) / constants.MAPL_CP)
            - ((constants.MAPL_ALHS * qi0_above) / constants.MAPL_CP)
        ) / exnmid0_above
        thl0_below = (
            t0_below
            - ((constants.MAPL_ALHL * ql0_below) / constants.MAPL_CP)
            - ((constants.MAPL_ALHS * qi0_below) / constants.MAPL_CP)
        ) / exnmid0_below

        # Compute slopes of environmental variables in each layer
        ssthl0 = slope(
            k_idx, thl0, thl0_above, thl0_below, pmid0, pmid0_above, pmid0_below
        )
        ssqt0 = slope(k_idx, qt0, qt0_above, qt0_below, pmid0, pmid0_above, pmid0_below)
        ssu0 = slope(k_idx, u0, u0_above, u0_below, pmid0, pmid0_above, pmid0_below)
        ssv0 = slope(k_idx, v0, v0_above, v0_below, pmid0, pmid0_above, pmid0_below)

        if dotransport == 1.0:
            n = 0
            while n < ncnst:
                sstr0[0, 0, 0][n] = slope(
                    k_idx,
                    tr0[0, 0, 0][n],
                    tr0[0, 0, 1][n],
                    tr0[0, 0, -1][n],
                    pmid0,
                    pmid0_above,
                    pmid0_below,
                )
                n += 1

    with computation(PARALLEL), interval(-1, None):
        # Compute interval environmental variables
        pmid0 = pmid0_in[0, 0, -1]
        pmid0_above = pmid0_in
        pmid0_below = pmid0_in[0, 0, -2]
        u0 = u0_in[0, 0, -1]
        u0_above = u0_in
        u0_below = u0_in[0, 0, -2]
        v0 = v0_in[0, 0, -1]
        v0_above = v0_in
        v0_below = v0_in[0, 0, -2]
        qv0 = qv0_in[0, 0, -1]
        qv0_above = qv0_in
        qv0_below = qv0_in[0, 0, -2]
        ql0 = ql0_in[0, 0, -1]
        ql0_above = ql0_in
        ql0_below = ql0_in[0, 0, -2]
        qi0 = qi0_in[0, 0, -1]
        qi0_above = qi0_in
        qi0_below = qi0_in[0, 0, -2]

        exnmid0 = exnmid0_in[0, 0, -1]
        exnmid0_above = exnmid0_in
        exnmid0_below = exnmid0_in[0, 0, -2]
        t0 = th0_in[0, 0, -1] * exnmid0
        t0_above = th0_in * exnmid0_above
        t0_below = th0_in[0, 0, -2] * exnmid0_below
        qt0: f32 = f32(qv0) + f32(ql0) + f32(qi0)
        qt0_above = qv0_above + ql0_above + qi0_above
        qt0_below = qv0_below + ql0_below + qi0_below
        thl0 = (
            t0
            - ((constants.MAPL_ALHL * ql0) / constants.MAPL_CP)
            - ((constants.MAPL_ALHS * qi0) / constants.MAPL_CP)
        ) / exnmid0
        thl0_above = (
            t0_above
            - ((constants.MAPL_ALHL * ql0_above) / constants.MAPL_CP)
            - ((constants.MAPL_ALHS * qi0_above) / constants.MAPL_CP)
        ) / exnmid0_above
        thl0_below = (
            t0_below
            - ((constants.MAPL_ALHL * ql0_below) / constants.MAPL_CP)
            - ((constants.MAPL_ALHS * qi0_below) / constants.MAPL_CP)
        ) / exnmid0_below

        if dotransport == 1.0:
            n = 0
            # Loop over tracers
            while n < ncnst:
                tr0[0, 0, 0][n] = tr0_inout[0, 0, 0][n]
                n += 1

        # Compute slopes of environmental variables in each layer
        ssthl0 = slope(
            k_idx, thl0, thl0_above, thl0_below, pmid0, pmid0_above, pmid0_below
        )
        ssqt0 = slope(k_idx, qt0, qt0_above, qt0_below, pmid0, pmid0_above, pmid0_below)
        ssu0 = slope(k_idx, u0, u0_above, u0_below, pmid0, pmid0_above, pmid0_below)
        ssv0 = slope(k_idx, v0, v0_above, v0_below, pmid0, pmid0_above, pmid0_below)

        # Calculate slope for each tracer by hand
        if dotransport == 1.0:
            n = 0
            while n < ncnst:
                sstr0[0, 0, 0][n] = slope(
                    k_idx,
                    tr0[0, 0, -1][n],
                    tr0[0, 0, 0][n],
                    tr0[0, 0, -2][n],
                    pmid0,
                    pmid0_above,
                    pmid0_below,
                )
                n += 1

    with computation(PARALLEL), interval(...):
        # Compute interval environmental variables
        pmid0 = pmid0_in
        exnmid0 = exnmid0_in
        qv0 = qv0_in
        ql0 = ql0_in
        qi0 = qi0_in
        qt0: f32 = f32(qv0) + f32(ql0) + f32(qi0)
        t0 = th0_in * exnmid0
        s0 = constants.MAPL_GRAV * zmid0 + constants.MAPL_CP * t0
        thl0 = (
            t0
            - constants.MAPL_LATENT_HEAT_VAPORIZATION * ql0 / constants.MAPL_CP
            - constants.MAPL_LATENT_HEAT_SUBLIMATION * qi0 / constants.MAPL_CP
        ) / exnmid0
        thvl0 = (1.0 + zvir * qt0) * thl0
        pifc0 = pifc0_in

        # Compute thv0 and thvl0 at top/bottom interfaces in each layer
        thl0bot = thl0 + ssthl0 * (pifc0 - pmid0)
        qt0bot = qt0 + ssqt0 * (pifc0 - pmid0)

        thj, qvj, qlj, qij, qse, id_check = conden(pifc0, thl0bot, qt0bot, ese, esx)

    with computation(FORWARD), interval(...):
        if id_check == 1:
            id_exit = True
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
        if id_exit == False:
            thv0bot = thj * (1.0 + zvir * qvj - qlj - qij)
            thvl0bot = thl0bot * (1.0 + zvir * qt0bot)

            thl0top = thl0 + ssthl0 * (pifc0_in[0, 0, 1] - pmid0)
            qt0top = qt0 + ssqt0 * (pifc0_in[0, 0, 1] - pmid0)

    with computation(PARALLEL), interval(0, -1):
        if id_exit == False:
            thj, qvj, qlj, qij, qse, id_check = conden(
                pifc0_in[0, 0, 1], thl0top, qt0top, ese, esx
            )
    with computation(FORWARD), interval(0, -1):
        if id_exit == False:
            if id_check == 1:
                id_exit = True
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
        if id_exit == False:
            thv0top = thj * (1.0 + zvir * qvj - qlj - qij)
            thvl0top = thl0top * (1.0 + zvir * qt0top)

    with computation(PARALLEL), interval(-1, None):
        if id_exit == False:
            thv0top = thv0bot
            thvl0top = thvl0bot

    with computation(PARALLEL), interval(...):
        if id_exit == False:
            # Save input and related environmental thermodynamic variables
            # for use at "iter_cin=2" when "del_CIN >= 0"
            qv0_o = qv0
            ql0_o = ql0
            qi0_o = qi0
            t0_o = t0
            s0_o = s0
            u0_o = u0
            v0_o = v0
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

            if dotransport == 1.0:
                n = 0
                while n < ncnst:
                    tr0_o[0, 0, 0][n] = tr0[0, 0, 0][n]
                    sstr0_o[0, 0, 0][n] = sstr0[0, 0, 0][n]
                    n += 1

    with computation(FORWARD), interval(...):
        if id_exit == False:
            # Initialize output variables at each grid point
            # umf[0, 0, 1] = 0.0
            dcm = 0.0
            # emf[0, 0, 1] = 0.0
            # slflx[0, 0, 1] = 0.0
            # qtflx[0, 0, 1] = 0.0
            # uflx[0, 0, 1] = 0.0
            # vflx[0, 0, 1] = 0.0
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
            # qldet       = 0.0
            # qidet       = 0.0
            qc_l = 0.0
            qc_i = 0.0
            cnt = f32(k0)
            cnb = 0.0
            qtten = 0.0
            slten = 0.0
            ufrc = 0.0

            # thlu[0, 0, 1] = constants.MAPL_UNDEF
            # qtu[0, 0, 1] = constants.MAPL_UNDEF
            # uu[0, 0, 1] = constants.MAPL_UNDEF
            # vu[0, 0, 1] = constants.MAPL_UNDEF
            # wu[0, 0, 1] = constants.MAPL_UNDEF
            # thvu[0, 0, 1] = constants.MAPL_UNDEF
            # thlu_emf[0, 0, 1] = constants.MAPL_UNDEF
            # qtu_emf[0, 0, 1] = constants.MAPL_UNDEF
            # uu_emf[0, 0, 1] = constants.MAPL_UNDEF
            # vu_emf[0, 0, 1] = constants.MAPL_UNDEF

            ufrcinvbase = 0.0
            ufrclcl = 0.0
            winvbase = 0.0
            wlcl = 0.0
            emfkbup = 0.0
            cbmflimit = 0.0

            # uemf[0, 0, 1] = 0.0
            comsub = 0.0
            qlten_sink = 0.0
            qiten_sink = 0.0
            qlten_det = 0.0
            qiten_det = 0.0
            # nlten_sink   = 0.0
            # iten_sink    = 0.0

            if dotransport == 1.0:
                m = 0
                while m < ncnst:
                    trflx[0, 0, 1][m] = 0.0
                    trten[0, 0, 0][m] = 0.0
                    tru[0, 0, 1][m] = 0.0
                    tru_emf[0, 0, 1][m] = 0.0
                    m += 1


def implicit_cin_closure(
    # Inputs:
    iteration: Int,
    windsrcavg: Int,
    qtsrchgt: Float,
    thlsrc_fac: Float,
    qtsrc_fac: Float,
    # rbuoy: Float,
    # epsvarw: Float,
    # use_CINcin: Int,
    # mumin1: Float,
    # rmaxfrac: Float,
    # PGFc: Float,
    cush: FloatFieldIJ,
    kpbl_in: IntFieldIJ,
    k0: Int,
    # dt: Float,
    pifc0: FloatField,
    # zifc0: FloatField,
    # exnifc0: FloatField,
    tke_in: FloatField,
    # rkfre: FloatFieldIJ,
    u0: FloatField,
    v0: FloatField,
    thvl0: FloatField,
    thvl0bot: FloatField,
    thvl0top: FloatField,
    zmid0: FloatField,
    qt0: FloatField,
    t0: FloatField,
    qv0: FloatField,
    shfx: FloatFieldIJ,
    evap: FloatFieldIJ,
    pmid0: FloatField,
    dotransport: Float,
    ncnst: Int,
    zvir: Float,
    tr0: FloatField_NTracers,
    ssu0: FloatField,
    ssv0: FloatField,
    ssthl0: FloatField,
    ssqt0: FloatField,
    # sstr0: FloatField_NTracers,
    thl0: FloatField,
    thv0bot: FloatField,
    thv0top: FloatField,
    # qv0_o: FloatField,
    # ql0_o: FloatField,
    # qi0_o: FloatField,
    # t0_o: FloatField,
    # s0_o: FloatField,
    # u0_o: FloatField,
    # v0_o: FloatField,
    # qt0_o: FloatField,
    # thl0_o: FloatField,
    # thvl0_o: FloatField,
    # ssthl0_o: FloatField,
    # ssqt0_o: FloatField,
    # thv0bot_o: FloatField,
    # thv0top_o: FloatField,
    # thvl0bot_o: FloatField,
    # thvl0top_o: FloatField,
    # ssu0_o: FloatField,
    # ssv0_o: FloatField,
    # tr0_o: FloatField,
    # sstr0_o: FloatField,
    # dp0: FloatField,
    ese: FloatField_Extra_Dim,
    esx: FloatField_Extra_Dim,
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
    # qldet_out: FloatField,
    # qidet_out: FloatField,
    qtflx_out: FloatField,
    slflx_out: FloatField,
    uflx_out: FloatField,
    vflx_out: FloatField,
    fer_out: FloatField,
    fdr_out: FloatField,
    id_exit: BoolField,
    # Outputs for testing:
    kinv: IntField,
    thvlavg: FloatField,
    tkeavg: FloatField,
    uavg: FloatField,
    vavg: FloatField,
    thvlmin: FloatField,
    qtavg: FloatField,
    dpi: FloatFieldIJ,
    thlsrc: FloatField,
    usrc: FloatField,
    vsrc: FloatField,
    trsrc: FloatField_NTracers,
    plcl: FloatField,
    klcl: IntField,
    thl0lcl: FloatField,
    qt0lcl: FloatField,
    thv0lcl: FloatField,
    plfc: FloatFieldIJ,
    cin: FloatField,
    thvubot: FloatField,
    thvutop: FloatField,
    thvlsrc: FloatField,
    thj2: FloatField,
    qvj2: FloatField,
    qlj2: FloatField,
    qij2: FloatField,
    qse2: FloatField,
    qtsrc: FloatField,
    test_var1: FloatField,
):
    """
    This 'iteration' loop is for implicit CIN closure

    It is important to note that this iterative cin loop is located at the outer
    shell of the code. Thus, source air properties can also be changed during the
    iterative cin calculation, because cumulus convection induces non-zero fluxes
    even at interfaces below PBL top height through 'fluxbelowinv' stencil.
    """
    with computation(FORWARD), interval(...):
        if id_exit == False:
            cush = -1.0
            qtavg = 0.0

            """
            Cumulus scale height
            In contrast to the premitive code, cumulus scale height is iteratively
            calculated at each time step, and at each iterative cin step.
            It is not clear whether I should locate below two lines within or out
            of the iterative cin loop.
            """
            tscaleh = cush

            """
            Find PBL top height interface index, 'kinv-1' where 'kinv' is the layer
            index with PBLH in it. When PBLH is exactly at interface, 'kinv' is the
            layer index having PBLH as a lower interface.
            In the previous code, I set the lower limit of 'kinv' by 2  in order to
            be consistent with the other parts of the code. However in the modified
            code, I allowed 'kinv' to be 1 & if 'kinv = 1', I just exit the program
            without performing cumulus convection. This new approach seems to be
            more reasonable: if PBL height is within 'kinv=1' layer, surface is STL
            interface (bflxs <= 0) and interface just above the surface should be
            either non-turbulent (Ri>0.19) or stably turbulent (0<=Ri<0.19 but this
            interface is identified as a base external interface of upperlying CL.
            Thus, when 'kinv=1', PBL scheme guarantees 'bflxs <= 0'.  For this case
            it is reasonable to assume that cumulus convection does not happen.
            When these is SBCL, PBL height from the PBL scheme is likely to be very
            close at 'kinv-1' interface, but not exactly, since 'zi' information is
            changed between two model time steps. In order to ensure correct identi
            fication of 'kinv' for general case including SBCL, I imposed an offset
            of 5 [m] in the below 'kinv' finding block.
            """
            # Invert kpbl index
            if kpbl_in > i64(k0 / 2):
                kinv = k0 - kpbl_in + 1
            else:
                kinv = 5

            if kinv <= i64(1):
                id_exit = True

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
                # print('------- UW ShCu: Exit, kinv<=1')

            if id_exit == False:
                if kinv >= i64(k0 / 4):
                    id_exit = True

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
                    # print('------- UW ShCu: Exit, kinv>k0/4')

                # kinv = min(kinv,7)
                # From here, it must be '7 >= kinv >= 2'.

    with computation(FORWARD), interval(0, 1):
        if id_exit == False:
            thvlmin = 1000.0
            thvlmin = min(
                thvlmin,
                min(thvl0bot, thvl0top),
            )
    with computation(FORWARD), interval(1, None):
        if id_exit == False:
            if i64(THIS_K) <= i64(kinv - 2):
                thvlmin = min(
                    thvlmin[0, 0, -1],
                    min(thvl0bot, thvl0top),
                )

    with computation(FORWARD), interval(...):
        if id_exit == False:

            """
            Find PBL averaged tke ('tkeavg') and minimum 'thvl' ('thvlmin') in the PBL
            In the current code, 'tkeavg' is obtained by averaging all interfacial TKE
            within the PBL. However, in order to be conceptually consistent with   PBL
            scheme, 'tkeavg' should be calculated by considering surface buoyancy flux.
            If surface buoyancy flux is positive ( bflxs >0 ), surface interfacial TKE
            should be included in calculating 'tkeavg', while if bflxs <= 0,   surface
            interfacial TKE should not be included in calculating 'tkeavg'.   I should
            modify the code when 'bflxs' is available as an input of cumulus scheme.
            'thvlmin' is a minimum 'thvl' within PBL obtained by comparing top &  base
            interface values of 'thvl' in each layers within the PBL.
            """

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

            # # weighted average over lowest 20mb
            # # dpsum = 0.
            # # if k <= kinv:
            # #   dpi = max(0.,(2e3+pmid0-pifc0.at(K=0))/2e3)
            # #   qtavg  = qtavg  + dpi*qt0
            # #   dpsum = dpsum + dpi
            # # qtavg   = qtavg/dpsum

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

            """
            Find characteristics of cumulus source air: qtsrc,thlsrc,usrc,vsrc
            Note that 'thlsrc' was concocted using 'thvlsrc' and 'qtsrc'.
            'qtsrc' is defined as the lowest layer mid-point value;   'thlsrc'
            is from 'qtsrc' and 'thvlmin=thvlsrc'; 'usrc' & 'vsrc' are defined
            as the values just below the PBL top interface.
            """
            if windsrcavg == 1:
                zrho = pifc0.at(K=0) / (
                    287.04 * (t0.at(K=0) * (1.0 + 0.608 * qv0.at(K=0)))
                )
                buoyflx = (
                    -shfx / constants.MAPL_CP - 0.608 * t0.at(K=0) * evap
                ) / zrho  # K m s-1
                # delzg = (zifc0.at(K=1)-zifc0.at(K=0))*constants.MAPL_GRAV
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
                    # qtsrc   = qt0.at(K=1) + qpert_out
                    # thvlsrc = thvlavg + tpert_out*(1.0+zvir*qtsrc) #/exnmid0.at(K=1)
                    thvlsrc = thvlmin + tpert_out * (1.0 + zvir * qtsrc)  # /exnmid0(1)
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

            if dotransport == 1.0:
                n = 0
                while n < ncnst:
                    trsrc[0, 0, 0][n] = tr0.at(K=0, ddim=[n])
                    n += 1

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
            """
            if pifc0.at(K=0) < 70000 or pifc0.at(K=0) > 115000.0:
                id_exit = True
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

            if id_exit == False:
                if qtsrc > 0.1 or qtsrc < 1e-8:
                    id_exit = True
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

            if id_exit == False:
                if thlsrc > 400.0 or thlsrc < 100.0:
                    id_exit = True
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

            if id_exit == False:
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

                klcl = max(1, klcl)

                if plcl < 60000.0:
                    id_exit = True
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

                if id_exit == False:
                    """
                    Calculate environmental virtual potential temperature at LCL,
                    'thv0lcl' which is solely used in the 'cin' calculation. Note
                    that 'thv0lcl' is calculated first by calculating  'thl0lcl'
                    and 'qt0lcl' at the LCL, and performing 'conden' afterward,
                    in fully consistent with the other parts of the code.
                    """
                    thl0lcl = thl0.at(K=klcl - 1) + ssthl0.at(K=klcl - 1) * (
                        plcl - pmid0.at(K=klcl - 1)
                    )
                    qt0lcl = qt0.at(K=klcl - 1) + ssqt0.at(K=klcl - 1) * (
                        plcl - pmid0.at(K=klcl - 1)
                    )
                    thj2, qvj2, qlj2, qij2, qse2, id_check = conden(
                        plcl, thl0lcl, qt0lcl, ese, esx
                    )

                    if id_check == 1:
                        id_exit = True
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

    with computation(FORWARD), interval(1, None):

        if id_exit == False:
            thv0lcl = thj2 * (1.0 + zvir * qvj2 - qlj2 - qij2)

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

                In the below code, 'klfc' is the layer index containing 'LFC' similarto
                'kinv' and 'klcl'.
            """

            cin = 0.0
            cinlcl = 0.0
            plfc = 0.0
            klfc = k0

            """
                #Case 1. LCL height is higher than PBL interface ( 'pLCL <=ps0(kinv-1)' )
            """

            # thvubot = thvlsrc

            stop35 = False
            klev = kinv
            if (klcl) >= kinv and stop35 == False and id_exit == False:
                while klev < k0 and stop35 == False and id_exit == False:
                    if klev < (klcl) and stop35 == False and id_exit == False:
                        if id_exit == False:
                            if THIS_K == klev - 1:
                                thvubot = thvlsrc
                                thvutop = thvlsrc

                                cin = cin[0, 0, -1] + single_cin(
                                    pifc0.at(K=klev - 1),
                                    thv0bot.at(K=klev - 1),
                                    pifc0.at(K=klev),
                                    thv0top.at(K=klev - 1),
                                    thvubot,
                                    thvutop,
                                )

                    elif klev == (klcl) and stop35 == False and id_exit == False:
                        if id_exit == False:
                            if THIS_K == klev - 1:
                                # ----- Bottom to LCL
                                thvubot = thvlsrc
                                thvutop = thvlsrc
                                cin = cin[0, 0, -1] + single_cin(
                                    pifc0.at(K=klev - 1),
                                    thv0bot.at(K=klev - 1),
                                    plcl,
                                    thv0lcl,
                                    thvubot,
                                    thvutop,
                                )
                                cinlcl = max(cin, 0.0)
                                cin = cinlcl

                                # ----- LCL to Top
                                thvubot = thvlsrc

                                thj3, qvj3, qlj3, qij3, qse3, id_check = conden(
                                    pifc0.at(K=klev), thlsrc, qtsrc, ese, esx
                                )

                                if id_check == 1:
                                    id_exit = True
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

                                if id_exit == False:
                                    thvutop = thj3 * (1.0 + zvir * qvj3 - qlj3 - qij3)

                                    plfc, cin = getbuoy(
                                        plcl,
                                        thv0lcl,
                                        pifc0.at(K=klev),
                                        thv0top.at(K=klev - 1),
                                        thvubot,
                                        thvutop,
                                        cin,
                                        plfc,
                                    )

                                    if plfc > 0.0:
                                        klfc = klev
                                        cin = 0.0
                                        cinlcl = 0.0
                                        plfc = 0.0
                                        thvubot = 0.0
                                        thvutop = 0.0
                                        thj3 = 0.0
                                        qvj3 = 0.0
                                        qij3 = 0.0
                                        qlj3 = 0.0
                                        qse3 = 0.0
                                        id_check = 0.0
                                        stop35 = True

                    else:
                        if id_exit == False and stop35 == False:
                            if THIS_K == klev - 1:
                                thvubot = thvutop.at(K=klcl - 1)
                                thj3, qvj3, qlj3, qij3, qse3, id_check = conden(
                                    pifc0.at(K=klev),
                                    thlsrc,
                                    qtsrc,
                                    ese,
                                    esx,
                                )

                    klev += 1

                    #             if id_check == 1:
                    #                 id_exit = True
                    #                 umf_out[0, 0, 1] = 0.0
                    #                 dcm_out = 0.0
                    #                 qvten_out = 0.0
                    #                 qlten_out = 0.0
                    #                 qiten_out = 0.0
                    #                 sten_out = 0.0
                    #                 uten_out = 0.0
                    #                 vten_out = 0.0
                    #                 qrten_out = 0.0
                    #                 qsten_out = 0.0
                    #                 cufrc_out = 0.0
                    #                 cush_inout = -1.0
                    #                 qldet_out = 0.0
                    #                 qidet_out = 0.0
                    #                 qtflx_out[0, 0, 1] = 0.0
                    #                 slflx_out[0, 0, 1] = 0.0
                    #                 uflx_out[0, 0, 1] = 0.0
                    #                 vflx_out[0, 0, 1] = 0.0
                    #                 fer_out = constants.MAPL_UNDEF
                    #                 fdr_out = constants.MAPL_UNDEF

                    #             if id_exit == False:
                    #                 thvutop = thj3 * (1.0 + zvir * qvj3 - qlj3 - qij3)

                    #             #     plfc, cin = getbuoy(
                    #             #         pifc0.at(K=klev - 1),
                    #             #         thv0bot.at(K=klev - 1),
                    #             #         pifc0.at(K=klev),
                    #             #         thv0top.at(K=klev - 1),
                    #             #         thvubot,
                    #             #         thvutop,
                    #             #         cin,
                    #             #         plfc,
                    #             #     )

                    #             #     if plfc > 0.0:
                    #             #         klfc = klev
                    #             #         # cin = 0.0
                    #             #         # plfc = 0.0
                    #             #         stop35 = True
                    # klev += 1

    # else:
    #     """
    #     #Case 2. LCL height is lower than PBL interface ( 'pLCL > ps0(kinv-1)')
    #     """
    #     if id_exit == False and stop35 == False:
    #         cinlcl = 0.0
    #         while klev <= (k0 - 1) and stop35 == False and id_exit == False:
    #             kbelow = klev - 1
    #             thj2, qvj2, qlj2, qij2, qse2, id_check = conden(
    #                 pifc0.at(K=klev - 1), thlsrc, qtsrc, ese, esx
    #             )
    #             if id_check == 1:
    #                 id_exit = True
    #                 umf_out[0, 0, 1] = 0.0
    #                 dcm_out = 0.0
    #                 qvten_out = 0.0
    #                 qlten_out = 0.0
    #                 qiten_out = 0.0
    #                 sten_out = 0.0
    #                 uten_out = 0.0
    #                 vten_out = 0.0
    #                 qrten_out = 0.0
    #                 qsten_out = 0.0
    #                 cufrc_out = 0.0
    #                 cush_inout = -1.0
    #                 qldet_out = 0.0
    #                 qidet_out = 0.0
    #                 qtflx_out[0, 0, 1] = 0.0
    #                 slflx_out[0, 0, 1] = 0.0
    #                 uflx_out[0, 0, 1] = 0.0
    #                 vflx_out[0, 0, 1] = 0.0
    #                 fer_out = constants.MAPL_UNDEF
    #                 fdr_out = constants.MAPL_UNDEF

    #                 if id_exit == False:
    #                     thvubot = thj2 * (1.0 + zvir * qvj2 - qlj2 - qij2)
    #                     thj2, qvj2, qlj2, qij2, qse2, id_check = conden(
    #                         pifc0.at(K=klev),
    #                         thlsrc,
    #                         qtsrc,
    #                         ese,
    #                         esx,
    #                     )

    #                     if id_check == 1:
    #                         id_exit = True
    #                         umf_out[0, 0, 1] = 0.0
    #                         dcm_out = 0.0
    #                         qvten_out = 0.0
    #                         qlten_out = 0.0
    #                         qiten_out = 0.0
    #                         sten_out = 0.0
    #                         uten_out = 0.0
    #                         vten_out = 0.0
    #                         qrten_out = 0.0
    #                         qsten_out = 0.0
    #                         cufrc_out = 0.0
    #                         cush_inout = -1.0
    #                         qldet_out = 0.0
    #                         qidet_out = 0.0
    #                         qtflx_out[0, 0, 1] = 0.0
    #                         slflx_out[0, 0, 1] = 0.0
    #                         uflx_out[0, 0, 1] = 0.0
    #                         vflx_out[0, 0, 1] = 0.0
    #                         fer_out = constants.MAPL_UNDEF
    #                         fdr_out = constants.MAPL_UNDEF

    #                         if id_exit == False:
    #                             thvutop = thj2 * (
    #                                 1.0 + zvir * qvj2 - qlj2 - qij2
    #                             )

    #                             plfc, cin = getbuoy(
    #                                 pifc0.at(K=kbelow),
    #                                 thv0bot.at(K=klev),
    #                                 pifc0.at(K=klev),
    #                                 thv0top.at(K=klev),
    #                                 thvubot,
    #                                 thvutop,
    #                                 plfc,
    #                                 cin,
    #                             )

    #                             if plfc > 0.0:
    #                                 klfc = klev
    #                                 stop35 = True

    #             klev += 1  # End of CIN case selection

    # cin = max(0.0, cin)
    # # cin = max(cin,0.04*(lts-18.))   # kludge to reduce UW in StCu regions


'''
        if klfc >= k0:
            klfc = k0
            id_exit = True
            # go to 333

        """
        In order to calculate implicit 'cin' (or 'cinlcl'), save the initially 
        calculated 'cin' and 'cinlcl', and other related variables. These will 
        be restored after calculating implicit CIN. 
        """
        if iteration == 1:
            cin_i = cin
            cinlcl_i = cinlcl
            ke = rbuoy / (rkfre * tkeavg + epsvarw)
            kinv_o = kinv
            klcl_o = klcl
            klfc_o = klfc
            plcl_o = plcl
            plfc_o = plfc
            tkeavg_o = tkeavg
            thvlmin_o = thvlmin
            qtsrc_o = qtsrc
            thvlsrc_o = thvlsrc
            thlsrc_o = thlsrc
            usrc_o = usrc
            vsrc_o = vsrc
            thv0lcl_o = thv0lcl

            if dotransport == 1.0:
                n = 0
                while n < ncnst:
                    trsrc_o[0, 0, 0][n] = trsrc[0, 0, 0][n]
                    n += 1

    with computation(FORWARD), interval(...):
        """
        Modification : If I impose w = max(0.1, w) up to the top interface of
                     klfc, I should only use cinlfc.  That is, if I want to
                     use cinlcl, I should not impose w = max(0.1, w).
                     Using cinlcl is equivalent to treating only 'saturated'
                     moist convection. Note that in this sense, I should keep
                     the functionality of both cinlfc and cinlcl.
                     However, the treatment of penetrative entrainment level becomes
                     ambiguous if I choose 'cinlcl'. Thus, the best option is to use
                     'cinlfc'.

        Calculate implicit 'cin' by averaging initial and final cins.    Note that
        implicit CIN is adopted only when cumulus convection stabilized the system,
        i.e., only when 'del_CIN >0'. If 'del_CIN<=0', just use explicit CIN. Note
        also that since 'cinlcl' is set to zero whenever LCL is below the PBL top,
        (see above CIN calculation part), the use of 'implicit CIN=cinlcl'  is not
        good. Thus, when using implicit CIN, always try to only use 'implicit CIN=
        cin', not 'implicit CIN=cinlcl'. However, both 'CIN=cin' and 'CIN=cinlcl'
        are good when using explicit CIN.
        """
        if iteration != 1:
            cin_f = cin
            cinlcl_f = cinlcl
            if use_CINcin == 1:
                del_CIN = cin_f - cin_i
            else:
                del_CIN = cinlcl_f - cinlcl_i

            if del_CIN > 0.0:
                """
                Calculate implicit 'cin' and 'cinlcl'. Note that when we chose
                to use 'implicit CIN = cin', choose 'cinlcl = cinlcl_i' below:
                because iterative CIN only aims to obtain implicit CIN,  once
                we obtained 'implicit CIN=cin', it is good to use the original
                profiles information for all the other variables after that.
                Note 'cinlcl' will be explicitly used in calculating  'wlcl' &
                'ufrclcl' after calculating 'winv' & 'ufrcinv'  at the PBL top
                interface later, after calculating 'cbmf'.
                """
                alpha = compute_alpha(del_CIN, ke)
                cin = cin_i + alpha * del_CIN
                if use_CINcin == 1:
                    cinlcl = cinlcl_i
                # else:
                # cinlcl = cinlcl_i + alpha * del_cinlcl

                """
                Restore the original values from the previous 'iter_cin' step (1) 
                to compute correct tendencies for (n+1) time step by implicit CIN
                """

                kinv = kinv_o
                klcl = klcl_o
                klfc = klfc_o
                plcl = plcl_o
                plfc = plfc_o
                tkeavg = tkeavg_o
                thvlmin = thvlmin_o
                qtsrc = qtsrc_o
                thvlsrc = thvlsrc_o
                thlsrc = thlsrc_o
                usrc = usrc_o
                vsrc = vsrc_o
                thv0lcl = thv0lcl_o

                if dotransport == 1.0:
                    n = 0
                    while n < ncnst:
                        trsrc[0, 0, 0][n] = trsrc_o[0, 0, 0][n]
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

                if dotransport == 1.0:
                    n = 0
                    while n < ncnst:
                        tr0[0, 0, 0][n] = tr0_o[0, 0, 0][n]
                        sstr0[0, 0, 0][n] = sstr0_o[0, 0, 0][n]
                        n += 1

                """
                Initialize all fluxes, tendencies, and other variables  
                in association with cumulus convection.   
                """

                umf[0, 0, 1] = 0.0
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

                if dotransport == 1.0:
                    n = 0
                    while n < ncnst:
                        trflx[0, 0, 1][n] = 0.0
                        trten[0, 0, 0][n] = 0.0
                        tru[0, 0, 1][n] = 0.0
                        tru_emf[0, 0, 1][n] = 0.0
                        n += 1

                # Below are diagnostic output variables for detailed
                # analysis of cumulus scheme.
                ufrcinvbase = 0.0
                ufrclcl = 0.0
                winvbase = 0.0
                wlcl = 0.0
                emfkbup = 0.0
                cbmflimit = 0.0

            else:  # When 'del_CIN < 0', use explicit CIN instead of implicit CIN.

                # Identifier showing whether explicit or implicit CIN is used
                ind_delcin = 1.0

                """
                Restore original output values of "iter_cin = 1" and exit 
                """

                umf_out[0, 0, 1] = umf_s[0, 0, 1]
                if THIS_K >= 0 and THIS_K <= (kinv - 1):
                    kbelow = kinv - 1
                    umf_out = umf_s.at(K=kbelow) * zifc0 / zifc0.at(K=kbelow)

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

                # Below are diagnostic output variables for detailed analysis of cumulus scheme.
                # The order of vertical index is reversed for this internal diagnostic output.
                fer_out = fer_s
                fdr_out = fdr_s

                id_exit = False
                # go to 333

    with computation(PARALLEL), interval(...):
        """
        Define a release level, 'prel' and release layer, 'krel'.
        'prel' is the lowest level from which buoyancy sorting occurs, and
        'krel' is the layer index containing 'prel' in it, similar to  the
        previous definitions of 'kinv', 'klcl', and 'klfc'.    In order to
        ensure that only PBL scheme works within the PBL,  if LCL is below
        PBL top height, then 'krel = kinv', while if LCL is above  PBL top
        height, then 'krel = klcl'.   Note however that regardless of  the
        definition of 'krel', cumulus convection induces fluxes within PBL
        through 'fluxbelowinv'.  We can make cumulus convection start from
        any level, even within the PBL by appropriately defining 'krel'  &
        'prel' here. Then it must be accompanied by appropriate definition
        of source air properties, CIN, and re-setting of 'fluxbelowinv', &
        many other stuffs.
        Note that even when 'prel' is located above the PBL top height, we
        still have cumulus convection between PBL top height and 'prel':
        we simply assume that no lateral mixing occurs in this range.
        """
        if klcl < kinv:
            krel = kinv
            krel_below = krel - 1
            prel = pifc0.at(K=krel_below)
            thv0rel = thv0bot.at(K=krel_below)
        else:
            krel = klcl
            prel = plcl
            thv0rel = thv0lcl

        """
       ! Calculate cumulus base mass flux ('cbmf'), fractional area ('ufrcinv'), and !
       ! and mean vertical velocity (winv) of cumulus updraft at PBL top interface.  !
       ! Also, calculate updraft fractional area (ufrclcl) and vertical velocity  at !
       ! the LCL (wlcl). When LCL is below PBLH, cinlcl = 0 and 'ufrclcl = ufrcinv', !
       ! and 'wlcl = winv.                                                           !
       ! Only updrafts strong enough to overcome CIN can rise over PBL top interface.! 
       ! Thus,  in order to calculate cumulus mass flux at PBL top interface, 'cbmf',!
       ! we need to know 'CIN' ( the strength of potential energy barrier ) and      !
       ! 'sigmaw' ( a standard deviation of updraft vertical velocity at the PBL top !
       ! interface, a measure of turbulentce strength in the PBL ).   Naturally, the !
       ! ratio of these two variables, 'mu' - normalized CIN by TKE- is key variable !
       ! controlling 'cbmf'.  If 'mu' becomes large, only small fraction of updrafts !
       ! with very strong TKE can rise over the PBL - both 'cbmf' and 'ufrc' becomes !
       ! small, but 'winv' becomes large ( this can be easily understood by PDF of w !
       ! at PBL top ).  If 'mu' becomes small, lots of updraft can rise over the PBL !
       ! top - both 'cbmf' and 'ufrc' becomes large, but 'winv' becomes small. Thus, !
       ! all of the key variables associated with cumulus convection  at the PBL top !
       ! - 'cbmf', 'ufrc', 'winv' where 'cbmf = rho*ufrc*winv' - are a unique functi !
       ! ons of 'mu', normalized CIN. Although these are uniquely determined by 'mu',! 
       ! we usually impose two comstraints on 'cbmf' and 'ufrc': (1) because we will !
       ! simply assume that subsidence warming and drying of 'kinv-1' layer in assoc !
       ! iation with 'cbmf' at PBL top interface is confined only in 'kinv-1' layer, !
       ! cbmf must not be larger than the mass within the 'kinv-1' layer. Otherwise, !
       ! instability will occur due to the breaking of stability con. If we consider !
       ! semi-Lagrangian vertical advection scheme and explicitly consider the exten !
       ! t of vertical movement of each layer in association with cumulus mass flux, !
       ! we don't need to impose this constraint. However,  using a  semi-Lagrangian !
       ! scheme is a future research subject. Note that this constraint should be ap !
       ! plied for all interfaces above PBL top as well as PBL top interface.   As a !
       ! result, this 'cbmf' constraint impose a 'lower' limit on mu - 'mumin0'. (2) !
       ! in order for mass flux parameterization - rho*(w'a')= M*(a_c-a_e) - to   be !
       ! valid, cumulus updraft fractional area should be much smaller than 1.    In !
       ! current code, we impose 'rmaxfrac = 0.1 ~ 0.2'   through the whole vertical !
       ! layers where cumulus convection occurs. At the PBL top interface,  the same !
       ! constraint is made by imposing another lower 'lower' limit on mu, 'mumin1'. !
       ! After that, also limit 'ufrclcl' to be smaller than 'rmaxfrac' by 'mumin2'. !
       ! --------------------------------------------------------------------------- !
       
       ! --------------------------------------------------------------------------- !
       ! Calculate normalized CIN, 'mu' satisfying all the three constraints imposed !
       ! on 'cbmf'('mumin0'), 'ufrc' at the PBL top - 'ufrcinv' - ( by 'mumin1' from !
       ! a parameter sentence), and 'ufrc' at the LCL - 'ufrclcl' ( by 'mumin2').    !
       ! Note that 'cbmf' does not change between PBL top and LCL  because we assume !
       ! that buoyancy sorting does not occur when cumulus updraft is unsaturated.   !
       ! ----------------------------------------------------------------------------!
       """

        if use_CINcin == 1:
            wcrit = sqrt(2.0 * cin * rbuoy)
        else:
            wcrit = sqrt(2.0 * cinlcl * rbuoy)
        sigmaw = sqrt(rkfre * tkeavg + epsvarw)
        mu = wcrit / sigmaw / 1.4142
        if mu >= 3.0:
            id_exit = True
            # go to 333
        kbelow = kinv - 1
        rho0inv = pifc0.at(K=kbelow) / (
            constants.MAPL_RDRY * thv0top.at(K=kbelow) * exnifc0.at(K=kbelow)
        )
        cbmf = (rho0inv * sigmaw / 2.5066) * exp(-(mu**2))

        # 1. 'cbmf' constraint
        cbmflimit = 0.9 * dp0.at(K=kbelow) / constants.MAPL_GRAV / dt
        mumin0 = 0.0
        if cbmf > cbmflimit:
            mumin0 = sqrt(-log(2.5066 * cbmflimit / rho0inv / sigmaw))
        # 2. 'ufrcinv' constraint
        mu = max(max(mu, mumin0), mumin1)
        # 3. 'ufrclcl' constraint
        mulcl = sqrt(2.0 * cinlcl * rbuoy) / 1.4142 / sigmaw
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
            if mu == mumin2:
                limit_ufrc = 1.0
        if mu == mumin0:
            limit_cbmf = 1.0
        if mu == mumin1:
            limit_ufrc = 1.0

        """
        Calculate final ['cbmf','ufrcinv','winv'] at the PBL top interface. 
        Note that final 'cbmf' here is obtained in such that 'ufrcinv' and  
        'ufrclcl' are smaller than ufrcmax with no instability.             
        """

        cbmf = rkfre * (rho0inv * sigmaw / 2.5066) * exp(-(mu**2))
        winv = sigmaw * (2.0 / 2.5066) * exp(-(mu**2)) / erfc(mu)
        ufrcinv = cbmf / winv / rho0inv

        """
        Calculate ['ufrclcl','wlcl'] at the LCL. When LCL is below PBL top, 
        it automatically becomes 'ufrclcl = ufrcinv' & 'wlcl = winv', since 
        it was already set to 'cinlcl=0' if LCL is below PBL top interface. 
        Note 'cbmf' at the PBL top is the same as 'cbmf' at the LCL.  Note  
        also that final 'cbmf' here is obtained in such that 'ufrcinv' and  
        'ufrclcl' are smaller than ufrcmax and there is no instability.     
        By construction, it must be 'wlcl > 0' but for assurance, I checked 
        this again in the below block. If 'ufrclcl < 0.1%', just exit.
        """
        wtw = winv * winv - 2.0 * cinlcl * rbuoy
        if wtw <= 0.0:
            id_exit = True
            # go to 333
        wlcl = sqrt(wtw)
        ufrclcl = cbmf / wlcl / rho0inv
        wrel = wlcl
        if ufrclcl <= 0.0001:
            id_exit = True
            # go to 333
        if THIS_K == (krel - 1):
            ufrc = ufrclcl

        # Below is just diagnostic output for detailed analysis of cumulus scheme
        ufrcinvbase = ufrcinv
        winvbase = winv
        if THIS_K >= (kinv - 1) and THIS_K <= (krel - 1):
            umf = cbmf
            wu = winv

        """
        Define updraft properties at the level where buoyancy sorting starts to be 
        happening, i.e., by definition, at 'prel' level within the release layer.  
        Because no lateral entrainment occurs upto 'prel', conservative scalars of  
        cumulus updraft at release level is same as those of source air.  However,  
        horizontal momentums of source air are modified by horizontal PGF forcings  
        from PBL top interface to 'prel'.  For this case, we should add additional 
        horizontal momentum from PBL top interface to 'prel' as will be done below 
        to 'usrc' and 'vsrc'. Note that below cumulus updraft properties - umf, wu,
        thlu, qtu, thvu, uu, vu - are defined all interfaces not at the layer mid- 
        point. From the index notation of cumulus scheme, wu(k) is the cumulus up- 
        draft vertical velocity at the top interface of k layer.                   
        Diabatic horizontal momentum forcing should be treated as a kind of 'body' 
        forcing without actual mass exchange between convective updraft and        
        environment, but still taking horizontal momentum from the environment to  
        the convective updrafts. Thus, diabatic convective momentum transport      
        vertically redistributes environmental horizontal momentum. 
        """
        if THIS_K == (krel - 1):
            emf = 0.0
            umf = cbmf
            wu = wrel
            thlu = thlsrc
            qtu = qtsrc
        thj, qvj, qlj, qij, qse, id_check = conden(prel, thlsrc, qtsrc, ese, esx)
        if id_check == 1:
            id_exit = True
            # go to 333
        if THIS_K == (krel - 1):
            thvu = thj * (1.0 + zvir * qvj - qlj - qij)

        uplus = 0.0
        vplus = 0.0
        if krel == kinv:
            kbelow = kinv - 1
            uplus = PGFc * ssu0.at(K=kinv) * (prel - pifc0.at(K=kbelow))
            vplus = PGFc * ssv0.at(K=kinv) * (prel - pifc0.at(K=kbelow))
        else:
            if THIS_K >= kinv and THIS_K <= max(krel - 1, kinv):
                kinv_bel = kinv - 1
                uplus = uplus + PGFc * ssu0.at(K=kinv) * (
                    pifc0.at(K=kinv) - pifc0.at(K=kinv_bel)
                )
                vplus = vplus + PGFc * ssv0.at(K=kinv) * (
                    pifc0.at(K=kinv) - pifc0.at(K=kinv_bel)
                )
            kbelow = krel - 1
            uplus = uplus + PGFc * ssu0.at(K=krel) * (prel - pifc0.at(K=kbelow))
            vplus = vplus + PGFc * ssv0.at(K=krel) * (prel - pifc0.at(K=kbelow))

        if THIS_K == (krel - 1):
            uu = usrc + uplus
            vu = vsrc + vplus

        if dotransport == 1.0:
            n = 0
            while n < ncnst:
                # REVISIT THIS!!!
                # tru[0, 0, krel - 1][n] = trsrc[0, 0, 0][n]
                n += 1

        """
        Define environmental properties at the level where buoyancy sorting occurs 
        ('pe', normally, layer midpoint except in the 'krel' layer). In the 'krel' 
        layer where buoyancy sorting starts to occur, however, 'pe' is defined     
        differently because LCL is regarded as lower interface for mixing purpose. 
        """
        pe = 0.5 * (prel + pifc0.at(K=krel))
        qsat_pe = 0.5 * (prel + pifc0.at(K=krel))
        dpe = prel - pifc0.at(K=krel)
        exne = exnerfn(pe)
        thvebot = thv0rel
        thle = thl0.at(K=krel) + ssthl0.at(K=krel) * (pe - pmid0.at(K=krel))
        qte = qt0.at(K=krel) + ssqt0.at(K=krel) * (pe - pmid0.at(K=krel))
        ue = u0.at(K=krel) + ssu0.at(K=krel) * (pe - pmid0.at(K=krel))
        ve = v0.at(K=krel) + ssv0.at(K=krel) * (pe - pmid0.at(K=krel))
        if dotransport == 1.0:
            n = 0
            while n < ncnst:
                tre[0, 0, 0][n] = tr0.at(K=krel, ddim=[n]) + sstr0.at(
                    K=krel, ddim=[n]
                ) * (pe - pmid0.at(K=krel))
                n += 1
'''


# def buoyancy_sorting_mixing(
#     id_exit: BoolField,
#     tscaleh: FloatFieldIJ,
#     krel: IntField,
#     wlcl: FloatFieldIJ,
#     prel: FloatFieldIJ,
#     pifc0: FloatField,
#     thv0rel: FloatFieldIJ,
#     thl0: FloatField,
#     ssthl0: FloatField,
#     pmid0: FloatField,
#     qt0: FloatField,
#     ssqt0: FloatField,
#     u0: FloatField,
#     ssu0: FloatField,
#     v0: FloatField,
#     ssv0: FloatField,
#     dotransport: Float,
#     ncnst: Int,
#     tr0: FloatField_NTracers,
#     sstr0: FloatField_NTracers,
#     k0: Int,
#     thlu: FloatFieldIJ,
#     qtu: FloatFieldIJ,
#     wu: FloatFieldIJ,
#     ese: FloatField_Extra_Dim,
#     esx: FloatField_Extra_Dim,
#     qsat_pe: FloatFieldIJ,
#     umf_out: FloatField,
#     dcm_out: FloatField,
#     qvten_out: FloatField,
#     qlten_out: FloatField,
#     qiten_out: FloatField,
#     sten_out: FloatField,
#     uten_out: FloatField,
#     vten_out: FloatField,
#     qrten_out: FloatField,
#     qsten_out: FloatField,
#     cufrc_out: FloatField,
#     cush_inout: FloatFieldIJ,
#     qldet_out: FloatField,
#     qidet_out: FloatField,
#     qtflx_out: FloatField,
#     slflx_out: FloatField,
#     uflx_out: FloatField,
#     vflx_out: FloatField,
#     fer_out: FloatField,
#     fdr_out: FloatField,
#     rbuoy: Float,
#     zifc0: FloatField,
#     zmid0: FloatField,
#     dp0: FloatField,
#     dt: Float,
#     umf: FloatField,
#     niter_xc: Int,
#     criqc: Float,
#     rle: Float,
#     mixscale: Float,
#     rkm: Float,
#     detrhgt: Float,
#     cridist_opt: Int,
#     PGFc: Float,
#     tru: FloatField_NTracers,
#     tre: FloatField_NTracers,
#     exnifc0: FloatField,
#     thv0top: FloatField,
#     thv0bot: FloatField,
#     rmaxfrac: Float,
#     exnmid0: FloatField,
#     emf: FloatField,
#     tru_emf: FloatField_NTracers,
#     rdrag: Float,
#     use_self_detrain: Int,
#     use_cumpenent: Int,
#     rpen: Float,
#     qtsrc: FloatFieldIJ,
#     kinv: IntField,
#     cbmf: FloatField,
#     thlsrc: FloatFieldIJ,
#     usrc: FloatFieldIJ,
#     vsrc: FloatFieldIJ,
#     trsrc: FloatField_NTracers,
#     trflx: FloatField_NTracers,
#     use_momenflx: Int,
#     ql0: FloatField,
#     qi0: FloatField,
#     frc_rasn: Float,
#     qv0: FloatField,
#     s0: FloatField,
#     t0: FloatField,
#     iteration: Int,
#     tr0_s: FloatField_NTracers,
#     trten: FloatField_NTracers,
#     umf_s: FloatField,
#     slflx_s: FloatField,
#     qtflx_s: FloatField,
#     uflx_s: FloatField,
#     vflx_s: FloatField,
#     cin: FloatFieldIJ,
#     cinlcl: FloatFieldIJ,
#     ufrc_s: FloatField,
#     ufrclcl: FloatFieldIJ,
# ):
#     with computation(FORWARD), interval(...):
#         # Define cumulus scale height.
#         # Cumulus scale height is defined as the maximum height cumulus can reach.
#         # In case of premitive code, cumulus scale height ('cush')  at the current
#         # time step was assumed to be the same as 'cush' of previous time step.
#         # However, I directly calculated cush at each time step using an iterative
#         # method. Note that within the cumulus scheme, 'cush' information is  used
#         # only at two places during buoyancy-sorting process:
#         # (1) Even negatively buoyancy mixtures with strong vertical velocity
#         #      enough to rise up to 'rle*scaleh' (rle = 0.1) from pe are entrained
#         #      into cumulus updraft,
#         # (2) The amount of mass that is involved in buoyancy-sorting mixing
#         #       process at pe is rei(k) = rkm/scaleh/rho*g [Pa-1]
#         # In terms of (1), I think critical stopping distance might be replaced by
#         # layer thickness. In future, we will use rei(k) = (0.5*rkm/z0(k)/rho/g).
#         # In the premitive code,  'scaleh' was largely responsible for the jumping
#         # variation of precipitation amount.

#         if id_exit == False:
#             scaleh = tscaleh
#             if tscaleh <= 0.0:
#                 scaleh = 1000

#             # Save time : Set iter_scaleh = 1. This will automatically use 'cush' from the previous
#             # time step at the first implicit iteration. At the second implicit iteration, it will
#             # use the updated 'cush' by the first implicit cin. So, this updating has an effect of
#             # doing one iteration for cush calculation, which is good. So, only this setting of
#             # 'iter_scaleh = 1' is sufficient-enough to save computation time.

#             iter_scaleh = 1

#             # Initialization of 'kbup' and 'kpen'
#             # 'kbup' is the top-most layer in which cloud buoyancy is positive
#             # both at the top and bottom interface of the layer. 'kpen' is the
#             # layer upto which cumulus panetrates ,i.e., cumulus w at the base
#             # interface is positive, but becomes negative at the top interface.
#             # Here, we initialize 'kbup' and 'kpen'. These initializations are
#             # not trivial but important, expecially   in calculating turbulent
#             # fluxes without confliction among several physics as explained in
#             # detail in the part of turbulent fluxes calculation later.   Note
#             # that regardless of whether 'kbup' and 'kpen' are updated or  not
#             # during updraft motion,  penetrative entrainments are dumped down
#             # across the top interface of 'kbup' later.      More specifically,
#             # penetrative entrainment heat and moisture fluxes are  calculated
#             # from the top interface of 'kbup' layer  to the base interface of
#             # 'kpen' layer. Because of this, initialization of 'kbup' & 'kpen'
#             # influence the convection system when there are not updated.  The
#             # below initialization of 'kbup = krel' assures  that  penetrative
#             # entrainment fluxes always occur at interfaces above the PBL  top
#             # interfaces (i.e., only at interfaces k >=kinv ), which seems  to
#             # be attractable considering that the most correct fluxes  at  the
#             # PBL top interface can be ontained from the 'fluxbelowinv'  using
#             # reconstructed PBL height.

#             kbup = krel
#             kpen = krel

#             # Since 'wtw' is continuously updated during vertical motion,
#             # I need below initialization command within this 'iter_scaleh'
#             # do loop. Similarily, I need initializations of environmental
#             # properties at 'krel' layer as below.

#             wtw = wlcl * wlcl
#             pe = 0.5 * (prel + pifc0.at(K=krel))
#             dpe = prel - pifc0.at(K=krel)
#             exne = exnerfn(pe)
#             thvebot = thv0rel
#             thle = thl0.at(K=krel) + ssthl0.at(K=krel) * (pe - pmid0.at(K=krel))
#             qte = qt0.at(K=krel) + ssqt0.at(K=krel) * (pe - pmid0.at(K=krel))
#             ue = u0.at(K=krel) + ssu0.at(K=krel) * (pe - pmid0.at(K=krel))
#             ve = v0.at(K=krel) + ssv0.at(K=krel) * (pe - pmid0.at(K=krel))

#             if dotransport == 1.0:
#                 n = 0
#                 # Loop over tracers
#                 while n < ncnst:
#                     tre[0, 0, 0][n] = tr0.at(K=krel, ddim=[n]) + sstr0.at(
#                         K=krel, ddim=[n]
#                     ) * (pe - pmid0.at(K=krel, ddim=[n]))
#                     n += 1

#             # Cumulus rises upward from 'prel' ( or base interface of  'krel' layer )
#             # until updraft vertical velocity becomes zero.
#             # Buoyancy sorting is performed via two stages. (1) Using cumulus updraft
#             # properties at the base interface of each layer,perform buoyancy sorting
#             # at the layer mid-point, 'pe',  and update cumulus properties at the top
#             # interface, and then  (2) by averaging updated cumulus properties at the
#             # top interface and cumulus properties at the base interface,   calculate
#             # cumulus updraft properties at pe that will be used  in buoyancy sorting
#             # mixing - thlue, qtue and, wue.  Using this averaged properties, perform
#             # buoyancy sorting again at pe, and re-calculate fer(k) and fdr(k). Using
#             # this recalculated fer(k) and fdr(k),  finally calculate cumulus updraft
#             # properties at the top interface - thlu, qtu, thvu, uu, vu. In the below,
#             # 'iter_xc = 1' performs the first stage, while 'iter_xc= 2' performs the
#             # second stage. We can increase the number of iterations, 'nter_xc'.as we
#             # want, but a sample test indicated that about 3 - 5 iterations  produced
#             # satisfactory converent solution. Finally, identify 'kbup' and 'kpen'.
#             k = krel
#             while k <= (k0 - 1):
#                 km1 = k - 1

#                 thlue = thlu.at(K=k)
#                 qtue = qtu.at(K=k)
#                 wue = wu.at(K=k)
#                 wtwb = wtw

#                 iter_xc = 1
#                 while iter_xc <= niter_xc:
#                     wtw = wu.at(K=km1) * wu.at(K=km1)

#                     # Calculate environmental and cumulus saturation 'excess' at 'pe'.
#                     # Note that in order to calculate saturation excess, we should use
#                     # liquid water temperature instead of temperature  as the argument
#                     # of "qsat". But note normal argument of "qsat" is temperature.

#                     thj, qvj, qlj, qij, qse, id_check = conden(pe, thle, qte, ese, esx)
#                     if id_check == 1:
#                         id_exit = True
#                         umf_out[0, 0, 1] = 0.0
#                         dcm_out = 0.0
#                         qvten_out = 0.0
#                         qlten_out = 0.0
#                         qiten_out = 0.0
#                         sten_out = 0.0
#                         uten_out = 0.0
#                         vten_out = 0.0
#                         qrten_out = 0.0
#                         qsten_out = 0.0
#                         cufrc_out = 0.0
#                         cush_inout = -1.0
#                         qldet_out = 0.0
#                         qidet_out = 0.0
#                         qtflx_out[0, 0, 1] = 0.0
#                         slflx_out[0, 0, 1] = 0.0
#                         uflx_out[0, 0, 1] = 0.0
#                         vflx_out[0, 0, 1] = 0.0
#                         fer_out = constants.MAPL_UNDEF
#                         fdr_out = constants.MAPL_UNDEF

#                     if id_exit == False:
#                         thv0j = thj * (1.0 + zvir * qvj - qlj - qij)
#                         rhomid0j = pe / (constants.MAPL_RDRY * thv0j * exne)
#                         qsat_arg = thle * exne
#                         qs, _ = QSat_Float(ese, esx, qsat_arg, qsat_pe / 100.0)
#                         excess0 = qte - qs

#                         thj, qvj, qlj, qij, qse, id_check = conden(
#                             pe, thlue, qtue, ese, esx
#                         )
#                         if id_check == 1:
#                             id_exit = True
#                             umf_out[0, 0, 1] = 0.0
#                             dcm_out = 0.0
#                             qvten_out = 0.0
#                             qlten_out = 0.0
#                             qiten_out = 0.0
#                             sten_out = 0.0
#                             uten_out = 0.0
#                             vten_out = 0.0
#                             qrten_out = 0.0
#                             qsten_out = 0.0
#                             cufrc_out = 0.0
#                             cush_inout = -1.0
#                             qldet_out = 0.0
#                             qidet_out = 0.0
#                             qtflx_out[0, 0, 1] = 0.0
#                             slflx_out[0, 0, 1] = 0.0
#                             uflx_out[0, 0, 1] = 0.0
#                             vflx_out[0, 0, 1] = 0.0
#                             fer_out = constants.MAPL_UNDEF
#                             fdr_out = constants.MAPL_UNDEF

#                     if id_exit == False:
#                         # Detrain excessive condensate larger than 'criqc' from the cumulus
#                         # updraft before performing buoyancy sorting. All I should to do is
#                         # to update 'thlue' &  'que' here. Below modification is completely
#                         # compatible with the other part of the code since 'thule' & 'qtue'
#                         # are used only for buoyancy sorting. I found that as long as I use
#                         # 'niter_xc >= 2',  detraining excessive condensate before buoyancy
#                         # sorting has negligible influence on the buoyancy sorting results.

#                         if (qlj + qij) > criqc:  # DONT FORGET TO SERIALIZE criqc !!!
#                             exql = ((qlj + qij) - criqc) * qlj / (qlj + qij)
#                             exqi = ((qlj + qij) - criqc) * qij / (qlj + qij)
#                             qtue = qtue - exql - exqi
#                             thlue = (
#                                 thlue
#                                 + (
#                                     constants.MAPL_LATENT_HEAT_VAPORIZATION
#                                     / constants.MAPL_CP
#                                     / exne
#                                 )
#                                 * exql
#                                 + (
#                                     constants.MAPL_LATENT_HEAT_SUBLIMATION
#                                     / constants.MAPL_CP
#                                     / exne
#                                 )
#                                 * exqi
#                             )

#                         thj, qvj, qlj, qij, qse, id_check = conden(
#                             pe, thlue, qtue, ese, esx
#                         )
#                         if id_check == 1:
#                             id_exit = True
#                             umf_out[0, 0, 1] = 0.0
#                             dcm_out = 0.0
#                             qvten_out = 0.0
#                             qlten_out = 0.0
#                             qiten_out = 0.0
#                             sten_out = 0.0
#                             uten_out = 0.0
#                             vten_out = 0.0
#                             qrten_out = 0.0
#                             qsten_out = 0.0
#                             cufrc_out = 0.0
#                             cush_inout = -1.0
#                             qldet_out = 0.0
#                             qidet_out = 0.0
#                             qtflx_out[0, 0, 1] = 0.0
#                             slflx_out[0, 0, 1] = 0.0
#                             uflx_out[0, 0, 1] = 0.0
#                             vflx_out[0, 0, 1] = 0.0
#                             fer_out = constants.MAPL_UNDEF
#                             fdr_out = constants.MAPL_UNDEF

#                     if id_exit == False:
#                         thvj = thj * (1.0 + zvir * qvj - qlj - qij)
#                         tj = (
#                             thj * exne
#                         )  # This 'tj' is used for computing thermo. coeffs. below
#                         qsat_arg = thlue * exne
#                         qs = QSat_Float(ese, esx, qsat_arg, qsat_pe / 100.0)
#                         excessu = qtue - qs

#                         # Calculate critical mixing fraction, 'xc'. Mixture with mixing ratio
#                         # smaller than 'xc' will be entrained into cumulus updraft.  Both the
#                         # saturated updrafts with 'positive buoyancy' or 'negative buoyancy +
#                         # strong vertical velocity enough to rise certain threshold distance'
#                         # are kept into the updraft in the below program. If the core updraft
#                         # is unsaturated, we can set 'xc = 0' and let the cumulus  convection
#                         # still works or we may exit.
#                         # Current below code does not entrain unsaturated mixture. However it
#                         # should be modified such that it also entrain unsaturated mixture.

#                         # cridis : Critical stopping distance for buoyancy sorting purpose.
#                         #        scaleh is only used here.

#                         if cridist_opt == 0:
#                             cridis = rle * scaleh  # Original code
#                         else:
#                             cridis = rle * (
#                                 zifc0.at(K=k) - zifc0.at(K=k - 1)
#                             )  # New code

#                         # Buoyancy Sorting
#                         # Case 1 : When both cumulus and env. are unsaturated or saturated.
#                         xsat = 0.0

#                         if (excessu <= 0.0 and excess0 <= 0.0) or (
#                             excessu >= 0.0 and excess0 >= 0.0
#                         ):
#                             xc = min(
#                                 1.0,
#                                 max(
#                                     0.0,
#                                     1.0
#                                     - 2.0
#                                     * rbuoy
#                                     * constants.MAPL_GRAV
#                                     * cridis
#                                     / wue**2.0
#                                     * (1.0 - thvj / thv0j),
#                                 ),
#                             )
#                             aquad = 0.0
#                             bquad = 0.0
#                             cquad = 0.0
#                             if excessu > 0.0:
#                                 xsat = 1.0
#                             else:
#                                 xsat = 0.0
#                         else:
#                             # Case 2 : When either cumulus or env. is saturated. !
#                             xsat = excessu / (excessu - excess0)
#                             thlxsat = thlue + xsat * (thle - thlue)
#                             qtxsat = qtue + xsat * (qte - qtue)
#                             thj, qvj, qlj, qij, qse, id_check = conden(
#                                 pe, thlxsat, qtxsat, ese, esx
#                             )
#                             if id_check == 1:
#                                 id_exit = True
#                                 umf_out[0, 0, 1] = 0.0
#                                 dcm_out = 0.0
#                                 qvten_out = 0.0
#                                 qlten_out = 0.0
#                                 qiten_out = 0.0
#                                 sten_out = 0.0
#                                 uten_out = 0.0
#                                 vten_out = 0.0
#                                 qrten_out = 0.0
#                                 qsten_out = 0.0
#                                 cufrc_out = 0.0
#                                 cush_inout = -1.0
#                                 qldet_out = 0.0
#                                 qidet_out = 0.0
#                                 qtflx_out[0, 0, 1] = 0.0
#                                 slflx_out[0, 0, 1] = 0.0
#                                 uflx_out[0, 0, 1] = 0.0
#                                 vflx_out[0, 0, 1] = 0.0
#                                 fer_out = constants.MAPL_UNDEF
#                                 fdr_out = constants.MAPL_UNDEF

#                             if id_exit == False:
#                                 thvxsat = thj * (1.0 + zvir * qvj - qlj - qij)

#                                 # kk=1 : Cumulus segment, kk=2 : Environment segment
#                                 kk = 1
#                                 while kk <= 2:
#                                     if xsat == 1.0:
#                                         xsat = 1.0 + 1e-6
#                                     if kk == 1:
#                                         thv_x0 = thvj
#                                         thv_x1 = (1.0 - 1.0 / xsat) * thvj + (
#                                             1.0 / xsat
#                                         ) * thvxsat
#                                     else:
#                                         thv_x1 = thv0j
#                                         thv_x0 = (xsat / (xsat - 1.0)) * thv0j + (
#                                             1.0 / (1.0 - xsat)
#                                         ) * thvxsat

#                                     aquad = wue**2
#                                     bquad = (
#                                         2.0
#                                         * rbuoy
#                                         * constants.MAPL_GRAV
#                                         * cridis
#                                         * (thv_x1 - thv_x0)
#                                         / thv0j
#                                         - 2.0 * wue**2
#                                     )
#                                     cquad = (
#                                         2.0
#                                         * rbuoy
#                                         * constants.MAPL_GRAV
#                                         * cridis
#                                         * (thv_x0 - thv0j)
#                                         / thv0j
#                                         + wue**2
#                                     )
#                                     if kk == 1:
#                                         if (bquad**2 - 4.0 * aquad * cquad) >= 0.0:
#                                             xs1, xs2, status = roots(
#                                                 aquad, bquad, cquad
#                                             )
#                                             x_cu = min(
#                                                 1.0, max(0.0, min(xsat, min(xs1, xs2)))
#                                             )
#                                         else:
#                                             x_cu = xsat

#                                     else:
#                                         if (bquad**2 - 4.0 * aquad * cquad) >= 0.0:
#                                             xs1, xs2, status = roots(
#                                                 aquad, bquad, cquad
#                                             )
#                                             x_en = min(
#                                                 1.0, max(0.0, max(xsat, min(xs1, xs2)))
#                                             )
#                                         else:
#                                             x_en = 1.0

#                                     kk += 1

#                                 if x_cu == xsat:
#                                     xc = max(x_cu, x_en)
#                                 else:
#                                     xc = x_cu

#                                 # Compute fractional lateral entrainment & detrainment rate in each layers.
#                                 # The unit of rei(k), fer(k), and fdr(k) is [Pa-1].  Alternative choice of
#                                 # 'rei(k)' is also shown below, where coefficient 0.5 was from approximate
#                                 # tuning against the BOMEX case.
#                                 # In order to prevent the onset of instability in association with cumulus
#                                 # induced subsidence advection, cumulus mass flux at the top interface  in
#                                 # any layer should be smaller than ( 90% of ) total mass within that layer.
#                                 # I imposed limits on 'rei(k)' as below,  in such that stability condition
#                                 # is always satisfied.
#                                 # Below limiter of 'rei(k)' becomes negative for some cases, causing error.
#                                 # So, for the time being, I came back to the original limiter.
#                                 ee2 = xc**2
#                                 ud2 = 1.0 - 2.0 * xc + xc**2
#                                 if min(scaleh, mixscale) != 0.0:
#                                     if THIS_K == k:
#                                         rei = (
#                                             (
#                                                 rkm
#                                                 + max(
#                                                     0.0,
#                                                     (zmid0.at(K=k) - detrhgt) / 200.0,
#                                                 )
#                                             )
#                                             / min(scaleh, mixscale)
#                                             / constants.MAPL_GRAV
#                                             / rhomid0j
#                                         )
#                                 else:
#                                     if THIS_K == k:
#                                         rei = (
#                                             0.5
#                                             * rkm
#                                             / zmid0.at(K=k)
#                                             / constants.MAPL_GRAV
#                                             / rhomid0j
#                                         )

#                                 if xc > 0.5:
#                                     if THIS_K == k:
#                                         rei = min(
#                                             rei(k),
#                                             0.9
#                                             * log(
#                                                 dp0.at(K=k)
#                                                 / constants.MAPL_GRAV
#                                                 / dt
#                                                 / umf.at(K=km1)
#                                                 + 1.0
#                                             )
#                                             / dpe
#                                             / (2.0 * xc - 1.0),
#                                         )
#                                 if THIS_K == k:
#                                     fer = rei.at(K=k) * ee2
#                                     fdr = rei.at(K=k) * ud2
#                                     xco = xc

#                                 # Iteration Start due to 'maxufrc' constraint
#                                 # Calculate cumulus updraft mass flux and penetrative entrainment mass flux.
#                                 # Note that  non-zero penetrative entrainment mass flux will be asigned only
#                                 # to interfaces from the top interface of 'kbup' layer to the base interface
#                                 # of 'kpen' layer as will be shown later.
#                                 if THIS_K == k:
#                                     umf[0, 0, 1] = umf.at(K=km1) * exp(
#                                         dpe * (fer.at(K=k) - fdr.at(K=k))
#                                     )
#                                     emf[0, 0, 1] = 0.0

#                                     dcm = (
#                                         0.5
#                                         * (umf.at(K=k) + umf.at(K=km1))
#                                         * rei.at(K=k)
#                                         * dpe
#                                         * min(1.0, max(0.0, xsat - xc))
#                                     )
#                                     # dcm(k) = min(1.,max(0.,xsat-xc))

#                                     # Compute cumulus updraft properties at the top interface.
#                                     # Also use Tayler expansion in order to treat limiting case

#                                     if fer.at(K=k) * dpe < 1.0e-4:
#                                         thlu = (
#                                             thlu.at(K=km1)
#                                             + (
#                                                 thle
#                                                 + ssthl0.at(K=k) * dpe / 2.0
#                                                 - thlu.at(K=km1)
#                                             )
#                                             * fer.at(K=k)
#                                             * dpe
#                                         )
#                                         qtu = (
#                                             qtu.at(K=km1)
#                                             + (
#                                                 qte
#                                                 + ssqt0.at(K=k) * dpe / 2.0
#                                                 - qtu.at(K=km1)
#                                             )
#                                             * fer.at(K=k)
#                                             * dpe
#                                         )
#                                         uu = (
#                                             uu.at(K=km1)
#                                             + (
#                                                 ue
#                                                 + ssu0.at(K=k) * dpe / 2.0
#                                                 - uu.at(K=km1)
#                                             )
#                                             * fer.at(K=k)
#                                             * dpe
#                                             - PGFc * ssu0.at(K=k) * dpe
#                                         )
#                                         vu = (
#                                             vu.at(K=km1)
#                                             + (
#                                                 ve
#                                                 + ssv0.at(K=k) * dpe / 2.0
#                                                 - vu.at(K=km1)
#                                             )
#                                             * fer.at(K=k)
#                                             * dpe
#                                             - PGFc * ssv0.at(K=k) * dpe
#                                         )
#                                         if dotransport == 1.0:
#                                             n = 0
#                                             while n < ncnst:
#                                                 tru[0, 0, 0][n] = (
#                                                     tru.at(K=km1, ddim=[n])
#                                                     + (
#                                                         tre[0, 0, 0][n]
#                                                         + sstr0.at(K=k, ddim=[n])
#                                                         * dpe
#                                                         / 2.0
#                                                         - tru.at(K=km1, ddim=[n])
#                                                     )
#                                                     * fer.at(K=k)
#                                                     * dpe
#                                                 )
#                                                 n += 1

#                                     else:
#                                         thlu = (
#                                             thle
#                                             + ssthl0.at(K=k) / fer.at(K=k)
#                                             - ssthl0.at(K=k) * dpe / 2.0
#                                         ) - (
#                                             thle
#                                             + ssthl0.at(K=k) * dpe / 2.0
#                                             - thlu.at(K=km1)
#                                             + ssthl0.at(K=k) / fer.at(K=k)
#                                         ) * exp(
#                                             -fer.at(K=k) * dpe
#                                         )
#                                         qtu = (
#                                             qte
#                                             + ssqt0.at(K=k) / fer.at(K=k)
#                                             - ssqt0.at(K=k) * dpe / 2.0
#                                         ) - (
#                                             qte
#                                             + ssqt0.at(K=k) * dpe / 2.0
#                                             - qtu.at(K=km1)
#                                             + ssqt0.at(K=k) / fer.at(K=k)
#                                         ) * exp(
#                                             -fer.at(K=k) * dpe
#                                         )
#                                         uu = (
#                                             ue
#                                             + (1.0 - PGFc) * ssu0.at(K=k) / fer.at(K=k)
#                                             - ssu0.at(K=k) * dpe / 2.0
#                                         ) - (
#                                             ue
#                                             + ssu0.at(K=k) * dpe / 2.0
#                                             - uu.at(K=km1)
#                                             + (1.0 - PGFc) * ssu0.at(K=k) / fer.at(K=k)
#                                         ) * exp(
#                                             -fer.at(K=k) * dpe
#                                         )
#                                         vu = (
#                                             ve
#                                             + (1.0 - PGFc) * ssv0.at(K=k) / fer.at(K=k)
#                                             - ssv0.at(K=k) * dpe / 2.0
#                                         ) - (
#                                             ve
#                                             + ssv0.at(K=k) * dpe / 2.0
#                                             - vu.at(K=km1)
#                                             + (1.0 - PGFc) * ssv0.at(K=k) / fer.at(K=k)
#                                         ) * exp(
#                                             -fer.at(K=k) * dpe
#                                         )
#                                         if dotransport == 1.0:
#                                             n = 0
#                                             while n < ncnst:
#                                                 tru[0, 0, 0][n] = (
#                                                     tre[0, 0, 0][n]
#                                                     + sstr0.at(K=k, ddim=[n])
#                                                     / fer.at(K=k)
#                                                     - sstr0.at(K=k, ddim=[n])
#                                                     * dpe
#                                                     / 2.0
#                                                 ) - (
#                                                     tre[0, 0, 0][n]
#                                                     + sstr0.at(K=k, ddim=[n])
#                                                     * dpe
#                                                     / 2.0
#                                                     - tru.at(K=km1, ddim=[n])
#                                                     + sstr0.at(K=k, ddim=[n])
#                                                     / fer.at(K=k)
#                                                 ) * exp(
#                                                     -fer.at(K=k) * dpe
#                                                 )
#                                                 n += 1

#                                 # Expel some of cloud water and ice from cumulus  updraft at the top
#                                 # interface.  Note that this is not 'detrainment' term  but a 'sink'
#                                 # term of cumulus updraft qt ( or one part of 'source' term of  mean
#                                 # environmental qt ). At this stage, as the most simplest choice, if
#                                 # condensate amount within cumulus updraft is larger than a critical
#                                 # value, 'criqc', expels the surplus condensate from cumulus updraft
#                                 # to the environment. A certain fraction ( e.g., 'frc_sus' ) of this
#                                 # expelled condesnate will be in a form that can be suspended in the
#                                 # layer k where it was formed, while the other fraction, '1-frc_sus'
#                                 # will be in a form of precipitatble (e.g.,can potentially fall down
#                                 # across the base interface of layer k ). In turn we should describe
#                                 # subsequent falling of precipitable condensate ('1-frc_sus') across
#                                 # the base interface of the layer k, &  evaporation of precipitating
#                                 # water in the below layer k-1 and associated evaporative cooling of
#                                 # the later, k-1, and falling of 'non-evaporated precipitating water
#                                 # ( which was initially formed in layer k ) and a newly-formed preci
#                                 # pitable water in the layer, k-1', across the base interface of the
#                                 # lower layer k-1.  Cloud microphysics should correctly describe all
#                                 # of these process.  In a near future, I should significantly modify
#                                 # this cloud microphysics, including precipitation-induced downdraft
#                                 # also.

#                                 thj, qvj, qlj, qij, qse, id_check = conden(
#                                     pifc0.at(K=k), thlu.at(K=k), qtu.at(K=k), ese, esx
#                                 )
#                                 if id_check == 1:
#                                     id_exit = True
#                                     umf_out[0, 0, 1] = 0.0
#                                     dcm_out = 0.0
#                                     qvten_out = 0.0
#                                     qlten_out = 0.0
#                                     qiten_out = 0.0
#                                     sten_out = 0.0
#                                     uten_out = 0.0
#                                     vten_out = 0.0
#                                     qrten_out = 0.0
#                                     qsten_out = 0.0
#                                     cufrc_out = 0.0
#                                     cush_inout = -1.0
#                                     qldet_out = 0.0
#                                     qidet_out = 0.0
#                                     qtflx_out[0, 0, 1] = 0.0
#                                     slflx_out[0, 0, 1] = 0.0
#                                     uflx_out[0, 0, 1] = 0.0
#                                     vflx_out[0, 0, 1] = 0.0
#                                     fer_out = constants.MAPL_UNDEF
#                                     fdr_out = constants.MAPL_UNDEF

#                             if id_exit == False:
#                                 if (qlj + qij) > criqc:
#                                     exql = ((qlj + qij) - criqc) * qlj / (qlj + qij)
#                                     exqi = ((qlj + qij) - criqc) * qij / (qlj + qij)

#                                     # It is very important to re-update 'qtu' and 'thlu'  at the upper
#                                     # interface after expelling condensate from cumulus updraft at the
#                                     # top interface of the layer. As mentioned above, this is a 'sink'
#                                     # of cumulus qt (or equivalently, a 'source' of environmentasl qt),
#                                     # not a regular convective'detrainment'.
#                                     if THIS_K == k:
#                                         qtu = qtu.at(K=k) - exql - exqi
#                                         thlu = (
#                                             thlu.at(K=k)
#                                             + (
#                                                 constants.MAPL_LATENT_HEAT_VAPORIZATION
#                                                 / exnifc0.at(K=k)
#                                                 / constants.MAPL_CP
#                                             )
#                                             * exql
#                                             + (
#                                                 constants.MAPL_LATENT_HEAT_SUBLIMATION
#                                                 / exnifc0.at(K=k)
#                                                 / constants.MAPL_CP
#                                             )
#                                             * exqi
#                                         )

#                                         # Expelled cloud condensate into the environment from the updraft.
#                                         # After all the calculation later, 'dwten' and 'diten' will have a
#                                         # unit of [ kg/kg/s ], because it is a tendency of qt. Restoration
#                                         # of 'dwten' and 'diten' to this correct unit through  multiplying
#                                         # 'umf(k)*g/dp0(k)' will be performed later after finally updating
#                                         # 'umf' using a 'rmaxfrac' constraint near the end of this updraft
#                                         # buoyancy sorting loop.

#                                         dwten = exql
#                                         diten = exqi
#                                 else:
#                                     if THIS_K == k:
#                                         dwten = 0.0
#                                         diten = 0.0

#                                 # Update 'thvu(k)' after detraining condensate from cumulus updraft.
#                                 thj, qvj, qlj, qij, qse, id_check = conden(
#                                     pifc0.at(K=k), thlu.at(K=k), qtu.at(K=k), ese, esx
#                                 )
#                                 if id_check == 1:
#                                     id_exit = True
#                                     umf_out[0, 0, 1] = 0.0
#                                     dcm_out = 0.0
#                                     qvten_out = 0.0
#                                     qlten_out = 0.0
#                                     qiten_out = 0.0
#                                     sten_out = 0.0
#                                     uten_out = 0.0
#                                     vten_out = 0.0
#                                     qrten_out = 0.0
#                                     qsten_out = 0.0
#                                     cufrc_out = 0.0
#                                     cush_inout = -1.0
#                                     qldet_out = 0.0
#                                     qidet_out = 0.0
#                                     qtflx_out[0, 0, 1] = 0.0
#                                     slflx_out[0, 0, 1] = 0.0
#                                     uflx_out[0, 0, 1] = 0.0
#                                     vflx_out[0, 0, 1] = 0.0
#                                     fer_out = constants.MAPL_UNDEF
#                                     fdr_out = constants.MAPL_UNDEF

#                                 if id_exit == False:
#                                     if THIS_K == k:
#                                         thvu = thj * (1.0 + zvir * qvj - qlj - qij)

#                                     # Calculate updraft vertical velocity at the upper interface.
#                                     # In order to calculate 'wtw' at the upper interface, we use
#                                     # 'wtw' at the lower interface. Note  'wtw'  is continuously
#                                     # updated as cumulus updraft rises.

#                                     bogbot = rbuoy * (
#                                         thvu.at(K=km1) / thvebot - 1.0
#                                     )  # Cloud buoyancy at base interface
#                                     bogtop = rbuoy * (
#                                         thvu.at(K=k) / thv0top.at(K=k) - 1.0
#                                     )  # Cloud buoyancy at top  interface

#                                     delbog = bogtop - bogbot
#                                     drage = fer.at(K=k) * (1.0 + rdrag)
#                                     expfac = exp(-2.0 * drage * dpe)

#                                     wtwb = wtw
#                                     if drage * dpe > 1.0e-3:
#                                         wtw = wtw * expfac + (
#                                             delbog
#                                             + (1.0 - expfac)
#                                             * (bogbot + delbog / (-2.0 * drage * dpe))
#                                         ) / (rhomid0j * drage)
#                                     else:
#                                         wtw = wtw + dpe * (bogbot + bogtop) / rhomid0j

#                                     # Force the plume rise at least to klfc of the undiluted plume.
#                                     # Because even the below is not complete, I decided not to include this.

#                                     # if( k .le. klfc ) then
#                                     #      wtw = max( 1.e-2, wtw )
#                                     # endif

#                                     # Repeat 'iter_xc' iteration loop until 'iter_xc = niter_xc'.
#                                     # Also treat the case even when wtw < 0 at the 'kpen' interface.

#                                     if wtw > 0.0:
#                                         thlue = 0.5 * (thlu.at(K=km1) + thlu.at(K=k))
#                                         qtue = 0.5 * (qtu.at(K=km1) + qtu.at(K=k))
#                                         wue = 0.5 * sqrt(max(wtwb + wtw, 0.0))
#                                     else:
#                                         iter_xc = (
#                                             niter_xc + 1
#                                         )  # Break out of iter_xc loop

#                                     # end iter_xc loop
#                     iter_xc += 1

#                     # Add the contribution of self-detrainment  to vertical variations of cumulus
#                     # updraft mass flux. The reason why we are trying to include self-detrainment
#                     # is as follows.  In current scheme,  vertical variation of updraft mass flux
#                     # is not fully consistent with the vertical variation of updraft vertical w.
#                     # For example, within a given layer, let's assume that  cumulus w is positive
#                     # at the base interface, while negative at the top interface. This means that
#                     # cumulus updraft cannot reach to the top interface of the layer. However,
#                     # cumulus updraft mass flux at the top interface is not zero according to the
#                     # vertical tendency equation of cumulus mass flux.   Ideally, cumulus updraft
#                     # mass flux at the top interface should be zero for this case. In order to
#                     # assures that cumulus updraft mass flux goes to zero when cumulus updraft
#                     # vertical velocity goes to zero, we are imposing self-detrainment term as
#                     # below by considering layer-mean cloud buoyancy and cumulus updraft vertical
#                     # velocity square at the top interface. Use of auto-detrainment term will  be
#                     # determined by setting 'use_self_detrain=.true.' in the parameter sentence.

#                     if use_self_detrain == 1:
#                         autodet = min(
#                             0.5
#                             * constants.MAPL_GRAV
#                             * (bogbot + bogtop)
#                             / (max(wtw, 0.0) + 1.0e-4),
#                             0.0,
#                         )
#                         if THIS_K == k:
#                             umf[0, 0, 1] = umf.at(K=k) * exp(
#                                 0.637 * (dpe / rhomid0j / constants.MAPL_GRAV) * autodet
#                             )

#                     if umf.at(K=k) == 0.0:
#                         wtw = -1.0

#                     # 'kbup' is the upper most layer in which cloud buoyancy  is positive
#                     # both at the base and top interface.  'kpen' is the upper most layer
#                     # up to cumulus can reach. Usually, 'kpen' is located higher than the
#                     # 'kbup'. Note we initialized these by 'kbup = krel' & 'kpen = krel'.
#                     # As explained before, it is possible that only 'kpen' is updated,
#                     # while 'kbup' keeps its initialization value. For this case, current
#                     # scheme will simply turns-off penetrative entrainment fluxes and use
#                     # normal buoyancy-sorting fluxes for 'kbup <= k <= kpen-1' interfaces,
#                     # in order to describe shallow continental cumulus convection.

#                     # if( bogbot .gt. 0. .and. bogtop .gt. 0. ) then
#                     # if( bogtop .gt. 0. ) then
#                     if bogtop > 0.0 and wtw > 0.0:
#                         kbup = k

#                     if wtw <= 0.0:
#                         kpen = k
#                         k = k0 + 1  # Break out of kloop

#                     if k <= (k0 - 1):
#                         if THIS_K == k:
#                             wu = sqrt(wtw)

#                         if wu.at(K=k) > 100.0:
#                             id_exit = True
#                             umf_out[0, 0, 1] = 0.0
#                             dcm_out = 0.0
#                             qvten_out = 0.0
#                             qlten_out = 0.0
#                             qiten_out = 0.0
#                             sten_out = 0.0
#                             uten_out = 0.0
#                             vten_out = 0.0
#                             qrten_out = 0.0
#                             qsten_out = 0.0
#                             cufrc_out = 0.0
#                             cush_inout = -1.0
#                             qldet_out = 0.0
#                             qidet_out = 0.0
#                             qtflx_out[0, 0, 1] = 0.0
#                             slflx_out[0, 0, 1] = 0.0
#                             uflx_out[0, 0, 1] = 0.0
#                             vflx_out[0, 0, 1] = 0.0
#                             fer_out = constants.MAPL_UNDEF
#                             fdr_out = constants.MAPL_UNDEF

#                         if id_exit == False:
#                             # Iteration end due to 'rmaxfrac' constraint

#                             # Calculate updraft fractional area at the upper interface and set upper
#                             # limit to 'ufrc' by 'rmaxfrac'. In order to keep the consistency  among
#                             # ['ufrc','umf','wu (or wtw)'], if ufrc is limited by 'rmaxfrac', either
#                             # 'umf' or 'wu' should be changed. Although both 'umf' and 'wu (wtw)' at
#                             # the current upper interface are used for updating 'umf' & 'wu'  at the
#                             # next upper interface, 'umf' is a passive variable not influencing  the
#                             # buoyancy sorting process in contrast to 'wtw'. This is a reason why we
#                             # adjusted 'umf' instead of 'wtw'. In turn we updated 'fdr' here instead
#                             # of 'fer',  which guarantees  that all previously updated thermodynamic
#                             # variables at the upper interface before applying 'rmaxfrac' constraint
#                             # are already internally consistent,  even though 'ufrc'  is  limited by
#                             # 'rmaxfrac'. Thus, we don't need to go through interation loop again.If
#                             # If we update 'fer' however, we should go through above iteration loop.

#                             rhoifc0j = pifc0.at(K=k) / (
#                                 constants.MAPL_RDRY
#                                 * 0.5
#                                 * (thv0bot.at(K=k + 1) + thv0top.at(K=k))
#                                 * exnifc0.at(K=k)
#                             )
#                             if THIS_K == k:
#                                 ufrc = umf.at(K=k) / (rhoifc0j * wu.at(K=k))
#                                 if ufrc.at(K=k) > rmaxfrac:
#                                     ufrc = rmaxfrac
#                                     umf = rmaxfrac * rhoifc0j * wu.at(K=k)
#                                     fdr = (
#                                         fer.at(K=k)
#                                         - log(umf.at(K=k) / umf.at(K=km1)) / dpe
#                                     )

#                             # Update environmental properties for at the mid-point of next
#                             # upper layer for use in buoyancy sorting.

#                             pe = pmid0.at(K=k + 1)
#                             dpe = dp0.at(K=k + 1)
#                             exne = exnmid0.at(K=k + 1)
#                             thvebot = thv0bot.at(K=k + 1)
#                             thle = thl0.at(K=k + 1)
#                             qte = qt0.at(K=k + 1)
#                             ue = u0.at(K=k + 1)
#                             ve = v0.at(K=k + 1)
#                             if dotransport == 1.0:
#                                 n = 0
#                                 while n < ncnst:
#                                     tre[0, 0, 0][n] = tr0.at(K=k + 1, ddim=[n])
#                                     n += 1

#                 k += 1  # end k loop

#             if id_exit == False:
#                 # Up to this point, we finished all of buoyancy sorting processes from the 'krel'
#                 # layer to 'kpen' layer: at the top interface of individual layers, we calculated
#                 # updraft and penetrative mass fluxes [ umf(k) & emf(k) = 0 ], updraft fractional
#                 # area [ ufrc(k) ],  updraft vertical velocity [ wu(k) ],  updraft  thermodynamic
#                 # variables [thlu(k),qtu(k),uu(k),vu(k),thvu(k)]. In the layer,we also calculated
#                 # fractional entrainment-detrainment rate [ fer(k), fdr(k) ], and detrainment ten
#                 # dency of water and ice from cumulus updraft [ dwten(k), diten(k) ]. In addition,
#                 # we updated and identified 'krel' and 'kpen' layer index, if any.  In the 'kpen'
#                 # layer, we calculated everything mentioned above except the 'wu(k)' and 'ufrc(k)'
#                 # since a real value of updraft vertical velocity is not defined at the kpen  top
#                 # interface (note 'ufrc' at the top interface of layer is calculated from 'umf(k)'
#                 # and 'wu(k)'). As mentioned before, special treatment is required when 'kbup' is
#                 # not updated and so 'kbup = krel'.

#                 # During the 'iter_scaleh' iteration loop, non-physical ( with non-zero values )
#                 # values can remain in the variable arrays above (also 'including' in case of wu
#                 # and ufrc at the top interface) the 'kpen' layer. This can happen when the kpen
#                 # layer index identified from the 'iter_scaleh = 1' iteration loop is located at
#                 # above the kpen layer index identified from   'iter_scaleh = 3' iteration loop.
#                 # Thus, in the following calculations, we should only use the values in each
#                 # variables only up to finally identified 'kpen' layer & 'kpen' interface except
#                 # 'wu' and 'ufrc' at the top interface of 'kpen' layer.    Note that in order to
#                 # prevent any problems due to these non-physical values, I re-initialized    the
#                 # values of [ umf(kpen:k0), emf(kpen:k0), dwten(kpen+1:k0), diten(kpen+1:k0),!
#                 # fer(kpen:k0), fdr(kpen+1:k0), ufrc(kpen:k0) ] to be zero after 'iter_scaleh'!
#                 # do loop.

#                 # Calculate 'ppen( < 0 )', updraft penetrative distance from the lower interface
#                 # of 'kpen' layer. Note that bogbot & bogtop at the 'kpen' layer either when fer
#                 # is zero or non-zero was already calculated above.
#                 # It seems that below qudarature solving formula is valid only when bogbot < 0.
#                 # Below solving equation is clearly wrong ! I should revise this !

#                 if drage == 0.0:
#                     aquad = (bogtop - bogbot) / (
#                         pifc0.at(K=kpen) - pifc0.at(K=kpen - 1)
#                     )
#                     bquad = 2.0 * bogbot
#                     cquad = -wu.at(K=kpen - 1) ** 2 * rhomid0j
#                     xc1, xc2, status = roots(aquad, bquad, cquad)
#                     if status == 0:
#                         if xc1 <= 0.0 and xc2 <= 0.0:
#                             ppen = max(xc1, xc2)
#                             ppen = min(0.0, max(-dp0.at(K=kpen), ppen))
#                         elif xc1 > 0.0 and xc2 > 0.0:
#                             ppen = -1 * dp0.at(K=kpen)
#                         else:
#                             ppen = min(xc1, xc2)
#                             ppen = min(0.0, max(-dp0.at(K=kpen), ppen))

#                     else:
#                         ppen = -1 * dp0(kpen)

#                 else:
#                     ppen = compute_ppen(
#                         wtwb, drage, bogbot, bogtop, rhomid0j, dp0.at(K=kpen)
#                     )

#                 # Re-calculate the amount of expelled condensate from cloud updraft
#                 # at the cumulus top. This is necessary for refined calculations of
#                 # bulk cloud microphysics at the cumulus top. Note that ppen < 0.
#                 # In the below, I explicitly calculate 'thlu_top' & 'qtu_top' by
#                 # using non-zero 'fer(kpen)'.

#                 if fer.at(K=kpen) * (-ppen) < 1.0e-4:
#                     thlu_top = thlu.at(K=kpen - 1) + (
#                         thl0.at(K=kpen)
#                         + ssthl0.at(K=kpen) * (-ppen) / 2.0
#                         - thlu.at(K=kpen - 1)
#                     ) * fer.at(K=kpen) * (-ppen)
#                     qtu_top = qtu.at(K=kpen - 1) + (
#                         qt0.at(K=kpen)
#                         + ssqt0.at(K=kpen) * (-ppen) / 2.0
#                         - qtu.at(K=kpen - 1)
#                     ) * fer.at(K=kpen) * (-ppen)
#                 else:
#                     thlu_top = (
#                         thl0.at(K=kpen)
#                         + ssthl0.at(K=kpen) / fer.at(K=kpen)
#                         - ssthl0.at(K=kpen) * (-ppen) / 2.0
#                     ) - (
#                         thl0(kpen)
#                         + ssthl0(kpen) * (-ppen) / 2.0
#                         - thlu(kpen - 1)
#                         + ssthl0(kpen) / fer(kpen)
#                     ) * exp(
#                         -fer(kpen) * (-ppen)
#                     )
#                     qtu_top = (
#                         qt0.at(K=kpen)
#                         + ssqt0.at(K=kpen) / fer.at(K=kpen)
#                         - ssqt0.at(K=kpen) * (-ppen) / 2.0
#                     ) - (
#                         qt0.at(K=kpen)
#                         + ssqt0.at(K=kpen) * (-ppen) / 2.0
#                         - qtu.at(K=kpen - 1)
#                         + ssqt0.at(K=kpen) / fer.at(K=kpen)
#                     ) * exp(
#                         -fer.at(K=kpen) * (-ppen)
#                     )

#                 thj, qvj, qlj, qij, qse, id_check = conden(
#                     pifc0.at(K=kpen - 1) + ppen, thlu_top, qtu_top, ese, esx
#                 )
#                 if id_check == 1:
#                     id_exit = True
#                     umf_out[0, 0, 1] = 0.0
#                     dcm_out = 0.0
#                     qvten_out = 0.0
#                     qlten_out = 0.0
#                     qiten_out = 0.0
#                     sten_out = 0.0
#                     uten_out = 0.0
#                     vten_out = 0.0
#                     qrten_out = 0.0
#                     qsten_out = 0.0
#                     cufrc_out = 0.0
#                     cush_inout = -1.0
#                     qldet_out = 0.0
#                     qidet_out = 0.0
#                     qtflx_out[0, 0, 1] = 0.0
#                     slflx_out[0, 0, 1] = 0.0
#                     uflx_out[0, 0, 1] = 0.0
#                     vflx_out[0, 0, 1] = 0.0
#                     fer_out = constants.MAPL_UNDEF
#                     fdr_out = constants.MAPL_UNDEF

#             if id_exit == False:
#                 p00 = 1e5
#                 rovcp = constants.MAPL_RGAS / constants.MAPL_CP

#                 exntop = ((pifc0.at(K=kpen - 1) + ppen) / p00) ** rovcp
#                 if (qlj + qij) > criqc:
#                     if THIS_K == kpen:
#                         dwten = ((qlj + qij) - criqc) * qlj / (qlj + qij)
#                         diten = ((qlj + qij) - criqc) * qij / (qlj + qij)
#                     qtu_top = qtu_top - dwten.at(K=kpen) - diten.at(K=kpen)
#                     thlu_top = (
#                         thlu_top
#                         + (
#                             constants.MAPL_LATENT_HEAT_VAPORIZATION
#                             / constants.MAPL_CP
#                             / exntop
#                         )
#                         * dwten(kpen)
#                         + (
#                             constants.MAPL_LATENT_HEAT_SUBLIMATION
#                             / constants.MAPL_CP
#                             / exntop
#                         )
#                         * diten(kpen)
#                     )
#                 else:
#                     if THIS_K == kpen:
#                         dwten = 0.0
#                         diten = 0.0

#                 # Calculate cumulus scale height as the top height that cumulus can reach.
#                 rhoifc0j = pifc0.at(K=kpen - 1) / (
#                     constants.MAPL_RDRY
#                     * 0.5
#                     * (thv0bot.at(K=kpen) + thv0top.at(K=kpen - 1))
#                     * exnifc0.at(K=kpen - 1)
#                 )
#                 cush = zifc0.at(K=kpen - 1) - ppen / rhoifc0j / constants.MAPL_GRAV
#                 scaleh = cush

#                 # The 'forcedCu' is logical identifier saying whether cumulus updraft
#                 # overcome the buoyancy barrier just above the PBL top. If it is true,
#                 # cumulus did not overcome the barrier -  this is a shallow convection
#                 # with negative cloud buoyancy, mimicking  shallow continental cumulus
#                 # convection. Depending on 'forcedCu' parameter, treatment of heat  &
#                 # moisture fluxes at the entraining interfaces, 'kbup <= k < kpen - 1'
#                 # will be set up in a different ways, as will be shown later.

#                 if kbup == krel:
#                     forcedCu = True
#                 else:
#                     forcedCu = False

#                 # Filtering of unerasonable cumulus adjustment here.  This is a very
#                 # important process which should be done cautiously. Various ways of
#                 # filtering are possible depending on cases mainly using the indices
#                 # of key layers - 'klcl','kinv','krel','klfc','kbup','kpen'. At this
#                 # stage, the followings are all possible : 'kinv >= 2', 'klcl >= 1',
#                 # 'krel >= kinv', 'kbup >= krel', 'kpen >= krel'. I must design this
#                 # filtering very cautiously, in such that none of  realistic cumulus
#                 # convection is arbitrarily turned-off. Potentially, I might turn-off
#                 # cumulus convection if layer-mean 'ql > 0' in the 'kinv-1' layer,in
#                 # order to suppress cumulus convection growing, based at the Sc top.
#                 # This is one of potential future modifications. Note that ppen < 0.

#                 cldhgt = pifc0.at(K=kpen - 1) + ppen
#                 if forcedCu == True:
#                     id_exit = True
#                     umf_out[0, 0, 1] = 0.0
#                     dcm_out = 0.0
#                     qvten_out = 0.0
#                     qlten_out = 0.0
#                     qiten_out = 0.0
#                     sten_out = 0.0
#                     uten_out = 0.0
#                     vten_out = 0.0
#                     qrten_out = 0.0
#                     qsten_out = 0.0
#                     cufrc_out = 0.0
#                     cush_inout = -1.0
#                     qldet_out = 0.0
#                     qidet_out = 0.0
#                     qtflx_out[0, 0, 1] = 0.0
#                     slflx_out[0, 0, 1] = 0.0
#                     uflx_out[0, 0, 1] = 0.0
#                     vflx_out[0, 0, 1] = 0.0
#                     fer_out = constants.MAPL_UNDEF
#                     fdr_out = constants.MAPL_UNDEF

#             if id_exit == False:
#                 # Re-initializing some key variables above the 'kpen' layer in order to suppress
#                 # the influence of non-physical values above 'kpen', in association with the use
#                 # of 'iter_scaleh' loop. Note that umf, emf,  ufrc are defined at the interfaces
#                 # (0:k0), while 'dwten','diten', 'fer', 'fdr' are defined at layer mid-points.
#                 # Initialization of 'fer' and 'fdr' is for correct writing purpose of diagnostic
#                 # output. Note that we set umf(kpen)=emf(kpen)=ufrc(kpen)=0, in consistent  with
#                 # wtw < 0  at the top interface of 'kpen' layer. However, we still have non-zero
#                 # expelled cloud condensate in the 'kpen' layer.
#                 if THIS_K >= kpen and THIS_K <= k0:
#                     umf[0, 0, 1] = 0.0
#                     emf[0, 0, 1] = 0.0
#                     ufrc[0, 0, 1] = 0.0
#                 if THIS_K >= kpen + 1 and THIS_K < k0:
#                     dwten = 0.0
#                     diten = 0.0
#                     fer = 0.0
#                     fdr = 0.0
#                     xco = 0.0

#                 # Calculate downward penetrative entrainment mass flux, 'emf(k) < 0',  and
#                 # thermodynamic properties of penetratively entrained airs at   entraining
#                 # interfaces. emf(k) is defined from the top interface of the  layer  kbup
#                 # to the bottom interface of the layer 'kpen'. Note even when  kbup = krel,
#                 # i.e.,even when 'kbup' was not updated in the above buoyancy  sorting  do
#                 # loop (i.e., 'kbup' remains as the initialization value),   below do loop
#                 # of penetrative entrainment flux can be performed without  any conceptual
#                 # or logical problems, because we have already computed all  the variables
#                 # necessary for performing below penetrative entrainment block.
#                 # In the below 'do' loop, 'k' is an interface index at which non-zero 'emf'
#                 # (penetrative entrainment mass flux) is calculated. Since cumulus updraft
#                 # is negatively buoyant in the layers between the top interface of 'kbup'
#                 # layer (interface index, kbup) and the top interface of 'kpen' layer, the
#                 # fractional lateral entrainment, fer(k) within these layers will be close
#                 # to zero - so it is likely that only strong lateral detrainment occurs in
#                 # thses layers. Under this situation,we can easily calculate the amount of
#                 # detrainment cumulus air into these negatively buoyanct layers by  simply
#                 # comparing cumulus updraft mass fluxes between the base and top interface
#                 # of each layer: emf(k) = emf(k-1)*exp(-fdr(k)*dp0(k))
#                 #                       ~ emf(k-1)*(1-rei(k)*dp0(k))
#                 #                emf(k-1)-emf(k) ~ emf(k-1)*rei(k)*dp0(k)
#                 # Current code assumes that about 'rpen~10' times of these detrained  mass
#                 # are penetratively re-entrained down into the 'k-1' interface. And all of
#                 # these detrained masses are finally dumped down into the top interface of
#                 # 'kbup' layer. Thus, the amount of penetratively entrained air across the
#                 # top interface of 'kbup' layer with 'rpen~10' becomes too large.
#                 # Note that this penetrative entrainment part can be completely turned-off
#                 # and we can simply use normal buoyancy-sorting involved turbulent  fluxes
#                 # by modifying 'penetrative entrainment fluxes' part below.

#                 # Calculate entrainment mass flux and conservative scalars of entraining
#                 # free air at interfaces of 'kbup <= k < kpen - 1'

#                 kk = 0
#                 while kk >= 0 and kk <= k0:
#                     if THIS_K == kk:
#                         thlu_emf = thlu.at(K=kk)
#                         qtu_emf = qtu.at(K=kk)
#                         uu_emf = uu.at(K=kk)
#                         vu_emf = vu.at(K=kk)
#                         if dotransport == 1.0:
#                             n = 0
#                             while n < ncnst:
#                                 tru_emf[0, 0, 0][n] = tru.at(K=kk, ddim=[n])
#                                 n += 1
#                     kk += 1

#     with computation(BACKWARD), interval(...):
#         if id_exit == False:
#             k = kpen - 1
#             while k <= kbup:  # Here, 'k' is an interface index at which
#                 rhoifc0j = pifc0.at(K=k) / (
#                     constants.MAPL_RDRY
#                     * 0.5
#                     * (thv0bot.at(K=k + 1) + thv0top.at(K=k))
#                     * exnifc0.at(K=k)
#                 )

#                 if k == (kpen - 1):
#                     # Note that 'ppen' has already been calculated in the above 'iter_scaleh'
#                     # loop assuming zero lateral entrainmentin the layer 'kpen'.

#                     # Calculate returning mass flux, emf ( < 0 )
#                     # Current penetrative entrainment rate with 'rpen~10' is too large and
#                     # future refinement is necessary including the definition of 'thl','qt'
#                     # of penetratively entrained air.  Penetratively entrained airs across
#                     # the 'kpen-1' interface is assumed to have the properties of the base
#                     # interface of 'kpen' layer. Note that 'emf ~ - umf/ufrc = - w * rho'.
#                     # Thus, below limit sets an upper limit of |emf| to be ~ 10cm/s, which
#                     # is very loose constraint. Here, I used more restricted constraint on
#                     # the limit of emf, assuming 'emf' cannot exceed a net mass within the
#                     # layer above the interface. Similar to the case of warming and drying
#                     # due to cumulus updraft induced compensating subsidence,  penetrative
#                     # entrainment induces compensating upwelling -     in order to prevent
#                     # numerical instability in association with compensating upwelling, we
#                     # should similarily limit the amount of penetrative entrainment at the
#                     # interface by the amount of masses within the layer just above the
#                     # penetratively entraining interface.

#                     if THIS_K == k:
#                         emf = max(
#                             max(
#                                 umf.at(K=k) * ppen * rei.at(K=kpen) * rpen,
#                                 -0.1 * rhoifc0j,
#                             ),
#                             -0.9 * dp0.at(K=kpen) / constants.MAPL_GRAV / dt,
#                         )
#                         thlu_emf = thl0.at(K=kpen) + ssthl0.at(K=kpen) * (
#                             pifc0.at(K=k) - pmid0.at(K=kpen)
#                         )
#                         qtu_emf = qt0.at(K=kpen) + ssqt0.at(K=kpen) * (
#                             pifc0.at(K=k) - pmid0.at(K=kpen)
#                         )
#                         uu_emf = u0.at(K=kpen) + ssu0.at(K=kpen) * (
#                             pifc0.at(K=k) - pmid0.at(K=kpen)
#                         )
#                         vu_emf = v0.at(K=kpen) + ssv0.at(K=kpen) * (
#                             pifc0.at(K=k) - pmid0.at(K=kpen)
#                         )
#                         if dotransport == 1.0:
#                             n = 0
#                             while n < ncnst:
#                                 tru_emf[0, 0, 0][n] = tr0.at(
#                                     K=kpen, ddim=[n]
#                                 ) + sstr0.at(K=kpen, ddim=[n]) * (
#                                     pifc0.at(K=k) - pmid0.at(K=kpen)
#                                 )
#                                 n += 1
#                 else:

#                     # Note we are coming down from the higher interfaces to the lower interfaces.
#                     # Also note that 'emf < 0'. So, below operation is a summing not subtracting.
#                     # In order to ensure numerical stability, I imposed a modified correct limit
#                     # of '-0.9*dp0(k+1)/g/dt' on emf(k).

#                     if (
#                         use_cumpenent == 1
#                     ):  # Original Cumulative Penetrative Entrainment
#                         if THIS_K == k:
#                             emf = max(
#                                 max(
#                                     emf.at(K=k + 1)
#                                     - umf.at(K=k)
#                                     * dp0.at(K=k + 1)
#                                     * rei.at(K=k + 1)
#                                     * rpen,
#                                     -0.1 * rhoifc0j,
#                                 ),
#                                 -0.9 * dp0.at(K=k + 1) / constants.MAPL_GRAV / dt,
#                             )
#                             if abs(emf.at(K=k)) > abs(emf.at(K=k + 1)):
#                                 thlu_emf = (
#                                     thlu_emf.at(K=k + 1) * emf.at(K=k + 1)
#                                     + thl0.at(K=k + 1) * (emf.at(K=k) - emf.at(K=k + 1))
#                                 ) / emf.at(K=k)
#                                 qtu_emf = (
#                                     qtu_emf.at(K=k + 1) * emf.at(K=k + 1)
#                                     + qt0.at(K=k + 1) * (emf.at(K=k) - emf.at(K=k + 1))
#                                 ) / emf.at(K=k)
#                                 uu_emf = (
#                                     uu_emf.at(K=k + 1) * emf.at(K=k + 1)
#                                     + u0.at(K=k + 1) * (emf.at(K=k) - emf.at(K=k + 1))
#                                 ) / emf.at(K=k)
#                                 vu_emf = (
#                                     vu_emf.at(K=k + 1) * emf.at(K=k + 1)
#                                     + v0.at(K=k + 1) * (emf.at(K=k) - emf.at(K=k + 1))
#                                 ) / emf.at(K=k)
#                                 if dotransport == 1.0:
#                                     n = 0
#                                     while n < ncnst:
#                                         tru_emf[0, 0, 0][n] = (
#                                             tru_emf.at(K=k + 1, ddim=[n])
#                                             * emf.at(K=k + 1)
#                                             + tr0.at(K=k + 1, ddim=[n])
#                                             * (emf.at(K=k) - emf.at(K=k + 1))
#                                         ) / emf.at(K=k)
#                                         n += 1
#                             else:
#                                 thlu_emf = thl0.at(K=k + 1)
#                                 qtu_emf = qt0.at(K=k + 1)
#                                 uu_emf = u0.at(K=k + 1)
#                                 vu_emf = v0.at(K=k + 1)
#                                 if dotransport == 1.0:
#                                     n = 0
#                                     while n < ncnst:
#                                         tru_emf[0, 0, 0][n] = tr0.at(K=k + 1, ddim=[n])
#                                         n += 1

#                     else:  # Alternative Non-Cumulative Penetrative Entrainment
#                         if THIS_K == k:
#                             emf = max(
#                                 max(
#                                     -umf(k) * dp0(k + 1) * rei(k + 1) * rpen,
#                                     -0.1 * rhoifc0j,
#                                 ),
#                                 -0.9 * dp0(k + 1) / g / dt,
#                             )
#                             thlu_emf = thl0.at(K=k + 1)
#                             qtu_emf = qt0.at(K=k + 1)
#                             uu_emf = u0.at(K=k + 1)
#                             vu_emf = v0.at(K=k + 1)
#                             if dotransport == 1.0:
#                                 n = 0
#                                 while n < ncnst:
#                                     tru_emf[0, 0, 0][n] = tr0.at(K=k + 1, ddim=[n])
#                                     n += 1

#                 # In this GCM modeling framework,  all what we should do is to calculate  heat
#                 # and moisture fluxes at the given geometrically-fixed height interfaces -  we
#                 # don't need to worry about movement of material height surface in association
#                 # with compensating subsidence or unwelling, in contrast to the bulk modeling.
#                 # In this geometrically fixed height coordinate system, heat and moisture flux
#                 # at the geometrically fixed height handle everything - a movement of material
#                 # surface is implicitly treated automatically. Note that in terms of turbulent
#                 # heat and moisture fluxes at model interfaces, both the cumulus updraft  mass
#                 # flux and penetratively entraining mass flux play the same role -both of them
#                 # warms and dries the 'kbup' layer, cools and moistens the 'kpen' layer,   and
#                 # cools and moistens any intervening layers between 'kbup' and 'kpen' layers.
#                 # It is important to note these identical roles on turbulent heat and moisture
#                 # fluxes of 'umf' and 'emf'.
#                 # When 'kbup' is a stratocumulus-topped PBL top interface,  increase of 'rpen'
#                 # is likely to strongly diffuse stratocumulus top interface,  resulting in the
#                 # reduction of cloud fraction. In this sense, the 'kbup' interface has a  very
#                 # important meaning and role : across the 'kbup' interface, strong penetrative
#                 # entrainment occurs, thus any sharp gradient properties across that interface
#                 # are easily diffused through strong mass exchange. Thus, an initialization of
#                 # 'kbup' (and also 'kpen') should be done very cautiously as mentioned before.
#                 # In order to prevent this stron diffusion for the shallow cumulus convection
#                 # based at the Sc top, it seems to be good to initialize 'kbup = krel', rather
#                 # that 'kbup = krel-1'.

#                 k += 1

#     with computation(FORWARD), interval(...):
#         # Compute turbulent heat, moisture, momentum flux at all interfaces
#         # 1. PBL fluxes :  0 <= k <= kinv - 1
#         # All the information necessary to reconstruct PBL
#         # height are passed to 'fluxbelowinv'.
#         if id_exit == False:
#             xsrc = qtsrc
#             xmean = qt0.at(K=kinv)
#             xtop = qt0.at(K=kinv + 1) + ssqt0.at(K=kinv + 1) * (
#                 pifc0.at(K=kinv) - pmid0.at(K=kinv + 1)
#             )
#             xbot = qt0.at(K=kinv - 1) + ssqt0.at(K=kinv - 1) * (
#                 pifc0.at(K=kinv - 1) - pmid0.at(K=kinv - 1)
#             )

#             xflx[0, 0, 1] = 0.0

#     with computation(PARALLEL), interval(...):
#         if id_exit == False:
#             xflx = 0.0
#             k_below = kinv - 1
#             dp = pifc0.at(K=k_below) - pifc0.at(K=kinv)

#             # Compute reconstructed inversion height
#             xtop_ori = xtop
#             xbot_ori = xbot
#             rcbmf = (
#                 cbmf * constants.MAPL_GRAV * dt
#             ) / dp  # Can be larger than 1 : 'OK'

#             if xbot >= xtop:
#                 rpeff = (xmean - xtop) / max(1.0e-20, xbot - xtop)
#             else:
#                 rpeff = (xmean - xtop) / min(-1.0e-20, xbot - xtop)

#             rpeff = min(max(0.0, rpeff), 1.0)  # As of this, 0<= rpeff <= 1
#             if rpeff == 0.0 or rpeff == 1.0:
#                 xbot = xmean
#                 xtop = xmean

#             # Below two commented-out lines are the old code replacing the above 'if' block.
#             # if(rpeff.eq.1) xbot = xmean
#             # if(rpeff.eq.0) xtop = xmean

#             rr = rpeff / rcbmf
#             pinv = pifc0.at(K=k_below) - rpeff * dp  # "pinv" before detraining mass
#             pinv_eff = (
#                 pifc0.at(K=k_below) + (rcbmf - rpeff) * dp
#             )  # Effective "pinv" after detraining mass

#             # Compute turbulent fluxes.
#             # Below two cases exactly converges at 'kinv-1' interface when rr = 1.
#             if THIS_K <= k_below:
#                 xflx = (
#                     cbmf
#                     * (xsrc - xbot)
#                     * (pifc0.at(K=0) - pifc0)
#                     / (pifc0.at(K=0) - pinv)
#                 )
#             if THIS_K == k_below:
#                 if rr <= 1:
#                     xflx = xflx - (1.0 - rr) * cbmf * (xtop_ori - xbot_ori)

#     with computation(FORWARD), interval(...):
#         if id_exit == False:
#             if THIS_K <= (kinv - 1):
#                 qtflx = xflx

#             xsrc = thlsrc
#             xmean = thl0.at(K=kinv)
#             xtop = thl0.at(K=kinv + 1) + ssthl0.at(K=kinv + 1) * (
#                 pifc0.at(K=kinv) - pmid0.at(K=kinv + 1)
#             )
#             xbot = thl0.at(K=kinv - 1) + ssthl0.at(K=kinv - 1) * (
#                 pifc0.at(K=kinv - 1) - pmid0.at(K=kinv - 1)
#             )
#             xflx[0, 0, 1] = 0.0

#     with computation(PARALLEL), interval(...):
#         if id_exit == False:
#             xflx = 0.0
#             k_below = kinv - 1
#             dp = pifc0.at(K=k_below) - pifc0.at(K=kinv)

#             # Compute reconstructed inversion height
#             xtop_ori = xtop
#             xbot_ori = xbot
#             rcbmf = (
#                 cbmf * constants.MAPL_GRAV * dt
#             ) / dp  # Can be larger than 1 : 'OK'

#             if xbot >= xtop:
#                 rpeff = (xmean - xtop) / max(1.0e-20, xbot - xtop)
#             else:
#                 rpeff = (xmean - xtop) / min(-1.0e-20, xbot - xtop)

#             rpeff = min(max(0.0, rpeff), 1.0)  # As of this, 0<= rpeff <= 1
#             if rpeff == 0.0 or rpeff == 1.0:
#                 xbot = xmean
#                 xtop = xmean

#             # Below two commented-out lines are the old code replacing the above 'if' block.
#             # if(rpeff.eq.1) xbot = xmean
#             # if(rpeff.eq.0) xtop = xmean

#             rr = rpeff / rcbmf
#             pinv = pifc0.at(K=k_below) - rpeff * dp  # "pinv" before detraining mass
#             pinv_eff = (
#                 pifc0.at(K=k_below) + (rcbmf - rpeff) * dp
#             )  # Effective "pinv" after detraining mass

#             # Compute turbulent fluxes.
#             # Below two cases exactly converges at 'kinv-1' interface when rr = 1.
#             if THIS_K <= k_below:
#                 xflx = (
#                     cbmf
#                     * (xsrc - xbot)
#                     * (pifc0.at(K=0) - pifc0)
#                     / (pifc0.at(K=0) - pinv)
#                 )
#             if THIS_K == k_below:
#                 if rr <= 1:
#                     xflx = xflx - (1.0 - rr) * cbmf * (xtop_ori - xbot_ori)

#     with computation(FORWARD), interval(...):
#         if id_exit == False:
#             if THIS_K <= (kinv - 1):
#                 slflx = constants.MAPL_CP * exnifc0 * xflx

#             xsrc = usrc
#             xmean = u0.at(K=kinv)
#             xtop = u0.at(K=kinv + 1) + ssu0.at(K=kinv + 1) * (
#                 pifc0.at(K=kinv) - pmid0.at(K=kinv + 1)
#             )
#             xbot = u0.at(K=kinv - 1) + ssu0.at(K=kinv - 1) * (
#                 pifc0.at(K=kinv - 1) - pmid0.at(K=kinv - 1)
#             )

#             xflx[0, 0, 1] = 0.0

#     with computation(PARALLEL), interval(...):
#         if id_exit == False:
#             xflx = 0.0
#             k_below = kinv - 1
#             dp = pifc0.at(K=k_below) - pifc0.at(K=kinv)

#             # Compute reconstructed inversion height
#             xtop_ori = xtop
#             xbot_ori = xbot
#             rcbmf = (
#                 cbmf * constants.MAPL_GRAV * dt
#             ) / dp  # Can be larger than 1 : 'OK'

#             if xbot >= xtop:
#                 rpeff = (xmean - xtop) / max(1.0e-20, xbot - xtop)
#             else:
#                 rpeff = (xmean - xtop) / min(-1.0e-20, xbot - xtop)

#             rpeff = min(max(0.0, rpeff), 1.0)  # As of this, 0<= rpeff <= 1
#             if rpeff == 0.0 or rpeff == 1.0:
#                 xbot = xmean
#                 xtop = xmean

#             # Below two commented-out lines are the old code replacing the above 'if' block.
#             # if(rpeff.eq.1) xbot = xmean
#             # if(rpeff.eq.0) xtop = xmean

#             rr = rpeff / rcbmf
#             pinv = pifc0.at(K=k_below) - rpeff * dp  # "pinv" before detraining mass
#             pinv_eff = (
#                 pifc0.at(K=k_below) + (rcbmf - rpeff) * dp
#             )  # Effective "pinv" after detraining mass

#             # Compute turbulent fluxes.
#             # Below two cases exactly converges at 'kinv-1' interface when rr = 1.
#             if THIS_K <= k_below:
#                 xflx = (
#                     cbmf
#                     * (xsrc - xbot)
#                     * (pifc0.at(K=0) - pifc0)
#                     / (pifc0.at(K=0) - pinv)
#                 )
#             if THIS_K == k_below:
#                 if rr <= 1:
#                     xflx = xflx - (1.0 - rr) * cbmf * (xtop_ori - xbot_ori)

#     with computation(FORWARD), interval(...):
#         if id_exit == False:
#             if THIS_K <= (kinv - 1):
#                 uflx = xflx

#             xsrc = vsrc
#             xmean = v0.at(K=kinv)
#             xtop = v0.at(K=kinv + 1) + ssv0.at(K=kinv + 1) * (
#                 pifc0.at(K=kinv) - pmid0.at(K=kinv + 1)
#             )
#             xbot = v0.at(K=kinv - 1) + ssv0.at(K=kinv - 1) * (
#                 pifc0.at(K=kinv - 1) - pmid0.at(K=kinv - 1)
#             )

#     with computation(PARALLEL), interval(...):
#         if id_exit == False:
#             xflx = 0.0
#             k_below = kinv - 1
#             dp = pifc0.at(K=k_below) - pifc0.at(K=kinv)

#             # Compute reconstructed inversion height
#             xtop_ori = xtop
#             xbot_ori = xbot
#             rcbmf = (
#                 cbmf * constants.MAPL_GRAV * dt
#             ) / dp  # Can be larger than 1 : 'OK'

#             if xbot >= xtop:
#                 rpeff = (xmean - xtop) / max(1.0e-20, xbot - xtop)
#             else:
#                 rpeff = (xmean - xtop) / min(-1.0e-20, xbot - xtop)

#             rpeff = min(max(0.0, rpeff), 1.0)  # As of this, 0<= rpeff <= 1
#             if rpeff == 0.0 or rpeff == 1.0:
#                 xbot = xmean
#                 xtop = xmean

#             # Below two commented-out lines are the old code replacing the above 'if' block.
#             # if(rpeff.eq.1) xbot = xmean
#             # if(rpeff.eq.0) xtop = xmean

#             rr = rpeff / rcbmf
#             pinv = pifc0.at(K=k_below) - rpeff * dp  # "pinv" before detraining mass
#             pinv_eff = (
#                 pifc0.at(K=k_below) + (rcbmf - rpeff) * dp
#             )  # Effective "pinv" after detraining mass

#             # Compute turbulent fluxes.
#             # Below two cases exactly converges at 'kinv-1' interface when rr = 1.
#             if THIS_K <= k_below:
#                 xflx = (
#                     cbmf
#                     * (xsrc - xbot)
#                     * (pifc0.at(K=0) - pifc0)
#                     / (pifc0.at(K=0) - pinv)
#                 )
#             if THIS_K == k_below:
#                 if rr <= 1:
#                     xflx = xflx - (1.0 - rr) * cbmf * (xtop_ori - xbot_ori)

#     with computation(FORWARD), interval(...):
#         if id_exit == False:
#             if THIS_K <= (kinv - 1):
#                 vflx = xflx

#             if dotransport == 1.0:
#                 n = 0
#                 while n < ncnst:
#                     xsrc = trsrc[0, 0, 0][n]
#                     xmean = tr0.at(K=kinv, ddim=[n])
#                     xtop = tr0.at(K=kinv + 1, ddim=[n]) + sstr0.at(
#                         K=kinv + 1, ddim=[n]
#                     ) * (pifc0.at(K=kinv) - pmid0.at(K=kinv + 1))
#                     xbot = tr0.at(K=kinv - 1, ddim=[n]) + sstr0.at(
#                         K=kinv - 1, ddim=[n]
#                     ) * (pifc0.at(K=kinv - 1) - pmid0.at(K=kinv - 1))
#                     # REVISIT THIS CALL TO FLUXBELOWINV
#                     # call fluxbelowinv( cbmf, pifc0(0:k0), k0, kinv, dt, xsrc, xmean, xtop, xbot, xflx )
#                     if THIS_K <= (kinv - 1):
#                         trflx[0, 0, 0][n] = xflx
#                     n += 1

#             # 2. Non-buoyancy sorting fluxes : kinv <= k <= krel - 1
#             # Note that when 'krel = kinv', below block is never executed
#             # as in a desirable, expected way ( but I must check  if this
#             # is the case ). The non-buoyancy sorting fluxes are computed
#             # only when 'krel > kinv'.

#             uplus = 0.0
#             vplus = 0.0
#             k = kinv
#             while k <= (krel - 1):
#                 kp1 = k + 1
#                 if THIS_K == k:
#                     qtflx = cbmf * (
#                         qtsrc
#                         - (
#                             qt0.at(K=kp1)
#                             + ssqt0.at(K=kp1) * (pifc0.at(K=k) - pmid0.at(K=kp1))
#                         )
#                     )
#                     slflx = (
#                         cbmf
#                         * (
#                             thlsrc
#                             - (
#                                 thl0.at(K=kp1)
#                                 + ssthl0.at(K=kp1) * (pifc0.at(K=k) - pmid0.at(K=kp1))
#                             )
#                         )
#                         * constants.MAPL_CP
#                         * exnifc0.at(K=k)
#                     )
#                 uplus = uplus + PGFc * ssu0.at(K=k) * (
#                     pifc0.at(K=k) - pifc0.at(K=k - 1)
#                 )
#                 vplus = vplus + PGFc * ssv0.at(K=k) * (
#                     pifc0.at(K=k) - pifc0.at(K=k - 1)
#                 )
#                 if THIS_K == k:
#                     uflx = cbmf * (
#                         usrc
#                         + uplus
#                         - (
#                             u0.at(K=kp1)
#                             + ssu0.at(K=kp1) * (pifc0.at(K=k) - pmid0.at(K=kp1))
#                         )
#                     )
#                     vflx = cbmf * (
#                         vsrc
#                         + vplus
#                         - (
#                             v0.at(K=kp1)
#                             + ssv0.at(K=kp1) * (pifc0.at(K=k) - pmid0.at(K=kp1))
#                         )
#                     )
#                 if dotransport == 1.0:
#                     n = 0
#                     while n < ncnst:
#                         if THIS_K == k:
#                             trflx[0, 0, 0][n] = cbmf * (
#                                 trsrc[0, 0, 0][n]
#                                 - (
#                                     tr0.at(K=kp1, ddim=[n])
#                                     + sstr0.at(K=kp1, ddim=[n])
#                                     * (pifc0.at(K=k) - pmid0.at(K=kp1))
#                                 )
#                             )
#                         n += 1
#                 k += 1

#             # 3. Buoyancy sorting fluxes : krel <= k <= kbup - 1
#             # In case that 'kbup = krel - 1 ' ( or even in case 'kbup = krel' ),
#             # buoyancy sorting fluxes are not calculated, which is consistent,
#             # desirable feature.

#             k = krel
#             while k <= (kbup - 1):
#                 if THIS_K == k:
#                     kp1 = k + 1
#                     slflx = (
#                         constants.MAPL_CP
#                         * exnifc0.at(K=k)
#                         * umf.at(K=k)
#                         * (
#                             thlu.at(K=k)
#                             - (
#                                 thl0.at(K=kp1)
#                                 + ssthl0.at(K=kp1) * (pifc0.at(K=k) - pmid0.at(K=kp1))
#                             )
#                         )
#                     )
#                     qtflx = umf.at(K=k) * (
#                         qtu.at(K=k)
#                         - (
#                             qt0.at(K=kp1)
#                             + ssqt0.at(K=kp1) * (pifc0.at(K=k) - pmid0.at(K=kp1))
#                         )
#                     )
#                     uflx = umf.at(K=k) * (
#                         uu.at(K=k)
#                         - (
#                             u0.at(K=kp1)
#                             + ssu0.at(K=kp1) * (pifc0.at(K=k) - pmid0.at(K=kp1))
#                         )
#                     )
#                     vflx = umf.at(K=k) * (
#                         vu.at(K=k)
#                         - (
#                             v0.at(K=kp1)
#                             + ssv0.at(K=kp1) * (pifc0.at(K=k) - pmid0.at(K=kp1))
#                         )
#                     )
#                     if dotransport == 1.0:
#                         n = 0
#                         while n < ncnst:
#                             trflx[0, 0, 0][n] = umf.at(K=k) * (
#                                 tru.at(K=k, ddim=[n])
#                                 - (
#                                     tr0.at(K=kp1, ddim=[n])
#                                     + sstr0.at(K=kp1, ddim=[n])
#                                     * (pifc0.at(K=k) - pmid0.at(K=kp1))
#                                 )
#                             )
#                             n += 1

#             # 4. Penetrative entrainment fluxes : kbup <= k <= kpen - 1
#             # The only confliction that can happen is when 'kbup = kinv-1'. For this
#             # case, turbulent flux at kinv-1 is calculated  both from 'fluxbelowinv'
#             # and here as penetrative entrainment fluxes.  Since penetrative flux is
#             # calculated later, flux at 'kinv - 1 ' will be that of penetrative flux.
#             # However, turbulent flux calculated at 'kinv - 1' from penetrative entr.
#             # is less attractable,  since more reasonable turbulent flux at 'kinv-1'
#             # should be obtained from 'fluxbelowinv', by considering  re-constructed
#             # inversion base height. This conflicting problem can be solved if we can
#             # initialize 'kbup = krel', instead of kbup = krel - 1. This choice seems
#             # to be more reasonable since it is not conflicted with 'fluxbelowinv' in
#             # calculating fluxes at 'kinv - 1' ( for this case, flux at 'kinv-1' is
#             # always from 'fluxbelowinv' ), and flux at 'krel-1' is calculated from
#             # the non-buoyancy sorting flux without being competed with penetrative
#             # entrainment fluxes. Even when we use normal cumulus flux instead of
#             # penetrative entrainment fluxes at 'kbup <= k <= kpen-1' interfaces,
#             # the initialization of kbup=krel perfectly works without any conceptual
#             # confliction. Thus it seems to be much better to choose 'kbup = krel'
#             # initialization of 'kbup', which is current choice.
#             # Note that below formula uses conventional updraft cumulus fluxes for
#             # shallow cumulus which did not overcome the first buoyancy barrier above
#             # PBL top while uses penetrative entrainment fluxes for the other cases
#             # 'kbup <= k <= kpen-1' interfaces. Depending on cases, however, I can
#             # selelct different choice.
#             k = kbup
#             while k <= (kpen - 1):
#                 if THIS_K == k:
#                     kp1 = k + 1
#                     slflx = (
#                         constants.MAPL_CP
#                         * exnifc0.at(K=k)
#                         * emf.at(K=k)
#                         * (
#                             thlu_emf.at(K=k)
#                             - (
#                                 thl0.at(K=k)
#                                 + ssthl0.at(K=k) * (pifc0.at(K=k) - pmid0.at(K=k))
#                             )
#                         )
#                     )
#                     qtflx = emf.at(K=k) * (
#                         qtu_emf.at(K=k)
#                         - (
#                             qt0.at(K=k)
#                             + ssqt0.at(K=k) * (pifc0.at(K=k) - pmid0.at(K=k))
#                         )
#                     )
#                     uflx = emf.at(K=k) * (
#                         uu_emf.at(K=k)
#                         - (u0.at(K=k) + ssu0.at(K=k) * (pifc0.at(K=k) - pmid0.at(K=k)))
#                     )
#                     vflx = emf.at(K=k) * (
#                         vu_emf.at(K=k)
#                         - (v0.at(K=k) + ssv0.at(K=k) * (pifc0.at(K=k) - pmid0.at(K=k)))
#                     )
#                     if dotransport == 1.0:
#                         n = 0
#                         while n < ncnst:
#                             trflx[0, 0, 0][n] = emf.at(K=k) * (
#                                 tru_emf.at(K=k, ddim=[n])
#                                 - (
#                                     tr0.at(K=k, ddim=[n])
#                                     + sstr0.at(K=k, ddim=[n])
#                                     * (pifc0.at(K=k) - pmid0.at(K=k))
#                                 )
#                             )
#                 k += 1

#             # Turn-off cumulus momentum flux as an option
#             if use_momenflx == 0:
#                 uflx[0, 0, 1] = 0.0
#                 vflx[0, 0, 1] = 0.0

#             # Condensate tendency by compensating subsidence/upwelling
#             uemf[0, 0, 1] = 0.0

#             k = 0
#             while k <= (kinv - 2):  # Assume linear updraft mass flux within the PBL.
#                 if THIS_K == k:
#                     uemf = (
#                         cbmf
#                         * (pifc0.at(K=0) - pifc0.at(K=k))
#                         / (pifc0.at(K=0) - pifc0.at(K=kinv - 1))
#                     )
#                 k += 1

#             if THIS_K >= (kinv - 1) and THIS_K <= (krel - 1):
#                 uemf = cbmf
#             if THIS_K >= krel and THIS_K <= (kbup - 1):
#                 uemf = umf
#             if THIS_K >= kbup and THIS_K <= (kpen - 1):
#                 uemf = emf  # Only use penetrative entrainment flux consistently.

#             if THIS_K > 0:
#                 comsub = 0.0

#             k = 0
#             while k <= kpen:
#                 if THIS_K == k:
#                     comsub = 0.5 * (
#                         uemf.at(K=k) + uemf.at(K=k - 1)
#                     )  # comsub defined on interfaces
#                 k += 1

#             k = 0
#             while k <= kpen:
#                 if id_exit == False:
#                     if comsub.at(K=k) >= 0.0:
#                         if k == k0:
#                             thlten_sub = 0.0
#                             qtten_sub = 0.0
#                             qlten_sub = 0.0
#                             qiten_sub = 0.0
#                             nlten_sub = 0.0
#                             niten_sub = 0.0
#                         else:
#                             thlten_sub = (
#                                 constants.MAPL_GRAV
#                                 * comsub.at(K=k)
#                                 * (thl0.at(K=k + 1) - thl0.at(K=k))
#                                 / (pmid0.at(K=k) - pmid0.at(K=k + 1))
#                             )
#                             qtten_sub = (
#                                 constants.MAPL_GRAV
#                                 * comsub.at(K=k)
#                                 * (qt0.at(K=k + 1) - qt0.at(K=k))
#                                 / (pmid0.at(K=k) - pmid0.at(K=k + 1))
#                             )
#                             qlten_sub = (
#                                 constants.MAPL_GRAV
#                                 * comsub.at(K=k)
#                                 * (ql0.at(K=k + 1) - ql0.at(K=k))
#                                 / (pmid0.at(K=k) - pmid0.at(K=k + 1))
#                             )
#                             qiten_sub = (
#                                 constants.MAPL_GRAV
#                                 * comsub.at(K=k)
#                                 * (qi0.at(K=k + 1) - qi0.at(K=k))
#                                 / (pmid0.at(K=k) - pmid0.at(K=k + 1))
#                             )
#                             # nlten_sub  = g * comsub(k) * (  tr0(k+1,ixnumliq) -  tr0(k,ixnumliq) ) / ( pmid0(k) - pmid0(k+1) )
#                             # niten_sub  = g * comsub(k) * (  tr0(k+1,ixnumice) -  tr0(k,ixnumice) ) / ( pmid0(k) - pmid0(k+1) )

#                     else:
#                         if k == 1:
#                             thlten_sub = 0.0
#                             qtten_sub = 0.0
#                             qlten_sub = 0.0
#                             qiten_sub = 0.0
#                             nlten_sub = 0.0
#                             niten_sub = 0.0
#                         else:
#                             thlten_sub = (
#                                 constants.MAPL_GRAV
#                                 * comsub.at(K=k)
#                                 * (thl0.at(K=k) - thl0.at(K=k - 1))
#                                 / (pmid0.at(K=k - 1) - pmid0.at(K=k))
#                             )
#                             qtten_sub = (
#                                 constants.MAPL_GRAV
#                                 * comsub.at(K=k)
#                                 * (qt0.at(K=k) - qt0.at(K=k - 1))
#                                 / (pmid0.at(K=k - 1) - pmid0.at(K=k))
#                             )
#                             qlten_sub = (
#                                 constants.MAPL_GRAV
#                                 * comsub.at(K=k)
#                                 * (ql0.at(K=k) - ql0.at(K=k - 1))
#                                 / (pmid0.at(K=k - 1) - pmid0.at(K=k))
#                             )
#                             qiten_sub = (
#                                 constants.MAPL_GRAV
#                                 * comsub.at(K=k)
#                                 * (qi0.at(K=k) - qi0.at(K=k - 1))
#                                 / (pmid0.at(K=k - 1) - pmid0.at(K=k))
#                             )
#                             # nlten_sub  = g * comsub(k) * (  tr0(k,ixnumliq) -  tr0(k-1,ixnumliq) ) / ( pmid0(k-1) - pmid0(k) )
#                             # niten_sub  = g * comsub(k) * (  tr0(k,ixnumice) -  tr0(k-1,ixnumice) ) / ( pmid0(k-1) - pmid0(k) )

#                     thl_prog = thl0.at(K=k) + thlten_sub * dt
#                     qt_prog = max(qt0.at(K=k) + qtten_sub * dt, 1.0e-12)
#                     thj, qvj, qlj, qij, qse, id_check = conden(
#                         pmid0.at(K=k), thl_prog, qt_prog, ese, esx
#                     )
#                     if id_check == 1:
#                         id_exit = True
#                         umf_out[0, 0, 1] = 0.0
#                         dcm_out = 0.0
#                         qvten_out = 0.0
#                         qlten_out = 0.0
#                         qiten_out = 0.0
#                         sten_out = 0.0
#                         uten_out = 0.0
#                         vten_out = 0.0
#                         qrten_out = 0.0
#                         qsten_out = 0.0
#                         cufrc_out = 0.0
#                         cush_inout = -1.0
#                         qldet_out = 0.0
#                         qidet_out = 0.0
#                         qtflx_out[0, 0, 1] = 0.0
#                         slflx_out[0, 0, 1] = 0.0
#                         uflx_out[0, 0, 1] = 0.0
#                         vflx_out[0, 0, 1] = 0.0
#                         fer_out = constants.MAPL_UNDEF
#                         fdr_out = constants.MAPL_UNDEF

#                 if id_exit == False:
#                     if THIS_K == k:
#                         qlten_sink = max(
#                             qlten_sub, -ql0.at(K=k) / dt
#                         )  # For consistency with prognostic macrophysics scheme
#                         qiten_sink = max(
#                             qiten_sub, -qi0.at(K=k) / dt
#                         )  # For consistency with prognostic macrophysics scheme

#                 k += 1

#     with computation(FORWARD), interval(...):
#         if id_exit == False:
#             # Calculate convective tendencies at each layer
#             # Momentum tendency
#             k = 0
#             while k <= kpen:
#                 km1 = k - 1
#                 if THIS_K == k:
#                     uten = (
#                         (uflx.at(K=km1) - uflx.at(K=k))
#                         * constants.MAPL_GRAV
#                         / dp0.at(K=k)
#                     )
#                     vten = (
#                         (vflx.at(K=km1) - vflx.at(K=k))
#                         * constants.MAPL_GRAV
#                         / dp0.at(K=k)
#                     )
#                     uf = u0.at(K=k) + uten.at(K=k) * dt
#                     vf = v0.at(K=k) + vten.at(K=k) * dt
#                 k += 1

#             # Tendencies of thermodynamic variables.
#             # This part requires a careful treatment of bulk cloud microphysics.
#             # Relocations of 'precipitable condensates' either into the surface
#             # or into the tendency of 'krel' layer will be performed just after
#             # finishing the below 'do-loop'.

#             k = 0
#             while k <= kpen:
#                 # Compute 'slten', 'qtten', 'qvten', 'qlten', 'qiten', and 'sten'

#                 # Key assumptions made in this 'cumulus scheme' are :
#                 # 1. Cumulus updraft expels condensate into the environment at the top interface
#                 #    of each layer. Note that in addition to this expel process ('source' term),
#                 #    cumulus updraft can modify layer mean condensate through normal detrainment
#                 #    forcing or compensating subsidence.
#                 # 2. Expelled water can be either 'sustaining' or 'precipitating' condensate. By
#                 #    definition, 'suataining condensate' will remain in the layer where it was
#                 #    formed, while 'precipitating condensate' will fall across the base of the
#                 #    layer where it was formed.
#                 # 3. All precipitating condensates are assumed to fall into the release layer or
#                 #    ground as soon as it was formed without being evaporated during the falling
#                 #    process down to the desinated layer ( either release layer of surface ).

#                 # 'dwten(k)','diten(k)' : Production rate of condensate  within the layer k
#                 #      [ kg/kg/s ]        by the expels of condensate from cumulus updraft.
#                 # It is important to note that in terms of moisture tendency equation, this
#                 # is a 'source' term of enviromental 'qt'.  More importantly,  these source
#                 # are already counted in the turbulent heat and moisture fluxes we computed
#                 # until now, assuming all the expelled condensate remain in the layer where
#                 # it was formed. Thus, in calculation of 'qtten' and 'slten' below, we MUST
#                 # NOT add or subtract these terms explicitly in order not to double or miss
#                 # count, unless some expelled condensates fall down out of the layer.  Note
#                 # this falling-down process ( i.e., precipitation process ) and  associated
#                 # 'qtten' and 'slten' and production of surface precipitation flux  will be
#                 # treated later in 'zm_conv_evap' in 'convect_shallow_tend' subroutine.
#                 # In below, we are converting expelled cloud condensate into correct unit.
#                 # I found that below use of '0.5 * (umf(k-1) + umf(k))' causes conservation
#                 # errors at some columns in global simulation. So, I returned to originals.
#                 # This will cause no precipitation flux at 'kpen' layer since umf(kpen)=0.
#                 km1 = k - 1
#                 if THIS_K == k:
#                     dwten = (
#                         dwten.at(K=k)
#                         * 0.5
#                         * (umf.at(K=k - 1) + umf.at(K=k))
#                         * constants.MAPL_GRAV
#                         / dp0.at(K=k)
#                     )  # [ kg/kg/s ]
#                     diten = (
#                         diten.at(K=k)
#                         * 0.5
#                         * (umf.at(K=k - 1) + umf.at(K=k))
#                         * constants.MAPL_GRAV
#                         / dp0.at(K=k)
#                     )  # [ kg/kg/s ]

#                     # 'qrten(k)','qsten(k)' : Production rate of rain and snow within the layer k
#                     # [ kg/kg/s ]         by cumulus expels of condensates to the environment.
#                     qrten = frc_rasn * dwten.at(K=k)
#                     qsten = frc_rasn * diten.at(K=k)

#                 # 'slten(k)','qtten(k)'
#                 # Note that 'slflx(k)' and 'qtflx(k)' we have calculated already included
#                 # all the contributions of (1) expels of condensate (dwten(k), diten(k)),
#                 # (2) mass detrainment ( delta * umf * ( qtu - qt ) ), & (3) compensating
#                 # subsidence ( M * dqt / dz ). Thus 'slflx(k)' and 'qtflx(k)' we computed
#                 # is a hybrid turbulent flux containing one part of 'source' term - expel
#                 # of condensate. In order to calculate 'slten' and 'qtten', we should add
#                 # additional 'source' term, if any. If the expelled condensate falls down
#                 # across the base of the layer, it will be another sink (negative source)
#                 # term.  Note also that we included frictional heating terms in the below
#                 # alculation of 'slten'.
#                 if THIS_K == k:
#                     slten = (
#                         (slflx.at(K=km1) - slflx.at(K=k))
#                         * constants.MAPL_GRAV
#                         / dp0.at(K=k)
#                     )
#                 if k == 1:
#                     if THIS_K == k:
#                         slten = slten.at(K=k) - constants.MAPL_GRAV / 4.0 / dp0.at(
#                             K=k
#                         ) * (
#                             uflx.at(K=k)
#                             * (
#                                 uf.at(K=k + 1)
#                                 - uf.at(K=k)
#                                 + u0.at(K=k + 1)
#                                 - u0.at(K=k)
#                             )
#                             + vflx.at(K=k)
#                             * (
#                                 vf.at(K=k + 1)
#                                 - vf.at(K=k)
#                                 + v0.at(K=k + 1)
#                                 - v0.at(K=k)
#                             )
#                         )
#                 elif k >= 2 and k <= (kpen - 1):
#                     if THIS_K == k:
#                         slten = slten.at(K=k) - constants.MAPL_GRAV / 4.0 / dp0.at(
#                             K=k
#                         ) * (
#                             uflx.at(K=k)
#                             * (
#                                 uf.at(K=k + 1)
#                                 - uf.at(K=k)
#                                 + u0.at(K=k + 1)
#                                 - u0.at(K=k)
#                             )
#                             + uflx.at(K=k - 1)
#                             * (
#                                 uf.at(K=k)
#                                 - uf.at(K=k - 1)
#                                 + u0.at(K=k)
#                                 - u0.at(K=k - 1)
#                             )
#                             + vflx.at(K=k)
#                             * (
#                                 vf.at(K=k + 1)
#                                 - vf.at(K=k)
#                                 + v0.at(K=k + 1)
#                                 - v0.at(K=k)
#                             )
#                             + vflx.at(K=k - 1)
#                             * (
#                                 vf.at(K=k)
#                                 - vf.at(K=k - 1)
#                                 + v0.at(K=k)
#                                 - v0.at(K=k - 1)
#                             )
#                         )
#                 elif k == kpen:
#                     if THIS_K == k:
#                         slten = slten.at(K=k) - constants.MAPL_GRAV / 4.0 / dp0.at(
#                             K=k
#                         ) * (
#                             uflx.at(K=k - 1)
#                             * (
#                                 uf.at(K=k)
#                                 - uf.at(K=k - 1)
#                                 + u0.at(K=k)
#                                 - u0.at(K=k - 1)
#                             )
#                             + vflx.at(K=k - 1)
#                             * (
#                                 vf.at(K=k)
#                                 - vf.at(K=k - 1)
#                                 + v0.at(K=k)
#                                 - v0.at(K=k - 1)
#                             )
#                         )
#                 if THIS_K == k:
#                     qtten = (
#                         (qtflx.at(K=km1) - qtflx.at(K=k))
#                         * constants.MAPL_GRAV
#                         / dp0.at(K=k)
#                     )

#                 # Compute condensate tendency, including reserved condensate
#                 # We assume that eventual detachment and detrainment occurs in kbup layer  due
#                 # to downdraft buoyancy sorting. In the layer above the kbup, only penetrative
#                 # entrainment exists. Penetrative entrained air is assumed not to contain any
#                 # condensate.

#                 # Compute in-cumulus condensate at the layer mid-point.
#                 if k < krel or k > kpen:
#                     qlu_mid = 0.0
#                     qiu_mid = 0.0
#                     qlj = 0.0
#                     qij = 0.0
#                 elif k == krel:
#                     thj, qvj, qlj, qij, qse, id_check = conden(
#                         prel, thlu.at(K=krel - 1), qtu.at(K=krel - 1), ese, esx
#                     )
#                     if id_check == 1:
#                         id_exit = True
#                         umf_out[0, 0, 1] = 0.0
#                         dcm_out = 0.0
#                         qvten_out = 0.0
#                         qlten_out = 0.0
#                         qiten_out = 0.0
#                         sten_out = 0.0
#                         uten_out = 0.0
#                         vten_out = 0.0
#                         qrten_out = 0.0
#                         qsten_out = 0.0
#                         cufrc_out = 0.0
#                         cush_inout = -1.0
#                         qldet_out = 0.0
#                         qidet_out = 0.0
#                         qtflx_out[0, 0, 1] = 0.0
#                         slflx_out[0, 0, 1] = 0.0
#                         uflx_out[0, 0, 1] = 0.0
#                         vflx_out[0, 0, 1] = 0.0
#                         fer_out = constants.MAPL_UNDEF
#                         fdr_out = constants.MAPL_UNDEF
#                     if id_exit == False:
#                         qlubelow = qlj
#                         qiubelow = qij
#                         thj, qvj, qlj, qij, qse, id_check = conden(
#                             pifc0.at(K=k), thlu.at(K=k), qtu.at(K=k), ese, esx
#                         )
#                         if id_check == 1:
#                             id_exit = True
#                             umf_out[0, 0, 1] = 0.0
#                             dcm_out = 0.0
#                             qvten_out = 0.0
#                             qlten_out = 0.0
#                             qiten_out = 0.0
#                             sten_out = 0.0
#                             uten_out = 0.0
#                             vten_out = 0.0
#                             qrten_out = 0.0
#                             qsten_out = 0.0
#                             cufrc_out = 0.0
#                             cush_inout = -1.0
#                             qldet_out = 0.0
#                             qidet_out = 0.0
#                             qtflx_out[0, 0, 1] = 0.0
#                             slflx_out[0, 0, 1] = 0.0
#                             uflx_out[0, 0, 1] = 0.0
#                             vflx_out[0, 0, 1] = 0.0
#                             fer_out = constants.MAPL_UNDEF
#                             fdr_out = constants.MAPL_UNDEF
#                     if id_exit == False:
#                         qlu_mid = (
#                             0.5
#                             * (qlubelow + qlj)
#                             * (prel - pifc0.at(K=k))
#                             / (pifc0.at(K=k - 1) - pifc0.at(K=k))
#                         )
#                         qiu_mid = (
#                             0.5
#                             * (qiubelow + qij)
#                             * (prel - pifc0.at(K=k))
#                             / (pifc0.at(K=k - 1) - pifc0.at(K=k))
#                         )
#                 elif k == kpen:
#                     if id_exit == False:
#                         thj, qvj, qlj, qij, qse, id_check = conden(
#                             pifc0.at(K=k - 1) + ppen, thlu_top, qtu_top, ese, esx
#                         )
#                         if id_check == 1:
#                             id_exit = True
#                             umf_out[0, 0, 1] = 0.0
#                             dcm_out = 0.0
#                             qvten_out = 0.0
#                             qlten_out = 0.0
#                             qiten_out = 0.0
#                             sten_out = 0.0
#                             uten_out = 0.0
#                             vten_out = 0.0
#                             qrten_out = 0.0
#                             qsten_out = 0.0
#                             cufrc_out = 0.0
#                             cush_inout = -1.0
#                             qldet_out = 0.0
#                             qidet_out = 0.0
#                             qtflx_out[0, 0, 1] = 0.0
#                             slflx_out[0, 0, 1] = 0.0
#                             uflx_out[0, 0, 1] = 0.0
#                             vflx_out[0, 0, 1] = 0.0
#                             fer_out = constants.MAPL_UNDEF
#                             fdr_out = constants.MAPL_UNDEF
#                     if id_exit == False:
#                         qlu_mid = (
#                             0.5
#                             * (qlubelow + qlj)
#                             * (-ppen)
#                             / (pifc0.at(K=k - 1) - pifc0.at(K=k))
#                         )
#                         qiu_mid = (
#                             0.5
#                             * (qiubelow + qij)
#                             * (-ppen)
#                             / (pifc0.at(K=k - 1) - pifc0.at(K=k))
#                         )
#                         qlu_top = qlj
#                         qiu_top = qij
#                 else:
#                     if id_exit == False:
#                         thj, qvj, qlj, qij, qse, id_check = conden(
#                             pifc0.at(K=k), thlu.at(K=k), qtu.at(K=k), ese, esx
#                         )
#                         if id_check == 1:
#                             id_exit = True
#                             umf_out[0, 0, 1] = 0.0
#                             dcm_out = 0.0
#                             qvten_out = 0.0
#                             qlten_out = 0.0
#                             qiten_out = 0.0
#                             sten_out = 0.0
#                             uten_out = 0.0
#                             vten_out = 0.0
#                             qrten_out = 0.0
#                             qsten_out = 0.0
#                             cufrc_out = 0.0
#                             cush_inout = -1.0
#                             qldet_out = 0.0
#                             qidet_out = 0.0
#                             qtflx_out[0, 0, 1] = 0.0
#                             slflx_out[0, 0, 1] = 0.0
#                             uflx_out[0, 0, 1] = 0.0
#                             vflx_out[0, 0, 1] = 0.0
#                             fer_out = constants.MAPL_UNDEF
#                             fdr_out = constants.MAPL_UNDEF
#                     if id_exit == False:
#                         qlu_mid = 0.5 * (qlubelow + qlj)
#                         qiu_mid = 0.5 * (qiubelow + qij)

#                 if id_exit == False:
#                     qlubelow = qlj
#                     qiubelow = qij

#                 # 1. Non-precipitating portion of expelled condensate
#                 if THIS_K == k:
#                     qc_l = (1.0 - frc_rasn) * dwten.at(K=k)  # [ kg/kg/s ]
#                     qc_i = (1.0 - frc_rasn) * diten.at(K=k)  # [ kg/kg/s ]

#                     # 2. Detrained Condensate
#                     if k <= kbup:
#                         if THIS_K == k:
#                             qc_l = (
#                                 qc_l.at(K=k)
#                                 + constants.MAPL_GRAV
#                                 * 0.5
#                                 * (umf.at(K=k - 1) + umf.at(K=k))
#                                 * fdr.at(K=k)
#                                 * qlu_mid
#                             )  # [ kg/kg/s ]
#                             qc_i = (
#                                 qc_i.at(K=k)
#                                 + constants.MAPL_GRAV
#                                 * 0.5
#                                 * (umf.at(K=k - 1) + umf.at(K=k))
#                                 * fdr.at(K=k)
#                                 * qiu_mid
#                             )  # [ kg/kg/s ]
#                         qc_lm = (
#                             -constants.MAPL_GRAV
#                             * 0.5
#                             * (umf.at(K=k - 1) + umf.at(K=k))
#                             * fdr.at(K=k)
#                             * ql0.at(K=k)
#                         )
#                         qc_im = (
#                             -constants.MAPL_GRAV
#                             * 0.5
#                             * (umf.at(K=k - 1) + umf.at(K=k))
#                             * fdr.at(K=k)
#                             * qi0.at(K=k)
#                         )
#                         # Below 'nc_lm', 'nc_im' should be used only when frc_rasn = 1.
#                         # nc_lm   =         - g * 0.5 * ( umf(k-1) + umf(k) ) * fdr(k) * tr0(k,ixnumliq)
#                         # nc_im   =         - g * 0.5 * ( umf(k-1) + umf(k) ) * fdr(k) * tr0(k,ixnumice)
#                     else:
#                         qc_lm = 0.0
#                         qc_im = 0.0
#                         nc_lm = 0.0
#                         nc_im = 0.0

#                     # 3. Detached Updraft
#                     if k == kbup:
#                         if THIS_K == k:
#                             qc_l = qc_l.at(K=k) + constants.MAPL_GRAV * umf.at(
#                                 K=k
#                             ) * qlj / (
#                                 pifc0.at(K=k - 1) - pifc0.at(K=k)
#                             )  # [ kg/kg/s ]
#                             qc_i = qc_i.at(K=k) + constants.MAPL_GRAV * umf.at(
#                                 K=k
#                             ) * qij / (
#                                 pifc0.at(K=k - 1) - pifc0.at(K=k)
#                             )  # [ kg/kg/s ]
#                         qc_lm = qc_lm - constants.MAPL_GRAV * umf.at(K=k) * ql0.at(
#                             K=k
#                         ) / (
#                             pifc0.at(K=k - 1) - pifc0.at(K=k)
#                         )  # [ kg/kg/s ]
#                         qc_im = qc_im - constants.MAPL_GRAV * umf.at(K=k) * qi0.at(
#                             K=k
#                         ) / (
#                             pifc0.at(K=k - 1) - pifc0.at(K=k)
#                         )  # [ kg/kg/s ]
#                         # nc_lm   = nc_lm   - g * umf(k) * tr0(k,ixnumliq)  / ( pifc0(k-1) - pifc0(k) ) ! [ kg/kg/s ]
#                         # nc_im   = nc_im   - g * umf(k) * tr0(k,ixnumice)  / ( pifc0(k-1) - pifc0(k) ) ! [ kg/kg/s ]

#                     # 4. Cumulative Penetrative entrainment detrained in the 'kbup' layer
#                     # Explicitly compute the properties detrained penetrative entrained airs in k = kbup layer.

#                 if k == kbup:
#                     thj, qvj, ql_emf_kbup, qi_emf_kbup, qse, id_check = conden(
#                         pmid0.at(K=k),
#                         thlu_emf.at(K=k),
#                         qtu_emf.at(K=k),
#                         ese,
#                         esx,
#                     )
#                     if id_check == 1:
#                         id_exit = True
#                         umf_out[0, 0, 1] = 0.0
#                         dcm_out = 0.0
#                         qvten_out = 0.0
#                         qlten_out = 0.0
#                         qiten_out = 0.0
#                         sten_out = 0.0
#                         uten_out = 0.0
#                         vten_out = 0.0
#                         qrten_out = 0.0
#                         qsten_out = 0.0
#                         cufrc_out = 0.0
#                         cush_inout = -1.0
#                         qldet_out = 0.0
#                         qidet_out = 0.0
#                         qtflx_out[0, 0, 1] = 0.0
#                         slflx_out[0, 0, 1] = 0.0
#                         uflx_out[0, 0, 1] = 0.0
#                         vflx_out[0, 0, 1] = 0.0
#                         fer_out = constants.MAPL_UNDEF
#                         fdr_out = constants.MAPL_UNDEF

#                     if id_exit == False:
#                         if ql_emf_kbup < 0.0:
#                             nl_emf_kbup = 0.0

#                         if qi_emf_kbup < 0.0:
#                             ni_emf_kbup = 0.0

#                         qc_lm = qc_lm - constants.MAPL_GRAV * emf.at(K=k) * (
#                             ql_emf_kbup - ql0.at(K=k)
#                         ) / (
#                             pifc0.at(K=k - 1) - pifc0.at(K=k)
#                         )  # [ kg/kg/s ]
#                         qc_im = qc_im - constants.MAPL_GRAV * emf.at(K=k) * (
#                             qi_emf_kbup - qi0.at(K=k)
#                         ) / (
#                             pifc0.at(K=k - 1) - pifc0.at(K=k)
#                         )  # [ kg/kg/s ]
#                         # nc_lm   = nc_lm   - g * emf(k) * ( nl_emf_kbup - tr0(k,ixnumliq) ) / ( pifc0(k-1) - pifc0(k) ) ! [ kg/kg/s ]
#                         # nc_im   = nc_im   - g * emf(k) * ( ni_emf_kbup - tr0(k,ixnumice) ) / ( pifc0(k-1) - pifc0(k) ) ! [ kg/kg/s ]
#                 if id_exit == False:
#                     if THIS_K == k:
#                         qlten_det = qc_l.at(K=k) + qc_lm
#                         qiten_det = qc_i.at(K=k) + qc_im

#                     if ((qc_lm + qlten_sink.at(K=k)) * dt + ql0.at(K=k)) < 0.0:
#                         totsink = qc_lm + qlten_sink.at(K=k)
#                         if totsink != 0.0:
#                             qc_lm = -(ql0.at(K=k) / dt) * qc_lm / totsink
#                             if THIS_K == k:
#                                 qlten_sink = (
#                                     -(ql0.at(K=k) / dt) * qlten_sink.at(K=k) / totsink
#                                 )
#                                 qlten_det = qc_l.at(K=k) + qc_lm

#                     if ((qc_im + qiten_sink.at(K=k)) * dt + qi0.at(K=k)) < 0.0:
#                         totsink = qc_im + qiten_sink.at(K=k)
#                         if totsink != 0.0:
#                             qc_im = -(qi0.at(K=k) / dt) * qc_im / totsink
#                             if THIS_K == k:
#                                 qiten_sink = (
#                                     -(qi0.at(K=k) / dt) * qiten_sink.at(K=k) / totsink
#                                 )
#                                 qiten_det = qc_i.at(K=k) + qc_im

#                     if THIS_K == k:
#                         qlten = qrten.at(K=k) + qlten_sink.at(K=k) + qlten_det.at(K=k)
#                         qiten = qsten.at(K=k) + qiten_sink.at(K=k) + qiten_det.at(K=k)

#                         qvten = qtten.at(K=k) - qlten.at(K=k) - qiten.at(K=k)
#                         sten = (
#                             slten.at(K=k)
#                             + constants.MAPL_LATENT_HEAT_VAPORIZATION * qlten.at(K=k)
#                             + constants.MAPL_LATENT_HEAT_SUBLIMATION * qiten.at(K=k)
#                         )

#                         qc = qc_l.at(K=k) + qc_i.at(K=k)

#                         qlten = qlten.at(K=k) - qrten.at(K=k)
#                         qiten = qiten.at(K=k) - qsten.at(K=k)
#                         qtten = qlten.at(K=k) + qiten.at(K=k) + qvten.at(K=k)

#                         slten = (
#                             sten.at(K=k)
#                             - constants.MAPL_LATENT_HEAT_VAPORIZATION * qlten.at(K=k)
#                             - constants.MAPL_LATENT_HEAT_SUBLIMATION * qiten.at(K=k)
#                         )
#                         slten = (
#                             slten.at(K=k)
#                             + constants.MAPL_LATENT_HEAT_VAPORIZATION * qrten.at(K=k)
#                             + constants.MAPL_LATENT_HEAT_SUBLIMATION * qsten.at(K=k)
#                         )
#                         sten = (
#                             slten.at(K=k)
#                             + constants.MAPL_LATENT_HEAT_VAPORIZATION * qlten.at(K=k)
#                             + constants.MAPL_LATENT_HEAT_SUBLIMATION * qiten.at(K=k)
#                         )

#                 k += 1

#     with computation(FORWARD), interval(...):
#         if id_exit == False:
#             # Prevent the onset-of negative condensate at the next time step
#             # Potentially, this block can be moved just in front of the above
#             # block.

#             # Modification : I should check whether this 'positive_moisture_single' routine is
#             #                consistent with the one used in UW PBL and cloud macrophysics schemes.
#             # Modification : Below may overestimate resulting 'ql, qi' if we use the new 'qc_l', 'qc_i'
#             #                in combination with the original computation of qlten, qiten. However,
#             #                if we use new 'qlten,qiten', there is no problem.

#             qv0_star = qv0 + qvten * dt
#             ql0_star = ql0 + qlten * dt
#             qi0_star = qi0 + qiten * dt
#             s0_star = s0 + sten * dt

#     with computation(BACKWARD), interval(...):
#         if id_exit == False:
#             ixcldice = 1
#             ixcldliq = 2
#             dql: f64 = max(f64(0.0), f64(1.0) * qmin.at(K=ixcldliq) - ql)
#             dqi: f64 = max(f64(0.0), f64(1.0) * qmin.at(K=ixcldice) - qi)
#             qlten = qlten + dql / dt
#             qiten = qiten + dqi / dt
#             qvten = qvten - (dql + dqi) / dt
#             sten = (
#                 sten
#                 + constants.MAPL_LATENT_HEAT_VAPORIZATION * (dql / dt)
#                 + constants.MAPL_LATENT_HEAT_SUBLIMATION * (dqi / dt)
#             )
#             ql = ql + dql
#             qi = qi + dqi
#             qv = qv - dql - dqi
#             s = (
#                 s
#                 + constants.MAPL_LATENT_HEAT_VAPORIZATION * dql
#                 + constants.MAPL_LATENT_HEAT_SUBLIMATION * dqi
#             )
#             dqv = max(0.0, 1.0 * qmin.at(K=0) - qv)
#             qvten = qvten + dqv / dt
#             qv = qv + dqv

#     with computation(BACKWARD), interval(1, None):
#         if id_exit == False:
#             qv[0, 0, -1] = qv[0, 0, -1] - dqv * dp / dp[0, 0, -1]
#             qvten[0, 0, -1] = qvten[0, 0, -1] - dqv * dp / dp[0, 0, -1] / dt

#     with computation(BACKWARD), interval(...):
#         if id_exit == False:
#             qv = max(qv, qmin.at(K=0))
#             ql = max(ql, qmin.at(K=ixcldliq))
#             qi = max(qi, qmin.at(K=ixcldice))

#     with computation(PARALLEL), interval(...):
#         if id_exit == False:
#             # Extra moisture used to satisfy 'qv(i,1)=qvmin' is proportionally
#             # extracted from all the layers that has 'qv > 2*qvmin'. This fully
#             # preserves column moisture.
#             if dqv > f64(1.0e-20):
#                 sum: f64 = 0.0
#                 if THIS_K <= k0:
#                     if qv > f64(2.0) * qmin.at(K=0):
#                         sum = sum + qv * dp
#                 aa: f64 = dqv * dp.at(K=1) / max(f64(1.0e-20), sum)
#                 if aa < f64(0.5):
#                     if THIS_K <= k0:
#                         if qv > f64(2.0) * qmin.at(K=0):
#                             dum: f64 = aa * qv
#                             qv = qv - dum
#                             qvten = qvten - dum / dt
#                 # else:
#                 # if scverbose == True:
#                 # print('Full positive_moisture is impossible in uwshcu')

#     with computation(FORWARD), interval(...):
#         if id_exit == False:
#             qtten = qvten + qlten + qiten
#             slten = (
#                 sten
#                 - constants.MAPL_LATENT_HEAT_VAPORIZATION * qlten
#                 - constants.MAPL_LATENT_HEAT_SUBLIMATION * qiten
#             )

#             # Tendencies of tracers

#             # NEED TO LOOK AT THIS SECTION AGAIN !!!
#             # if dotransport == 1.0:
#             #     do m = 1, ncnst

#             # trmin = 0. !qmin(m)

#             # trflx_d(0:k0) = 0.
#             # trflx_u(0:k0) = 0.
#             # do k = 1, k0-1
#             #         pdelx = dp0(k)

#             #         km1 = k - 1
#             #         dum = ( tr0(k,m) - trmin ) *  pdelx / g / dt + trflx(km1,m) - trflx(k,m) + trflx_d(km1)
#             #         trflx_d(k) = min( 0., dum )

#             # do k = k0, 2, -1
#             #         pdelx = dp0(k)
#             #         km1 = k - 1
#             #         dum = ( tr0(k,m) - trmin ) * pdelx / g / dt + trflx(km1,m) - trflx(k,m) + &
#             #                                                 trflx_d(km1) - trflx_d(k) - trflx_u(k)
#             #         trflx_u(km1) = max( 0., -dum )

#             # do k = 1, k0
#             #         pdelx = dp0(k)

#             #     km1 = k - 1
#             #     trten(k,m) = ( trflx(km1,m) - trflx(k,m) + &
#             #                     trflx_d(km1) - trflx_d(k) + &
#             #                     trflx_u(km1) - trflx_u(k) ) * g / pdelx

#     with computation(FORWARD), interval(...):
#         if id_exit == False:
#             # Compute default diagnostic outputs
#             # Note that since 'qtu(krel-1:kpen-1)' & 'thlu(krel-1:kpen-1)' has
#             # been adjusted after detraining cloud condensate into environment
#             # during cumulus updraft motion,  below calculations will  exactly
#             # reproduce in-cloud properties as shown in the output analysis.

#             thj, qvj, qlj, qij, qse, id_check = conden(
#                 prel, thlu(krel - 1), qtu(krel - 1), ese, esx
#             )
#             if id_check == 1.0:
#                 id_exit = True
#                 umf_out[0, 0, 1] = 0.0
#                 dcm_out = 0.0
#                 qvten_out = 0.0
#                 qlten_out = 0.0
#                 qiten_out = 0.0
#                 sten_out = 0.0
#                 uten_out = 0.0
#                 vten_out = 0.0
#                 qrten_out = 0.0
#                 qsten_out = 0.0
#                 cufrc_out = 0.0
#                 cush_inout = -1.0
#                 qldet_out = 0.0
#                 qidet_out = 0.0
#                 qtflx_out[0, 0, 1] = 0.0
#                 slflx_out[0, 0, 1] = 0.0
#                 uflx_out[0, 0, 1] = 0.0
#                 vflx_out[0, 0, 1] = 0.0
#                 fer_out = constants.MAPL_UNDEF
#                 fdr_out = constants.MAPL_UNDEF

#         if id_exit == False:
#             qcubelow = qlj + qij
#             qlubelow = qlj
#             qiubelow = qij
#             rcwp = 0.0
#             rlwp = 0.0
#             riwp = 0.0

#             k = krel
#             while k <= kpen:
#                 if id_exit == False:
#                     # Calculate cumulus condensate at the upper interface of each layer.
#                     # Note 'ppen < 0' and at 'k=kpen' layer, I used 'thlu_top'&'qtu_top'
#                     # which explicitly considered zero or non-zero 'fer(kpen)'.

#                     if k == kpen:
#                         thj, qvj, qlj, qij, qse, id_check = conden(
#                             pifc0.at(K=k - 1) + ppen, thlu_top, qtu_top, ese, esx
#                         )
#                     else:
#                         thj, qvj, qlj, qij, qse, id_check = conden(
#                             pifc0.at(K=k), thlu.at(K=k), qtu.at(K=k), ese, esx
#                         )

#                     if id_check == 1:
#                         id_exit = True
#                         umf_out[0, 0, 1] = 0.0
#                         dcm_out = 0.0
#                         qvten_out = 0.0
#                         qlten_out = 0.0
#                         qiten_out = 0.0
#                         sten_out = 0.0
#                         uten_out = 0.0
#                         vten_out = 0.0
#                         qrten_out = 0.0
#                         qsten_out = 0.0
#                         cufrc_out = 0.0
#                         cush_inout = -1.0
#                         qldet_out = 0.0
#                         qidet_out = 0.0
#                         qtflx_out[0, 0, 1] = 0.0
#                         slflx_out[0, 0, 1] = 0.0
#                         uflx_out[0, 0, 1] = 0.0
#                         vflx_out[0, 0, 1] = 0.0
#                         fer_out = constants.MAPL_UNDEF
#                         fdr_out = constants.MAPL_UNDEF

#                     if id_exit == False:
#                         # Calculate in-cloud mean LWC ( qlu(k) ), IWC ( qiu(k) ),  & layer
#                         # mean cumulus fraction ( cufrc(k) ),  vertically-integrated layer
#                         # mean LWP and IWP. Expel some of in-cloud condensate at the upper
#                         # interface if it is largr than criqc. Note cumulus cloud fraction
#                         # is assumed to be twice of core updraft fractional area. Thus LWP
#                         # and IWP will be twice of actual value coming from our scheme.
#                         if THIS_K == k:
#                             qcu = 0.5 * (qcubelow + qlj + qij)
#                             qlu = 0.5 * (qlubelow + qlj)
#                             qiu = 0.5 * (qiubelow + qij)
#                             cufrc = ufrc.at(K=k - 1) + ufrc.at(K=k)

#                         if k == krel:
#                             if THIS_K == k:
#                                 cufrc = (
#                                     (ufrclcl + ufrc.at(K=k))
#                                     * (prel - pifc0.at(K=k))
#                                     / (pifc0.at(K=k - 1) - pifc0.at(K=k))
#                                 )
#                         elif k == kpen:
#                             if THIS_K == k:
#                                 cufrc = (
#                                     (ufrc.at(K=k - 1) + 0.0)
#                                     * (-ppen)
#                                     / (pifc0.at(K=k - 1) - pifc0.at(K=k))
#                                 )
#                                 if (qlj + qij) > criqc:
#                                     qcu = 0.5 * (qcubelow + criqc)
#                                     qlu = 0.5 * (qlubelow + criqc * qlj / (qlj + qij))
#                                     qiu = 0.5 * (qiubelow + criqc * qij / (qlj + qij))
#                         rcwp = rcwp + (qlu.at(K=k) + qiu.at(K=k)) * (
#                             pifc0.at(K=k - 1) - pifc0.at(K=k)
#                         ) / constants.MAPL_GRAV * cufrc.at(K=k)
#                         rlwp = rlwp + qlu.at(K=k) * (
#                             pifc0.at(K=k - 1) - pifc0.at(K=k)
#                         ) / constants.MAPL_GRAV * cufrc.at(K=k)
#                         riwp = riwp + qiu.at(K=k) * (
#                             pifc0.at(K=k - 1) - pifc0.at(K=k)
#                         ) / constants.MAPL_GRAV * cufrc.at(K=k)
#                         qcubelow = qlj + qij
#                         qlubelow = qlj
#                         qiubelow = qij

#                 k += 1

#                 if id_exit == False:
#                     # Cloud top and base interface indices
#                     cnt = f32(kpen)
#                     cnb = f32(krel - 1)
#                     iter_cin = 2

#                     # End of formal calculation. Below blocks are for implicit CIN calculations
#                     # with re-initialization and save variables at iter_cin = 1.

#                     # Adjust the original input profiles for implicit CIN calculation
#                     if iteration != iter_cin:

#                         # Save the output from "iter_cin = 1"
#                         # These output will be writed-out if "iter_cin = 1" was not performed
#                         # for some reasons.

#                         qv0_s = qv0 + qvten * dt
#                         ql0_s = ql0 + qlten * dt
#                         qi0_s = qi0 + qiten * dt
#                         s0_s = s0 + sten * dt
#                         u0_s = u0 + uten * dt
#                         v0_s = v0 + vten * dt

#                         t0_s = t0 + sten * dt / constants.MAPL_CP

#                         if dotransport == 1.0:
#                             n = 0
#                             while n < ncnst:
#                                 tr0_s[0, 0, 0][n] = (
#                                     tr0[0, 0, 0][n] + trten[0, 0, 0][n] * dt
#                                 )
#                                 n += 1

#                         umf_s[0, 0, 1] = umf[0, 0, 1]
#                         dcm_s = dcm
#                         qvten_s = qvten
#                         qlten_s = qlten
#                         qiten_s = qiten
#                         sten_s = sten
#                         uten_s = uten
#                         vten_s = vten
#                         qrten_s = qrten
#                         qsten_s = qsten
#                         cush_s = cush
#                         cufrc_s = cufrc
#                         slflx_s[0, 0, 1] = slflx[0, 0, 1]
#                         qtflx_s[0, 0, 1] = qtflx[0, 0, 1]
#                         uflx_s[0, 0, 1] = uflx[0, 0, 1]
#                         vflx_s[0, 0, 1] = vflx[0, 0, 1]
#                         qcu_s = qcu
#                         qlu_s = qlu
#                         qiu_s = qiu
#                         fer_s = fer
#                         fdr_s = fdr
#                         xc_s = xco
#                         cin_s = cin
#                         cinlcl_s = cinlcl
#                         cbmf_s = cbmf
#                         qc_s = qc
#                         qldet_s = qlten_det
#                         qidet_s = qiten_det
#                         qlsub_s = qlten_sink
#                         qisub_s = qiten_sink

#                         ufrc_s[0, 0, 1] = ufrc[0, 0, 1]

#                         # Recalculate environmental variables for new cin calculation at "iter_cin = 2"
#                         # using the updated state variables. Perform only for variables necessary  for
#                         # the new cin calculation.

#                         qv0 = qv0_s
#                         ql0 = ql0_s
#                         qi0 = qi0_s
#                         s0 = s0_s
#                         t0 = t0_s

#                         qt0 = qv0 + ql0 + qi0
#                         thl0 = (
#                             t0
#                             - constants.MAPL_LATENT_HEAT_VAPORIZATION
#                             * ql0
#                             / constants.MAPL_CP
#                             - constants.MAPL_LATENT_HEAT_SUBLIMATION
#                             * qi0
#                             / constants.MAPL_CP
#                         ) / exnmid0
#                         thvl0 = (1.0 + zvir * qt0) * thl0

#                         # NEED TO REVISIT THIS !!!
#                         ssthl0 = slope(k0, thl0, pmid0)
#                         ssqt0 = slope(k0, qt0, pmid0)
#                         ssu0 = slope(k0, u0, pmid0)
#                         ssv0 = slope(k0, v0, pmid0)
#                         if dotransport == 1.0:
#                             n = 0
#                             while n < ncnst:
#                                 sstr0[0, 0, 0][n] = slope(k0, tr0[0, 0, 0][n], pmid0)
#                                 n += 1

#                         thl0bot = thl0.at(K=k) + ssthl0.at(K=k) * (
#                             pifc0.at(K=k - 1) - pmid0.at(K=k)
#                         )
#                         qt0bot = qt0.at(K=k) + ssqt0.at(K=k) * (
#                             pifc0.at(K=k - 1) - pmid0.at(K=k)
#                         )
#                         thj, qvj, qlj, qij, qse, id_check = conden(
#                             pifc0.at(K=k - 1), thl0bot, qt0bot, ese, esx
#                         )
#                         if id_check == 1:
#                             id_exit = True
#                             umf_out[0, 0, 1] = 0.0
#                             dcm_out = 0.0
#                             qvten_out = 0.0
#                             qlten_out = 0.0
#                             qiten_out = 0.0
#                             sten_out = 0.0
#                             uten_out = 0.0
#                             vten_out = 0.0
#                             qrten_out = 0.0
#                             qsten_out = 0.0
#                             cufrc_out = 0.0
#                             cush_inout = -1.0
#                             qldet_out = 0.0
#                             qidet_out = 0.0
#                             qtflx_out[0, 0, 1] = 0.0
#                             slflx_out[0, 0, 1] = 0.0
#                             uflx_out[0, 0, 1] = 0.0
#                             vflx_out[0, 0, 1] = 0.0
#                             fer_out = constants.MAPL_UNDEF
#                             fdr_out = constants.MAPL_UNDEF

#                         if id_exit == False:
#                             if THIS_K == k:
#                                 thv0bot = thj * (1.0 + zvir * qvj - qlj - qij)
#                                 thvl0bot = thl0bot * (1.0 + zvir * qt0bot)

#                             thl0top = thl0.at(K=k) + ssthl0.at(K=k) * (
#                                 pifc0.at(K=k) - pmid0.at(K=k)
#                             )
#                             qt0top = qt0.at(K=k) + ssqt0.at(K=k) * (
#                                 pifc0.at(K=k) - pmid0.at(K=k)
#                             )
#                             thj, qvj, qlj, qij, qse, id_check = conden(
#                                 pifc0.at(K=k), thl0top, qt0top, ese, esx
#                             )
#                             if id_check == 1:
#                                 id_exit = True
#                                 umf_out[0, 0, 1] = 0.0
#                                 dcm_out = 0.0
#                                 qvten_out = 0.0
#                                 qlten_out = 0.0
#                                 qiten_out = 0.0
#                                 sten_out = 0.0
#                                 uten_out = 0.0
#                                 vten_out = 0.0
#                                 qrten_out = 0.0
#                                 qsten_out = 0.0
#                                 cufrc_out = 0.0
#                                 cush_inout = -1.0
#                                 qldet_out = 0.0
#                                 qidet_out = 0.0
#                                 qtflx_out[0, 0, 1] = 0.0
#                                 slflx_out[0, 0, 1] = 0.0
#                                 uflx_out[0, 0, 1] = 0.0
#                                 vflx_out[0, 0, 1] = 0.0
#                                 fer_out = constants.MAPL_UNDEF
#                                 fdr_out = constants.MAPL_UNDEF

#                             if id_exit == False:
#                                 if THIS_K == k:
#                                     thv0top = thj * (1.0 + zvir * qvj - qlj - qij)
#                                     thvl0top = thl0top * (1.0 + zvir * qt0top)

#                             # End of iter loop


# def update_output_variables(
#     umf: FloatField,
#     zifc0: FloatField,
#     kinv: IntField,
#     dcm: FloatField,
#     qvten: FloatField,
#     qlten: FloatField,
#     qiten: FloatField,
#     sten: FloatField,
#     uten: FloatField,
#     vten: FloatField,
#     qrten: FloatField,
#     qsten: FloatField,
#     cufrc: FloatField,
#     cush: FloatFieldIJ,
#     qlten_det: FloatField,
#     qiten_det: FloatField,
#     qlten_sink: FloatField,
#     qiten_sink: FloatField,
#     rdrop: Float,
#     qtflx: FloatField,
#     slflx: FloatField,
#     uflx: FloatField,
#     vflx: FloatField,
#     dotransport: Float,
#     tr0_inout: FloatField,
#     trten: FloatField,
#     dt: Int,
#     fer: FloatField,
#     fdr: FloatField,
#     kpen: IntField,
# ):
#     with computation(FORWARD), interval(...):
#         if id_exit == False:
#             umf_out[0, 0, 1] = umf[0, 0, 1]

#             if THIS_K <= kinv - 1:
#                 umf_out[0, 0, 1] = umf.at(K=kinv - 1) * zifc0 / zifc0.at(K=kinv - 1)

#             dcm_out = dcm
#             qvten_out = qvten
#             qlten_out = qlten
#             qiten_out = qiten
#             sten_out = sten
#             uten_out = uten
#             vten_out = vten
#             qrten_out = qrten
#             qsten_out = qsten
#             cufrc_out = cufrc
#             cush_inout = cush
#             qldet_out = qlten_det
#             qidet_out = qiten_det
#             qlsub_out = qlten_sink
#             qisub_out = qiten_sink
#             ndrop_out = qlten_det / (4188.787 * rdrop**3)

#             nice_out = qiten_det / (3.0e-10)  # /crystal mass
#             qtflx_out[0, 0, 1] = qtflx
#             slflx_out[0, 0, 1] = slflx
#             uflx_out[0, 0, 1] = uflx
#             vflx_out[0, 0, 1] = vflx

#             if dotransport == 1.0:
#                 n = 0
#                 while n < ncnst:
#                     tr0_inout[0, 0, 0][n] = (
#                         tr0_inout[0, 0, 0][n] + trten[0, 0, 0][n] * dt
#                     )
#                     n += 1

#             # Below are specific diagnostic output for detailed
#             # analysis of cumulus scheme

#             if THIS_K < kpen:
#                 fer_out = fer
#                 fdr_out = fdr


class ComputeUwshcu:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        formulation: SaturationFormulation = SaturationFormulation.Staars,
    ) -> None:
        """
        Initialize the ComputeUwshcu class.

        Parameters:
        stencil_factory (StencilFactory): Factory for creating stencil computations.
        quantity_factory (QuantityFactory): Factory for creating quantities.
        formulation: Saturation Formulation used for QSat
        """

        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory
        grid_indexing = stencil_factory.grid_indexing

        self._compute_thermodynamic_variables = self.stencil_factory.from_dims_halo(
            func=compute_thermodynamic_variables,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._implicit_cin_closure = self.stencil_factory.from_dims_halo(
            func=implicit_cin_closure,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        # self._buoyancy_sorting_mixing = self.stencil_factory.from_dims_halo(
        #     func=buoyancy_sorting_mixing,
        #     compute_dims=[X_DIM, Y_DIM, Z_DIM],
        # )

        # self._update_output_variables = self.stencil_factory.from_dims_halo(
        #     func=update_output_variables,
        #     compute_dims=[X_DIM, Y_DIM, Z_DIM],
        # )

        self._k_idx = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        for k in range(0, self._k_idx.view[:].shape[2]):
            self._k_idx.view[:, :, k] = k

        self.id_exit = self.quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_DIM], "n/a", dtype=bool
        )

    @staticmethod
    def make_ntracers_quantity_factory(ijk_quantity_factory: QuantityFactory):
        ntracers_quantity_factory = copy.deepcopy(ijk_quantity_factory)
        ntracers_quantity_factory.set_extra_dim_lengths(
            **{
                "ntracers": constants.NCNST,
            }
        )
        return ntracers_quantity_factory

    def __call__(
        self,
        windsrcavg: Int,
        qtsrchgt: Float,
        qtsrc_fac: Float,
        thlsrc_fac: Float,
        frc_rasn: Float,
        rbuoy: Float,
        epsvarw: Float,
        use_CINcin: Int,
        mumin1: Float,
        rmaxfrac: Float,
        PGFc: Float,
        # niter_xc: Int,
        # criqc: Float,
        # rle: Float,
        # mixscale: Float,
        # rkm: Float,
        # detrhgt: Float,
        # cridist_opt: Int,
        # rdrag: Float,
        # use_self_detrain: Int,
        # use_cumpenent: Int,
        # rpen: Float,
        # use_momenflx: Int,
        k0: Int,
        dt: Float,
        ncnst: Int,
        pifc0_in: FloatField,
        zifc0_in: FloatField,
        exnifc0_in: FloatField,
        pmid0_in: FloatField,
        zmid0_in: FloatField,
        zmid0: FloatField,
        exnmid0_in: FloatField,
        dp0_in: FloatField,
        u0_in: FloatField,
        v0_in: FloatField,
        qv0_in: FloatField,
        ql0_in: FloatField,
        qi0_in: FloatField,
        th0_in: FloatField,
        tr0_inout: FloatField_NTracers,
        kpbl_in: IntFieldIJ,
        frland_in: FloatFieldIJ,
        tke_in: FloatField,
        rkfre: FloatFieldIJ,
        cush_inout: FloatFieldIJ,
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
        fer_out: FloatField,
        fdr_out: FloatField,
        # qldet_out: FloatField,
        # qidet_out: FloatField,
        # qlsub_out: FloatField,
        # qisub_out: FloatField,
        # ndrop_out: FloatField,
        # nice_out: FloatField,
        shfx: FloatFieldIJ,
        evap: FloatFieldIJ,
        # cnvtr: FloatFieldIJ,
        # tpert_out: FloatFieldIJ,
        # qpert_out: FloatFieldIJ,
        qtflx_out: FloatField,
        slflx_out: FloatField,
        uflx_out: FloatField,
        vflx_out: FloatField,
        dotransport: Int,
        # Outputs for testing
        ssthl0: FloatField,
        ssqt0: FloatField,
        ssu0: FloatField,
        ssv0: FloatField,
        sstr0: FloatField_NTracers,
        tr0: FloatField_NTracers,
        qt0: FloatField,
        thj: FloatField,
        qij: FloatField,
        qlj: FloatField,
        qvj: FloatField,
        qse: FloatField,
        id_check: IntField,
        thv0top: FloatField,
        thvl0top: FloatField,
        thvl0: FloatField,
        thl0top: FloatField,
        qt0top: FloatField,
        thvl0bot: FloatField,
        # thl0bot: FloatFieldIJ,
        tr0_o: FloatField_NTracers,
        sstr0_o: FloatField_NTracers,
        trflx: FloatField_NTracers,
        trten: FloatField_NTracers,
        tru: FloatField_NTracers,
        tru_emf: FloatField_NTracers,
        kinv: IntField,
        thvlavg: FloatField,
        tkeavg: FloatField,
        uavg: FloatField,
        vavg: FloatField,
        thvlmin: FloatField,
        qtavg: FloatField,
        dpi: FloatFieldIJ,
        t0: FloatField,
        qv0: FloatField,
        pmid0: FloatField,
        thl0: FloatField,
        thlsrc: FloatField,
        usrc: FloatField,
        vsrc: FloatField,
        trsrc: FloatField_NTracers,
        qtsrc: FloatField,
        plcl: FloatField,
        klcl: IntField,
        thl0lcl: FloatField,
        qt0lcl: FloatField,
        thv0lcl: FloatField,
        thv0bot: FloatField,
        plfc: FloatFieldIJ,
        cin: FloatField,
        thvubot: FloatField,
        thvutop: FloatField,
        thvlsrc: FloatField,
        thj2: FloatField,
        qvj2: FloatField,
        qlj2: FloatField,
        qij2: FloatField,
        qse2: FloatField,
        test_var1: FloatField,
        formulation: SaturationFormulation = SaturationFormulation.Staars,
    ):
        self.qsat = QSat(
            self.stencil_factory,
            self.quantity_factory,
            formulation=formulation,
        )

        self.id_exit.view[:, :, :] = False

        self._compute_thermodynamic_variables(
            k0=k0,
            dt=dt,
            ncnst=ncnst,
            pifc0_in=pifc0_in,
            zifc0_in=zifc0_in,
            exnifc0_in=exnifc0_in,
            pmid0_in=pmid0_in,
            zmid0_in=zmid0_in,
            exnmid0_in=exnmid0_in,
            dp0_in=dp0_in,
            u0_in=u0_in,
            v0_in=v0_in,
            qv0_in=qv0_in,
            ql0_in=ql0_in,
            qi0_in=qi0_in,
            th0_in=th0_in,
            tr0_inout=tr0_inout,
            kpbl_in=kpbl_in,
            frland_in=frland_in,
            # tke_in=tke_in,
            rkfre=rkfre,
            cush_inout=cush_inout,
            umf_out=umf_out,
            dcm_out=dcm_out,
            qvten_out=qvten_out,
            qlten_out=qvten_out,
            qiten_out=qiten_out,
            sten_out=sten_out,
            uten_out=uten_out,
            vten_out=vten_out,
            qrten_out=qrten_out,
            qsten_out=qsten_out,
            cufrc_out=cufrc_out,
            fer_out=fer_out,
            fdr_out=fdr_out,
            # qldet_out=qldet_out,
            # qidet_out=qidet_out,
            # qlsub_out=qlsub_out,
            # qisub_out=qisub_out,
            # ndrop_out=ndrop_out,
            # nice_out=nice_out,
            shfx=shfx,
            evap=evap,
            # cush=cush,
            # cnvtr=cnvtr,
            # tpert_out=tpert_out,
            # qpert_out=qpert_out,
            qtflx_out=qtflx_out,
            slflx_out=slflx_out,
            uflx_out=uflx_out,
            vflx_out=vflx_out,
            # zifc0=zifc0,
            # exnifc0=exnifc0,
            # pifc0=pifc0,
            # tke=tke_in,
            # u0=u0,
            # v0=v0,
            thvl0=thvl0,
            thv0bot=thv0bot,
            thvl0bot=thvl0bot,
            thl0top=thl0top,
            qt0top=qt0top,
            # thl0bot=thl0bot,
            # thvl0top=thvl0top,
            zmid0=zmid0,
            qt0=qt0,
            t0=t0,
            qv0=qv0,
            pmid0=pmid0,
            # zvir=zvir,
            tr0=tr0,
            ssu0=ssu0,
            ssv0=ssv0,
            ssthl0=ssthl0,
            ssqt0=ssqt0,
            sstr0=sstr0,
            thl0=thl0,
            # qv0_o=qv0_o,
            # ql0_o=ql0_o,
            # qi0_o=qi0_o,
            # t0_o=t0_o,
            # s0_o=s0_o,
            # u0_o=u0_o,
            # v0_o=v0_o,
            # qt0_o=qt0_o,
            # thl0_o=thl0_o,
            # thvl0_o=thvl0_o,
            # ssthl0_o=ssthl0_o,
            # ssqt0_o=ssqt0_o,
            # thv0bot_o=thv0bot_o,
            # thv0top_o=thv0_top_o,
            # thvl0bot_o=thvl0bot_o,
            # thvl0top_o=thvl0top_o,
            # ssu0_o=ssu0_o,
            # ssv0_o=ssv0_o,
            # tr0_o=tr0_o,
            # sstr0_o=sstr0_o,
            # dp0=dp0,
            thj=thj,
            qij=qij,
            qlj=qlj,
            qvj=qvj,
            qse=qse,
            id_check=id_check,
            thv0top=thv0top,
            thvl0top=thvl0top,
            tr0_o=tr0_o,
            sstr0_o=sstr0_o,
            trflx=trflx,
            trten=trten,
            tru=tru,
            tru_emf=tru_emf,
            dotransport=dotransport,
            k_idx=self._k_idx,
            ese=self.qsat.ese,
            esx=self.qsat.esx,
            id_exit=self.id_exit,
        )

        if (self.id_exit.view[:, :, :] == True).any():
            raise NotImplementedError(
                "Expected id_exit == False, got id_exit == True! "
            )

        iteration = 1
        iter_cin = 1  # Change to 2
        while iteration <= iter_cin:
            self._implicit_cin_closure(
                # Inputs:
                iteration=iteration,
                windsrcavg=windsrcavg,
                qtsrchgt=qtsrchgt,
                qtsrc_fac=qtsrc_fac,
                thlsrc_fac=thlsrc_fac,
                # rbuoy=rbuoy,
                # epsvarw=epsvarw,
                # use_CINcin=use_CINcin,
                # mumin1=mumin1,
                # rmaxfrac=rmaxfrac,
                # PGFc=PGFc,
                cush=cush_inout,
                kpbl_in=kpbl_in,
                k0=k0,
                # dt=dt,
                pifc0=pifc0_in,
                # zifc0=zifc0,
                # exnifc0=exnifc0,
                # rkfre=rkfre,
                tke_in=tke_in,
                u0=u0_in,
                v0=v0_in,
                thvl0=thvl0,
                thvl0bot=thvl0bot,
                thvl0top=thvl0top,
                zmid0=zmid0,
                qt0=qt0,
                t0=t0,
                qv0=qv0,
                shfx=shfx,
                evap=evap,
                pmid0=pmid0,
                dotransport=dotransport,
                ncnst=ncnst,
                zvir=zvir,
                tr0=tr0,
                ssu0=ssu0,
                ssv0=ssv0,
                ssthl0=ssthl0,
                ssqt0=ssqt0,
                # sstr0=sstr0,
                thl0=thl0,
                thv0bot=thv0bot,
                thv0top=thv0top,
                # qv0_o=qv0_o,
                # ql0_o=ql0_o,
                # qi0_o=qi0_o,
                # t0_o=t0_o,
                # s0_o=s0_o,
                # u0_o=u0_o,
                # v0_o=v0_o,
                # qt0_o=qt0_o,
                # thl0_o=thl0_o,
                # thvl0_o=thvl0_o,
                # ssthl0_o=ssthl0_o,
                # ssqt0_o=ssqt0_o,
                # thv0bot_o=thv0bot_o,
                # thv0top_o=thv0_top_o,
                # thvl0bot_o=thvl0bot_o,
                # thvl0top_o=thvl0top_o,
                # ssu0_o=ssu0_o,
                # ssv0_o=ssv0_o,
                # tr0_o=tr0_o,
                # sstr0_o=sstr0_o,
                # dp0=dp0,
                umf_out=umf_out,
                dcm_out=dcm_out,
                qvten_out=qvten_out,
                qlten_out=qlten_out,
                qiten_out=qiten_out,
                sten_out=sten_out,
                uten_out=uten_out,
                vten_out=vten_out,
                qrten_out=qrten_out,
                qsten_out=qsten_out,
                cufrc_out=cufrc_out,
                cush_inout=cush_inout,
                # qldet_out=qldet_out,
                # qidet_out=qidet_out,
                qtflx_out=qtflx_out,
                slflx_out=slflx_out,
                uflx_out=uflx_out,
                vflx_out=vflx_out,
                fer_out=fer_out,
                fdr_out=fdr_out,
                ese=self.qsat.ese,
                esx=self.qsat.esx,
                id_exit=self.id_exit,
                # Outputs for testing
                kinv=kinv,
                thvlavg=thvlavg,
                tkeavg=tkeavg,
                uavg=uavg,
                vavg=vavg,
                thvlmin=thvlmin,
                qtavg=qtavg,
                dpi=dpi,
                thlsrc=thlsrc,
                usrc=usrc,
                vsrc=vsrc,
                trsrc=trsrc,
                qtsrc=qtsrc,
                plcl=plcl,
                klcl=klcl,
                thl0lcl=thl0lcl,
                qt0lcl=qt0lcl,
                thv0lcl=thv0lcl,
                plfc=plfc,
                cin=cin,
                thvubot=thvubot,
                thvutop=thvutop,
                thvlsrc=thvlsrc,
                thj2=thj2,
                qvj2=qvj2,
                qlj2=qlj2,
                qij2=qij2,
                qse2=qse2,
                test_var1=test_var1,
            )

            # self._buoyancy_sorting_mixing(
            #     id_exit=id_exit,
            #     tscaleh=tscaleh,
            #     krel=krel,
            #     wlcl=wlcl,
            #     prel=prel,
            #     pifc0=pifc0,
            #     thv0rel=thv0rel,
            #     thl0=thl0,
            #     ssthl0=ssthl0,
            #     pmid0=pmid0,
            #     qt0=qt0,
            #     ssqt0=ssqt0,
            #     u0=u0,
            #     ssu0=ssu0,
            #     v0=v0,
            #     ssv0=ssv0,
            #     dotransport=dotransport,
            #     ncnst=ncnst,
            #     tr0=tr0,
            #     sstr0=sstr0,
            #     k0=k0,
            #     thlu=thlu,
            #     qtu=qtu,
            #     wu=wu,
            #     ese=self.qsat.ese,
            #     esx=self.qsat.esx,
            #     qsat_pe=qsat_pe,
            #     umf_out=umf_out,
            #     dcm_out=dcm_out,
            #     qvten_out=qvten_out,
            #     qlten_out=qlten_out,
            #     qiten_out=qiten_out,
            #     sten_out=sten_out,
            #     uten_out=uten_out,
            #     vten_out=vten_out,
            #     qrten_out=qrten_out,
            #     qsten_out=qsten_out,
            #     cufrc_out=cufrc_out,
            #     cush_inout=cush_inout,
            #     qldet_out=qldet_out,
            #     qidet_out=qidet_out,
            #     qtflx_out=qtflx_out,
            #     slflx_out= slflx_out,
            #     uflx_out=uflx_out,
            #     vflx_out=vflx_out,
            #     fer_out=fer_out,
            #     fdr_out=fdr_out,
            #     rbuoy=rbuoy,
            #     zifc0=zifc0,
            #     zmid0=zmid0,
            #     dp0=dp0,
            #     dt=dt,
            #     umf=umf,
            #     niter_xc=niter_xc,
            #     criqc=criqc,
            #     rle=rle,
            #     mixscale=mixscale,
            #     rkm=rkm,
            #     detrhgt=detrhgt,
            #     cridist_opt=cridist_opt,
            #     PGFc=PGFc,
            #     tru=tru,
            #     tre=tre,
            #     exnifc0=exnifc0,
            #     thv0top=thv0top,
            #     thv0bot=thv0bot,
            #     rmaxfrac=rmaxfrac,
            #     exnmid0=exnmid0,
            #     emf=emf,
            #     tru_emf=tru_emf,
            #     rdrag=rdrag,
            #     use_self_detrain=use_self_detrain,
            #     use_cumpenent=use_cumpenent,
            #     rpen=rpen,
            #     qtsrc=qtsrc,
            #     kinv=kinv,
            #     cbmf=cbmf,
            #     thlsrc=thlsrc,
            #     usrc=usrc,
            #     vsrc=vsrc,
            #     trsrc=trsrc,
            #     trflx=trflx,
            #     use_momenflx=use_momenflx,
            #     ql0=ql0,
            #     qi0=qi0,
            #     frc_rasn=frc_rasn,
            #     qv0=qv0,
            #     s0=s0,
            #     t0=t0,
            #     iteration=iteration,
            #     tr0_s=tr0_s,
            #     trten=trten,
            #     umf_s=umf_s,
            #     slflx_s=slflx_s,
            #     qtflx_s=qtflx_s,
            #     uflx_s=uflx_s,
            #     vflx_s=vflx_s,
            #     cin=cin,
            #     cinlcl=cinlcl,
            #     ufrc_s=ufrc_s,
            #     ufrclcl=ufrclcl)

            # if (self.id_exit.view[:, :, :] == True).any():
            #     raise NotImplementedError(
            #         "Expected id_exit == False, got id_exit == True! " "Exit UWSHCU! "
            #     )

            iteration += 1

        # self._update_output_variables(
        #     umf=umf,
        #     zifc0=zifc0,
        #     kinv=kinv,
        #     dcm=dcm,
        #     qvten=qvten,
        #     qlten=qlten,
        #     qiten=qiten,
        #     sten=sten,
        #     uten=uten,
        #     vten=vten,
        #     qrten=qrten,
        #     qsten=qsten,
        #     cufrc=cufrc,
        #     cush=cush,
        #     qlten_det=qlten_det,
        #     qiten_det=qiten_det,
        #     qlten_sink=qlten_sink,
        #     qiten_sink=qiten_sink,
        #     rdrop=rdrop,
        #     qtflx=qtflx,
        #     slflx=slflx,
        #     uflx=uflx,
        #     vflx=vflx,
        #     dotransport=dotransport,
        #     tr0_inout=tr0_inout,
        #     trten=trten,
        #     dt=dt,
        #     fer=fer,
        #     fdr=fdr,
        #     kpen=kpen,
        # )
