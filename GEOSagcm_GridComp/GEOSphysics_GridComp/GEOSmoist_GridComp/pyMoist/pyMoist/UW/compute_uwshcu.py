import copy
import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import (
    computation,
    interval,
    PARALLEL,
    FORWARD,
    exp,
    sqrt,
    log,
    erfc,
    THIS_K,
    f32,
    i64,
    i32,
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
from pyMoist.saturation.qsat import QSat, FloatField_Extra_Dim
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
    exnerfn,
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
    # thvl0top: FloatField,
    zmid0: FloatField,
    qt0: FloatField,
    # t0: FloatField,
    # qv0: FloatField,
    # pmid0: FloatField,
    # zvir: Float,
    tr0: FloatField_NTracers,
    ssu0: FloatField,
    ssv0: FloatField,
    ssthl0: FloatField,
    ssqt0: FloatField,
    sstr0: FloatField_NTracers,
    # thl0: FloatField,
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
        tke = tke_in[0, 0, 1]
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
        qt0 = qv0 + ql0 + qi0
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
        qt0 = qv0 + ql0 + qi0
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
        qt0 = qv0 + ql0 + qi0
        t0 = th0_in * exnmid0
        s0 = constants.MAPL_GRAV * zmid0 + constants.MAPL_CP * t0
        thl0 = (
            t0
            - constants.MAPL_ALHL * ql0 / constants.MAPL_CP
            - constants.MAPL_ALHS * qi0 / constants.MAPL_CP
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
    tke: FloatField,
    # rkfre: FloatFieldIJ,
    u0: FloatField,
    v0: FloatField,
    thvl0: FloatField,
    thvl0bot: FloatField,
    thvl0top: FloatField,
    zmid0: FloatField,
    qt0: FloatField,
    # t0: FloatField,
    # qv0: FloatField,
    shfx: FloatFieldIJ,
    evap: FloatFieldIJ,
    # pmid0: FloatField,
    dotransport: Int,
    ncnst: Int,
    zvir: Float,
    tr0: FloatField_NTracers,
    ssu0: FloatField,
    ssv0: FloatField,
    ssthl0: FloatField,
    ssqt0: FloatField,
    # sstr0: FloatField_NTracers,
    # thl0: FloatField,
    # thv0bot: FloatField,
    # thv0top: FloatField,
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
    kinv: IntFieldIJ,
    thvlavg: FloatFieldIJ,
    tkeavg: FloatFieldIJ,
    uavg: FloatFieldIJ,
    vavg: FloatFieldIJ,
    thvlmin: FloatFieldIJ,
    qtavg: FloatFieldIJ,
    dpi: FloatFieldIJ,
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
            tkeavg = 0.0
            qtavg = 0.0
            uavg = 0.0
            vavg = 0.0

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

    with computation(FORWARD), interval(...):

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
        if id_exit == False:
            dpsum = 0.0
            thvlmin = 1000.0
            thvlavg = 0.0

            lev = 0
            while lev < kinv:
                kbelow = lev - 1
                dpi = pifc0.at(K=kbelow) - pifc0.at(K=lev)
                dpsum = dpsum + dpi
                tkeavg = tkeavg + dpi * tke.at(K=lev)
                uavg = uavg + dpi * u0.at(K=lev)
                vavg = vavg + dpi * v0.at(K=lev)
                thvlavg = thvlavg + dpi * thvl0.at(K=lev)
                if lev != kinv:
                    thvlmin = min(thvlmin, min(thvl0bot.at(K=lev), thvl0top.at(K=lev)))
                lev += 1

            tkeavg = tkeavg / dpsum
            uavg = uavg / dpsum
            vavg = vavg / dpsum
            thvlavg = thvlavg / dpsum

            # weighted average over lowest 20mb
            # dpsum = 0.
            # if k <= kinv:
            #   dpi = max(0.,(2e3+pmid0-pifc0.at(K=0))/2e3)
            #   qtavg  = qtavg  + dpi*qt0
            #   dpsum = dpsum + dpi
            # qtavg   = qtavg/dpsum

            # Interpolate qt to specified height
            lev = 0
            while zmid0.at(K=lev) < qtsrchgt:
                lev += 1
            if lev > 1:
                kbelow = lev - 1
                qtavg = qt0.at(K=kbelow) * (zmid0.at(K=lev) - qtsrchgt) + qt0.at(
                    K=lev
                ) * (qtsrchgt - zmid0.at(K=kbelow))
                qtavg = qtavg / (zmid0.at(K=lev) - zmid0.at(K=kbelow))
            else:
                qtavg = qt0.at(K=0)


'''
            """
            Find characteristics of cumulus source air: qtsrc,thlsrc,usrc,vsrc 
            Note that 'thlsrc' was concocted using 'thvlsrc' and 'qtsrc'.      
            'qtsrc' is defined as the lowest layer mid-point value;   'thlsrc' 
            is from 'qtsrc' and 'thvlmin=thvlsrc'; 'usrc' & 'vsrc' are defined 
            as the values just below the PBL top interface. 
            """
            if windsrcavg == 1:
                zrho = pifc0.at(K=0) / (
                    287.04 * (t0.at(K=1) * (1.0 + 0.608 * qv0.at(K=1)))
                )
                buoyflx = (
                    -shfx / constants.MAPL_CP - 0.608 * t0(1) * evap
                ) / zrho  # K m s-1
                # delzg = (zifc0.at(K=1)-zifc0.at(K=0))*constants.MAPL_GRAV
                delzg = (50.0) * constants.MAPL_GRAV  # assume 50m surface scale
                wstar = max(0.0, 0.001 - 0.41 * buoyflx * delzg / t0.at(K=1))  # m3 s-3
                qpert_out = 0.0
                tpert_out = 0.0
                if wstar > 0.001:
                    wstar = 1.0 * wstar**0.3333
                    tpert_out = (
                        thlsrc_fac * shfx / (zrho * wstar * constants.MAPL_CP)
                    )  # K
                    qpert_out = qtsrc_fac * evap / (zrho * wstar)  # kg kg-1
                    qpert_out = max(
                        min(qpert_out, 0.02 * qt0.at(K=1)), 0.0
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
                    qtsrc = qt0.at(K=1)
                    thvlsrc = thvlmin
                    thlsrc = thvlsrc / (1.0 + zvir * qtsrc)
                    kbelow = kinv - 1
                    usrc = u0.at(K=kbelow) + ssu0.at(K=kbelow) * (
                        pifc0.at(K=kbelow) - pmid0.at(K=kbelow)
                    )
                    vsrc = v0.at(K=kbelow) + ssv0.at(K=kbelow) * (
                        pifc0.at(K=kbelow) - pmid0.at(K=kbelow)
                    )

            if dotransport == 1.0:
                n = 0
                while n < ncnst:
                    trsrc[0, 0, 0][n] = tr0.at(K=1, ddim=[n])
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

    
    with computation(FORWARD), interval(...):
        k = 0
        while k <= k0:
            if pifc0.at(K=k) < plcl:
                klcl = k
                k = k0 + 1  # Break out of loop
            else:
                klcl = k0
            k += 1

        klcl = max(1, klcl)

        if plcl < 60000.0:
            id_exit = True
            # go to 333

        """ 
        Calculate environmental virtual potential temperature at LCL, 
        'thv0lcl' which is solely used in the 'cin' calculation. Note  
        that 'thv0lcl' is calculated first by calculating  'thl0lcl'  
        and 'qt0lcl' at the LCL, and performing 'conden' afterward,   
        in fully consistent with the other parts of the code. 
        """

        thl0lcl = thl0.at(K=klcl) + ssthl0.at(K=klcl) * (plcl - pmid0.at(K=klcl))
        qt0lcl = qt0.at(K=klcl) + ssqt0.at(K=klcl) * (plcl - pmid0.at(K=klcl))
        thj, qvj, qlj, qij, qse, id_check = conden(plcl, thl0lcl, qt0lcl, ese, esx)
        if id_check == 1:
            id_exit = True
            # go to 333
        thv0lcl = thj * (1.0 + zvir * qvj - qlj - qij)


    with computation(PARALLEL), interval(...):
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
        stop35 = False
        if klcl >= kinv and stop35 == False:
            while kinv <= (k0 - 1) and stop35 == False:
                kbelow = kinv - 1
                if kinv < klcl:
                    thvubot = thvlsrc
                    thvutop = thvlsrc
                    cin = cin + single_cin(
                        pifc0.at(K=kbelow),
                        thv0bot.at(K=kinv),
                        pifc0.at(K=kinv),
                        thv0top.at(K=kinv),
                        thvubot,
                        thvutop,
                    )
                elif kinv == klcl:
                    # ----- Bottom to LCL
                    thvubot = thvlsrc
                    thvutop = thvlsrc
                    cin = cin + single_cin(
                        pifc0.at(K=kbelow),
                        thv0bot.at(K=kinv),
                        plcl,
                        thv0lcl,
                        thvubot,
                        thvutop,
                    )
                    cinlcl = max(cin, 0.0)
                    cin = cinlcl
                    # ----- LCL to Top
                    thvubot = thvlsrc
                    thj, qvj, qlj, qij, qse, id_check = conden(
                        pifc0.at(K=kinv), thlsrc, qtsrc, ese, esx
                    )
                    if id_check == 1:
                        id_exit = True
                        # go to 333
                    thvutop = thj * (1.0 + zvir * qvj - qlj - qij)

                    if thvubot > thv0lcl and thvutop > thv0top.at(K=kinv):
                        plfc, _ = getbuoy(
                            plcl,
                            thv0lcl,
                            pifc0.at(K=kinv),
                            thv0top.at(K=kinv),
                            thvubot,
                            thvutop,
                            cin,
                            plfc,
                        )
                    elif thvubot < thv0lcl and thvutop < thv0top.at(K=kinv):
                        _, cin = getbuoy(
                            plcl,
                            thv0lcl,
                            pifc0.at(K=kinv),
                            thv0top.at(K=kinv),
                            thvubot,
                            thvutop,
                            cin,
                            plfc,
                        )
                    elif thvubot > thv0lcl and thvutop < thv0top.at(K=kinv):
                        _, cin = getbuoy(
                            plcl,
                            thv0lcl,
                            pifc0.at(K=kinv),
                            thv0top.at(K=kinv),
                            thvubot,
                            thvutop,
                            cin,
                            plfc,
                        )
                    else:
                        plfc, cin = getbuoy(
                            plcl,
                            thv0lcl,
                            pifc0.at(K=kinv),
                            thv0top.at(K=kinv),
                            thvubot,
                            thvutop,
                            cin,
                            plfc,
                        )

                    if plfc > 0.0:
                        klfc = kinv
                        stop35 = True
                else:
                    thvubot = thvutop
                    thj, qvj, qlj, qij, qse, id_check = conden(
                        pifc0.at(K=kinv), thlsrc, qtsrc, ese, esx
                    )
                    if id_check == 1:
                        id_exit = True
                        # go to 333
                    thvutop = thj * (1.0 + zvir * qvj - qlj - qij)

                    if thvubot > thv0bot.at(K=kinv) and thvutop > thv0top.at(K=kinv):
                        plfc, _ = getbuoy(
                            pifc0.at(K=kbelow),
                            thv0bot.at(K=kinv),
                            pifc0.at(K=kinv),
                            thv0top.at(K=kinv),
                            thvubot,
                            thvutop,
                            cin,
                            plfc,
                        )
                    elif thvubot < thv0bot.at(K=kinv) and thvutop < thv0top.at(K=kinv):
                        _, cin = getbuoy(
                            pifc0.at(K=kbelow),
                            thv0bot.at(K=kinv),
                            pifc0.at(K=kinv),
                            thv0top.at(K=kinv),
                            thvubot,
                            thvutop,
                            cin,
                            plfc,
                        )
                    elif thvubot > thv0bot.at(K=kinv) and thvutop < thv0top.at(K=kinv):
                        _, cin = getbuoy(
                            pifc0.at(K=kbelow),
                            thv0bot.at(K=kinv),
                            pifc0.at(K=kinv),
                            thv0top.at(K=kinv),
                            thvubot,
                            thvutop,
                            cin,
                            plfc,
                        )
                    else:
                        plfc, cin = getbuoy(
                            pifc0.at(K=kbelow),
                            thv0bot.at(K=kinv),
                            pifc0.at(K=kinv),
                            thv0top.at(K=kinv),
                            thvubot,
                            thvutop,
                            cin,
                            plfc,
                        )

                    if plfc > 0.0:
                        klfc = kinv
                        stop35 = True
                kinv += 1

        else:
            """
            #Case 2. LCL height is lower than PBL interface ( 'pLCL > ps0(kinv-1)')
            """
            cinlcl = 0.0
            while kinv <= (k0 - 1) and stop35 == False:
                kbelow = kinv - 1
                thj, qvj, qlj, qij, qse, id_check = conden(
                    pifc0.at(K=kbelow), thlsrc, qtsrc, ese, esx
                )
                if id_check == 1:
                    id_exit = True
                    # go to 333
                thvubot = thj * (1.0 + zvir * qvj - qlj - qij)
                thj, qvj, qlj, qij, qse, id_check = conden(
                    pifc0.at(K=kinv), thlsrc, qtsrc, ese, esx
                )
                if id_check == 1:
                    id_exit = True
                    # go to 333
                thvutop = thj * (1.0 + zvir * qvj - qlj - qij)

                if thvubot > thv0bot.at(K=kinv) and thvutop > thv0top.at(K=kinv):
                    plfc, _ = getbuoy(
                        pifc0.at(K=kbelow),
                        thv0bot.at(K=kinv),
                        pifc0.at(K=kinv),
                        thv0top.at(K=kinv),
                        thvubot,
                        thvutop,
                        plfc,
                        cin,
                    )
                elif thvubot < thv0bot.at(K=kinv) and thvutop < thv0top.at(K=kinv):
                    _, cin = getbuoy(
                        pifc0.at(K=kbelow),
                        thv0bot.at(K=kinv),
                        pifc0.at(K=kinv),
                        thv0top.at(K=kinv),
                        thvubot,
                        thvutop,
                        plfc,
                        cin,
                    )
                elif thvubot > thv0bot.at(K=kinv) and thvutop < thv0top.at(K=kinv):
                    _, cin = getbuoy(
                        pifc0.at(K=kbelow),
                        thv0bot.at(K=kinv),
                        pifc0.at(K=kinv),
                        thv0top.at(K=kinv),
                        thvubot,
                        thvutop,
                        plfc,
                        cin,
                    )
                else:
                    plfc, cin = getbuoy(
                        pifc0.at(K=kbelow),
                        thv0bot.at(K=kinv),
                        pifc0.at(K=kinv),
                        thv0top.at(K=kinv),
                        thvubot,
                        thvutop,
                        plfc,
                        cin,
                    )

                if plfc > 0.0:
                    klfc = kinv
                    stop35 = True

                kinv += 1  # End of CIN case selection

        cin = max(0.0, cin)
        # cin = max(cin,0.04*(lts-18.))   # kludge to reduce UW in StCu regions

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
        thvl0bot: FloatField,
        tr0_o: FloatField_NTracers,
        sstr0_o: FloatField_NTracers,
        trflx: FloatField_NTracers,
        trten: FloatField_NTracers,
        tru: FloatField_NTracers,
        tru_emf: FloatField_NTracers,
        kinv: IntFieldIJ,
        thvlavg: FloatFieldIJ,
        tkeavg: FloatFieldIJ,
        uavg: FloatFieldIJ,
        vavg: FloatFieldIJ,
        thvlmin: FloatFieldIJ,
        qtavg: FloatFieldIJ,
        dpi: FloatFieldIJ,
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
            tke_in=tke_in,
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
            # tke=tke,
            # u0=u0,
            # v0=v0,
            thvl0=thvl0,
            thvl0bot=thvl0bot,
            # thvl0top=thvl0top,
            zmid0=zmid0,
            qt0=qt0,
            # t0=t0,
            # qv0=qv0,
            # pmid0=pmid0,
            # zvir=zvir,
            tr0=tr0,
            ssu0=ssu0,
            ssv0=ssv0,
            ssthl0=ssthl0,
            ssqt0=ssqt0,
            sstr0=sstr0,
            # thl0=thl0,
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

        # if (self.id_exit.view[:, :, :] == True).any():
        #     raise NotImplementedError(
        #         "Expected id_exit == False, got id_exit == True! "
        #         "Exit UWSHCU! "
        #         "This code has not been ported."
        #     )

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
                tke=tke_in,
                # rkfre=rkfre,
                u0=u0_in,
                v0=v0_in,
                thvl0=thvl0,
                thvl0bot=thvl0bot,
                thvl0top=thvl0top,
                zmid0=zmid0,
                qt0=qt0,
                # t0=t0,
                # qv0=qv0,
                shfx=shfx,
                evap=evap,
                # pmid0=pmid0,
                dotransport=dotransport,
                ncnst=ncnst,
                zvir=zvir,
                tr0=tr0,
                ssu0=ssu0,
                ssv0=ssv0,
                ssthl0=ssthl0,
                ssqt0=ssqt0,
                # sstr0=sstr0,
                # thl0=thl0,
                # thv0bot=thv0bot,
                # thv0top=thv0top,
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
            )

            # if (self.id_exit.view[:, :, :] == True).any():
            #     raise NotImplementedError(
            #         "Expected id_exit == False, got id_exit == True! "
            #         "Exit UWSHCU! "
            #         "This code has not been ported."
            #     )

            iteration += 1
