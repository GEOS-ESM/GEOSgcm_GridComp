from dataclasses import dataclass

from ndsl import Local, NDSLRuntime, QuantityFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Int


@dataclass
class UWLocals:
    ssthl0: Local
    ssqt0: Local
    ssu0: Local
    ssv0: Local
    thj: Local
    qlj: Local
    qvj: Local
    qse: Local
    qij: Local
    thv0top: Local
    thv0bot: Local
    thvl0top: Local
    dcm_out: Local
    qvten_out: Local
    qlten_out: Local
    qiten_out: Local
    sten_out: Local
    uten_out: Local
    vten_out: Local
    qrten_out: Local
    qsten_out: Local
    cufrc_out: Local
    fer_out: Local
    fdr_out: Local
    thvlavg: Local
    tkeavg: Local
    uavg: Local
    vavg: Local
    thvlmin: Local
    qtavg: Local
    zmid0: Local
    qt0: Local
    thvl0: Local
    thvl0bot: Local
    t0: Local
    qv0: Local
    pmid0: Local
    thl0: Local
    thlsrc: Local
    usrc: Local
    vsrc: Local
    plcl: Local
    thl0lcl: Local
    qt0lcl: Local
    thv0lcl: Local
    plfc: Local
    fer_outvar: Local
    fdr_outvar: Local
    cin: Local
    thvubot: Local
    thvutop: Local
    thvlsrc: Local
    thl0top: Local
    qt0top: Local
    qldet_outvar: Local
    qidet_outvar: Local
    qlsub_outvar: Local
    qisub_outvar: Local
    dcm_outvar: Local
    qvten_outvar: Local
    qlten_outvar: Local
    qiten_outvar: Local
    sten_outvar: Local
    uten_outvar: Local
    vten_outvar: Local
    qrten_outvar: Local
    qsten_outvar: Local
    cufrc_outvar: Local
    usrc_o: Local
    vsrc_o: Local
    thv0lcl_o: Local
    ql0_o: Local
    qi0_o: Local
    t0_o: Local
    s0_o: Local
    u0_o: Local
    v0_o: Local
    qt0_o: Local
    thl0_o: Local
    thvl0_o: Local
    ssthl0_o: Local
    ssqt0_o: Local
    thv0bot_o: Local
    thv0top_o: Local
    thvl0bot_o: Local
    thvl0top_o: Local
    ssu0_o: Local
    ssv0_o: Local
    dcm_s: Local
    qvten_s: Local
    qlten_s: Local
    qiten_s: Local
    sten_s: Local
    uten_s: Local
    vten_s: Local
    qrten_s: Local
    qsten_s: Local
    qldet_s: Local
    qidet_s: Local
    qlsub_s: Local
    qisub_s: Local
    cush_s: Local
    cufrc_s: Local
    fer_s: Local
    fdr_s: Local
    qtsrc_o: Local
    thvlsrc_o: Local
    thlsrc_o: Local
    qldet_out: Local
    qidet_out: Local
    qlsub_out: Local
    qisub_out: Local
    ndrop_out: Local
    nice_out: Local
    dcm: Local
    xco: Local
    qc: Local
    qlten_det: Local
    qiten_det: Local
    qv0_s: Local
    ql0_s: Local
    qi0_s: Local
    s0_s: Local
    t0_s: Local
    u0_s: Local
    v0_s: Local
    slten: Local
    qv0_o: Local
    plcl_o: Local
    plfc_o: Local
    tkeavg_o: Local
    thvlmin_o: Local
    ufrclcl: Local
    qcu: Local
    qlu: Local
    qiu: Local
    cufrc: Local
    qtsrc: Local
    uplus_3D: Local
    vplus_3D: Local
    prel: Local
    thv0rel: Local
    winv: Local
    cbmf: Local
    rho0inv: Local
    ufrcinv: Local
    wlcl: Local
    qsat_pe: Local
    thlue: Local
    qtue: Local
    wue: Local
    rei: Local
    fer: Local
    dwten: Local
    diten: Local
    ql0: Local
    qi0: Local
    uten: Local
    vten: Local
    uf: Local
    vf: Local
    dwten_temp: Local
    diten_temp: Local
    fdr: Local
    qlten_sink: Local
    qiten_sink: Local
    qrten: Local
    qsten: Local
    s0: Local
    qvten: Local
    qlten: Local
    sten: Local
    qiten: Local
    qmin: Local
    pmid0_in: Local
    u0_in: Local
    v0_in: Local
    zmid0_in: Local
    exnmid0_in: Local
    dp0_in: Local
    qv0_in: Local
    ql0_in: Local
    qi0_in: Local
    th0_in: Local
    cush_inout: Local
    dpi: Local
    thvlmin_IJ: Local
    wcrit: Local
    alpha: Local
    del_CIN: Local
    cin_IJ: Local
    plfc_IJ: Local
    cinlcl_IJ: Local
    pe: Local
    thle: Local
    qte: Local
    dpe: Local
    exne: Local
    thvebot: Local
    ue: Local
    ve: Local
    drage: Local
    bogbot: Local
    bogtop: Local
    rhomid0j: Local
    cush_inoutvar: Local
    uplus: Local
    vplus: Local
    cin_i: Local
    cinlcl_i: Local
    ke: Local
    thlu_top: Local
    qtu_top: Local
    cldhgt: Local
    qlubelow: Local
    qiubelow: Local
    qlj_2D: Local
    qij_2D: Local
    qcubelow: Local
    rcwp: Local
    rlwp: Local
    riwp: Local
    ppen: Local
    tscaleh: Local
    wtwb: Local
    cnvtrmax: Local
    qtu_emf: Local
    umf_out: Local
    qtflx_out: Local
    slflx_out: Local
    slflx: Local
    thlu_emf: Local
    uu_emf: Local
    vu_emf: Local
    uemf: Local
    uflx_out: Local
    vflx_out: Local
    ufrc: Local
    wu: Local
    emf: Local
    thlu: Local
    qtu: Local
    thvu: Local
    uu: Local
    vu: Local
    umf_zint: Local
    umf_outvar: Local
    qtflx_outvar: Local
    slflx_outvar: Local
    uflx_outvar: Local
    vflx_outvar: Local
    slflx_s: Local
    qtflx_s: Local
    uflx_s: Local
    vflx_s: Local
    qtflx: Local
    uflx: Local
    ufrc_s: Local
    xflx: Local
    vflx: Local
    umf_temp: Local
    umf_s: Local
    tke_in: Local
    pifc0_in: Local
    zifc0_in: Local
    exnifc0_in: Local
    kinv: Local
    klcl: Local
    klfc: Local
    kinv_o: Local
    klcl_o: Local
    klfc_o: Local
    kbup: Local
    krel: Local
    kpen: Local
    kbup_IJ: Local
    klfc_IJ: Local
    kpen_IJ: Local
    kpbl_in: Local
    tr0_temp: Local
    u0: Local
    v0: Local
    cinlcl: Local

    @classmethod
    def make(cls, runtime: NDSLRuntime, quantity_factory: QuantityFactory):
        # FloatFields

        ssthl0 = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        ssqt0 = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        ssu0 = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        ssv0 = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        thj = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qlj = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qvj = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qse = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qij = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        tr0_temp = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        thv0top = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        thv0bot = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        thvl0top = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        dcm_out = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qvten_out = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qlten_out = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qiten_out = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        sten_out = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        uten_out = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        vten_out = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qrten_out = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qsten_out = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        cufrc_out = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        fer_out = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        fdr_out = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        thvlavg = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        tkeavg = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        uavg = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        vavg = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        thvlmin = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qtavg = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        zmid0 = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qt0 = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        thvl0 = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        thvl0bot = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        t0 = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qv0 = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        pmid0 = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        thl0 = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        thlsrc = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        usrc = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        vsrc = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        plcl = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        thl0lcl = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qt0lcl = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        thv0lcl = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        plfc = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        fer_outvar = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        fdr_outvar = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        cin = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        thvubot = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        thvutop = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        thvlsrc = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        thl0top = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qt0top = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qldet_outvar = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qidet_outvar = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qlsub_outvar = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qisub_outvar = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        dcm_outvar = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qvten_outvar = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qlten_outvar = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qiten_outvar = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        sten_outvar = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        uten_outvar = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        vten_outvar = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qrten_outvar = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qsten_outvar = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        cufrc_outvar = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        usrc_o = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        vsrc_o = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        thv0lcl_o = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        ql0_o = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qi0_o = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        t0_o = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        s0_o = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        u0_o = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        v0_o = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qt0_o = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        thl0_o = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        thvl0_o = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        ssthl0_o = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        ssqt0_o = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        thv0bot_o = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        thv0top_o = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        thvl0bot_o = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        thvl0top_o = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        ssu0_o = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        ssv0_o = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        dcm_s = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qvten_s = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qlten_s = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qiten_s = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        sten_s = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        uten_s = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        vten_s = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qrten_s = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qsten_s = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qldet_s = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qidet_s = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qlsub_s = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qisub_s = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        cush_s = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        cufrc_s = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        fer_s = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        fdr_s = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qtsrc_o = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        thvlsrc_o = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        thlsrc_o = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qldet_out = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qidet_out = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qlsub_out = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qisub_out = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        ndrop_out = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        nice_out = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        dcm = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        xco = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qc = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qlten_det = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qiten_det = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qv0_s = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        ql0_s = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qi0_s = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        s0_s = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        t0_s = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        u0_s = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        v0_s = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        slten = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qv0_o = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        plcl_o = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        plfc_o = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        tkeavg_o = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        thvlmin_o = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        ufrclcl = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qcu = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qlu = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qiu = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        cufrc = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qtsrc = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        uplus_3D = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        vplus_3D = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        prel = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        thv0rel = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        winv = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        cbmf = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        rho0inv = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        ufrcinv = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        wlcl = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qsat_pe = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        thlue = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qtue = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        wue = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        rei = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        fer = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        dwten = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        diten = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        ql0 = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qi0 = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        uten = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        vten = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        uf = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        vf = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        dwten_temp = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        diten_temp = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        fdr = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qlten_sink = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qiten_sink = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qrten = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qsten = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        s0 = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qvten = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qlten = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        sten = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qiten = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qmin = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        pmid0_in = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        u0_in = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        v0_in = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        u0 = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        v0 = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        zmid0_in = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        exnmid0_in = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        dp0_in = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qv0_in = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        ql0_in = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        qi0_in = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        th0_in = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        cinlcl = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        # FloatFieldIJs
        cush_inout = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        dpi = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        thvlmin_IJ = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        wcrit = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        alpha = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        del_CIN = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        cin_IJ = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        plfc_IJ = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        cinlcl_IJ = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        pe = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        thle = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        qte = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        dpe = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        exne = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        thvebot = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        ue = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        ve = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        drage = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        bogbot = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        bogtop = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        rhomid0j = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        cush_inoutvar = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        uplus = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        vplus = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        cin_i = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        cinlcl_i = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        ke = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        thlu_top = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        qtu_top = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        cldhgt = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        qlubelow = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        qiubelow = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        qlj_2D = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        qij_2D = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        qcubelow = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        rcwp = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        rlwp = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        riwp = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        ppen = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        tscaleh = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        wtwb = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        cnvtrmax = runtime.make_local(quantity_factory, [X_DIM, Y_DIM])
        # Interface FloatFields
        qtu_emf = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        umf_out = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        qtflx_out = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        slflx_out = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        slflx = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        thlu_emf = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        uu_emf = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        vu_emf = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        uemf = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        uflx_out = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        vflx_out = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        ufrc = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        wu = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        emf = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        thlu = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        qtu = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        uu = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        vu = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        umf_zint = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        thvu = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        umf_outvar = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        qtflx_outvar = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        slflx_outvar = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        uflx_outvar = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        vflx_outvar = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        slflx_s = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        qtflx_s = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        uflx_s = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        vflx_s = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        qtflx = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        uflx = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        ufrc_s = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        xflx = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        vflx = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        umf_temp = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        umf_s = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        tke_in = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        pifc0_in = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        zifc0_in = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        exnifc0_in = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        # IntFields
        kinv = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM], dtype=Int)
        klcl = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM], dtype=Int)
        klfc = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM], dtype=Int)
        kinv_o = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM], dtype=Int)
        klcl_o = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM], dtype=Int)
        klfc_o = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM], dtype=Int)
        kbup = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM], dtype=Int)
        krel = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM], dtype=Int)
        kpen = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM], dtype=Int)

        # IntFieldIJs
        kbup_IJ = runtime.make_local(quantity_factory, [X_DIM, Y_DIM], dtype=Int)
        klfc_IJ = runtime.make_local(quantity_factory, [X_DIM, Y_DIM], dtype=Int)
        kpen_IJ = runtime.make_local(quantity_factory, [X_DIM, Y_DIM], dtype=Int)
        kpbl_in = runtime.make_local(quantity_factory, [X_DIM, Y_DIM], dtype=Int)

        return cls(
            ssthl0,
            ssqt0,
            ssu0,
            ssv0,
            thj,
            qlj,
            qvj,
            qse,
            qij,
            thv0top,
            thv0bot,
            thvl0top,
            dcm_out,
            qvten_out,
            qlten_out,
            qiten_out,
            sten_out,
            uten_out,
            vten_out,
            qrten_out,
            qsten_out,
            cufrc_out,
            fer_out,
            fdr_out,
            thvlavg,
            tkeavg,
            uavg,
            vavg,
            thvlmin,
            qtavg,
            zmid0,
            qt0,
            thvl0,
            thvl0bot,
            t0,
            qv0,
            pmid0,
            thl0,
            thlsrc,
            usrc,
            vsrc,
            plcl,
            thl0lcl,
            qt0lcl,
            thv0lcl,
            plfc,
            fer_outvar,
            fdr_outvar,
            cin,
            thvubot,
            thvutop,
            thvlsrc,
            thl0top,
            qt0top,
            qldet_outvar,
            qidet_outvar,
            qlsub_outvar,
            qisub_outvar,
            dcm_outvar,
            qvten_outvar,
            qlten_outvar,
            qiten_outvar,
            sten_outvar,
            uten_outvar,
            vten_outvar,
            qrten_outvar,
            qsten_outvar,
            cufrc_outvar,
            usrc_o,
            vsrc_o,
            thv0lcl_o,
            ql0_o,
            qi0_o,
            t0_o,
            s0_o,
            u0_o,
            v0_o,
            qt0_o,
            thl0_o,
            thvl0_o,
            ssthl0_o,
            ssqt0_o,
            thv0bot_o,
            thv0top_o,
            thvl0bot_o,
            thvl0top_o,
            ssu0_o,
            ssv0_o,
            dcm_s,
            qvten_s,
            qlten_s,
            qiten_s,
            sten_s,
            uten_s,
            vten_s,
            qrten_s,
            qsten_s,
            qldet_s,
            qidet_s,
            qlsub_s,
            qisub_s,
            cush_s,
            cufrc_s,
            fer_s,
            fdr_s,
            qtsrc_o,
            thvlsrc_o,
            thlsrc_o,
            qldet_out,
            qidet_out,
            qlsub_out,
            qisub_out,
            ndrop_out,
            nice_out,
            dcm,
            xco,
            qc,
            qlten_det,
            qiten_det,
            qv0_s,
            ql0_s,
            qi0_s,
            s0_s,
            t0_s,
            u0_s,
            v0_s,
            slten,
            qv0_o,
            plcl_o,
            plfc_o,
            tkeavg_o,
            thvlmin_o,
            ufrclcl,
            qcu,
            qlu,
            qiu,
            cufrc,
            qtsrc,
            uplus_3D,
            vplus_3D,
            prel,
            thv0rel,
            winv,
            cbmf,
            rho0inv,
            ufrcinv,
            wlcl,
            qsat_pe,
            thlue,
            qtue,
            wue,
            rei,
            fer,
            dwten,
            diten,
            ql0,
            qi0,
            uten,
            vten,
            uf,
            vf,
            dwten_temp,
            diten_temp,
            fdr,
            qlten_sink,
            qiten_sink,
            qrten,
            qsten,
            s0,
            qvten,
            qlten,
            sten,
            qiten,
            qmin,
            pmid0_in,
            u0_in,
            v0_in,
            zmid0_in,
            exnmid0_in,
            dp0_in,
            qv0_in,
            ql0_in,
            qi0_in,
            th0_in,
            cush_inout,
            dpi,
            thvlmin_IJ,
            wcrit,
            alpha,
            del_CIN,
            cin_IJ,
            plfc_IJ,
            cinlcl_IJ,
            pe,
            thle,
            qte,
            dpe,
            exne,
            thvebot,
            ue,
            ve,
            drage,
            bogbot,
            bogtop,
            rhomid0j,
            cush_inoutvar,
            uplus,
            vplus,
            cin_i,
            cinlcl_i,
            ke,
            thlu_top,
            qtu_top,
            cldhgt,
            qlubelow,
            qiubelow,
            qlj_2D,
            qij_2D,
            qcubelow,
            rcwp,
            rlwp,
            riwp,
            ppen,
            tscaleh,
            wtwb,
            cnvtrmax,
            qtu_emf,
            umf_out,
            qtflx_out,
            slflx_out,
            slflx,
            thlu_emf,
            uu_emf,
            vu_emf,
            uemf,
            uflx_out,
            vflx_out,
            ufrc,
            wu,
            emf,
            thlu,
            qtu,
            thvu,
            uu,
            vu,
            umf_zint,
            umf_outvar,
            qtflx_outvar,
            slflx_outvar,
            uflx_outvar,
            vflx_outvar,
            slflx_s,
            qtflx_s,
            uflx_s,
            vflx_s,
            qtflx,
            uflx,
            ufrc_s,
            xflx,
            vflx,
            umf_temp,
            umf_s,
            tke_in,
            pifc0_in,
            zifc0_in,
            exnifc0_in,
            kinv,
            klcl,
            klfc,
            kinv_o,
            klcl_o,
            klfc_o,
            kbup,
            krel,
            kpen,
            kbup_IJ,
            klfc_IJ,
            kpen_IJ,
            kpbl_in,
            tr0_temp,
            u0,
            v0,
            cinlcl,
        )
