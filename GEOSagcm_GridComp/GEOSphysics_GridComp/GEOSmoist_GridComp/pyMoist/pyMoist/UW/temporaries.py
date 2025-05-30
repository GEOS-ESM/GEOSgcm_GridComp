from dataclasses import dataclass

from ndsl import Quantity, QuantityFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Int


@dataclass
class Temporaries:
    ssthl0: Quantity
    ssqt0: Quantity
    ssu0: Quantity
    ssv0: Quantity
    thj: Quantity
    qlj: Quantity
    qvj: Quantity
    qse: Quantity
    qij: Quantity
    thv0top: Quantity
    thv0bot: Quantity
    thvl0top: Quantity
    dcm_out: Quantity
    qvten_out: Quantity
    qlten_out: Quantity
    qiten_out: Quantity
    sten_out: Quantity
    uten_out: Quantity
    vten_out: Quantity
    qrten_out: Quantity
    qsten_out: Quantity
    cufrc_out: Quantity
    fer_out: Quantity
    fdr_out: Quantity
    thvlavg: Quantity
    tkeavg: Quantity
    uavg: Quantity
    vavg: Quantity
    thvlmin: Quantity
    qtavg: Quantity
    zmid0: Quantity
    qt0: Quantity
    thvl0: Quantity
    thvl0bot: Quantity
    t0: Quantity
    qv0: Quantity
    pmid0: Quantity
    thl0: Quantity
    thlsrc: Quantity
    usrc: Quantity
    vsrc: Quantity
    plcl: Quantity
    thl0lcl: Quantity
    qt0lcl: Quantity
    thv0lcl: Quantity
    plfc: Quantity
    fer_outvar: Quantity
    fdr_outvar: Quantity
    cin: Quantity
    thvubot: Quantity
    thvutop: Quantity
    thvlsrc: Quantity
    thl0top: Quantity
    qt0top: Quantity
    qldet_outvar: Quantity
    qidet_outvar: Quantity
    qlsub_outvar: Quantity
    qisub_outvar: Quantity
    dcm_outvar: Quantity
    qvten_outvar: Quantity
    qlten_outvar: Quantity
    qiten_outvar: Quantity
    sten_outvar: Quantity
    uten_outvar: Quantity
    vten_outvar: Quantity
    qrten_outvar: Quantity
    qsten_outvar: Quantity
    cufrc_outvar: Quantity
    usrc_o: Quantity
    vsrc_o: Quantity
    thv0lcl_o: Quantity
    ql0_o: Quantity
    qi0_o: Quantity
    t0_o: Quantity
    s0_o: Quantity
    u0_o: Quantity
    v0_o: Quantity
    qt0_o: Quantity
    thl0_o: Quantity
    thvl0_o: Quantity
    ssthl0_o: Quantity
    ssqt0_o: Quantity
    thv0bot_o: Quantity
    thv0top_o: Quantity
    thvl0bot_o: Quantity
    thvl0top_o: Quantity
    ssu0_o: Quantity
    ssv0_o: Quantity
    dcm_s: Quantity
    qvten_s: Quantity
    qlten_s: Quantity
    qiten_s: Quantity
    sten_s: Quantity
    uten_s: Quantity
    vten_s: Quantity
    qrten_s: Quantity
    qsten_s: Quantity
    qldet_s: Quantity
    qidet_s: Quantity
    qlsub_s: Quantity
    qisub_s: Quantity
    cush_s: Quantity
    cufrc_s: Quantity
    fer_s: Quantity
    fdr_s: Quantity
    qtsrc_o: Quantity
    thvlsrc_o: Quantity
    thlsrc_o: Quantity
    qldet_out: Quantity
    qidet_out: Quantity
    qlsub_out: Quantity
    qisub_out: Quantity
    ndrop_out: Quantity
    nice_out: Quantity
    dcm: Quantity
    xco: Quantity
    qc: Quantity
    qlten_det: Quantity
    qiten_det: Quantity
    qv0_s: Quantity
    ql0_s: Quantity
    qi0_s: Quantity
    s0_s: Quantity
    t0_s: Quantity
    slten: Quantity
    qv0_o: Quantity
    plcl_o: Quantity
    plfc_o: Quantity
    tkeavg_o: Quantity
    thvlmin_o: Quantity
    ufrclcl: Quantity
    qcu: Quantity
    qlu: Quantity
    qiu: Quantity
    cufrc: Quantity
    qtsrc: Quantity
    uplus_3D: Quantity
    vplus_3D: Quantity
    prel: Quantity
    thv0rel: Quantity
    winv: Quantity
    cbmf: Quantity
    rho0inv: Quantity
    ufrcinv: Quantity
    wlcl: Quantity
    qsat_pe: Quantity
    thlue: Quantity
    qtue: Quantity
    wue: Quantity
    rei: Quantity
    fer: Quantity
    dwten: Quantity
    diten: Quantity
    ql0: Quantity
    qi0: Quantity
    uten: Quantity
    vten: Quantity
    uf: Quantity
    vf: Quantity
    dwten_temp: Quantity
    diten_temp: Quantity
    fdr: Quantity
    qlten_sink: Quantity
    qiten_sink: Quantity
    qrten: Quantity
    qsten: Quantity
    s0: Quantity
    qvten: Quantity
    qlten: Quantity
    sten: Quantity
    qiten: Quantity
    qmin: Quantity
    pmid0_in: Quantity
    u0_in: Quantity
    v0_in: Quantity
    zmid0_in: Quantity
    exnmid0_in: Quantity
    dp0_in: Quantity
    qv0_in: Quantity
    ql0_in: Quantity
    qi0_in: Quantity
    th0_in: Quantity
    cush_inout: Quantity
    dpi: Quantity
    thvlmin_IJ: Quantity
    wcrit: Quantity
    alpha: Quantity
    del_CIN: Quantity
    cin_IJ: Quantity
    plfc_IJ: Quantity
    cinlcl_IJ: Quantity
    pe: Quantity
    thle: Quantity
    qte: Quantity
    dpe: Quantity
    exne: Quantity
    thvebot: Quantity
    ue: Quantity
    ve: Quantity
    drage: Quantity
    bogbot: Quantity
    bogtop: Quantity
    rhomid0j: Quantity
    cush_inoutvar: Quantity
    uplus: Quantity
    vplus: Quantity
    cin_i: Quantity
    cinlcl_i: Quantity
    ke: Quantity
    thlu_top: Quantity
    qtu_top: Quantity
    cldhgt: Quantity
    qlubelow: Quantity
    qiubelow: Quantity
    qlj_2D: Quantity
    qij_2D: Quantity
    qcubelow: Quantity
    rcwp: Quantity
    rlwp: Quantity
    riwp: Quantity
    ppen: Quantity
    tscaleh: Quantity
    wtwb: Quantity
    cnvtrmax: Quantity
    qtu_emf: Quantity
    umf_out: Quantity
    qtflx_out: Quantity
    slflx_out: Quantity
    slflx: Quantity
    thlu_emf: Quantity
    uu_emf: Quantity
    vu_emf: Quantity
    uemf: Quantity
    uflx_out: Quantity
    vflx_out: Quantity
    ufrc: Quantity
    wu: Quantity
    emf: Quantity
    thlu: Quantity
    qtu: Quantity
    thvu: Quantity
    uu: Quantity
    vu: Quantity
    umf_zint: Quantity
    umf_outvar: Quantity
    qtflx_outvar: Quantity
    slflx_outvar: Quantity
    uflx_outvar: Quantity
    vflx_outvar: Quantity
    slflx_s: Quantity
    qtflx_s: Quantity
    uflx_s: Quantity
    vflx_s: Quantity
    qtflx: Quantity
    uflx: Quantity
    ufrc_s: Quantity
    xflx: Quantity
    vflx: Quantity
    umf_temp: Quantity
    umf_s: Quantity
    tke_in: Quantity
    pifc0_in: Quantity
    zifc0_in: Quantity
    exnifc0_in: Quantity
    kinv: Quantity
    klcl: Quantity
    klfc: Quantity
    kinv_o: Quantity
    klcl_o: Quantity
    klfc_o: Quantity
    kbup: Quantity
    krel: Quantity
    kpen: Quantity
    kbup_IJ: Quantity
    klfc_IJ: Quantity
    kpen_IJ: Quantity
    kpbl_in: Quantity
    tr0_temp: Quantity
    u0: Quantity
    v0: Quantity
    cinlcl: Quantity

    @classmethod
    def make(cls, quantity_factory: QuantityFactory):
        # FloatFields
        ssthl0 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        ssqt0 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        ssu0 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        ssv0 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        thj = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qlj = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qvj = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qse = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qij = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        tr0_temp = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        thv0top = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        thv0bot = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        thvl0top = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dcm_out = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qvten_out = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qlten_out = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qiten_out = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        sten_out = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        uten_out = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        vten_out = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qrten_out = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qsten_out = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        cufrc_out = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        fer_out = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        fdr_out = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        thvlavg = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        tkeavg = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        uavg = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        vavg = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        thvlmin = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qtavg = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        zmid0 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qt0 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        thvl0 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        thvl0bot = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        t0 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qv0 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        pmid0 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        thl0 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        thlsrc = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        usrc = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        vsrc = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        plcl = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        thl0lcl = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qt0lcl = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        thv0lcl = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        plfc = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        fer_outvar = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        fdr_outvar = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        cin = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        thvubot = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        thvutop = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        thvlsrc = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        thl0top = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qt0top = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qldet_outvar = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qidet_outvar = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qlsub_outvar = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qisub_outvar = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dcm_outvar = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qvten_outvar = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qlten_outvar = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qiten_outvar = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        sten_outvar = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        uten_outvar = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        vten_outvar = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qrten_outvar = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qsten_outvar = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        cufrc_outvar = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        usrc_o = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        vsrc_o = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        thv0lcl_o = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        ql0_o = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qi0_o = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        t0_o = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        s0_o = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        u0_o = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        v0_o = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qt0_o = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        thl0_o = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        thvl0_o = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        ssthl0_o = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        ssqt0_o = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        thv0bot_o = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        thv0top_o = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        thvl0bot_o = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        thvl0top_o = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        ssu0_o = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        ssv0_o = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dcm_s = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qvten_s = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qlten_s = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qiten_s = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        sten_s = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        uten_s = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        vten_s = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qrten_s = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qsten_s = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qldet_s = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qidet_s = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qlsub_s = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qisub_s = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        cush_s = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        cufrc_s = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        fer_s = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        fdr_s = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qtsrc_o = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        thvlsrc_o = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        thlsrc_o = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qldet_out = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qidet_out = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qlsub_out = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qisub_out = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        ndrop_out = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        nice_out = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dcm = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        xco = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qc = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qlten_det = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qiten_det = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qv0_s = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        ql0_s = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qi0_s = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        s0_s = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        t0_s = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        slten = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qv0_o = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        plcl_o = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        plfc_o = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        tkeavg_o = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        thvlmin_o = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        ufrclcl = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qcu = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qlu = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qiu = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        cufrc = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qtsrc = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        uplus_3D = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        vplus_3D = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        prel = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        thv0rel = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        winv = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        cbmf = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        rho0inv = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        ufrcinv = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        wlcl = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qsat_pe = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        thlue = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qtue = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        wue = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        rei = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        fer = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dwten = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        diten = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        ql0 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qi0 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        uten = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        vten = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        uf = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        vf = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dwten_temp = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        diten_temp = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        fdr = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qlten_sink = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qiten_sink = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qrten = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qsten = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        s0 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qvten = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qlten = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        sten = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qiten = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qmin = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        pmid0_in = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        u0_in = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        v0_in = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        u0 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        v0 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        zmid0_in = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        exnmid0_in = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dp0_in = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qv0_in = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        ql0_in = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qi0_in = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        th0_in = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        cinlcl = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        # FloatFieldIJs
        cush_inout = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        dpi = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        thvlmin_IJ = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        wcrit = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        alpha = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        del_CIN = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        cin_IJ = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        plfc_IJ = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        cinlcl_IJ = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        pe = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        thle = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        qte = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        dpe = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        exne = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        thvebot = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        ue = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        ve = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        drage = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        bogbot = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        bogtop = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        rhomid0j = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        cush_inoutvar = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        uplus = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        vplus = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        cin_i = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        cinlcl_i = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        ke = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        thlu_top = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        qtu_top = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        cldhgt = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        qlubelow = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        qiubelow = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        qlj_2D = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        qij_2D = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        qcubelow = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        rcwp = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        rlwp = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        riwp = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        ppen = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        tscaleh = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        wtwb = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        cnvtrmax = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        # Interface FloatFields
        qtu_emf = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        umf_out = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        qtflx_out = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        slflx_out = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        slflx = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        thlu_emf = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        uu_emf = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        vu_emf = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        uemf = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        uflx_out = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        vflx_out = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        ufrc = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        wu = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        emf = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        thlu = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        qtu = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        uu = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        vu = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        umf_zint = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        thvu = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        umf_outvar = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        qtflx_outvar = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        slflx_outvar = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        uflx_outvar = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        vflx_outvar = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        slflx_s = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        qtflx_s = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        uflx_s = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        vflx_s = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        qtflx = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        uflx = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        ufrc_s = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        xflx = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        vflx = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        umf_temp = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        umf_s = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        tke_in = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        pifc0_in = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        zifc0_in = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        exnifc0_in = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        # IntFields
        kinv = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int)
        klcl = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int)
        klfc = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int)
        kinv_o = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int)
        klcl_o = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int)
        klfc_o = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int)
        kbup = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int)
        krel = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int)
        kpen = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int)

        # IntFieldIJs
        kbup_IJ = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a", dtype=Int)
        klfc_IJ = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a", dtype=Int)
        kpen_IJ = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a", dtype=Int)
        kpbl_in = quantity_factory.zeros([X_DIM, Y_DIM], units="n/a", dtype=Int)

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
