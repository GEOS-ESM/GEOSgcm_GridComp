from dataclasses import dataclass

from ndsl import Quantity, QuantityFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Int


@dataclass
class GF2020DriverTemporaries:
    do_this_column: Quantity
    ierr: Quantity
    jmin: Quantity
    klcl: Quantity
    k22: Quantity
    kbcon: Quantity
    ktop: Quantity
    kstabi: Quantity
    kstabm: Quantity
    xmb: Quantity
    cprr: Quantity
    edt: Quantity
    pwav: Quantity
    sigma: Quantity
    pcup: Quantity
    entr: Quantity
    up_massentr: Quantity
    up_massdetr: Quantity
    dd_massentr: Quantity
    dd_massdetr: Quantity
    zup: Quantity
    zdn: Quantity
    prup: Quantity
    prdn: Quantity
    clwup: Quantity
    tup: Quantity
    conv_cld_fr: Quantity
    SRC_T: Quantity
    SRC_Q: Quantity
    SRC_CI: Quantity
    SRC_U: Quantity
    SRC_V: Quantity
    SUB_MPQI: Quantity
    SUB_MPQL: Quantity
    SUB_MPCF: Quantity
    SRC_BUOY: Quantity
    REVSU_GF: Quantity
    PRFIL_GF: Quantity
    mpqi: Quantity
    mpql: Quantity
    mpcf: Quantity
    fixout_qv: Quantity
    revsu_gf_internal: Quantity
    prfil_gf_internal: Quantity
    outt: Quantity
    outu: Quantity
    outv: Quantity
    outq: Quantity
    outqc: Quantity
    outnice: Quantity
    outnliq: Quantity
    outbuoy: Quantity

    @classmethod
    def make(cls, quantity_factory: QuantityFactory):
        do_this_column = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        ierr = quantity_factory.zeros([X_DIM, Y_DIM, "maxiens"], "n/a", dtype=Int)
        jmin = quantity_factory.zeros([X_DIM, Y_DIM, "maxiens"], "n/a", dtype=Int)
        klcl = quantity_factory.zeros([X_DIM, Y_DIM, "maxiens"], "n/a", dtype=Int)
        k22 = quantity_factory.zeros([X_DIM, Y_DIM, "maxiens"], "n/a", dtype=Int)
        kbcon = quantity_factory.zeros([X_DIM, Y_DIM, "maxiens"], "n/a", dtype=Int)
        ktop = quantity_factory.zeros([X_DIM, Y_DIM, "maxiens"], "n/a", dtype=Int)
        kstabi = quantity_factory.zeros([X_DIM, Y_DIM, "maxiens"], "n/a", dtype=Int)
        kstabm = quantity_factory.zeros([X_DIM, Y_DIM, "maxiens"], "n/a", dtype=Int)
        xmb = quantity_factory.zeros([X_DIM, Y_DIM, "maxiens"], "n/a")
        cprr = quantity_factory.zeros([X_DIM, Y_DIM, "maxiens"], "n/a")
        edt = quantity_factory.zeros([X_DIM, Y_DIM, "maxiens"], "n/a")
        pwav = quantity_factory.zeros([X_DIM, Y_DIM, "maxiens"], "n/a")
        sigma = quantity_factory.zeros([X_DIM, Y_DIM, "maxiens"], "n/a")
        pcup = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM, "maxiens"], "n/a")
        entr = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM, "maxiens"], "n/a")
        up_massentr = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM, "maxiens"], "n/a")
        up_massdetr = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM, "maxiens"], "n/a")
        dd_massentr = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM, "maxiens"], "n/a")
        dd_massdetr = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM, "maxiens"], "n/a")
        zup = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM, "maxiens"], "n/a")
        zdn = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM, "maxiens"], "n/a")
        prup = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM, "maxiens"], "n/a")
        prdn = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM, "maxiens"], "n/a")
        clwup = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM, "maxiens"], "n/a")
        tup = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM, "maxiens"], "n/a")
        conv_cld_fr = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM, "maxiens"], "n/a")
        SRC_T = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        SRC_Q = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        SRC_CI = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        SRC_U = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        SRC_V = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        SUB_MPQI = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM, "nmp"], "n/a")
        SUB_MPQL = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM, "nmp"], "n/a")
        SUB_MPCF = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM, "nmp"], "n/a")
        SRC_BUOY = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        REVSU_GF = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        PRFIL_GF = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        mpqi = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM, "nmp"], "n/a")
        mpql = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM, "nmp"], "n/a")
        mpcf = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM, "nmp"], "n/a")

        ztexec = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        zqexec = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        fixout_qv = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        fixout_qv.field[:] = 1.0
        revsu_gf_internal = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        prfil_gf_internal = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        outt = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM, "maxiens"], "n/a")
        outu = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM, "maxiens"], "n/a")
        outv = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM, "maxiens"], "n/a")
        outq = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM, "maxiens"], "n/a")
        outqc = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM, "maxiens"], "n/a")
        outnice = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM, "maxiens"], "n/a")
        outnliq = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM, "maxiens"], "n/a")
        outbuoy = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM, "maxiens"], "n/a")
        psur = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        tsur = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        ter1l = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        xlons = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        xlats = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")

        return cls(
            do_this_column,
            ierr,
            jmin,
            klcl,
            k22,
            kbcon,
            ktop,
            kstabi,
            kstabm,
            xmb,
            cprr,
            edt,
            pwav,
            sigma,
            pcup,
            entr,
            up_massentr,
            up_massdetr,
            dd_massentr,
            dd_massdetr,
            zup,
            zdn,
            prup,
            prdn,
            clwup,
            tup,
            conv_cld_fr,
            SRC_T,
            SRC_Q,
            SRC_CI,
            SRC_U,
            SRC_V,
            SUB_MPQI,
            SUB_MPQL,
            SUB_MPCF,
            SRC_BUOY,
            REVSU_GF,
            PRFIL_GF,
            mpqi,
            mpql,
            mpcf,
            fixout_qv,
            revsu_gf_internal,
            prfil_gf_internal,
            outt,
            outu,
            outv,
            outq,
            outqc,
            outnice,
            outnliq,
            outbuoy,
        )


@dataclass
class GF2020DriverTemporaries:
    last_ierr: Quantity

    @classmethod
    def make(cls, quantity_factory: QuantityFactory):
        last_ierr = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")

        return cls(last_ierr)
