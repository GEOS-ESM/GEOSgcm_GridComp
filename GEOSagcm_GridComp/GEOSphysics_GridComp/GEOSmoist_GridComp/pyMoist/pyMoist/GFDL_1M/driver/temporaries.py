from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl import Quantity, QuantityFactory
from dataclasses import dataclass


@dataclass
class Temporaries:
    t1: Quantity
    dp1: Quantity
    omq: Quantity
    qv0: Quantity
    ql0: Quantity
    qr0: Quantity
    qi0: Quantity
    qs0: Quantity
    qg0: Quantity
    qa0: Quantity
    qv1: Quantity
    ql1: Quantity
    qr1: Quantity
    qi1: Quantity
    qs1: Quantity
    qg1: Quantity
    qa1: Quantity
    dz1: Quantity
    den: Quantity
    den1: Quantity
    denfac: Quantity
    p_dry: Quantity
    m1: Quantity
    u1: Quantity
    v1: Quantity
    w1: Quantity
    onemsig: Quantity
    ccn: Quantity
    c_praut: Quantity
    rh_limited: Quantity
    ze: Quantity
    zt: Quantity
    lhi: Quantity
    icpk: Quantity
    hold_data: Quantity
    vti: Quantity
    vts: Quantity
    vtg: Quantity
    vtr: Quantity
    m1_sol: Quantity
    m1_rain: Quantity
    rain1: Quantity
    graupel1: Quantity
    snow1: Quantity
    ice1: Quantity
    evap1: Quantity
    subl1: Quantity

    @classmethod
    def make(cls, quantity_factory: QuantityFactory):
        t1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dp1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        omq = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qv0 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        ql0 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qr0 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qi0 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qs0 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qg0 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qa0 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qv1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        ql1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qr1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qi1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qs1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qg1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qa1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dz1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        den = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        den1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        denfac = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        p_dry = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        m1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        u1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        v1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        w1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        onemsig = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        ccn = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        c_praut = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        rh_limited = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        ze = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a")
        zt = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a")
        lhi = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        icpk = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        hold_data = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        vti = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        vts = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        vtg = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        vtr = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        m1_sol = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        m1_rain = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        rain1 = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        graupel1 = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        snow1 = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        ice1 = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        evap1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        subl1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        return cls(
            t1,
            dp1,
            omq,
            qv0,
            ql0,
            qr0,
            qi0,
            qs0,
            qg0,
            qa0,
            qv1,
            ql1,
            qr1,
            qi1,
            qs1,
            qg1,
            qa1,
            dz1,
            den,
            den1,
            denfac,
            p_dry,
            m1,
            u1,
            v1,
            w1,
            onemsig,
            ccn,
            c_praut,
            rh_limited,
            ze,
            zt,
            lhi,
            icpk,
            hold_data,
            vti,
            vts,
            vtg,
            vtr,
            m1_sol,
            m1_rain,
            rain1,
            graupel1,
            snow1,
            ice1,
            evap1,
            subl1,
        )
