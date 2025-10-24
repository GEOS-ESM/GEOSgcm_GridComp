from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField, IntField, FloatFieldIJ, IntFieldIJ, Float
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
from gt4py.cartesian.gtscript import (
    PARALLEL,
    FORWARD,
    K,
    computation,
    interval,
    int32,
    abs,
)

beta3 = Float(-1.13)
alpha3 = Float(1.9)


def cup_dd_edt(
    # In
    ccn: FloatField,
    cumulus: IntField,
    edtmax: FloatField,
    edtmin: FloatField,
    kbcon: IntField,
    ktop: IntField,
    maxens2: IntField,
    p: FloatField,
    psum2: FloatField,
    psumh: FloatField,
    pw: FloatField,
    pwav: FloatField,
    pwev: FloatField,
    rho: FloatField,
    us: FloatField,
    vs: FloatField,
    z: FloatField,
    aeroevap: IntField,
    sdp: FloatFieldIJ,
    vshear: FloatFieldIJ,
    vws: FloatFieldIJ,
    pef: FloatFieldIJ,
    dp: FloatFieldIJ,
    # Out
    edt: FloatField,
    edtc: FloatField,
    ierr: IntField,
):

    with computation(FORWARD), interval(...):
        edt = 0.0
        vws = 0.0
        sdp = 0.0
        vshear = 0.0
        edtc = 0.0
    with computation(FORWARD), interval(...):
        if cumulus != cumulus_parameterization_constants.shallow:
            if ierr == 0:
                idx = kbcon - 1
                while idx >= kbcon - 1 and idx <= ktop - 1:
                    dp = p.at(K=idx) - p.at(K=idx + 1)
                    vws = vws + (
                        (
                            abs(
                                (us.at(K=idx + 1) - us.at(K=idx))
                                / (z.at(K=idx + 1) - z.at(K=idx))
                            )
                            + abs(
                                (vs.at(K=idx + 1) - vs.at(K=idx))
                                / (z.at(K=idx + 1) - z.at(K=idx))
                            )
                        )
                        * dp
                    )
                    sdp = sdp + dp
                    idx += 1
                vshear = 1.0e3 * vws / sdp

    with computation(FORWARD), interval(...):
        if cumulus != cumulus_parameterization_constants.shallow:
            if ierr == 0:
                pef = (
                    1.591
                    - 0.639 * vshear
                    + 0.0953 * (vshear * vshear)
                    - 0.00496 * (vshear * vshear * vshear)
                )
                pef = min(pef, 0.9)
                pef = max(pef, 0.1)
                zkbc = z.at(K=kbcon - 1) * 3.281e-3
                prezk = 0.02

                if zkbc > 3.0:
                    prezk = 0.96729352 + zkbc * (
                        -0.70034167
                        + zkbc
                        * (
                            0.162179896
                            + zkbc
                            * (-1.2569798e-2 + zkbc * (4.2772e-4 - zkbc * 5.44e-6))
                        )
                    )

                if zkbc > 25.0:
                    prezk = 2.4

                pefb = 1.0 / (1.0 + prezk)
                pefb = min(pefb, 0.9)
                pefb = max(pefb, 0.1)

                edt = 1.0 - 0.5 * (pefb + pef)

                # if aeroevap > 1:
                #     aeroadd = (cumulus_parameterization_constants.ccnclean**beta3) * (
                #         (psumh) ** (alpha3 - 1)
                #     )

                #     prop_c = 0.5 * (pefb + pef) / aeroadd
                #     aeroadd = (ccn**beta3) * ((psum2) ** (alpha3 - 1))

                #     aeroadd = prop_c * aeroadd
                #     pefc = aeroadd
                #     if pefc > 0.9:
                #         pefc = 0.9
                #     if pefc < 0.1:
                #         pefc = 0.1
                #     edt = 1.0 - pefc
                #     if aeroevap == 2:
                #         edt = 1.0 - 0.25 * (pefb + pef + 2.0 * pefc)

                einc = 0.2 * edt
                # if K <= maxens2 - 1:
                edtc = edt + float(int32(K) - 1) * einc

    with computation(PARALLEL), interval(...):
        if cumulus != cumulus_parameterization_constants.shallow:
            if ierr == 0:
                # if K <= maxens2 - 1:
                edtc = -edtc * pwav / pwev
                if edtc > edtmax:
                    edtc = edtmax
                if edtc < edtmin:
                    edtc = edtmin

                temp = edtc.at(K=0)
                edtc = temp


class CupDDEdt:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
    ) -> None:

        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory

        self._cup_dd_edt = self.stencil_factory.from_dims_halo(
            func=cup_dd_edt,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._sdp = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM], units="n/a"
        )
        self._vshear = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM], units="n/a"
        )
        self._vws = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM], units="n/a"
        )
        self._pef = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM], units="n/a"
        )
        self._dp = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM], units="n/a"
        )

    def __call__(
        self,
        # In
        ccn: FloatField,
        cumulus: IntField,
        edtmax: FloatField,
        edtmin: FloatField,
        kbcon: IntField,
        ktop: IntField,
        maxens2: IntField,
        p: FloatField,
        psum2: FloatField,
        psumh: FloatField,
        pw: FloatField,
        pwav: FloatField,
        pwev: FloatField,
        rho: FloatField,
        us: FloatField,
        vs: FloatField,
        z: FloatField,
        aeroevap: IntField,
        # Out
        edt: FloatField,
        edtc: FloatField,
        ierr: IntField,
    ):

        self._cup_dd_edt(
            # In
            ccn=ccn,
            cumulus=cumulus,
            edtmax=edtmax,
            edtmin=edtmin,
            kbcon=kbcon,
            ktop=ktop,
            maxens2=maxens2,
            p=p,
            psum2=psum2,
            psumh=psumh,
            pw=pw,
            pwav=pwav,
            pwev=pwev,
            rho=rho,
            us=us,
            vs=vs,
            z=z,
            aeroevap=aeroevap,
            sdp=self._sdp,
            vshear=self._vshear,
            vws=self._vws,
            pef=self._pef,
            dp=self._dp,
            # Out
            edt=edt,
            edtc=edtc,
            ierr=ierr,
        )
