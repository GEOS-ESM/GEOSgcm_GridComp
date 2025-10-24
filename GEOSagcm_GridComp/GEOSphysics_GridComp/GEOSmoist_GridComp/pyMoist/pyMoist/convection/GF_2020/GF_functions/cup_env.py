from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField, IntField, Float, Int
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
import pyMoist.constants as constants
from gt4py.cartesian.gtscript import (
    PARALLEL,
    computation,
    interval,
    int32,
    log,
    exp,
)
import gt4py.cartesian.gtscript as gtscript


@gtscript.function
def satvap(temp2: Float):
    # CAUTION: This function has not been verified!!!
    temp = temp2 - 273.155
    if temp < -20.0:
        toot = 273.16 / temp2
        toto = 1 / toot
        eilog = (
            -9.09718 * (toot - 1)
            - 3.56654 * (log(toot) / log(10.0))
            + 0.876793 * (1 - toto)
            + (log(6.1071) / log(10.0))
        )
        satvap = 10**eilog
    else:
        tsot = 373.16 / temp2
        ewlog = -7.90298 * (tsot - 1) + 5.02808 * (log(tsot) / log(10.0))
        ewlog2 = ewlog - 1.3816e-07 * (10 ** (11.344 * (1 - (1 / tsot))) - 1)
        ewlog3 = ewlog2 + 0.0081328 * (10 ** (-3.49149 * (tsot - 1)) - 1)
        ewlog4 = ewlog3 + (log(1013.246) / log(10.0))
        satvap = 10**ewlog4

    return satvap


# Parameters needed for satur_spec_hum
rd = 287.06
rv = 461.52
rtt = 273.16
retv = rv / rd - 1.0
r2es = 611.21 * rd / rv
r3les = 17.502
r3ies = 22.587
r4les = 32.19
r4ies = -0.7
rtwat = rtt
rtice = rtt - 23.0
rticecu = rtt - 23.0
rtwat_rtice_r = 1.0 / (rtwat - rtice)
rtwat_rticecu_r = 1.0 / (rtwat - rticecu)


@gtscript.function
def satur_spec_hum(
    pt: Float,
    press: Float,
):
    foealfcu = min(
        1.0, ((max(rticecu, min(rtwat, pt)) - rticecu) * rtwat_rticecu_r) ** 2
    )
    foeewmcu = r2es * (
        foealfcu * exp(r3les * (pt - rtt) / (pt - r4les))
        + (1.0 - foealfcu) * exp(r3ies * (pt - rtt) / (pt - r4ies))
    )

    zew = foeewmcu
    zqs = zew / (100.0 * press)
    if 1.0 - retv * zqs > 0.0:
        zcor = 1.0 / (1.0 - retv * zqs)
        pqsat = zqs * zcor
    else:
        pqsat = cumulus_parameterization_constants.max_qsat

    return pqsat


def cup_env(
    # In
    SATUR_CALC: IntField,
    itest: IntField,
    p: FloatField,
    psur: FloatField,
    q: FloatField,
    t: FloatField,
    z1: FloatField,
    # Out
    he: FloatField,
    hes: FloatField,
    ierr: IntField,
    qes: FloatField,
    z: FloatField,
):
    # CAUTION: This stencil has only been ported for the cases SATUR_CALC = 1
    # If this is not the case, an error will raised because that code has not been ported.
    with computation(PARALLEL), interval(...):
        he = 0.0
        hes = 0.0
        qes = 0.0
    with computation(PARALLEL), interval(0, -1):
        if SATUR_CALC == 0:
            if ierr == 0:
                e = satvap(t)
                qes = 0.622 * e / max(1.0e-8, (p - e))

                if qes <= 1.0e-08:
                    qes = 1.0e-08
                if qes > cumulus_parameterization_constants.max_qsat:
                    qes = cumulus_parameterization_constants.max_qsat
                if qes < q:
                    qes = q

                TV = t + 0.608 * q * t
        else:
            if ierr == 0:
                qes = satur_spec_hum(t, p)
                qes = min(
                    cumulus_parameterization_constants.max_qsat, max(1.0e-08, qes)
                )
                qes = max(qes, q)
                tv = t + 0.608 * q * t

    with computation(PARALLEL), interval(0, -1):
        if ierr == 0:
            if itest <= 0:
                he = (
                    constants.MAPL_GRAV * z
                    + cumulus_parameterization_constants.CP * t
                    + cumulus_parameterization_constants.xlv * q
                )

            hes = (
                constants.MAPL_GRAV * z
                + cumulus_parameterization_constants.CP * t
                + cumulus_parameterization_constants.xlv * qes
            )

            if he >= hes:
                he = hes


class CupEnv:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
    ) -> None:

        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory

        self._cup_env = self.stencil_factory.from_dims_halo(
            func=cup_env,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        # In
        SATUR_CALC: IntField,
        itest: IntField,
        p: FloatField,
        psur: FloatField,
        q: FloatField,
        t: FloatField,
        z1: FloatField,
        # Out
        he: FloatField,
        hes: FloatField,
        ierr: IntField,
        qes: FloatField,
        z: FloatField,
    ):

        if SATUR_CALC.view[:].all() != Int(1):
            raise NotImplementedError(
                f"Warning: This code has not been ported!! Expected SATUR_CALC = 1"
            )

        self._cup_env(
            # In
            SATUR_CALC=SATUR_CALC,
            itest=itest,
            p=p,
            psur=psur,
            q=q,
            t=t,
            z1=z1,
            # Out
            he=he,
            hes=hes,
            ierr=ierr,
            qes=qes,
            z=z,
        )
