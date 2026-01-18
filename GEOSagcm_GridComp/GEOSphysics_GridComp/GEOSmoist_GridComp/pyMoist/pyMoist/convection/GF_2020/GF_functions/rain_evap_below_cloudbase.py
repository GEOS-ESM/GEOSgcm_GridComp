from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField, FloatFieldIJ, IntField, Float

import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
import pyMoist.constants as constants
from gt4py.cartesian.gtscript import (
    sqrt,
    FORWARD,
    PARALLEL,
    BACKWARD,
    K,
    computation,
    interval,
)

# Parameters needed for rain_evap_below_cloudbase
alpha1 = Float(5.44e-4)
alpha2 = Float(5.09e-3)
alpha3 = Float(0.5777)
c_conv = Float(0.05)


def rain_evap_below_cloudbase(
    # In
    cumulus: IntField,
    edto: FloatField,
    ierr: IntField,
    kbcon: IntField,
    ktop: IntField,
    po_cup: FloatField,
    psur: FloatField,
    pwdo: FloatField,
    pwo: FloatField,
    qes_cup: FloatField,
    qo_cup: FloatField,
    xland: FloatField,
    xmb: FloatField,
    q_deficit: FloatFieldIJ,
    RH_cr: FloatFieldIJ,
    dp: FloatFieldIJ,
    # Out
    evap_bcb: FloatField,
    evap_flx: FloatField,
    outbuoy: FloatField,
    outq: FloatField,
    outt: FloatField,
    pre: FloatField,
    prec_flx: FloatField,
):
    with computation(FORWARD), interval(...):
        if cumulus == cumulus_parameterization_constants.SHALLOW:
            RH_cr_OCEAN = 1.0
            RH_cr_LAND = 1.0
            eff_c_conv = min(0.2, max(xmb, c_conv))
        else:
            RH_cr_OCEAN = 0.95
            RH_cr_LAND = 0.85
            eff_c_conv = c_conv

        prec_flx = 0.0
        evap_flx = 0.0
        tot_evap_bcb = 0.0

    with computation(FORWARD), interval(...):
        if ierr == 0:
            RH_cr = RH_cr_OCEAN * xland + RH_cr_LAND * (1.0 - xland)

    with computation(BACKWARD), interval(...):
        if ierr == 0:
            if K <= ktop - 1:
                dp = 100.0 * (po_cup - po_cup[0, 0, 1])

                if K <= kbcon - 1:
                    q_deficit = max(0.0, (RH_cr * qes_cup - qo_cup))

                    evap_bcb = (
                        eff_c_conv
                        * alpha1
                        * q_deficit
                        * (sqrt(po_cup / psur) / alpha2 * prec_flx[0, 0, 1] / eff_c_conv) ** alpha3
                    )

                    evap_bcb = evap_bcb * dp / constants.MAPL_GRAV

                else:

                    evap_bcb = 0.0

                prec_flx = prec_flx[0, 0, 1] - evap_bcb + xmb * (pwo + edto * pwdo)
                prec_flx = max(0.0, prec_flx)

                evap_flx = evap_flx[0, 0, 1] + evap_bcb - xmb * edto * pwdo
                evap_flx = max(0.0, evap_flx)

                tot_evap_bcb = tot_evap_bcb + evap_bcb

                del_q = evap_bcb * constants.MAPL_GRAV / dp
                del_t = (
                    -evap_bcb
                    * constants.MAPL_GRAV
                    / dp
                    * (cumulus_parameterization_constants.xlv / cumulus_parameterization_constants.CP)
                )

                outq = outq + del_q
                outt = outt + del_t
                outbuoy = (
                    outbuoy
                    + cumulus_parameterization_constants.CP * del_t
                    + cumulus_parameterization_constants.xlv * del_q
                )

    with computation(PARALLEL), interval(...):
        if ierr == 0:
            pre = pre - evap_bcb.at(K=0)


class RainEvapBelowCloudbase:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
    ) -> None:

        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory

        self._rain_evap_below_cloudbase = self.stencil_factory.from_dims_halo(
            func=rain_evap_below_cloudbase,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._q_deficit = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM], units="n/a")

        self._RH_cr = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM], units="n/a")

        self._dp = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM], units="n/a")

        self._evap_bcb = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM], units="n/a")

    def __call__(
        self,
        # In
        cumulus: IntField,
        edto: FloatField,
        ierr: IntField,
        kbcon: IntField,
        ktop: IntField,
        po_cup: FloatField,
        psur: FloatField,
        pwavo: FloatField,
        pwdo: FloatField,
        pwevo: FloatField,
        pwo: FloatField,
        qes_cup: FloatField,
        qo_cup: FloatField,
        t_cup: FloatField,
        xland: FloatField,
        xmb: FloatField,
        # Out
        evap_bcb: FloatField,
        evap_flx: FloatField,
        outbuoy: FloatField,
        outq: FloatField,
        outt: FloatField,
        pre: FloatField,
        prec_flx: FloatField,
    ):

        self._rain_evap_below_cloudbase(
            # In
            cumulus=cumulus,
            edto=edto,
            ierr=ierr,
            kbcon=kbcon,
            ktop=ktop,
            po_cup=po_cup,
            psur=psur,
            pwavo=pwavo,
            pwdo=pwdo,
            pwevo=pwevo,
            pwo=pwo,
            qes_cup=qes_cup,
            qo_cup=qo_cup,
            t_cup=t_cup,
            xland=xland,
            xmb=xmb,
            q_deficit=self._q_deficit,
            RH_cr=self._RH_cr,
            dp=self._dp,
            # Out
            evap_bcb=evap_bcb,
            evap_flx=evap_flx,
            outbuoy=outbuoy,
            outq=outq,
            outt=outt,
            pre=pre,
            prec_flx=prec_flx,
        )
