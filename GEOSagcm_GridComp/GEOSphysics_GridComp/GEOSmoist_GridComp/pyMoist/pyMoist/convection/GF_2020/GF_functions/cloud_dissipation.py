from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField, IntField, Float, Int

import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
import pyMoist.constants as constants
from gt4py.cartesian.gtscript import (
    BACKWARD,
    K,
    computation,
    interval,
)

# Parameters needed for cloud_dissipation
cloud_lifetime = Float(1800.0)
versionx = Int(2)


def cloud_dissipation(
    # In
    COUPL_MPHYSICS: IntField,
    dtime: FloatField,
    heso_cup: FloatField,
    ierr: IntField,
    kbcon: IntField,
    ktop: IntField,
    qeso_cup: FloatField,
    qo_cup: FloatField,
    rho_hydr: FloatField,
    sig: FloatField,
    tempco: FloatField,
    tn_cup: FloatField,
    use_cloud_dissipation: FloatField,
    vvel2d: FloatField,
    zo: FloatField,
    zuo: FloatField,
    xmb: FloatField,
    # Out
    qrco: FloatField,
    outqc: FloatField,
    outq: FloatField,
    outt: FloatField,
):
    with computation(BACKWARD), interval(...):
        if ierr == 0:
            if K <= ktop - 1 and K >= kbcon - 1:

                qrc_diss = max(0.0, qrco - outqc * dtime)

                frh = 0.0

                fractional_area = (xmb / sig) * zuo / (rho_hydr * vvel2d)

                outqc_diss = (qrc_diss * (1.0 - frh)) / cloud_lifetime

                if versionx == 1 or COUPL_MPHYSICS == 0:

                    outt_diss = -outqc_diss * (
                        cumulus_parameterization_constants.xlv
                        / cumulus_parameterization_constants.CP
                    )

                    qvx = qeso_cup
                    tempx = (
                        heso_cup
                        - constants.MAPL_GRAV * zo
                        - cumulus_parameterization_constants.xlv * qeso_cup
                    ) / cumulus_parameterization_constants.CP

                    outq_mix = (qvx - qo_cup) / cloud_lifetime

                    outt_mix = (tempx - tn_cup) / cloud_lifetime

                    del_q = (
                        (outqc_diss + outq_mix)
                        * use_cloud_dissipation
                        * fractional_area
                    )
                    del_t = (
                        (outt_diss + outt_mix) * use_cloud_dissipation * fractional_area
                    )

                    outq = outq + del_q
                    outt = outt + del_t

                else:

                    outqc = outqc + outqc_diss * fractional_area * use_cloud_dissipation

                qrco = max(
                    0.0,
                    qrco - outqc_diss * use_cloud_dissipation * fractional_area * dtime,
                )


class CloudDissipation:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
    ) -> None:

        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory

        self._cloud_dissipation = self.stencil_factory.from_dims_halo(
            func=cloud_dissipation,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        # In
        COUPL_MPHYSICS: IntField,
        dtime: FloatField,
        heso_cup: FloatField,
        ierr: IntField,
        kbcon: IntField,
        ktop: IntField,
        qeso_cup: FloatField,
        qo_cup: FloatField,
        rho_hydr: FloatField,
        sig: FloatField,
        tempco: FloatField,
        tn_cup: FloatField,
        use_cloud_dissipation: FloatField,
        vvel2d: FloatField,
        zo: FloatField,
        zuo: FloatField,
        xmb: FloatField,
        # Out
        qrco: FloatField,
        outqc: FloatField,
        outq: FloatField,
        outt: FloatField,
    ):
        if COUPL_MPHYSICS.view[:].all() == Int(0):
            raise NotImplementedError(
                f"Warning: This code has not been ported!! Expected COUPL_MPHYSICS > 0"
            )

        self._cloud_dissipation(
            # In
            COUPL_MPHYSICS=COUPL_MPHYSICS,
            dtime=dtime,
            heso_cup=heso_cup,
            ierr=ierr,
            kbcon=kbcon,
            ktop=ktop,
            qeso_cup=qeso_cup,
            qo_cup=qo_cup,
            rho_hydr=rho_hydr,
            sig=sig,
            tempco=tempco,
            tn_cup=tn_cup,
            use_cloud_dissipation=use_cloud_dissipation,
            vvel2d=vvel2d,
            zo=zo,
            zuo=zuo,
            xmb=xmb,
            # Out
            qrco=qrco,
            outqc=outqc,
            outq=outq,
            outt=outt,
        )
