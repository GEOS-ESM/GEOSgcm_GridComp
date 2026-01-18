from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField, IntField, FloatFieldIJ, IntFieldIJ
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
from gt4py.cartesian.gtscript import (
    PARALLEL,
    FORWARD,
    K,
    computation,
    interval,
)


def rates_up_pdf(
    # In
    OVERSHOOT: FloatField,
    entr_rate: FloatField,
    heo: FloatField,
    heso_cup: FloatField,
    hkbo: FloatField,
    kbcon: IntField,
    name: IntField,
    z_cup: FloatField,
    hcot: FloatField,
    ktop_flag1: FloatFieldIJ,
    ktop_flag2: FloatFieldIJ,
    # Out
    ierr: IntField,
    ktopIJ: IntFieldIJ,
    ktop: IntField,
):
    from __externals__ import k_end

    with computation(FORWARD), interval(...):
        delz_oversh = OVERSHOOT
        hcot = 0.0
        ktop_flag1 = 0.0
        ktop_flag2 = 0.0

    with computation(FORWARD), interval(...):
        if name != cumulus_parameterization_constants.SHALLOW:
            ktopIJ = k_end - 3

    with computation(PARALLEL), interval(...):
        if name != cumulus_parameterization_constants.SHALLOW:
            if ierr == 0:
                start_level = kbcon - 1
                if K <= start_level:
                    hcot = hkbo

    with computation(FORWARD), interval(1, None):
        if name != cumulus_parameterization_constants.SHALLOW:
            if ierr == 0:
                idx = start_level + 1
                while K >= idx and K <= k_end - 3:
                    dz = z_cup - z_cup[0, 0, -1]

                    hcot = (
                        (1.0 - 0.5 * entr_rate[0, 0, -1] * dz) * hcot[0, 0, -1]
                        + entr_rate[0, 0, -1] * dz * heo[0, 0, -1]
                    ) / (1.0 + 0.5 * entr_rate[0, 0, -1] * dz)
                    idx += 1

    with computation(FORWARD), interval(...):
        if name != cumulus_parameterization_constants.SHALLOW:
            if ierr == 0:
                idx = start_level + 1
                while K >= idx and K <= k_end - 3:
                    if hcot < heso_cup and ktop_flag1 == 0.0:
                        ktopIJ = K - 1
                        ktop_flag1 = 1.0
                    idx += 1

    with computation(PARALLEL), interval(...):
        if name != cumulus_parameterization_constants.SHALLOW:
            if ierr == 0:
                if ktopIJ <= kbcon:
                    ierr = 41

    with computation(FORWARD), interval(...):
        if name != cumulus_parameterization_constants.SHALLOW:
            if OVERSHOOT > 0.0 and ierr == 0:
                Z_overshoot = (1.0 + delz_oversh) * z_cup.at(K=ktopIJ)
                idx = ktopIJ
                while K >= idx and K <= k_end - 3:
                    if Z_overshoot < z_cup and ktop_flag2 == 0.0:
                        ktopIJ = min(K - 1, k_end - 3)
                        ktop_flag2 = 1.0
                    idx += 1

    with computation(PARALLEL), interval(...):
        if name != cumulus_parameterization_constants.SHALLOW:
            ktop = ktopIJ


class RatesUpPdf:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
    ) -> None:

        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory

        self._rates_up_pdf = self.stencil_factory.from_dims_halo(
            func=rates_up_pdf,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )
        self._hcot = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")

        self._ktop_flag1 = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM], units="n/a")
        self._ktop_flag2 = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM], units="n/a")

    def __call__(
        self,
        # In
        OVERSHOOT: FloatField,
        entr_rate: FloatField,
        heo: FloatField,
        heso_cup: FloatField,
        hkbo: FloatField,
        kbcon: IntField,
        name: IntField,
        z_cup: FloatField,
        # Out
        ierr: IntField,
        ktopIJ: IntFieldIJ,
        ktop: IntField,
    ):

        self._rates_up_pdf(
            # In
            OVERSHOOT=OVERSHOOT,
            entr_rate=entr_rate,
            heo=heo,
            heso_cup=heso_cup,
            hkbo=hkbo,
            kbcon=kbcon,
            name=name,
            z_cup=z_cup,
            hcot=self._hcot,
            ktop_flag1=self._ktop_flag1,
            ktop_flag2=self._ktop_flag2,
            ktopIJ=ktopIJ,
            # Out
            ierr=ierr,
            ktop=ktop,
        )
