from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import (
    FloatField,
    IntField,
)

import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
import pyMoist.constants as constants
from gt4py.cartesian.gtscript import (
    FORWARD,
    PARALLEL,
    K,
    computation,
    interval,
)


def cup_up_aa0(
    # In
    dby: FloatField,
    gamma_cup: FloatField,
    ktop: IntField,
    kbcon: IntField,
    k22: IntField,
    z_cup: FloatField,
    t_cup: FloatField,
    zu: FloatField,
    integ: IntField,
    integ_interval: IntField,
    # Out
    ierr: IntField,
    aa0: FloatField,
):
    from __externals__ import k_start

    with computation(PARALLEL), interval(...):
        if integ == 1:
            if integ_interval == cumulus_parameterization_constants.BL:
                kbeg = k_start
                kend = kbcon - 2
            elif integ_interval == cumulus_parameterization_constants.CIN:
                kbeg = k22 - 1
                kend = kbcon - 2

        else:
            kbeg = kbcon - 1
            kend = ktop - 1

    with computation(FORWARD), interval(0, -1):
        if ierr == 0:
            if K >= kbeg and K <= kend:
                dz = z_cup[0, 0, 1] - z_cup
                aa_1 = (
                    zu
                    * (
                        constants.MAPL_GRAV
                        / (cumulus_parameterization_constants.CP * t_cup)
                    )
                    * dby
                    / (1.0 + gamma_cup)
                )
                aa_2 = (
                    zu[0, 0, 1]
                    * (
                        constants.MAPL_GRAV
                        / (cumulus_parameterization_constants.CP * t_cup[0, 0, 1])
                    )
                    * dby[0, 0, 1]
                    / (1.0 + gamma_cup[0, 0, 1])
                )
                da = 0.5 * (aa_1 + aa_2) * dz
                aa0_below = aa0.at(K=K - 1)
                aa0 = aa0_below + da

    with computation(PARALLEL), interval(...):
        if ierr == 0:
            aa0_temp = aa0.at(K=kend)
            aa0 = aa0_temp


class CupUpAa0:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
    ) -> None:

        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory

        self._cup_up_aa0 = self.stencil_factory.from_dims_halo(
            func=cup_up_aa0,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        # In
        dby: FloatField,
        gamma_cup: FloatField,
        ktop: IntField,
        kbcon: IntField,
        k22: IntField,
        z_cup: FloatField,
        t_cup: FloatField,
        zu: FloatField,
        integ: IntField,
        integ_interval: IntField,
        # Out
        ierr: IntField,
        aa0: FloatField,
    ):

        self._cup_up_aa0(
            # In
            dby=dby,
            gamma_cup=gamma_cup,
            ktop=ktop,
            kbcon=kbcon,
            k22=k22,
            z_cup=z_cup,
            t_cup=t_cup,
            zu=zu,
            integ=integ,
            integ_interval=integ_interval,
            # Out
            ierr=ierr,
            aa0=aa0,
        )
