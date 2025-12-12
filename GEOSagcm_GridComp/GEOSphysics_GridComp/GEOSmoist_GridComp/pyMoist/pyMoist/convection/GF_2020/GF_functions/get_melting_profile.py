from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import (
    FloatField,
    FloatFieldIJ,
    IntField,
)
import pyMoist.constants as constants
from gt4py.cartesian.gtscript import (
    FORWARD,
    PARALLEL,
    K,
    computation,
    interval,
)


def get_melting_profile(
    # In
    MELT_GLAC: IntField,
    cumulus: IntField,
    edto: FloatField,
    ierr: IntField,
    melting_layer: FloatField,
    p_liq_ice: FloatField,
    po_cup: FloatField,
    pwdo: FloatField,
    pwo: FloatField,
    qrco: FloatField,
    tn_cup: FloatField,
    total_pwo_solid_phase: FloatFieldIJ,
    # Out
    melting: FloatField,
):
    from __externals__ import k_end

    with computation(FORWARD), interval(...):
        ktf = k_end - 1
        if MELT_GLAC == True and cumulus == 1:
            pwo_solid_phase = 0.0
            pwo_eff = 0.0
            melting = 0.0

            if ierr > 0:
                melting = 0.0

            total_pwo_solid_phase = 0.0

    with computation(FORWARD), interval(...):
        if MELT_GLAC == True and cumulus == 1:
            if K <= ktf - 1:
                if ierr == 0:
                    dp = 100.0 * (po_cup - po_cup[0, 0, 1])

                    pwo_eff = 0.5 * (pwo + pwo[0, 0, 1])

                    pwo_solid_phase = (1.0 - p_liq_ice) * pwo_eff

                    total_pwo_solid_phase = total_pwo_solid_phase + pwo_solid_phase * dp / constants.MAPL_GRAV

    with computation(PARALLEL), interval(...):
        if MELT_GLAC == True and cumulus == 1:
            if K <= ktf:
                if ierr == 0:
                    melting = melting_layer * (
                        total_pwo_solid_phase
                        / (100 * (po_cup.at(K=0) - po_cup.at(K=ktf)) / constants.MAPL_GRAV)
                    )
        else:
            melting = 0.0


class GetMeltingProfile:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
    ) -> None:

        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory

        self._get_melting_profile = self.stencil_factory.from_dims_halo(
            func=get_melting_profile,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self.total_pwo_solid_phase = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM], units="n/a"
        )

    def __call__(
        self,
        # In
        MELT_GLAC: IntField,
        cumulus: IntField,
        edto: FloatField,
        ierr: IntField,
        melting_layer: FloatField,
        p_liq_ice: FloatField,
        po_cup: FloatField,
        pwdo: FloatField,
        pwo: FloatField,
        qrco: FloatField,
        tn_cup: FloatField,
        # Out
        melting: FloatField,
    ):

        self._get_melting_profile(
            # In
            MELT_GLAC=MELT_GLAC,
            cumulus=cumulus,
            edto=edto,
            ierr=ierr,
            melting_layer=melting_layer,
            p_liq_ice=p_liq_ice,
            po_cup=po_cup,
            pwdo=pwdo,
            pwo=pwo,
            qrco=qrco,
            tn_cup=tn_cup,
            total_pwo_solid_phase=self.total_pwo_solid_phase,
            # Out
            melting=melting,
        )
