import gt4py
import numpy as np
from gt4py.cartesian.gtscript import PARALLEL, computation, interval, stencil

import pyMoist.radiation_coupling_constants as radconstants
from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField


def redist_clouds(
    CF: FloatField,
    QL: FloatField,
    QI: FloatField,
    CLCN: FloatField,
    CLLS: FloatField,
    QLCN: FloatField,
    QLLS: FloatField,
    QICN: FloatField,
    QILS: FloatField,
    QV: FloatField,
    TE: FloatField,
):

    with computation(PARALLEL), interval(...):
        # Constants from MAPL.h
        alhlbcp = radconstants.ALHLBCP
        alhsbcp = radconstants.ALHSBCP

        # Define FCN as a 3-d array
        FCN = CF

        # Fix cloud quants if too small
        if QL + QI < 1e-8:
            QV = QV + QL + QI
            TE = TE - alhlbcp * QL - alhsbcp * QI
            CF = 0.0
            QL = 0.0
            QI = 0.0

        if CF < 1e-5:
            QV = QV + QL + QI
            TE = TE - (alhlbcp * QL) - (alhsbcp * QI)
            CF = 0.0
            QL = 0.0
            QI = 0.0

        # Redistribute liquid CN/LS portions based on prior fractions
        FCN = 0.0
        if QLCN + QLLS > 0.0:
            FCN = min(max(QLCN / (QLCN + QLLS), 0.0), 1.0)

        # Put all new condensate into LS
        DQC = QL - (QLCN + QLLS)
        if DQC > 0.0:
            QLLS = QLLS + DQC
            DQC = 0.0

        # Any loss of condensate uses the FCN ratio
        QLCN = QLCN + DQC * FCN
        QLLS = QLLS + DQC * (1.0 - FCN)

        # Redistribute ice CN/LS portions based on prior fractions
        FCN = 0.0
        if QICN + QILS > 0.0:
            FCN = min(max(QICN / (QICN + QILS), 0.0), 1.0)

        # Put all new condensate into LS
        DQC = QI - (QICN + QILS)
        if DQC > 0.0:
            QILS = QILS + DQC
            DQC = 0.0

        # Any loss of condensate uses the FCN ratio
        QICN = QICN + DQC * FCN
        QILS = QILS + DQC * (1.0 - FCN)

        # Redistribute cloud-fraction CN/LS portions based on prior fractions
        FCN = 0.0
        if CLCN + CLLS > 0.0:
            FCN = min(max(CLCN / (CLCN + CLLS), 0.0), 1.0)

        # Put all new condensate into LS
        DQC = CF - (CLCN + CLLS)
        if DQC > 0.0:
            CLLS = CLLS + DQC
            DQC = 0.0

        # Any loss of condensate uses the FCN ratio
        CLCN = CLCN + DQC * FCN
        CLLS = CLLS + DQC * (1.0 - FCN)


class RedistributeClouds:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
    ) -> None:

        self._redist_clouds = stencil_factory.from_dims_halo(
            func=redist_clouds,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        RAD_CF: FloatField,
        RAD_QL: FloatField,
        RAD_QI: FloatField,
        CLCN: FloatField,
        CLLS: FloatField,
        QLCN: FloatField,
        QLLS: FloatField,
        QICN: FloatField,
        QILS: FloatField,
        RAD_QV: FloatField,
        T: FloatField,
    ):
        """
        Perform the redistribute clouds process.

        Parameters:
        RAD_CF (3D inout): Radiation cloud fraction.
        RAD_QL (3D inout): Radiation liquid cloud mixing ratio.
        RAD_QI (3D inout): Radiation ice cloud mixing ratio.
        CLCN (3D inout): Cloud fraction (anvil).
        CLLS (3D inout): Cloud fraction (large-scale).
        QLCN (3D inout): Liquid cloud mixing ratio (anvil).
        QLLS (3D inout): Liquid cloud mixing ratio (large-scale).
        QICN (3D inout): Ice cloud mixing ratio (anvil).
        QILS (3D inout): Ice cloud mixing ratio (large-scale).
        RAD_QV (3D inout): Radiation water vapor mixing ratio.
        T (3D inout): Temperature.
        """
        self._redist_clouds(
            CF=RAD_CF,
            QL=RAD_QL,
            QI=RAD_QI,
            CLCN=CLCN,
            CLLS=CLLS,
            QLCN=QLCN,
            QLLS=QLLS,
            QICN=QICN,
            QILS=QILS,
            QV=RAD_QV,
            TE=T,
        )
