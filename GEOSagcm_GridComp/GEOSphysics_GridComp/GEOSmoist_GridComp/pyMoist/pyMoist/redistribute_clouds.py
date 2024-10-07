from gt4py.cartesian.gtscript import (
    PARALLEL,
    computation,
    interval,
)
from ndsl import StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField

import pyMoist.radiation_coupling_constants as radconstants

# ruff: noqa: PLR0913


def redist_clouds(
    cf: FloatField,
    ql: FloatField,
    qi: FloatField,
    clcn: FloatField,
    clls: FloatField,
    qlcn: FloatField,
    qlls: FloatField,
    qinc: FloatField,
    qils: FloatField,
    qv: FloatField,
    te: FloatField,
):
    with computation(PARALLEL), interval(...):
        # Constants from MAPL.h
        alhlbcp = radconstants.ALHLBCP
        alhsbcp = radconstants.ALHSBCP

        # Define fcn as a 3-d array
        fcn = cf

        # Fix cloud quants if too small
        if ql + qi < 1e-8:
            qv = qv + ql + qi
            te = te - alhlbcp * ql - alhsbcp * qi
            cf = 0.0
            ql = 0.0
            qi = 0.0

        if cf < 1e-5:
            qv = qv + ql + qi
            te = te - (alhlbcp * ql) - (alhsbcp * qi)
            cf = 0.0
            ql = 0.0
            qi = 0.0

        # Redistribute liquid CN/LS portions based on prior fractions
        fcn = 0.0
        if qlcn + qlls > 0.0:
            fcn = min(max(qlcn / (qlcn + qlls), 0.0), 1.0)

        # Put all new condensate into LS
        dqc = ql - (qlcn + qlls)
        if dqc > 0.0:
            qlls = qlls + dqc
            dqc = 0.0

        # Any loss of condensate uses the FCN ratio
        qlcn = qlcn + dqc * fcn
        qlls = qlls + dqc * (1.0 - fcn)

        # Redistribute ice CN/LS portions based on prior fractions
        fcn = 0.0
        if qinc + qils > 0.0:
            fcn = min(max(qinc / (qinc + qils), 0.0), 1.0)

        # Put all new condensate into LS
        dqc = qi - (qinc + qils)
        if dqc > 0.0:
            qils = qils + dqc
            dqc = 0.0

        # Any loss of condensate uses the FCN ratio
        qinc = qinc + dqc * fcn
        qils = qils + dqc * (1.0 - fcn)

        # Redistribute cloud-fraction CN/LS portions based on prior fractions
        fcn = 0.0
        if clcn + clls > 0.0:
            fcn = min(max(clcn / (clcn + clls), 0.0), 1.0)

        # Put all new condensate into LS
        dqc = cf - (clcn + clls)
        if dqc > 0.0:
            clls = clls + dqc
            dqc = 0.0

        # Any loss of condensate uses the FCN ratio
        clcn = clcn + dqc * fcn
        clls = clls + dqc * (1.0 - fcn)


class RedistributeClouds:
    def __init__(
        self,
        stencil_factory: StencilFactory,
    ) -> None:
        self._redist_clouds = stencil_factory.from_dims_halo(
            func=redist_clouds,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        rad_cf: FloatField,
        rad_ql: FloatField,
        rad_qi: FloatField,
        clcn: FloatField,
        clls: FloatField,
        qlcn: FloatField,
        qlls: FloatField,
        qicn: FloatField,
        q_ils: FloatField,
        rad_qv: FloatField,
        t: FloatField,
    ):
        """
        Perform the redistribute clouds process.

        Parameters:
        rad_cf (3D in/out): Radiation cloud fraction.
        rad_ql (3D in/out): Radiation liquid cloud mixing ratio.
        rad_qi (3D in/out): Radiation ice cloud mixing ratio.
        clcn (3D in/out): Cloud fraction (anvil).
        clls (3D in/out): Cloud fraction (large-scale).
        qlcn (3D in/out): Liquid cloud mixing ratio (anvil).
        qlls (3D in/out): Liquid cloud mixing ratio (large-scale).
        qicn (3D in/out): Ice cloud mixing ratio (anvil).
        q_ils (3D in/out): Ice cloud mixing ratio (large-scale).
        rad_qv (3D in/out): Radiation water vapor mixing ratio.
        t (3D in/out): Temperature.
        """
        self._redist_clouds(
            cf=rad_cf,
            ql=rad_ql,
            qi=rad_qi,
            clcn=clcn,
            clls=clls,
            qlcn=qlcn,
            qlls=qlls,
            qicn=qicn,
            qils=q_ils,
            qv=rad_qv,
            te=t,
        )
