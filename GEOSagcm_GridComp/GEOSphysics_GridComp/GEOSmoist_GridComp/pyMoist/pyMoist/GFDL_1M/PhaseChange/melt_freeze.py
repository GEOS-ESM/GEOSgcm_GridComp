from gt4py.cartesian.gtscript import PARALLEL, computation, exp, interval

import pyMoist.constants as constants
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ
from pyMoist.shared_incloud_processes import ice_fraction


def melt_freeze(
    cnv_frc: FloatFieldIJ,
    srf_type: FloatFieldIJ,
    t: FloatField,
    qlcn: FloatField,
    qicn: FloatField,
):
    from __externals__ import DT_MOIST

    with computation(PARALLEL), interval(...):
        if t <= constants.MAPL_TICE:
            f_qi = ice_fraction(t, cnv_frc, srf_type)
            d_qil = qlcn * (1.0 - exp(-DT_MOIST * f_qi / constants.TAUFRZ))
            d_qil = max(0.0, d_qil)
            qicn = qicn + d_qil
            qlcn = qlcn - d_qil
            t = (
                t
                + (
                    constants.MAPL_LATENT_HEAT_SUBLIMATION
                    - constants.MAPL_LATENT_HEAT_VAPORIZATION
                )
                * d_qil
                / constants.MAPL_CP
            )
        else:
            d_qil = -qicn
            d_qil = min(0.0, d_qil)
            qicn = qicn + d_qil
            qlcn = qlcn - d_qil
            t = (
                t
                + (
                    constants.MAPL_LATENT_HEAT_SUBLIMATION
                    - constants.MAPL_LATENT_HEAT_VAPORIZATION
                )
                * d_qil
                / constants.MAPL_CP
            )
