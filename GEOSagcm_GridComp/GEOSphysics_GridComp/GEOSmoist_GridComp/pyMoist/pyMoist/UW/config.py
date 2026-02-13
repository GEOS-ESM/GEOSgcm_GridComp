from dataclasses import dataclass

from ndsl.dsl.typing import Float, Int


@dataclass
class UWConfiguration:
    NCNST: Int
    k0: Int
    windsrcavg: Int
    dotransport: Int
    qtsrchgt: Float
    qtsrc_fac: Float
    thlsrc_fac: Float
    frc_rasn: Float
    rbuoy: Float
    epsvarw: Float
    use_CINcin: Int
    mumin1: Float
    rmaxfrac: Float
    PGFc: Float
    dt: Float
    niter_xc: Int
    criqc: Float
    rle: Float
    cridist_opt: Int
    mixscale: Float
    rdrag: Float
    rkm: Float
    use_self_detrain: Int
    detrhgt: Float
    use_cumpenent: Int
    rpen: Float
    use_momenflx: Int
    rdrop: Float
    iter_cin: Int
    SCLM_SHALLOW: Float
