from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import cffi
import numpy as np


@dataclass
class MoistFlags:
    # Grid layout
    npx: int = 0
    npy: int = 0
    npz: int = 0
    layout_x: int = 1
    layout_y: int = 1
    n_tiles: int = 6
    # Aer Activation
    n_modes: int = 0
    # Magic number
    mn_123456789: int = 0


@dataclass
class GFDL1MFlags:
    # GFDL_1M driver configuration. Initial values are not true defaults.
    phys_hydrostatic: bool = False
    hydrostatic: bool = False
    dt_moist: np.float32 = 0
    mp_time: np.float32 = 0
    t_min: np.float32 = 0
    t_sub: np.float32 = 0
    tau_r2g: np.float32 = 0
    tau_smlt: np.float32 = 0
    tau_g2r: np.float32 = 0
    dw_land: np.float32 = 0
    dw_ocean: np.float32 = 0
    vi_fac: np.float32 = 0
    vr_fac: np.float32 = 0
    vs_fac: np.float32 = 0
    vg_fac: np.float32 = 0
    ql_mlt: np.float32 = 0
    do_qa: bool = False
    fix_negative: bool = False
    vi_max: np.float32 = 0
    vs_max: np.float32 = 0
    vg_max: np.float32 = 0
    vr_max: np.float32 = 0
    qs_mlt: np.float32 = 0
    qs0_crt: np.float32 = 0
    qi_gen: np.float32 = 0
    ql0_max: np.float32 = 0
    qi0_max: np.float32 = 0
    qi0_crt: np.float32 = 0
    qr0_crt: np.float32 = 0
    fast_sat_adj: bool = False
    rh_inc: np.float32 = 0
    rh_ins: np.float32 = 0
    rh_inr: np.float32 = 0
    const_vi: bool = False
    const_vs: bool = False
    const_vg: bool = False
    const_vr: bool = False
    use_ccn: bool = False
    rthreshu: np.float32 = 0
    rthreshs: np.float32 = 0
    ccn_l: np.float32 = 0
    ccn_o: np.float32 = 0
    qc_crt: np.float32 = 0
    tau_g2v: np.float32 = 0
    tau_v2g: np.float32 = 0
    tau_s2v: np.float32 = 0
    tau_v2s: np.float32 = 0
    tau_revp: np.float32 = 0
    tau_frz: np.float32 = 0
    do_bigg: bool = False
    do_evap: bool = False
    do_subl: bool = False
    sat_adj0: np.float32 = 0
    c_piacr: np.float32 = 0
    tau_imlt: np.float32 = 0
    tau_v2l: np.float32 = 0
    tau_l2v: np.float32 = 0
    tau_i2v: np.float32 = 0
    tau_i2s: np.float32 = 0
    tau_l2r: np.float32 = 0
    qi_lim: np.float32 = 0
    ql_gen: np.float32 = 0
    c_paut: np.float32 = 0
    c_psaci: np.float32 = 0
    c_pgacs: np.float32 = 0
    c_pgaci: np.float32 = 0
    z_slope_liq: bool = False
    z_slope_ice: bool = False
    prog_ccn: bool = False
    c_cracw: np.float32 = 0
    alin: np.float32 = 0
    clin: np.float32 = 0
    preciprad: bool = False
    cld_min: np.float32 = 0
    use_ppm: bool = False
    mono_prof: bool = False
    do_sedi_heat: bool = False
    sedi_transport: bool = False
    do_sedi_w: bool = False
    de_ice: bool = False
    icloud_f: np.int32 = 0
    irain_f: np.int32 = 0
    mp_print: bool = False
    # Magic number
    mn_123456789: int = 0


def _generic_config_bridge(
    py_flags: Any,
    f_flags: cffi.FFI.CData,
):
    keys = list(filter(lambda k: not k.startswith("__"), dir(type(py_flags))))
    for k in keys:
        if hasattr(f_flags, k):
            setattr(py_flags, k, getattr(f_flags, k))


def moist_flags_f_to_python(
    f_flags: cffi.FFI.CData,
) -> MoistFlags:
    if f_flags.mn_123456789 != 123456789:
        raise RuntimeError(
            "Magic number failed, pyMoist interface is broken on the python side"
        )

    py_flags = MoistFlags()
    _generic_config_bridge(py_flags, f_flags)
    return py_flags


def gfdl_1m_flags_f_to_python(
    f_flags: cffi.FFI.CData,
) -> GFDL1MFlags:
    if f_flags.mn_123456789 != 123456789:
        raise RuntimeError(
            "Magic number failed, pyMoist interface is broken on the python side"
        )

    py_flags = GFDL1MFlags()
    _generic_config_bridge(py_flags, f_flags)
    return py_flags
