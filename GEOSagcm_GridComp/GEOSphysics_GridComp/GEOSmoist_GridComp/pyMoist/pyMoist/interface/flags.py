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
    DT_MOIST: np.float32 = 0
    MP_TIME: np.float32 = 0
    T_MIN: np.float32 = 0
    T_SUB: np.float32 = 0
    TAU_R2G: np.float32 = 0
    TAU_SMLT: np.float32 = 0
    TAU_G2R: np.float32 = 0
    DW_LAND: np.float32 = 0
    DW_OCEAN: np.float32 = 0
    VI_FAC: np.float32 = 0
    VR_FAC: np.float32 = 0
    VS_FAC: np.float32 = 0
    VG_FAC: np.float32 = 0
    QL_MLT: np.float32 = 0
    DO_QA: bool = False
    FIX_NEGATIVE: bool = False
    VI_MAX: np.float32 = 0
    VS_MAX: np.float32 = 0
    VG_MAX: np.float32 = 0
    VR_MAX: np.float32 = 0
    QS_MLT: np.float32 = 0
    QS0_CRT: np.float32 = 0
    QI_GEN: np.float32 = 0
    QL0_MAX: np.float32 = 0
    QI0_MAX: np.float32 = 0
    QI0_CRT: np.float32 = 0
    QR0_CRT: np.float32 = 0
    FAST_SAT_ADJ: bool = False
    RH_INC: np.float32 = 0
    RH_INS: np.float32 = 0
    RH_INR: np.float32 = 0
    CONST_VI: bool = False
    CONST_VS: bool = False
    CONST_VG: bool = False
    CONST_VR: bool = False
    USE_CCN: bool = False
    RTHRESHU: np.float32 = 0
    RTHRESHS: np.float32 = 0
    CCN_L: np.float32 = 0
    CCN_O: np.float32 = 0
    QC_CRT: np.float32 = 0
    TAU_G2V: np.float32 = 0
    TAU_V2G: np.float32 = 0
    TAU_S2V: np.float32 = 0
    TAU_V2S: np.float32 = 0
    TAU_REVP: np.float32 = 0
    TAU_FRZ: np.float32 = 0
    DO_BIGG: bool = False
    DO_EVAP: bool = False
    DO_SUBL: bool = False
    SAT_ADJ0: np.float32 = 0
    C_PIACR: np.float32 = 0
    TAU_IMLT: np.float32 = 0
    TAU_V2L: np.float32 = 0
    TAU_L2V: np.float32 = 0
    TAU_I2V: np.float32 = 0
    TAU_I2S: np.float32 = 0
    TAU_L2R: np.float32 = 0
    QI_LIM: np.float32 = 0
    QL_GEN: np.float32 = 0
    C_PAUT: np.float32 = 0
    C_PSACI: np.float32 = 0
    C_PGACS: np.float32 = 0
    C_PGACI: np.float32 = 0
    Z_SLOPE_LIQ: bool = False
    Z_SLOPE_ICE: bool = False
    PROG_CCN: bool = False
    C_CRACW: np.float32 = 0
    ALIN: np.float32 = 0
    CLIN: np.float32 = 0
    PRECIPRAD: bool = False
    CLD_MIN: np.float32 = 0
    USE_PPM: bool = False
    MONO_PROF: bool = False
    DO_SEDI_HEAT: bool = False
    SEDI_TRANSPORT: bool = False
    DO_SEDI_W: bool = False
    DE_ICE: bool = False
    ICLOUD_F: np.int32 = 0
    IRAIN_F: np.int32 = 0
    MP_PRINT: bool = False
    USE_BERGERON: bool = False
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
        raise RuntimeError("Magic number failed, pyMoist interface is broken on the python side")

    py_flags = MoistFlags()
    _generic_config_bridge(py_flags, f_flags)
    return py_flags


def gfdl_1m_flags_f_to_python(
    f_flags: cffi.FFI.CData,
) -> GFDL1MFlags:
    if f_flags.mn_123456789 != 123456789:
        raise RuntimeError("Magic number failed, pyMoist interface is broken on the python side")

    py_flags = GFDL1MFlags()
    _generic_config_bridge(py_flags, f_flags)
    return py_flags
