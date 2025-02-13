from _cffi_backend import _CDataBase as CFFIObj  # type: ignore
import dataclasses
from pyMKIAU.f_py_conversion import FortranPythonConversion
from pyMKIAU.cuda_profiler import TimedCUDAProfiler
import numpy as np
from typing import Dict, List


@dataclasses.dataclass
class FPYOptions:
    npx: int = 0
    npy: int = 0
    npz: int = 0
    mn_123456789: int = 0


def options_fortran_to_python(
    f_options: CFFIObj,
) -> FPYOptions:
    if f_options.mn_123456789 != 123456789:  # type:ignore
        raise RuntimeError(
            "Magic number failed, pyMoist interface is broken on the python side"
        )

    py_flags = FPYOptions()
    keys = list(filter(lambda k: not k.startswith("__"), dir(type(py_flags))))
    for k in keys:
        if hasattr(f_options, k):
            setattr(py_flags, k, getattr(f_options, k))
    return py_flags


F_PY_MEMORY_CONV = None


def pyMKIAU_init():
    print("[pyMKIAU] Init called")


def pyMKIAU_run(
    f_options: CFFIObj,
    f_in_buffer: CFFIObj,
    f_out_buffer: CFFIObj,
):
    print("[pyMKIAU] Run called")
    options = options_fortran_to_python(f_options)
    print(f"[pyMKIAU] Options: {options}")

    # Dev Note: this should be doen better in it's own class
    #           and the `np` should be driven by the user code requirements
    #           for GPU or CPU memory
    global F_PY_MEMORY_CONV
    if F_PY_MEMORY_CONV is None:
        F_PY_MEMORY_CONV = FortranPythonConversion(
            options.npx,
            options.npy,
            options.npz,
            np,
        )

    # Move memory into a manipulable numpy array
    in_buffer = F_PY_MEMORY_CONV.fortran_to_python(f_in_buffer)
    out_buffer = F_PY_MEMORY_CONV.fortran_to_python(f_out_buffer)

    # Here goes math and dragons
    timings: Dict[str, List[float]] = {}
    with TimedCUDAProfiler("pyMKIAU bogus math", timings):
        out_buffer[:, :, :] = in_buffer[:, :, :] * 2

    print(f"[pyMKIAU] At 5,5,5 in python OUT is: {out_buffer[5,5,5]}")
    print(f"[pyMKIAU] Timers: {timings}")

    # Go back to fortran
    F_PY_MEMORY_CONV.python_to_fortran(out_buffer, f_out_buffer)
