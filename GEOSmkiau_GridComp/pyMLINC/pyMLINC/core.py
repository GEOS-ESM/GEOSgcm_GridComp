import numpy
import typing
import dataclasses
from _cffi_backend import _CDataBase as CFFIObj  # type: ignore

from pyMLINC.f_py_conversion import FortranPythonConversion
from pyMLINC.cuda_profiler import TimedCUDAProfiler


F_PY_MEMORY_CONV = None


def check_magic_number(magic_number: int):
    if magic_number != 123456789:  # type:ignore
        raise RuntimeError("Magic number failed on the Python side")


def pyMLINC_init(magic_number: int):
    check_magic_number(magic_number)
    print(f"[pyMLINC] init", flush=True)


def pyMLINC_run(
        # input
        xdim: int,
        ydim: int,
        zdim: int,
        u_f: CFFIObj,
        qv_f: CFFIObj,
        # output
        dtdt_f: CFFIObj,
        # LAST ARGUMENT - input
        magic_number: int
):
    check_magic_number(magic_number)
    global F_PY_MEMORY_CONV
    if F_PY_MEMORY_CONV is None:
        F_PY_MEMORY_CONV = FortranPythonConversion(xdim, ydim, zdim, numpy)

    # Move memory into a manipulable numpy array
    u = F_PY_MEMORY_CONV.fortran_to_python(u_f)
    qv = F_PY_MEMORY_CONV.fortran_to_python(qv_f)
    print("[pyMLINC] run - qv:", numpy.sum(qv), numpy.min(qv), numpy.max(qv))
    dtdt = F_PY_MEMORY_CONV.fortran_to_python(dtdt_f)

    # Here goes math and dragons
    timings: typing.Dict[str, typing.List[float]] = {}
    with TimedCUDAProfiler("pyMLINC bogus math", timings):
        dtdt[:, :, :] = qv[:, :, :] * 2

    print(f"[pyMLINC] run - dtdt:", numpy.sum(dtdt), numpy.min(dtdt), numpy.max(dtdt))
    print(f"[pyMLINC] run - timers: {timings}", flush=True)

    # Go back to fortran
    F_PY_MEMORY_CONV.python_to_fortran(dtdt, dtdt_f)
