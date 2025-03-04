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
        v_f: CFFIObj,
        t_f: CFFIObj,
        qv_f: CFFIObj,
        ql_f: CFFIObj,
        qi_f: CFFIObj,
        qr_f: CFFIObj,
        qs_f: CFFIObj,
        qg_f: CFFIObj,
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
    v = F_PY_MEMORY_CONV.fortran_to_python(v_f)
    t = F_PY_MEMORY_CONV.fortran_to_python(t_f)
    qv = F_PY_MEMORY_CONV.fortran_to_python(qv_f)
    ql = F_PY_MEMORY_CONV.fortran_to_python(ql_f)
    qi = F_PY_MEMORY_CONV.fortran_to_python(qi_f)
    qr = F_PY_MEMORY_CONV.fortran_to_python(qr_f)
    qs = F_PY_MEMORY_CONV.fortran_to_python(qs_f)
    qg = F_PY_MEMORY_CONV.fortran_to_python(qg_f)
    print("[pyMLINC] run - u:", numpy.sum(u), numpy.min(u), numpy.max(u))
    print("[pyMLINC] run - v:", numpy.sum(v), numpy.min(v), numpy.max(v))
    print("[pyMLINC] run - t:", numpy.sum(t), numpy.min(t), numpy.max(t))
    print("[pyMLINC] run - qv:", numpy.sum(qv), numpy.min(qv), numpy.max(qv))
    print("[pyMLINC] run - ql:", numpy.sum(ql), numpy.min(ql), numpy.max(ql))
    print("[pyMLINC] run - qi:", numpy.sum(qi), numpy.min(qi), numpy.max(qi))
    print("[pyMLINC] run - qr:", numpy.sum(qr), numpy.min(qr), numpy.max(qr))
    print("[pyMLINC] run - qs:", numpy.sum(qs), numpy.min(qs), numpy.max(qs))
    print("[pyMLINC] run - qg:", numpy.sum(qg), numpy.min(qg), numpy.max(qg))
    dtdt = F_PY_MEMORY_CONV.fortran_to_python(dtdt_f)

    # Here goes math and dragons
    timings: typing.Dict[str, typing.List[float]] = {}
    with TimedCUDAProfiler("pyMLINC bogus math", timings):
        dtdt[:, :, :] = qv[:, :, :] * 2

    print(f"[pyMLINC] run - dtdt:", numpy.sum(dtdt), numpy.min(dtdt), numpy.max(dtdt))
    print(f"[pyMLINC] run - timers: {timings}", flush=True)

    # Go back to fortran
    F_PY_MEMORY_CONV.python_to_fortran(dtdt, dtdt_f)
