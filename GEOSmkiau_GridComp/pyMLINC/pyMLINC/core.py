import numpy
import typing
import dataclasses
from _cffi_backend import _CDataBase as CFFIObj  # type: ignore

from pyMLINC.f_py_conversion import FortranPythonConversion
from pyMLINC.cuda_profiler import TimedCUDAProfiler
from geos_state_bias.processor import Processor


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
        ps_f: CFFIObj,
        # output
        dtdt_f: CFFIObj,
        # LAST ARGUMENT - input
        magic_number: int
):
    check_magic_number(magic_number)
    global F_PY_MEMORY_CONV
    if F_PY_MEMORY_CONV is None:
        F_PY_MEMORY_CONV = FortranPythonConversion(xdim, ydim, zdim, numpy)

    # Move memory into a manipulatable numpy array
    u = F_PY_MEMORY_CONV.fortran_to_python(u_f).transpose()
    v = F_PY_MEMORY_CONV.fortran_to_python(v_f).transpose()
    t = F_PY_MEMORY_CONV.fortran_to_python(t_f).transpose()
    qv = F_PY_MEMORY_CONV.fortran_to_python(qv_f).transpose()
    ql = F_PY_MEMORY_CONV.fortran_to_python(ql_f).transpose()
    qi = F_PY_MEMORY_CONV.fortran_to_python(qi_f).transpose()
    qr = F_PY_MEMORY_CONV.fortran_to_python(qr_f).transpose()
    qs = F_PY_MEMORY_CONV.fortran_to_python(qs_f).transpose()
    qg = F_PY_MEMORY_CONV.fortran_to_python(qg_f).transpose()
    ps = F_PY_MEMORY_CONV.fortran_to_python(ps_f, dim=[xdim, ydim]).transpose()
    print("[pyMLINC] Python - u:", u.shape)
    print("[pyMLINC] Python - u:", numpy.sum(u), numpy.min(u), numpy.max(u))
    print("[pyMLINC] Python - v:", numpy.sum(v), numpy.min(v), numpy.max(v))
    print("[pyMLINC] Python - t:", numpy.sum(t), numpy.min(t), numpy.max(t))
    print("[pyMLINC] Python - qv:", numpy.sum(qv), numpy.min(qv), numpy.max(qv))
    print("[pyMLINC] Python - ql:", numpy.sum(ql), numpy.min(ql), numpy.max(ql))
    print("[pyMLINC] Python - qi:", numpy.sum(qi), numpy.min(qi), numpy.max(qi))
    print("[pyMLINC] Python - qr:", numpy.sum(qr), numpy.min(qr), numpy.max(qr))
    print("[pyMLINC] Python - qs:", numpy.sum(qs), numpy.min(qs), numpy.max(qs))
    print("[pyMLINC] Python - qg:", numpy.sum(qg), numpy.min(qg), numpy.max(qg))
    print("[pyMLINC] Python - ps:", numpy.sum(ps), numpy.min(ps), numpy.max(ps))
    print("[pyMLINC] Python - flushing buffer.", flush=True)
    

    # Order of vvariables as defined in processor::__init__
    # U, V, T, QV, QI, QL, QG, QR, QS, PS
    arrays = [u, v, t, qv, qi, ql, qg, qr, qs, ps]
    ckpt_root_path = "/discover/nobackup/jli30/geos_state/SmaAt-UNet/geos_state_bias/checkpoints/batch_2"

    # Here goes math and dragons
    timings: typing.Dict[str, typing.List[float]] = {}
    with TimedCUDAProfiler("pyMLINC bogus math", timings):
        processor = Processor(ckpt_root_path, *arrays)
        dtdt = processor.predict()

    print(f"[pyMLINC] run - dtdt:", numpy.sum(dtdt), numpy.min(dtdt), numpy.max(dtdt))
    print(f"[pyMLINC] run - timers: {timings}", flush=True)

    # Output goes back to fortran
    F_PY_MEMORY_CONV.python_to_fortran(dtdt.transpose(), dtdt_f)
