import cffi
from mpi4py import MPI
from distutils.sysconfig import get_config_var


TMPFILEBASE = "pyMoist_interface_py"

ffi = cffi.FFI()

# MPI_Comm can be int or void*
if MPI._sizeof(MPI.Comm) == ffi.sizeof("int"):
    _mpi_comm_t = "int"
else:
    _mpi_comm_t = "void*"

source = """
from {} import ffi
from datetime import datetime
from mpi4py import MPI
from pyMoist.interface.python_bridge import (
    pyMoist_init,
    pyMoist_run_AerActivation,
    pyMoist_run_GFDL1M,
    pyMoist_finalize
)
import traceback

@ffi.def_extern()
def pymoist_interface_py_init(flags) -> int:

    try:
        pyMoist_init(flags)
    except Exception as err:
        print("Error in Python:")
        print(traceback.format_exc())
        return -1
    return 0

@ffi.def_extern()
def pymoist_interface_py_run_AerActivation(
    aero_dgn, aero_num, aero_hygroscopicity, aero_sigma,
    frland, nn_ocean, nn_land,
    t, plo,
    qicn, qils, qlcn, qlls,
    vvel, tke,
    nacti, nwfa, nactl) -> int:

    try:
        pyMoist_run_AerActivation(
            aero_dgn, aero_num, aero_hygroscopicity, aero_sigma,
            frland, nn_ocean, nn_land,
            t, plo,
            qicn, qils, qlcn, qlls,
            vvel, tke,
            nacti, nwfa, nactl)
    except Exception as err:
        print("Error in Python:")
        print(traceback.format_exc())
        return -1
    return 0

@ffi.def_extern()
def pymoist_interface_py_run_GFDL1M(
    dw_land, dw_ocean, PDFSHAPE, TURNRHCRIT_PARAM,
    DT_MOIST, CCW_EVAP_EFF, CCI_EVAP_EFF,
    LMELTFRZ,
    AREA, CNV_FRC, SRF_TYPE,
    KLCL,
    EIS, PLmb, PLEmb, NACTL, NACTI, QST,
    T, Q, QLCN, QICN, QLLS, QILS, CLLS, CLCN,
    SUBLC, EVAPC, RHX):

    try:
        pyMoist_run_GFDL1M(
            dw_land, dw_ocean, PDFSHAPE, TURNRHCRIT_PARAM,
            DT_MOIST, CCW_EVAP_EFF, CCI_EVAP_EFF,
            LMELTFRZ,
            AREA, CNV_FRC, SRF_TYPE,
            KLCL,
            EIS, PLmb, PLEmb, NACTL, NACTI, QST,
            T, Q, QLCN, QICN, QLLS, QILS, CLLS, CLCN,
            SUBLC, EVAPC, RHX
        )
    except Exception as err:
        print("Error in Python:")
        print(traceback.format_exc())
        return -1
    return 0

@ffi.def_extern()
def pymoist_interface_py_finalize() -> int:
    try:
        pyMoist_finalize()
    except Exception as err:
        print("Error in Python:")
        print(traceback.format_exc())
        return -1
    return 0

""".format(
    TMPFILEBASE
)

with open("moist.h") as f:
    data = "".join([line for line in f if not line.startswith("#")])
    data = data.replace("CFFI_DLLEXPORT", "")
    ffi.embedding_api(data)

ffi.set_source(
    TMPFILEBASE,
    '#include "moist.h"',
    library_dirs=[get_config_var("LIBDIR")],
)
ffi.embedding_init_code(source)
ffi.compile(target="lib" + TMPFILEBASE + ".*", verbose=3)
