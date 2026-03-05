from MAPL_PythonBridge.python2fortran import MAPLPyAPI
from ndsl.constants import Float
from typing import no_type_check
import cffi
import os
from MAPL_PythonBridge.types import FFI
import platform
import dataclasses
import numpy.typing as npt
from MAPL_PythonBridge import get_MAPLPy


@dataclasses.dataclass
class CNVTracers:
    Q: npt.NDArray
    fscav: npt.NDArray
    Vect_Hcts: npt.NDArray
    use_gcc_washout: npt.NDArray


def _fortran_to_numpy(
    maplpy: MAPLPyAPI,
    c_void_ptr,
    dtype: npt.DTypeLike,
    size: list[int],
) -> npt.NDArray:
    data_as_typed_ptr = maplpy.fpy_converter.cast(dtype, c_void_ptr)
    return maplpy.fpy_converter.fortran_to_python(data_as_typed_ptr, size)


class MoistWorkarounds:
    """"""

    def __init__(self, ffi: cffi.FFI) -> None:
        # We leverage an environment variable to know where to look for the bridge
        # library
        geos_dir = os.getenv("GEOSDIR", "Not found")
        if geos_dir == "Not found":
            raise RuntimeError(
                "[pyMoist.fortran.moist_workarounds] Libary loads require a GEOSDIR environment variable"
                "pointing to the install directory of GEOS."
            )
        self.ffi = ffi

        # FFI & C library setup
        # TODO: we should be using the out-of-line API system and link to .a
        # that way we don't have to rely on a shady .so load
        if platform.system() == "Linux":
            ext = "so"
        elif platform.system() == "Darwin":
            ext = "dylib"
        else:
            raise SystemError(f"MAPLPyish unavailable on {platform.system()}")

        self.libGEOSmoist_GridComp = self.ffi.dlopen(f"{geos_dir}/lib/libGEOSmoist_GridComp.{ext}")

        # We use CFFI ABI mode, so we need to describe each function cdef
        # to the system
        self.ffi.cdef("int   get_CNV_Tracers_SOA__size();")
        self.ffi.cdef("void* get_CNV_Tracers_SOA__Q();")
        self.ffi.cdef("void* get_CNV_Tracers_SOA__fscav();")
        self.ffi.cdef("void* get_CNV_Tracers_SOA__Vect_Hcts();")
        self.ffi.cdef("void* get_CNV_Tracers_SOA__use_gcc_washout();")

    def __del__(self):
        self.ffi.dlclose(self.libGEOSmoist_GridComp)

    @no_type_check
    def CNV_Tracers(self) -> CNVTracers:

        maplpy = get_MAPLPy()
        size_ = self.libGEOSmoist_GridComp.get_CNV_Tracers_SOA__size()

        return CNVTracers(
            Q=_fortran_to_numpy(
                maplpy, self.libGEOSmoist_GridComp.get_CNV_Tracers_SOA__Q(), Float, maplpy.grid_dims + [size_]
            ),
            fscav=_fortran_to_numpy(
                maplpy, self.libGEOSmoist_GridComp.get_CNV_Tracers_SOA__fscav(), Float, [size_]
            ),
            Vect_Hcts=_fortran_to_numpy(
                maplpy, self.libGEOSmoist_GridComp.get_CNV_Tracers_SOA__Vect_Hcts(), Float, [4, size_]
            ),
            use_gcc_washout=_fortran_to_numpy(
                maplpy, self.libGEOSmoist_GridComp.get_CNV_Tracers_SOA__use_gcc_washout(), bool, [size_]
            ),
        )


MOIST_WORKAROUNDS = MoistWorkarounds(FFI)
