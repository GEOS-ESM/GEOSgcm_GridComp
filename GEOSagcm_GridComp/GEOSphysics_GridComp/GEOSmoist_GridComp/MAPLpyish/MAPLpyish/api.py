import cffi
import os
from _cffi_backend import _CDataBase as CFFIObj
from typing import Any
import numpy as np


class MAPLBridge:
    def __init__(self) -> None:
        # We leverage an environment variable to know where to look for the bridge
        # library
        geos_dir = os.getenv("GEOSDIR", "Not found")
        if geos_dir == "Not found":
            raise RuntimeError(
                "[MAPLPyish] Libary loads require a GEOSDIR environment variable"
                "pointing to the install directory of GEOS."
            )

        # FFI & C library setup
        self.ffi = cffi.FFI()
        self.mapl_c_bridge = self.ffi.dlopen(f"{geos_dir}lib/libMAPLpyish.so")

        # We use CFFI ABI mode, so we need to describe each function cdef
        # to the system

        # ESMF_AttributeGet
        self.ffi.cdef(
            "int MAPLPy_ESMF_AttributeGet_1D_int(void* esmf_state_c_ptr, char* name_c_ptr, int name_len);"
        )

        # MAPLPy_ESMF_MethodExecute
        self.ffi.cdef(
            "void MAPLPy_ESMF_MethodExecute(void* esmf_state_c_ptr, char* label_c_ptr, int label_len);"
        )

        # MAPLpy_GetPointer_via_ESMFAttr
        self.ffi.cdef(
            "void* MAPLpy_GetPointer_via_ESMFAttr(void* esmf_state_c_ptr, char* name_c_ptr, int name_len);"
        )

        # MAPLpy_GetPointer
        self.ffi.cdef(
            "void* MAPLpy_GetPointer(void* esmf_state_c_ptr, char* name_c_ptr, int name_len, bool alloc);"
        )

        # MAPL_GetResource
        self.ffi.cdef(
            "void* MAPLpy_GetResource_Bool(void* esmf_state_c_ptr, char* name_c_ptr, int name_len, bool default);"
        )
        self.ffi.cdef(
            "void* MAPLpy_GetResource_Int(void* esmf_state_c_ptr, char* name_c_ptr, int name_len, bool default);"
        )
        self.ffi.cdef(
            "void* MAPLpy_GetResource_Float(void* esmf_state_c_ptr, char* name_c_ptr, int name_len, bool default);"
        )

        # ESMF_TimeIntervalGet
        self.ffi.cdef("void* MAPLpy_ESMF_TimeIntervalGet(void* esmf_time_state_c_ptr);")

    def __del__(self):
        self.ffi.dlclose(self.mapl_c_bridge)

    def ESMF_AttributeGet(self, state: CFFIObj, name: str) -> Any:
        # TODO: depending on value type, redirect to correct bridge function
        return self.mapl_c_bridge.MAPLPy_ESMF_AttributeGet_1D_int(  # type: ignore
            state,
            self.ffi.new("char[]", name.encode()),
            len(name),
        )

    def ESMF_MethodExecute(self, state: CFFIObj, label: str) -> Any:
        self.mapl_c_bridge.MAPLPy_ESMF_MethodExecute(  # type: ignore
            state,
            self.ffi.new("char[]", label.encode()),
            len(label),
        )

    def MAPL_GetPointer_via_ESMFAttr(
        self,
        state: CFFIObj,
        name: str,
    ) -> Any:
        # TODO: depending on value type, redirect to correct bridge function
        return self.mapl_c_bridge.MAPLpy_GetPointer_via_ESMFAttr(  # type: ignore
            state,
            self.ffi.new("char[]", name.encode()),
            len(name),
        )

    def MAPL_GetPointer(
        self,
        state: CFFIObj,
        name: str,
        dtype: type,
        default: bool = False,
    ) -> Any:
        if isinstance(dtype, int):
            return self.mapl_c_bridge.MAPLpy_GetResource_Int(  # type: ignore
                state, self.ffi.new("char[]", name.encode()), len(name), default
            )
        elif isinstance(dtype, float):
            return self.mapl_c_bridge.MAPLpy_GetResource_Float(  # type: ignore
                state, self.ffi.new("char[]", name.encode()), len(name), default
            )
        elif isinstance(dtype, bool):
            return self.mapl_c_bridge.MAPLpy_GetResource_Bool(  # type: ignore
                state, self.ffi.new("char[]", name.encode()), len(name), default
            )
        raise NotImplementedError(f"MAPL_GetResource for type {dtype} not implemented.")

    def ESMF_TimeIntervalGet(
        self,
        time_state: CFFIObj,
    ) -> np.float64:
        return self.mapl_c_bridge.MAPLpy_ESMF_TimeIntervalGet(time_state)
