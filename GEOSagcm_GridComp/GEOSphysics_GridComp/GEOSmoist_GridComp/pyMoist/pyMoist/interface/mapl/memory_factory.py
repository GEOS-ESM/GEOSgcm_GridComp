from __future__ import annotations
from MAPLpyish import MAPLBridge, MAPLState
from pyMoist.interface.f_py_conversion import FortranPythonConversion
from ndsl import QuantityFactory
import numpy.typing as npt
from _cffi_backend import _CDataBase as CFFIObj
import dataclasses
import numpy as np


class MAPLMemoryRepository:
    """A factory capable of accessing memory handlded by MAPL and formatting it
    for direct use in DSL application."""

    @dataclasses.dataclass
    class FortranMemory:
        pointer: CFFIObj
        """C-binded Fortran layout memory"""
        associated: bool
        """Is the pointer properly associated in Fortran"""
        shape: tuple[int, ...]
        """Shape of the array"""
        python_array: npt.NDArray
        """Synced python memory, check dirty"""
        dirty: bool = False
        """Does the memory need to be synced"""

    def __init__(self, state: MAPLState, quantity_factory: QuantityFactory) -> None:
        self._quantity_factory = quantity_factory
        self._state = state
        self._bridge = MAPLBridge()
        self._f_py_converter = FortranPythonConversion(
            self._quantity_factory.sizer.nx,
            self._quantity_factory.sizer.ny,
            self._quantity_factory.sizer.nz,
            self._quantity_factory._numpy,
        )
        self._fortran_pointers: dict[str, MAPLMemoryRepository.FortranMemory] = {}

    def register(
        self,
        name: str,
        dtype: npt.DTypeLike,
        dims: list[str],
        alloc: bool = False,
    ):
        """Register the fortran memory with the factorty"""
        # MAPL Fortran call retrieve memory as a void* - we will cast
        if len(dims) == 3:
            void_fptr = self._bridge.MAPL_GetPointer_3D(self._state, name, alloc=alloc)
            is_associated = self._bridge.associated_3d(self._state, name, alloc=alloc)
        elif len(dims) == 2:
            void_fptr = self._bridge.MAPL_GetPointer_2D(self._state, name, alloc=alloc)
            is_associated = self._bridge.associated_2d(self._state, name, alloc=alloc)
        else:
            raise NotImplementedError(
                f"Only 2D & 3D fields implemented, missing support for {len(dims)}D arrays."
            )
        self._fortran_pointers[name] = MAPLMemoryRepository.FortranMemory(
            pointer=self._f_py_converter.cast(dtype, void_fptr),
            associated=is_associated,
            shape=self._quantity_factory.sizer.get_extent(dims),
            python_array=np.empty((0)),
        )

    def get_from_fortran(
        self,
        name: str,
    ) -> npt.NDArray:
        """Retrieve the data from Fortran. Prefer using a MAPLManager."""
        try:
            fmem = self._fortran_pointers[name]
        except KeyError:
            raise KeyError(f"Pointer {name} was never registered.")
        if not fmem.associated:
            return

        fmem.python_array = self._f_py_converter.fortran_to_python(fmem.pointer, dim=list(fmem.shape))

        return fmem.python_array

    def send_to_fortran(
        self,
        name: str,
    ) -> None:
        """Move the data back to Fortran. Prefer using a MAPLManager."""
        try:
            fmem = self._fortran_pointers[name]

        except KeyError:
            raise KeyError(f"Pointer {name} was never registered.")
        if not fmem.associated:
            return

        self._f_py_converter.python_to_fortran(fmem.python_array, fmem.pointer)

    def associated(self, name: str) -> bool:
        try:
            fmem = self._fortran_pointers[name]
        except KeyError:
            raise KeyError(f"Pointer {name} was never registered.")
        return fmem.associated

    def get_resource(self, name: str, dtype: npt.DTypeLike, default) -> npt.DTypeLike:
        return self._bridge.MAPL_GetResource(self._state, name, dtype(default))


class MAPLManagedMemory:
    """Context manager capable of get/send the memory from/to Fortran.

    Usage:
        with MAPLManagedMemory(mapl_memory_repository) as mmm:
            mmm.Field[...] # access memory "Field" as a NDArray
            mmm.associated('Field') # see if Fortran associatied the memory

    """

    def __init__(self, mapl_factory: MAPLMemoryRepository) -> None:
        self._mapl_factory = mapl_factory
        self._local_memory: dict[str, npt.NDArray] = {}

    def __enter__(self) -> MAPLManagedMemory:
        return self

    def __exit__(self, _exc_type, _exc_value, _traceback) -> None:
        for array_name in self._local_memory:
            self._mapl_factory.send_to_fortran(array_name)

    def associated(self, name: str) -> bool:
        return self._mapl_factory.associated(name)

    def __getattr__(self, name) -> npt.NDArray:
        if name in self._local_memory.keys():
            return self._local_memory[name]
        array = self._mapl_factory.get_from_fortran(name)
        self._local_memory[name] = array
        return array
