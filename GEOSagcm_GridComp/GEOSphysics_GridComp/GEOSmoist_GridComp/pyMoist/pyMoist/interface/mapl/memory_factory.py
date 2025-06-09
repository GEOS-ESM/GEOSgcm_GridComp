from MAPLpyish import MAPLBridge, MAPLState
from pyMoist.interface.f_py_conversion import FortranPythonConversion
from ndsl import QuantityFactory
import numpy.typing as npt


class MAPLMemoryFactory:
    """A factory capable of accessing memory handlded by MAPL and formatting it
    for direct use in DSL application."""

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

    def get_array(
        self,
        name: str,
        dtype: npt.DTypeLike,
        dims: list[str],
        alloc: bool = False,
    ) -> npt.ArrayLike:
        # MAPL Fortran call retrieve memory as a void*
        fptr = self._bridge.MAPL_GetPointer(self._state, name, alloc=alloc)

        # Turn dims into int-sized list
        shape = self._quantity_factory.sizer.get_shape(dims)

        return self._f_py_converter.fortran_to_python(fptr, dim=list(shape))

    def get_resource(
        self, name: str, dtype: npt.DTypeLike, default: bool = False
    ) -> npt.DTypeLike:
        return self._bridge.MAPL_GetResource(self._state, name, dtype, default=default)
