from typing import Tuple

import numpy.typing as npt

from ndsl import State
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Float
from ndsl.optional_imports import cupy as cp
from pyMoist.fortran.build_helper import InterfaceTransferType
from pyMoist.fortran.memory_factory import MAPLMemoryRepository


class MAPLManagedState:
    """Manage a NDSL <> MAPL shared state by linking MAPL pointers to NDSL state fields"""

    def __init__(
        self,
        py_state: State,
        transfer_type: InterfaceTransferType,
    ) -> None:
        self._ndsl_state = py_state
        self._state_to_mapl_mapping: dict[str, Tuple[MAPLMemoryRepository, str]] = {}
        self._transfer_type = transfer_type

    def register(
        self,
        ndsl_field_name: str,
        mapl_field_name: str,
        mapl_state: MAPLMemoryRepository,
        dtype: npt.DTypeLike = Float,
        dims: list[str] = [X_DIM, Y_DIM, Z_DIM],
        alloc: bool = False,
    ):
        mapl_state.register(mapl_field_name, dtype=dtype, dims=dims, alloc=alloc)
        self._state_to_mapl_mapping[ndsl_field_name] = (mapl_state, mapl_field_name)

    def register_2D(
        self,
        ndsl_field_name: str,
        mapl_field_name: str,
        mapl_state: MAPLMemoryRepository,
        dtype: npt.DTypeLike = Float,
        alloc: bool = False,
    ):
        mapl_state.register(mapl_field_name, dtype=dtype, dims=[X_DIM, Y_DIM], alloc=alloc)
        self._state_to_mapl_mapping[ndsl_field_name] = (mapl_state, mapl_field_name)

    @property
    def ndsl_state(self):
        return self._ndsl_state

    def fortran_to_ndsl(self) -> None:
        """Copy all Fortran memory in Python"""

        def _pull_from_fortran(
            mapl_field_: str,
            mapl_state_: MAPLMemoryRepository,
            ndsl_field_: str,
            ndsl_state_: State,
        ):
            if "." in ndsl_field_:
                inner_dataclass = ndsl_field_.split(".")[0]
                _pull_from_fortran(
                    mapl_field_,
                    mapl_state_,
                    "".join(ndsl_field_.split(".")[1:]),
                    getattr(ndsl_state_, inner_dataclass),
                )
            else:
                mapl_array = mapl_state_.get_from_fortran(mapl_field_)
                if self._transfer_type == InterfaceTransferType.CPU_TO_GPU_TO_CPU:
                    cp.cuda.runtime.deviceSynchronize()
                if mapl_array is None:
                    setattr(ndsl_state_, ndsl_field_, None)
                elif self._transfer_type == InterfaceTransferType.CPU_COPY:
                    getattr(ndsl_state_, ndsl_field_).field[:] = mapl_array[:]
                elif self._transfer_type == InterfaceTransferType.CPU_MAP:
                    getattr(ndsl_state_, ndsl_field_).data = mapl_array
                else:
                    raise ValueError("Transfer type unknown for Fortran/NDSL")

        for ndsl_field, (mapl_state, mapl_field) in self._state_to_mapl_mapping.items():
            try:
                _pull_from_fortran(mapl_field, mapl_state, ndsl_field, self._ndsl_state)
            except ValueError as e:
                e.add_note(f"Mapping {ndsl_field} to {mapl_field}")
                raise e

    def ndsl_to_fortran(self) -> None:
        """Copy all Python memory back in Fortran"""

        # Skip sending back - we are mapped
        if self._transfer_type == InterfaceTransferType.CPU_MAP:
            return

        def _push_back_to_fortran(
            mapl_field_: str,
            mapl_state_: MAPLMemoryRepository,
            ndsl_field_: str,
            ndsl_state_: State,
        ):
            if "." in ndsl_field_:
                inner_dataclass = ndsl_field_.split(".")[0]
                _push_back_to_fortran(
                    mapl_field_,
                    mapl_state_,
                    "".join(ndsl_field_.split(".")[1:]),
                    ndsl_state_.__getattribute__(inner_dataclass),
                )
            else:
                mapl_array = mapl_state_.get_from_fortran(mapl_field_)
                if mapl_array is None:
                    pass
                elif self._transfer_type == InterfaceTransferType.CPU_COPY:
                    ndsl_array = getattr(ndsl_state_, ndsl_field_).field[:]
                    mapl_array[:] = ndsl_array[:]
                    mapl_state_.send_to_fortran(mapl_field_)
                elif self._transfer_type == InterfaceTransferType.CPU_MAP:
                    raise RuntimeError("Coding issue. We should never send back mapped data")
                else:
                    raise ValueError("Transfer type unknown for NDSL/Fortran")

        for ndsl_field, (mapl_state, mapl_field) in self._state_to_mapl_mapping.items():
            _push_back_to_fortran(mapl_field, mapl_state, ndsl_field, self._ndsl_state)
        if self._transfer_type == InterfaceTransferType.CPU_TO_GPU_TO_CPU:
            cp.cuda.runtime.deviceSynchronize()
