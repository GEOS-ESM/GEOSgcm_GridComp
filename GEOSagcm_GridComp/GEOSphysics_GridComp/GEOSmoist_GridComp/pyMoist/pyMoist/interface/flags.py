from dataclasses import dataclass

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import cffi


@dataclass
class MoistFlags:
    npx: int = 0
    npy: int = 0
    npz: int = 0
    layout_x: int = 1
    layout_y: int = 1
    n_tiles: int = 6
    # Aer Activation
    n_modes: int = 0
    # Magic number
    mn_123456789: int = 0


def _generic_config_bridge(
    py_flags: MoistFlags,
    fv_flags: "cffi.FFI.CData",
):
    keys = list(filter(lambda k: not k.startswith("__"), dir(type(py_flags))))
    for k in keys:
        if hasattr(fv_flags, k):
            setattr(py_flags, k, getattr(fv_flags, k))


def flags_fv_to_python(
    fv_flags: "cffi.FFI.CData",
) -> MoistFlags:
    if fv_flags.mn_123456789 != 123456789:
        raise RuntimeError(
            "Magic number failed, pyMoist interface is broken on the python side"
        )

    py_flags = MoistFlags()
    _generic_config_bridge(py_flags, fv_flags)
    py_flags.layout = (
        getattr(fv_flags, "layout_x"),
        getattr(fv_flags, "layout_y"),
    )
    return py_flags
