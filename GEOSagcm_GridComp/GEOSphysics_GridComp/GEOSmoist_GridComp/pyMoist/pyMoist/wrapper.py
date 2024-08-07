"""
Wraps pyMoist for GEOS interface use.
"""

from dataclasses import dataclass
from typing import Any

from pyMoist.aer_activation import AerActivation


@dataclass
class MoistFlags:
    npx: int
    npy: int
    npz: int
    layout_x: int
    layout_y: int
    n_tiles: int
    # Aer Activation
    n_modes: int
    # Magic number
    mn_123456789: int


class GEOSPyMoistWrapper:
    def __init__(self) -> None:
        pass

    def __call__(self, *args: Any, **kwds: Any) -> Any:
        pass
