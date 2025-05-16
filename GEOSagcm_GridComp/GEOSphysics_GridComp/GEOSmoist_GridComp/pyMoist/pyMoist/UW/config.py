from dataclasses import dataclass

from ndsl.dsl.typing import Int


@dataclass
class UWConfiguration:
    NCNST: Int
    k0: Int
    windsrcavg: Int
