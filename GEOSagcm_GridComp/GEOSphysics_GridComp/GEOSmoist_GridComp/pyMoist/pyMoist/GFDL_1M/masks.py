from dataclasses import dataclass

from ndsl import Quantity, QuantityFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Int


@dataclass
class Masks:
    boolean_2d_mask: Quantity

    @classmethod
    def make(cls, quantity_factory: QuantityFactory):
        boolean_2d_mask = quantity_factory.zeros([X_DIM, Y_DIM], "n/a", dtype=bool)
        return cls(boolean_2d_mask)
