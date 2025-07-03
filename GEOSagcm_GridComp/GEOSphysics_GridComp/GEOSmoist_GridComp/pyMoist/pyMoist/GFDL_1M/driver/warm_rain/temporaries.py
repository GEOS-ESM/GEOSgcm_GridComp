from dataclasses import dataclass

from ndsl import Quantity, QuantityFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM


@dataclass
class Temporaries:
    dm: Quantity

    @classmethod
    def make(cls, quantity_factory: QuantityFactory):
        dm = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        return cls(dm)
