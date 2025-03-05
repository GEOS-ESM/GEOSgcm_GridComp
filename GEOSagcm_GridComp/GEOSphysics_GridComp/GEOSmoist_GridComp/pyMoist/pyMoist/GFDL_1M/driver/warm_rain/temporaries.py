from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl import Quantity, QuantityFactory
from dataclasses import dataclass


@dataclass
class Temporaries:
    dm: Quantity

    @classmethod
    def make(cls, quantity_factory: QuantityFactory):
        dm = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        return cls(dm)
