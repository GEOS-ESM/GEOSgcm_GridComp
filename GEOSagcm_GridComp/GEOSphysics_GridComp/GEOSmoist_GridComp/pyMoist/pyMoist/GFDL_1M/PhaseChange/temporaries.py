from dataclasses import dataclass

from ndsl import Quantity, QuantityFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM


@dataclass
class Temporaries:
    temporary: Quantity
    minrhcrit: Quantity
    alpha: Quantity

    @classmethod
    def make(cls, quantity_factory: QuantityFactory):
        temporary = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        minrhcrit = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        alpha = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        return cls(temporary, minrhcrit, alpha)
