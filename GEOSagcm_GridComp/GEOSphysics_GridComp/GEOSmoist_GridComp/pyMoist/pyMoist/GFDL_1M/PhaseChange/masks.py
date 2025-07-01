from dataclasses import dataclass

from ndsl import Quantity, QuantityFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM


@dataclass
class Masks:
    k_mask: Quantity

    @classmethod
    def make(cls, quantity_factory: QuantityFactory):
        k_mask = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        for i in range(0, k_mask.view[:].shape[0]):
            for j in range(0, k_mask.view[:].shape[1]):
                for k in range(0, k_mask.view[:].shape[2]):
                    k_mask.view[i, j, k] = k + 1
        return cls(k_mask)
