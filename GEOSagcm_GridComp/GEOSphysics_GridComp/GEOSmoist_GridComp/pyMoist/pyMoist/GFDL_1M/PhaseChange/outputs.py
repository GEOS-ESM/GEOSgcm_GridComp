from dataclasses import dataclass

from ndsl import Quantity, QuantityFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM


@dataclass
class Outputs:
    rhx: Quantity
    evapc: Quantity
    sublc: Quantity

    @classmethod
    def make(cls, quantity_factory: QuantityFactory):
        # -----------------------------------------------------------------------
        # initialize precipitation outputs
        # -----------------------------------------------------------------------
        rhx = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        evapc = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        sublc = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        return cls(rhx, evapc, sublc)
