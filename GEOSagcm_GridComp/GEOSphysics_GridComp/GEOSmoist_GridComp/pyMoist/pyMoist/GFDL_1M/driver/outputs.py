from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl import Quantity, QuantityFactory
from dataclasses import dataclass


@dataclass
class Outputs:
    rain: Quantity
    snow: Quantity
    ice: Quantity
    graupel: Quantity
    m2_rain: Quantity
    m2_sol: Quantity
    revap: Quantity
    isubl: Quantity

    @classmethod
    def make(cls, quantity_factory: QuantityFactory):
        # -----------------------------------------------------------------------
        # initialize precipitation outputs
        # -----------------------------------------------------------------------
        rain = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        snow = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        ice = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        graupel = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        m2_rain = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        m2_sol = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        revap = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        isubl = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        return cls(
            rain,
            snow,
            ice,
            graupel,
            m2_rain,
            m2_sol,
            revap,
            isubl,
        )
