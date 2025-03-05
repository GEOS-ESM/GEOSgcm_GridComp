from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl import Quantity, QuantityFactory
from dataclasses import dataclass


@dataclass
class Temporaries:
    lhi: Quantity
    icpk: Quantity
    cvm: Quantity
    m1: Quantity
    dm: Quantity
    need_2d_temporaries_feature: Quantity
    need_double_k_loop_feature: Quantity

    @classmethod
    def make(cls, quantity_factory: QuantityFactory):
        lhi = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        icpk = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        cvm = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        m1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dm = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        # used as placeholders to note what features are needed. will be r
        need_2d_temporaries_feature = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        need_double_k_loop_feature = quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_DIM], "n/a"
        )
        return cls(
            lhi,
            icpk,
            cvm,
            m1,
            dm,
            need_2d_temporaries_feature,
            need_double_k_loop_feature,
        )
