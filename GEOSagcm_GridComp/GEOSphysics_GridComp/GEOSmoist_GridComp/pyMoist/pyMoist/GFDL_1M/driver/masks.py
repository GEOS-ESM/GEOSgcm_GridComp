from dataclasses import dataclass

from ndsl import Quantity, QuantityFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Int


@dataclass
class Masks:
    is_frozen: Quantity
    precip_fall: Quantity
    melting_mask_1: Quantity
    melting_mask_2: Quantity
    current_k_level: Quantity

    @classmethod
    def make(cls, quantity_factory: QuantityFactory):
        is_frozen = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a", dtype=bool)
        precip_fall = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        # TODO: temporary code requires mask used within a double k loop
        # will be removed once a proper feature for double k loop is introduces
        melting_mask_1 = quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_DIM], "n/a", dtype=bool
        )
        melting_mask_2 = quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_DIM], "n/a", dtype=bool
        )
        current_k_level = quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_DIM], "n/a", dtype=Int
        )
        for k in range(current_k_level.view[:].shape[2]):
            current_k_level.view[:, :, k] = k
        return cls(
            is_frozen, precip_fall, melting_mask_1, melting_mask_2, current_k_level
        )
