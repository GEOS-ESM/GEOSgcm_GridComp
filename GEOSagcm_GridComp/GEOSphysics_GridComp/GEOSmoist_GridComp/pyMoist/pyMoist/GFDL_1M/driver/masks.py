from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Int


class Masks:
    def __init__(self, quantity_factory):
        self.is_frozen = quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_DIM], "n/a", dtype=bool
        )
        self.precip_fall = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        # TODO: temporary code requires mask used within a double k loop
        # will be removed once a proper feature for double k loop is introduces
        self.melting_mask_1 = quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_DIM], "n/a", dtype=bool
        )
        self.melting_mask_2 = quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_DIM], "n/a", dtype=bool
        )
        self.current_k_level = quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_DIM], "n/a", dtype=Int
        )
        for k in range(self.current_k_level.view[:].shape[2]):
            self.current_k_level.view[:, :, k] = k
