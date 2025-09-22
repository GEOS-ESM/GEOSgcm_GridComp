import dataclasses

from ndsl import Quantity, QuantityFactory
from ndsl.constants import X_DIM, Y_DIM

from ndsl.dsl.dace.orchestration import dace_inhibitor


@dataclasses.dataclass
class Masks:
    boolean_2d_mask: Quantity

    @classmethod
    def make(cls, quantity_factory: QuantityFactory):
        boolean_2d_mask = quantity_factory.zeros([X_DIM, Y_DIM], "n/a", dtype=bool)
        return cls(boolean_2d_mask)

    @dace_inhibitor
    def zeros(self):
        for field in dataclasses.fields(Masks):
            getattr(self, field.name).data[:] = 0
