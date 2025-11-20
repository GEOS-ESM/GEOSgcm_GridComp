from dataclasses import dataclass

from ndsl import Local, NDSLRuntime, QuantityFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM


@dataclass
class TerminalFallLocals:
    lhi: Local
    icpk: Local
    cvm: Local
    m1: Local
    dm: Local
    need_2d_temporaries_feature: Local
    need_double_k_loop_feature: Local

    @classmethod
    def make(cls, runtime: NDSLRuntime, quantity_factory: QuantityFactory):
        lhi = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        icpk = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        cvm = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        m1 = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])
        dm = runtime.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM])

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
