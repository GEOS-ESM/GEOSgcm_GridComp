import numpy as np
from gt4py.cartesian.gtscript import i32
from ndsl.dsl.typing import Float
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
import pyMoist.GFDL_1M.GFDL_1M_driver.constants as constants
from pyMoist.GFDL_1M.GFDL_1M_driver.config import config
from ndsl.dsl.typing import Int


class temporaries:
    def __init__(self, quantity_factory):
        self.lhi = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.icpk = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.cvm = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.m1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.dm = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        # used as placeholders to note what features are needed. will be r
        self.need_2d_temporaries_feature = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        self.need_double_k_loop_feature = quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_DIM], "n/a"
        )
