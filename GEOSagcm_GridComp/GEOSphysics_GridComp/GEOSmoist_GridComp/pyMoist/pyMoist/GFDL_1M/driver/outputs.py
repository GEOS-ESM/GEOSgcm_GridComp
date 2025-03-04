from ndsl.constants import X_DIM, Y_DIM, Z_DIM


class Outputs:
    def __init__(self, quantity_factory):
        # -----------------------------------------------------------------------
        # initialize precipitation outputs
        # -----------------------------------------------------------------------

        self.rain = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        self.snow = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        self.ice = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        self.graupel = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        self.m2_rain = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.m2_sol = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.revap = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.isubl = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
