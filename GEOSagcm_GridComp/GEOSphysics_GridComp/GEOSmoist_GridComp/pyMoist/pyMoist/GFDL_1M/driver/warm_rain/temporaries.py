from ndsl.constants import X_DIM, Y_DIM, Z_DIM


class Temporaries:
    def __init__(self, quantity_factory):
        self.dm = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
