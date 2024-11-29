import copy
import gt4py.cartesian.gtscript as gtscript
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float, Int
from ndsl import QuantityFactory


# class test:
#     def __init__(self, number) -> None:
#         self.x = number


# if __name__ == "__main__":
#     out = test(6)
#     print(out.x)


class data_dim_factory:
    """
    Factory capable of adding an extra dimension of specified size to any of the primary
    dimensions (X, Y, Z)
    """

    def __init__(
        self,
        quantity_factory: QuantityFactory,
        data_type: int,
        dim_size: int,
        x_flag: bool = False,
        y_flag: bool = False,
        z_flag: bool = False,
        units: str = "n/a",
    ):
        """
        Function for creating FloatField with extra dimension.

        data_type specified the type of created Field: Float (1), Int (2), bool (3)
        dim_size specifies size of extra dimension

        Dimension flags specify which axis are returned on final quantity. The dimension upon which
        the added dimension will be build is required, all others are optional.
        """
        flags = [x_flag, y_flag, z_flag]
        self.check_input(flags, data_type)

        self.data_dim_quantity_factory = self.make_data_dim_quantity_factory(
            quantity_factory, dim_size
        )

        if data_type == 1:
            if x_flag:
                self.data_dim_field = self.data_dim_quantity_factory.zeros(
                    [X_DIM, "data_axis"], units, dtype=Float
                )
            if y_flag:
                self.data_dim_field = self.data_dim_quantity_factory.zeros(
                    [Y_DIM, "data_axis"], units, dtype=Float
                )
            if z_flag:
                self.data_dim_field = self.data_dim_quantity_factory.zeros(
                    [Z_DIM, "data_axis"], units, dtype=Float
                )

            if x_flag and y_flag:
                self.data_dim_field = self.data_dim_quantity_factory.zeros(
                    [X_DIM, Y_DIM, "data_axis"], units, dtype=Float
                )
            if x_flag and z_flag:
                self.data_dim_field = self.data_dim_quantity_factory.zeros(
                    [X_DIM, Z_DIM, "data_axis"], units, dtype=Float
                )
            if y_flag and z_flag:
                self.data_dim_field = self.data_dim_quantity_factory.zeros(
                    [Y_DIM, Z_DIM, "data_axis"], units, dtype=Float
                )
            if x_flag and y_flag and z_flag:
                self.data_dim_field = self.data_dim_quantity_factory.zeros(
                    [X_DIM, Y_DIM, Z_DIM, "data_axis"], units, dtype=Float
                )

        if data_type == 2:
            if x_flag:
                self.data_dim_field = self.data_dim_quantity_factory.zeros(
                    [X_DIM, "data_axis"], units, dtype=Int
                )
            if y_flag:
                self.data_dim_field = self.data_dim_quantity_factory.zeros(
                    [Y_DIM, "data_axis"], units, dtype=Int
                )
            if z_flag:
                self.data_dim_field = self.data_dim_quantity_factory.zeros(
                    [Z_DIM, "data_axis"], units, dtype=Int
                )

            if x_flag and y_flag:
                self.data_dim_field = self.data_dim_quantity_factory.zeros(
                    [X_DIM, Y_DIM, "data_axis"], units, dtype=Int
                )
            if x_flag and z_flag:
                self.data_dim_field = self.data_dim_quantity_factory.zeros(
                    [X_DIM, Z_DIM, "data_axis"], units, dtype=Int
                )
            if y_flag and z_flag:
                self.data_dim_field = self.data_dim_quantity_factory.zeros(
                    [Y_DIM, Z_DIM, "data_axis"], units, dtype=Int
                )
            if x_flag and y_flag and z_flag:
                self.data_dim_field = self.data_dim_quantity_factory.zeros(
                    [X_DIM, Y_DIM, Z_DIM, "data_axis"], units, dtype=Int
                )

    def check_input(self, flags, data_type):
        if any(flags) is not True:
            raise ValueError(
                "At least one dimension must be present for the created Field"
            )
        if data_type != 1 and data_type != 2 and data_type != 3:
            raise ValueError(
                "data_type must specify a data type: Float (1), Int (2), bool (3)."
            )
        if data_type == 3:
            raise ValueError(
                "Extra dimension boolean fields not yet implemented. \
                    Contact development team for further assistance."
            )

    @staticmethod
    def make_data_dim_quantity_factory(
        ijk_quantity_factory: QuantityFactory, size: Int
    ):
        data_dim_quantity_factory = copy.deepcopy(ijk_quantity_factory)
        data_dim_quantity_factory.set_extra_dim_lengths(
            **{
                "data_axis": size,
            }
        )
        return data_dim_quantity_factory


if __name__ == "__main__":
    test = data_dim_factory(QuantityFactory, 1, 100, z_flag=True)
    print(test.data_dim_field)
