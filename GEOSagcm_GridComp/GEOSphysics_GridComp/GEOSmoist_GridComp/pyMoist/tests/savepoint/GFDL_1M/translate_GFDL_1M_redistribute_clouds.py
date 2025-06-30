from ndsl import Namelist, Quantity, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.redistribute_clouds import RedistributeClouds


class TranslateGFDL_1M_redistribute_clouds(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "RAD_CF": grid.compute_dict(),
            "RAD_QL": grid.compute_dict(),
            "RAD_QI": grid.compute_dict(),
            "CLCN": grid.compute_dict(),
            "CLLS": grid.compute_dict(),
            "QLCN": grid.compute_dict(),
            "QLLS": grid.compute_dict(),
            "QICN": grid.compute_dict(),
            "QILS": grid.compute_dict(),
            "RAD_QV": grid.compute_dict(),
            "T": grid.compute_dict(),
        }

        # FloatField Outputs
        self.out_vars = self.in_vars["data_vars"].copy()

    def make_ijk_quantity(self, data, interface: bool = False) -> Quantity:
        if interface is True:
            quantity = self.quantity_factory.empty([X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a")
            quantity.view[:, :, :] = quantity.np.asarray(data[:, :, :])
            return quantity
        else:
            quantity = self.quantity_factory.empty([X_DIM, Y_DIM, Z_DIM], "n/a")
            quantity.view[:, :, :] = quantity.np.asarray(data[:, :, :])
            return quantity

    # Calculated Outputs
    def compute(self, inputs):
        cloud_fraction = self.make_ijk_quantity(inputs.pop("RAD_CF"))
        convective_cloud_fraction = self.make_ijk_quantity(inputs.pop("CLCN"))
        large_scale_cloud_fraction = self.make_ijk_quantity(inputs.pop("CLLS"))
        liquid = self.make_ijk_quantity(inputs.pop("RAD_QL"))
        convective_liquid = self.make_ijk_quantity(inputs.pop("QLCN"))
        large_scale_liquid = self.make_ijk_quantity(inputs.pop("QLLS"))
        ice = self.make_ijk_quantity(inputs.pop("RAD_QI"))
        convective_ice = self.make_ijk_quantity(inputs.pop("QICN"))
        large_scale_ice = self.make_ijk_quantity(inputs.pop("QILS"))
        vapor = self.make_ijk_quantity(inputs.pop("RAD_QV"))
        temperature = self.make_ijk_quantity(inputs.pop("T"))
        redistribute_clouds = RedistributeClouds(
            stencil_factory=self.stencil_factory,
        )
        redistribute_clouds(
            cloud_fraction=cloud_fraction,
            convective_cloud_fraction=convective_cloud_fraction,
            large_scale_cloud_fraction=large_scale_cloud_fraction,
            liquid=liquid,
            convective_liquid=convective_liquid,
            large_scale_liquid=large_scale_liquid,
            ice=ice,
            convective_ice=convective_ice,
            large_scale_ice=large_scale_ice,
            vapor=vapor,
            temperature=temperature,
        )
        return {
            "RAD_CF": cloud_fraction.field,
            "CLCN": convective_cloud_fraction.field,
            "CLLS": large_scale_cloud_fraction.field,
            "RAD_QL": liquid.field,
            "QLCN": convective_liquid.field,
            "QLLS": large_scale_liquid.field,
            "RAD_QI": ice.field,
            "QICN": convective_ice.field,
            "QILS": large_scale_ice.field,
            "RAD_QV": vapor.field,
            "T": temperature.field,
        }
