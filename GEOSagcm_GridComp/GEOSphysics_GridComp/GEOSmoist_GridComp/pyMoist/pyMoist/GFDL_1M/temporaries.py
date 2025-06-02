from dataclasses import dataclass

from ndsl import Quantity, QuantityFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM


@dataclass
class Temporaries:
    p_interface_mb: Quantity
    p_mb: Quantity
    edge_height_above_surface: Quantity
    layer_height_above_surface: Quantity
    layer_thickness: Quantity
    dp: Quantity
    mass: Quantity
    qsat: Quantity
    dqsat: Quantity
    u_unmodified: Quantity
    v_unmodified: Quantity
    k_lcl: Quantity
    lower_tropospheric_stability: Quantity
    estimated_inversion_strength: Quantity
    test_field: Quantity

    @classmethod
    def make(cls, quantity_factory: QuantityFactory):
        p_interface_mb = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a")
        p_mb = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        edge_height_above_surface = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a")
        layer_height_above_surface = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        layer_thickness = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dp = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        mass = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qsat = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dqsat = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        u_unmodified = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        v_unmodified = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        k_lcl = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        lower_tropospheric_stability = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        estimated_inversion_stregnth = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        test_field = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        return cls(
            p_interface_mb,
            p_mb,
            edge_height_above_surface,
            layer_height_above_surface,
            layer_thickness,
            dp,
            mass,
            qsat,
            dqsat,
            u_unmodified,
            v_unmodified,
            k_lcl,
            lower_tropospheric_stability,
            estimated_inversion_stregnth,
            test_field,
        )
