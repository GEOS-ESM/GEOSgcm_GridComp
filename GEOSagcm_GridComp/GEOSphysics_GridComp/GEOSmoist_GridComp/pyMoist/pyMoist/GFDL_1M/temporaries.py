from dataclasses import dataclass

from ndsl import Quantity, QuantityFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Int


@dataclass
class Temporaries:
    p_interface_mb: Quantity
    p_mb: Quantity
    edge_height_above_surface: Quantity
    layer_height_above_surface: Quantity
    layer_thickness: Quantity
    layer_thickness_negative: Quantity
    dp: Quantity
    mass: Quantity
    qsat: Quantity
    dqsat: Quantity
    u_unmodified: Quantity
    v_unmodified: Quantity
    k_lcl: Quantity
    temporary_3d: Quantity
    th: Quantity
    temporary_2d_1: Quantity
    temporary_2d_2: Quantity
    all_zeros_3d: Quantity
    th700: Quantity
    t700: Quantity
    z700: Quantity
    du_dt: Quantity
    dv_dt: Quantity
    dt_dt: Quantity
    dvapor_dt: Quantity
    dliquid_dt: Quantity
    dice_dt: Quantity
    dcloud_fraction_dt: Quantity
    drain_dt: Quantity
    dsnow_dt: Quantity
    dgraupel_dt: Quantity

    @classmethod
    def make(cls, quantity_factory: QuantityFactory):
        p_interface_mb = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a")
        p_mb = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        edge_height_above_surface = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a")
        layer_height_above_surface = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        layer_thickness = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        layer_thickness_negative = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dp = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        mass = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        qsat = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dqsat = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        u_unmodified = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        v_unmodified = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        k_lcl = quantity_factory.zeros([X_DIM, Y_DIM], "n/a", dtype=Int)
        temporary_3d = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        th = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        temporary_2d_1 = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        temporary_2d_2 = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        all_zeros_3d = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        th700 = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        t700 = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        z700 = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        du_dt = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dv_dt = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dt_dt = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dvapor_dt = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dliquid_dt = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dice_dt = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dcloud_fraction_dt = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        drain_dt = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dsnow_dt = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dgraupel_dt = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        return cls(
            p_interface_mb,
            p_mb,
            edge_height_above_surface,
            layer_height_above_surface,
            layer_thickness,
            layer_thickness_negative,
            dp,
            mass,
            qsat,
            dqsat,
            u_unmodified,
            v_unmodified,
            k_lcl,
            temporary_3d,
            th,
            temporary_2d_1,
            temporary_2d_2,
            all_zeros_3d,
            th700,
            t700,
            z700,
            du_dt,
            dv_dt,
            dt_dt,
            dvapor_dt,
            dliquid_dt,
            dice_dt,
            dcloud_fraction_dt,
            drain_dt,
            dsnow_dt,
            dgraupel_dt,
        )
