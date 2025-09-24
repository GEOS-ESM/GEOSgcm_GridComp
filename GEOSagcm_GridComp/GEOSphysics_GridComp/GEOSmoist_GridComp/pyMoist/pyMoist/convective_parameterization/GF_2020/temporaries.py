from dataclasses import dataclass

from ndsl import Quantity, QuantityFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Int


@dataclass
class Temporaries:
    # first used in Setup
    p: Quantity
    p_kappa: Quantity
    edge_height_above_surface: Quantity
    layer_height_above_surface: Quantity
    mass: Quantity
    th: Quantity
    tpwi: Quantity
    tpwi_star: Quantity
    seed_convection: Quantity
    modified_area: Quantity
    # first used in DriverInterface
    aot500: Quantity
    temp2m: Quantity
    sflux_r: Quantity
    sflux_t: Quantity
    topt: Quantity
    xland: Quantity
    dx2d: Quantity
    kpbl: Quantity
    temp: Quantity
    press: Quantity
    rvap: Quantity
    up: Quantity
    vp: Quantity
    wp: Quantity
    zt3d: Quantity
    zm3d: Quantity
    dm3d: Quantity
    khloc: Quantity
    curr_rvap: Quantity
    mp_ice: Quantity
    mp_liq: Quantity
    mp_cf: Quantity
    buoy_exc: Quantity
    DZ: Quantity
    AIR_DEN: Quantity
    entr3d: Quantity

    @classmethod
    def make(cls, quantity_factory: QuantityFactory):
        # first used in Setup
        p = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        p_kappa = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        edge_height_above_surface = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a")
        layer_height_above_surface = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        mass = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        th = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        tpwi = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        tpwi_star = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        seed_convection = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        modified_area = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        # first used in DriverInterface
        aot500 = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        temp2m = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        sflux_r = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        sflux_t = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        topt = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        xland = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        dx2d = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        kpbl = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        temp = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        press = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        rvap = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        up = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        vp = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        wp = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        zt3d = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        zm3d = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dm3d = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        khloc = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        curr_rvap = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        mp_ice = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM, "nmp"], "n/a")
        mp_liq = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM, "nmp"], "n/a")
        mp_cf = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM, "nmp"], "n/a")
        buoy_exc = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        DZ = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        AIR_DEN = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        entr3d = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        return cls(
            p,
            p_kappa,
            edge_height_above_surface,
            layer_height_above_surface,
            mass,
            th,
            tpwi,
            tpwi_star,
            seed_convection,
            modified_area,
            aot500,
            temp2m,
            sflux_r,
            sflux_t,
            topt,
            xland,
            dx2d,
            kpbl,
            temp,
            press,
            rvap,
            up,
            vp,
            wp,
            zt3d,
            zm3d,
            dm3d,
            khloc,
            curr_rvap,
            mp_ice,
            mp_liq,
            mp_cf,
            buoy_exc,
            DZ,
            AIR_DEN,
            entr3d,
        )
