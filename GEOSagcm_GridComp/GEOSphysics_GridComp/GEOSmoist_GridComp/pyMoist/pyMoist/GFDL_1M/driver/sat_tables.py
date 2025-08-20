import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import (
    FORWARD,
    PARALLEL,
    K,
    computation,
    exp,
    interval,
    log,
    log10,
)

from ndsl.boilerplate import get_factories_single_tile
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float, FloatField, Int
from pyMoist.GFDL_1M.driver.constants import constants
from pyMoist.shared_incloud_processes import ice_fraction


# Workaround to create a 1d off-grid axis that can be written to
GlobalTable_driver_qsat = gtscript.GlobalTable[(Float, (int(constants.LENGTH)))]


def qs_table_1(length: Int, table1: FloatField, esupc: FloatField):
    """
    Compute saturation water vapor pressure table 1
    three phase table

    reference Fortran: gfdl_cloud_microphys.F90: subroutine qs_table
    """
    with computation(PARALLEL), interval(0, 1600):
        # -----------------------------------------------------------------------
        # compute es over ice between - 160 deg c and 0 deg c.
        # -----------------------------------------------------------------------
        tem = constants.TMIN + constants.DELT * K
        fac0 = (tem - constants.T_ICE) / (tem * constants.T_ICE)
        fac1 = fac0 * constants.LI2
        fac2 = (constants.D2ICE * log(tem / constants.T_ICE) + fac1) / constants.RVGAS
        table1 = constants.E_00 * exp(fac2)

    with computation(PARALLEL), interval(0, 1421):
        # -----------------------------------------------------------------------
        # compute es over water between - 40 deg c and 102 deg c.
        # -----------------------------------------------------------------------
        tem = 233.16 + constants.DELT * K
        fac0 = (tem - constants.T_ICE) / (tem * constants.T_ICE)
        fac1 = fac0 * constants.LV0
        fac2 = (constants.DC_VAP * log(tem / constants.T_ICE) + fac1) / constants.RVGAS
        esh40 = constants.E_00 * exp(fac2)
        if K < 400:
            esupc = esh40

    with computation(PARALLEL), interval(400 + 1200, 1421 + 1200):
        table1 = esh40[0, 0, -1200]

    with computation(PARALLEL), interval(0 + 1200, 400 + 1200):
        # -----------------------------------------------------------------------
        # derive blended es over ice and supercooled
        # water between - 40 deg c and 0 deg c
        # -----------------------------------------------------------------------
        tem = 233.16 + constants.DELT * (K - 1200)
        # GEOS ! WMP impose CALIPSO ice polynomial from 0 C to -40 C
        wice = ice_fraction(tem, 0.0, 0.0)
        wh2o = 1.0 - wice
        table1 = wice * table1 + wh2o * esupc[0, 0, -1200]


def qs_table_2(length: Int, table2: FloatField):
    """
    Compute saturation water vapor pressure table 2
    one phase table

    reference Fortran: gfdl_cloud_microphys.F90: subroutine qs_tablew
    """
    with computation(PARALLEL), interval(...):
        tem = constants.TMIN + constants.DELT * K
        fac0 = (tem - constants.T_ICE) / (tem * constants.T_ICE)
        fac1 = fac0 * constants.LV0
        fac2 = (constants.DC_VAP * log(tem / constants.T_ICE) + fac1) / constants.RVGAS
        table2 = constants.E_00 * exp(fac2)


def qs_table_3(length: Int, table3: FloatField, table1: FloatField):
    """
    Compute saturation water vapor pressure table 3
    two phase table

    reference Fortran: gfdl_cloud_microphys.F90: subroutine qs_table2
    """
    with computation(PARALLEL), interval(...):
        tem0 = constants.TMIN + constants.DELT * K
        fac0 = (tem0 - constants.T_ICE) / (tem0 * constants.T_ICE)
        if K < 1600:
            # -----------------------------------------------------------------------
            # compute es over ice between - 160 deg c and 0 deg c.
            # -----------------------------------------------------------------------
            fac1 = fac0 * constants.LI2
            fac2 = (
                constants.D2ICE * log(tem0 / constants.T_ICE) + fac1
            ) / constants.RVGAS
        else:
            # -----------------------------------------------------------------------
            # compute es over water between 0 deg c and 102 deg c.
            # -----------------------------------------------------------------------
            fac1 = fac0 * constants.LV0
            fac2 = (
                constants.DC_VAP * log(tem0 / constants.T_ICE) + fac1
            ) / constants.RVGAS
        table3 = constants.E_00 * exp(fac2)

    with computation(FORWARD), interval(1599, 1601):
        # -----------------------------------------------------------------------
        # smoother around 0 deg c
        # -----------------------------------------------------------------------
        if K == 1599:
            t0 = 0.25 * (table3[0, 0, -1] + 2.0 * table1 + table3[0, 0, 1])
            t1 = 0.25 * (table3 + 2.0 * table1[0, 0, 1] + table3[0, 0, 2])
            table3 = t0
        elif K == 1600:
            table3 = t1  # type: ignore


def qs_table_4(length: Int, table4: FloatField, table1: FloatField):
    """
    Compute saturation water vapor pressure table 4
    two phase table with " - 2 c" as the transition point

    reference Fortran: gfdl_cloud_microphys.F90: subroutine qs_table3
    """
    with computation(PARALLEL), interval(...):
        tem = constants.TMIN + constants.DELT * K
        if K < 1580:  # change to - 2 c
            # -----------------------------------------------------------------------
            # compute es over ice between - 160 deg c and 0 deg c.
            # see smithsonian meteorological tables page 350.
            # -----------------------------------------------------------------------
            aa = -9.09718 * (constants.TABLE_ICE / tem - 1.0)
            b = -3.56654 * log10(constants.TABLE_ICE / tem)
            c = 0.876793 * (1.0 - tem / constants.TABLE_ICE)
            e = log10(constants.ESBASI)
            table4 = 0.1 * 10 ** (aa + b + c + e)
        else:
            # -----------------------------------------------------------------------
            # compute es over water between - 2 deg c and 102 deg c.
            # see smithsonian meteorological tables page 350.
            # -----------------------------------------------------------------------
            aa = -7.90298 * (constants.TBASW / tem - 1.0)
            b = 5.02808 * log10(constants.TBASW / tem)
            c = -1.3816e-7 * (10 ** ((1.0 - tem / constants.TBASW) * 11.344) - 1.0)
            d = 8.1328e-3 * (10 ** ((constants.TBASW / tem - 1.0) * (-3.49149)) - 1.0)
            e = log10(constants.ESBASW)
            table4 = 0.1 * 10 ** (aa + b + c + d + e)

    with computation(FORWARD), interval(1579, 1581):
        # -----------------------------------------------------------------------
        # smoother around - 2 deg c
        # -----------------------------------------------------------------------
        if K == 1579:
            t0 = 0.25 * (table4[0, 0, -1] + 2.0 * table1 + table4[0, 0, 1])
            t1 = 0.25 * (table4 + 2.0 * table1[0, 0, 1] + table4[0, 0, 2])
            table4 = t0
        elif K == 1580:
            table4 = t1  # type: ignore


def des_tables(
    length: Int,
    des1: FloatField,
    des2: FloatField,
    des3: FloatField,
    des4: FloatField,
    table1: FloatField,
    table2: FloatField,
    table3: FloatField,
    table4: FloatField,
):
    with computation(FORWARD), interval(...):
        des1 = max(0.0, table1[0, 0, 1] - table1)
        des2 = max(0.0, table2[0, 0, 1] - table2)
        des3 = max(0.0, table3[0, 0, 1] - table3)
        des4 = max(0.0, table4[0, 0, 1] - table4)

    with computation(PARALLEL), interval(-1, None):
        des1 = des1[0, 0, -1]
        des2 = des2[0, 0, -1]
        des3 = des3[0, 0, -1]
        des4 = des4[0, 0, -1]


class GFDL_driver_tables:
    """
    Initializes lookup tables for saturation water vapor pressure
    for the utility routines that are designed to return qs
    consistent with the assumptions in FV3.

    Reference Fortran: gfdl_cloud_microphys.F90: qsmith_init.py
    """

    def __init__(self, backend):
        table_compute_domain = (1, 1, constants.LENGTH)

        stencil_factory, quantity_factory = get_factories_single_tile(
            table_compute_domain[0],
            table_compute_domain[1],
            table_compute_domain[2],
            0,
            backend,
        )

        self._table1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._table2 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._table3 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._table4 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._des1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._des2 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._des3 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._des4 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._esupc = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        compute_qs_table_1 = stencil_factory.from_origin_domain(
            func=qs_table_1,
            origin=(0, 0, 0),
            domain=table_compute_domain,
        )
        compute_qs_table_2 = stencil_factory.from_origin_domain(
            func=qs_table_2,
            origin=(0, 0, 0),
            domain=table_compute_domain,
        )
        compute_qs_table_3 = stencil_factory.from_origin_domain(
            func=qs_table_3,
            origin=(0, 0, 0),
            domain=table_compute_domain,
        )
        compute_qs_table_4 = stencil_factory.from_origin_domain(
            func=qs_table_4,
            origin=(0, 0, 0),
            domain=table_compute_domain,
        )
        compute_des_tables = stencil_factory.from_origin_domain(
            func=des_tables,
            origin=(0, 0, 0),
            domain=table_compute_domain,
        )

        compute_qs_table_1(constants.LENGTH, self._table1, self._esupc)
        compute_qs_table_2(constants.LENGTH, self._table2)
        compute_qs_table_3(constants.LENGTH, self._table3, self._table1)
        compute_qs_table_4(constants.LENGTH, self._table4, self._table1)
        compute_des_tables(
            constants.LENGTH,
            self._des1,
            self._des2,
            self._des3,
            self._des4,
            self._table1,
            self._table2,
            self._table3,
            self._table4,
        )

        self.table1 = self._table1.view[0, 0, :]
        self.table2 = self._table2.view[0, 0, :]
        self.table3 = self._table3.view[0, 0, :]
        self.table4 = self._table4.view[0, 0, :]
        self.des1 = self._des1.view[0, 0, :]
        self.des2 = self._des2.view[0, 0, :]
        self.des3 = self._des3.view[0, 0, :]
        self.des4 = self._des4.view[0, 0, :]


# Table needs to be calculated only once
_cached_table = {
    "driver_qsat": None,
}


def get_tables(backend):
    if _cached_table["driver_qsat"] is None:
        _cached_table["driver_qsat"] = GFDL_driver_tables(backend)

    return _cached_table["driver_qsat"]
