import copy
import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import computation, interval, PARALLEL, exp, log, log10
from ndsl import QuantityFactory, StencilFactory
from ndsl.boilerplate import get_factories_single_tile
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float, Int, FloatField
import pyMoist.GFDL_1M.GFDL_1M_driver.GFDL_1M_driver_constants as driver_constants
from pyMoist.shared_incloud_processes import ice_fraction


length = 2621
_FloatField_data_dim = gtscript.Field[gtscript.IJK, (Float, (int(length)))]
GlobalTable_driver_qsat = gtscript.GlobalTable[(Float, (length))]


def qs_table_1(length: Int, table1: _FloatField_data_dim, esupc: _FloatField_data_dim):
    """
    compute saturation water vapor pressure table 1
    three phase table

    reference Fortran: gfdl_cloud_microphys.F90: subroutine qs_table
    """
    with computation(PARALLEL), interval(...):
        delt = 0.1
        tmin = driver_constants.table_ice - 160.0

        # -----------------------------------------------------------------------
        # compute es over ice between - 160 deg c and 0 deg c.
        # -----------------------------------------------------------------------

        i = 0
        while i < 1600:
            tem = tmin + delt * i
            fac0 = (tem - driver_constants.t_ice) / (tem * driver_constants.t_ice)
            fac1 = fac0 * driver_constants.li2
            fac2 = (
                driver_constants.d2ice * log(tem / driver_constants.t_ice) + fac1
            ) / driver_constants.rvgas
            table1[0, 0, 0][i] = driver_constants.e00 * exp(fac2)
            i = i + 1

        # -----------------------------------------------------------------------
        # compute es over water between - 40 deg c and 102 deg c.
        # -----------------------------------------------------------------------

        i = 0
        while i < 1421:
            tem = 233.16 + delt * i
            fac0 = (tem - driver_constants.t_ice) / (tem * driver_constants.t_ice)
            fac1 = fac0 * driver_constants.lv0
            fac2 = (
                driver_constants.dc_vap * log(tem / driver_constants.t_ice) + fac1
            ) / driver_constants.rvgas
            esh40 = driver_constants.e00 * exp(fac2)
            if i < 400:
                esupc[0, 0, 0][i] = esh40
            else:
                table1[0, 0, 0][i + 1200] = esh40
            i = i + 1

        # -----------------------------------------------------------------------
        # derive blended es over ice and supercooled water between - 40 deg c and 0 deg c
        # -----------------------------------------------------------------------

        i = 0
        while i < 400:
            tem = 233.16 + delt * i
            # GEOS ! WMP impose CALIPSO ice polynomial from 0 C to -40 C
            wice = ice_fraction(tem, 0.0, 0.0)
            wh2o = 1.0 - wice
            table1[0, 0, 0][i + 1200] = (
                wice * table1[0, 0, 0][i + 1200] + wh2o * esupc[0, 0, 0][i]
            )
            i = i + 1


def qs_table_2(length: Int, table2: _FloatField_data_dim):
    """
    compute saturation water vapor pressure table 2
    one phase table

    reference Fortran: gfdl_cloud_microphys.F90: subroutine qs_tablew
    """
    with computation(PARALLEL), interval(...):
        delt = 0.1
        tmin = driver_constants.table_ice - 160.0

        i = 0
        while i < length:
            tem = tmin + delt * i
            fac0 = (tem - driver_constants.t_ice) / (tem * driver_constants.t_ice)
            fac1 = fac0 * driver_constants.lv0
            fac2 = (
                driver_constants.dc_vap * log(tem / driver_constants.t_ice) + fac1
            ) / driver_constants.rvgas
            table2[0, 0, 0][i] = driver_constants.e00 * exp(fac2)
            i = i + 1


def qs_table_3(length: Int, table3: _FloatField_data_dim, table1: _FloatField_data_dim):
    """
    compute saturation water vapor pressure table 3
    two phase table

    reference Fortran: gfdl_cloud_microphys.F90: subroutine qs_table2
    """
    with computation(PARALLEL), interval(...):
        delt = 0.1
        tmin = driver_constants.table_ice - 160.0

        i = 0
        while i < length:
            tem0 = tmin + delt * i
            fac0 = (tem0 - driver_constants.t_ice) / (tem0 * driver_constants.t_ice)
            if i < 1600:
                # -----------------------------------------------------------------------
                # compute es over ice between - 160 deg c and 0 deg c.
                # -----------------------------------------------------------------------
                fac1 = fac0 * driver_constants.li2
                fac2 = (
                    driver_constants.d2ice * log(tem0 / driver_constants.t_ice) + fac1
                ) / driver_constants.rvgas
            else:
                # -----------------------------------------------------------------------
                # compute es over water between 0 deg c and 102 deg c.
                # -----------------------------------------------------------------------
                fac1 = fac0 * driver_constants.lv0
                fac2 = (
                    driver_constants.dc_vap * log(tem0 / driver_constants.t_ice) + fac1
                ) / driver_constants.rvgas
            table3[0, 0, 0][i] = driver_constants.e00 * exp(fac2)
            i = i + 1

        # -----------------------------------------------------------------------
        # smoother around 0 deg c
        # -----------------------------------------------------------------------

        i0 = 1599
        i1 = 1600
        tem0 = 0.25 * (
            table3[0, 0, 0][i0 - 1]
            + 2.0 * table1[0, 0, 0][i0]
            + table3[0, 0, 0][i0 + 1]
        )
        tem1 = 0.25 * (
            table3[0, 0, 0][i1 - 1]
            + 2.0 * table1[0, 0, 0][i1]
            + table3[0, 0, 0][i1 + 1]
        )
        table3[0, 0, 0][i0] = tem0
        table3[0, 0, 0][i1] = tem1


def qs_table_4(length: Int, table4: _FloatField_data_dim, table1: _FloatField_data_dim):
    """
    compute saturation water vapor pressure table 4
    two phase table with " - 2 c" as the transition point

    reference Fortran: gfdl_cloud_microphys.F90: subroutine qs_table3
    """
    with computation(PARALLEL), interval(...):
        delt = 0.1
        esbasw = 1013246.0
        tbasw = driver_constants.table_ice + 100.0
        esbasi = 6107.1
        tmin = driver_constants.table_ice - 160.0

        i = 0
        while i < length:
            tem = tmin + delt * i
            if i < 1580:  # change to - 2 c
                # -----------------------------------------------------------------------
                # compute es over ice between - 160 deg c and 0 deg c.
                # see smithsonian meteorological tables page 350.
                # -----------------------------------------------------------------------
                aa = -9.09718 * (driver_constants.table_ice / tem - 1.0)
                b = -3.56654 * log10(driver_constants.table_ice / tem)
                c = 0.876793 * (1.0 - tem / driver_constants.table_ice)
                e = log10(esbasi)
                table4[0, 0, 0][i] = 0.1 * 10 ** (aa + b + c + e)
            else:
                # -----------------------------------------------------------------------
                # compute es over water between - 2 deg c and 102 deg c.
                # see smithsonian meteorological tables page 350.
                # -----------------------------------------------------------------------
                aa = -7.90298 * (tbasw / tem - 1.0)
                b = 5.02808 * log10(tbasw / tem)
                c = -1.3816e-7 * (10 ** ((1.0 - tem / tbasw) * 11.344) - 1.0)
                d = 8.1328e-3 * (10 ** ((tbasw / tem - 1.0) * (-3.49149)) - 1.0)
                e = log10(esbasw)
                table4[0, 0, 0][i] = 0.1 * 10 ** (aa + b + c + d + e)
            i = i + 1

        # -----------------------------------------------------------------------
        # smoother around - 2 deg c
        # -----------------------------------------------------------------------

        i0 = 1579
        i1 = 1580
        tem0 = 0.25 * (
            table4[0, 0, 0][i0 - 1]
            + 2.0 * table1[0, 0, 0][i0]
            + table4[0, 0, 0][i0 + 1]
        )
        tem1 = 0.25 * (
            table4[0, 0, 0][i1 - 1]
            + 2.0 * table1[0, 0, 0][i1]
            + table4[0, 0, 0][i1 + 1]
        )
        table4[0, 0, 0][i0] = tem0
        table4[0, 0, 0][i1] = tem1


def des_tables(
    length: Int,
    des1: _FloatField_data_dim,
    des2: _FloatField_data_dim,
    des3: _FloatField_data_dim,
    des4: _FloatField_data_dim,
    table1: _FloatField_data_dim,
    table2: _FloatField_data_dim,
    table3: _FloatField_data_dim,
    table4: _FloatField_data_dim,
):
    with computation(PARALLEL), interval(...):
        i = 0
        while i < length - 1:
            des1[0, 0, 0][i] = max(0.0, table1[0, 0, 0][i + 1] - table1[0, 0, 0][i])
            des2[0, 0, 0][i] = max(0.0, table2[0, 0, 0][i + 1] - table2[0, 0, 0][i])
            des3[0, 0, 0][i] = max(0.0, table3[0, 0, 0][i + 1] - table3[0, 0, 0][i])
            des4[0, 0, 0][i] = max(0.0, table4[0, 0, 0][i + 1] - table4[0, 0, 0][i])
            i = i + 1

        des1[0, 0, 0][length - 1] = des1[0, 0, 0][length - 2]
        des2[0, 0, 0][length - 1] = des2[0, 0, 0][length - 2]
        des3[0, 0, 0][length - 1] = des3[0, 0, 0][length - 2]
        des4[0, 0, 0][length - 1] = des4[0, 0, 0][length - 2]


class GFDL_driver_tables:
    """
    Initializes lookup tables for saturation water vapor pressure for the utility routines
    that are designed to return qs consistent with the assumptions in FV3.

    Reference Fortran: gfdl_cloud_microphys.F90: qsmith_init.py
    """

    def __init__(self):
        # NOTE: add backend as an input to make sure the same backend is always being used,
        # and a warning that explains why we are using gt:cpu_ifirst for now

        qsat_domain = (1, 1, 1)

        stencil_factory, quantity_factory_data_dim = get_factories_single_tile(
            qsat_domain[0], qsat_domain[1], qsat_domain[2], 0, backend="gt:cpu_ifirst"
        )
        quantity_factory_data_dim.set_extra_dim_lengths(
            **{
                "table_axis": length,
            }
        )

        self._table1 = quantity_factory_data_dim.zeros(
            [X_DIM, Y_DIM, Z_DIM, "table_axis"], "n/a"
        )
        self._table2 = quantity_factory_data_dim.zeros(
            [X_DIM, Y_DIM, Z_DIM, "table_axis"], "n/a"
        )
        self._table3 = quantity_factory_data_dim.zeros(
            [X_DIM, Y_DIM, Z_DIM, "table_axis"], "n/a"
        )
        self._table4 = quantity_factory_data_dim.zeros(
            [X_DIM, Y_DIM, Z_DIM, "table_axis"], "n/a"
        )
        self._des1 = quantity_factory_data_dim.zeros(
            [X_DIM, Y_DIM, Z_DIM, "table_axis"], "n/a"
        )
        self._des2 = quantity_factory_data_dim.zeros(
            [X_DIM, Y_DIM, Z_DIM, "table_axis"], "n/a"
        )
        self._des3 = quantity_factory_data_dim.zeros(
            [X_DIM, Y_DIM, Z_DIM, "table_axis"], "n/a"
        )
        self._des4 = quantity_factory_data_dim.zeros(
            [X_DIM, Y_DIM, Z_DIM, "table_axis"], "n/a"
        )
        self._esupc = quantity_factory_data_dim.zeros(
            [X_DIM, Y_DIM, Z_DIM, "table_axis"], "n/a"
        )

        compute_qs_table_1 = stencil_factory.from_origin_domain(
            func=qs_table_1,
            origin=(0, 0, 0),
            domain=qsat_domain,
        )
        compute_qs_table_2 = stencil_factory.from_origin_domain(
            func=qs_table_2,
            origin=(0, 0, 0),
            domain=qsat_domain,
        )
        compute_qs_table_3 = stencil_factory.from_origin_domain(
            func=qs_table_3,
            origin=(0, 0, 0),
            domain=qsat_domain,
        )
        compute_qs_table_4 = stencil_factory.from_origin_domain(
            func=qs_table_4,
            origin=(0, 0, 0),
            domain=qsat_domain,
        )
        compute_des_tables = stencil_factory.from_origin_domain(
            func=des_tables,
            origin=(0, 0, 0),
            domain=qsat_domain,
        )

        compute_qs_table_1(length, self._table1, self._esupc)
        compute_qs_table_2(length, self._table2)
        compute_qs_table_3(length, self._table3, self._table1)
        compute_qs_table_4(length, self._table4, self._table1)
        compute_des_tables(
            length,
            self._des1,
            self._des2,
            self._des3,
            self._des4,
            self._table1,
            self._table2,
            self._table3,
            self._table4,
        )

    @property
    def table1(
        self,
    ):
        return self._table1.view[0, 0, 0, :]

    @property
    def table2(
        self,
    ):
        return self._table2.view[0, 0, 0, :]

    @property
    def table3(
        self,
    ):
        return self._table3.view[0, 0, 0, :]

    @property
    def table4(
        self,
    ):
        return self._table4.view[0, 0, 0, :]

    @property
    def des1(
        self,
    ):
        return self._des1.view[0, 0, 0, :]

    @property
    def des2(
        self,
    ):
        return self._des2.view[0, 0, 0, :]

    @property
    def des3(
        self,
    ):
        return self._des3.view[0, 0, 0, :]

    @property
    def des4(
        self,
    ):
        return self._des4.view[0, 0, 0, :]

    @staticmethod
    def make_data_dim_quantity_factory(quantity_factory: QuantityFactory, length):
        data_dim_quantity_factory = copy.deepcopy(quantity_factory)
        data_dim_quantity_factory.set_extra_dim_lengths(
            **{
                "table_axis": length,
            }
        )
        return data_dim_quantity_factory


# Table needs to be calculated only once
_cached_table = {
    "driver_qsat": None,
}


def get_tables():
    if _cached_table["driver_qsat"] is None:
        _cached_table["driver_qsat"] = GFDL_driver_tables()

    return _cached_table["driver_qsat"]