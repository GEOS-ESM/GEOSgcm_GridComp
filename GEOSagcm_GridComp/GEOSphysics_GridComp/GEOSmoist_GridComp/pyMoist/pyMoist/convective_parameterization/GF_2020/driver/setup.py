import copy
from ndsl import StencilFactory, QuantityFactory
from ndsl.dsl.gt4py import PARALLEL, interval, computation, FORWARD, sqrt, max, min, abs, floor
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField, FloatFieldIJ
from pyMoist.convective_parameterization.GF_2020.config import GF2020Config
import pyMoist.constants as constants
from pyMoist.convective_parameterization.GF_2020.temporaries import Temporaries
from pyMoist.convective_parameterization.GF_2020.state import MixingRatios
from pyMoist.saturation_tables.qsat_functions import saturation_specific_humidity
from pyMoist.field_types import GlobalTable_saturaion_tables
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.convective_parameterization.GF_2020.driver.temporaries import DriverTemporaries
from pyMoist.field_types import FloatField_maxiens, FloatField_nmp
from pyMoist.convective_parameterization.GF_2020.driver.constants import driver_constants


def setup(
    # outputs used later in the driver
    ztexec: FloatFieldIJ,
    zqexec: FloatFieldIJ,
    fixout_qv: FloatFieldIJ,
):
    """ """
    from __externals__ import C1, ADV_TRIGGER, AUTOCONV, DT, USE_TRACER_TRANSP

    with computation(PARALLEL), interval(...):
        if C1 > 0:
            USE_C1D = True
        if ADV_TRIGGER == 2:
            option_not_implemented_placeholder = True

    with computation(FORWARD), interval(0, 1):
        rtgt = 1.0

    with computation(PARALLEL), interval(...):
        if ADV_TRIGGER == 2:
            option_not_implemented_placeholder = True

        if AUTOCONV == 2:
            # TODO is 3d, swap to 2d when 2d temporaries are available
            ccn = max(100.0, (370.37 * (0.01 + max(0.0, aot500))) ** 1.555)
        else:
            ccn = 100.0

        # TODO this block: all are 3d, swap to 2d when 2d temporaries are available
        psur = sfc_press * 1.0e-2  # mbar
        tsur = temp2m
        ter11 = max(0.0, topt)
        xlons = lons * 180.0 / 3.14159
        xlats = lats * 180.0 / 3.14159

    with computation(PARALLEL), interval(0, -1):
        # heigths, current pressure, temp and water vapor mix ratio
        zo = zt * rtgt + topt
        po = press * 1.0e-2  # mbar
        temp_old = temp

        qv_old = rvap  # @ begin of the timestep
        qv_curr = curr_rvap  # current (after dynamics + physical processes called before GF)

        # air density, TKE and cloud liq water mixing ratio
        rhoi = 1.0e2 * po / (287.04 * temp_old * (1.0 + 0.608 * qv_old))
        tkeg = driver_constants.tkmin
        rcpg = 0.0

        # wind velocities
        us = u
        vs = v
        ws = w
        dm2d = dm
        omeg = -driver_constants.g * rhoi * w
        # buoyancy excess
        buoy_exc2d = buoy_exc
        # temp/water vapor modified only by advection
        temp_new_ADV = temp_old + (rth_advten) * DT
        qv_new_ADV = qv_old + (rqvften) * DT

        if APPLY_SUB_MP == 1:
            # microphysics ice and liq mixing ratio, and cloud fraction of the host model
            # (only subsidence is applied)
            mpqi[0, 0, 0][0] = mp_ice[0, 0, 0][0]  # kg/kg
            mpql[0, 0, 0][0] = mp_liq[0, 0, 0][0]  # kg/kg
            mpcf[0, 0, 0][0] = mp_cf[0, 0, 0][0]  # 1
            mpqi[0, 0, 0][1] = mp_ice[0, 0, 0][1]  # kg/kg
            mpql[0, 0, 0][1] = mp_liq[0, 0, 0][1]  # kg/kg
            mpcf[0, 0, 0][1] = mp_cf[0, 0, 0][1]  # 1

        if USE_TRACER_TRANSP == 1:
            option_not_implemented_placeholder = True

    with computation(PARALLEL), interval(...):
        # TODO is 3d, swap to 2d when 2d temporaries are available
        pbl = zo.at(K=kpbl) - topt

        # get execess T and Q for source air parcels
        pten = temp_old.at(K=0)
        pqen = qv_old.at(K=0)
        paph = 100.0 * psur
        zrho = paph / (287.04 * (temp_old(i, 1) * (1.0 + 0.608 * qv_old(i, 1))))
        # sensible and latent sfc fluxes for the heat-engine closure
        # TODO is 3d, swap to 2d when 2d temporaries are available
        h_sfc_flux = zrho * driver_constants.cp * sflux_t  # W/m^2
        # TODO is 3d, swap to 2d when 2d temporaries are available
        le_sfc_flux = zrho * driver_constants.xlv * sflux_r  # W/m^2
        # local le and h fluxes for W*
        pahfs = -sflux_t * zrho * 1004.64  # W/m^2
        pqhfl = -sflux_r  # kg/m^2/s
        # buoyancy flux (h+le)
        zkhvfl = (pahfs / 1004.64 + 0.608 * pten * pqhfl) / zrho  # K m s-1
        # depth of 1st model layer
        # (zo(1)-top is ~ 1/2 of the depth of 1st model layer, => mult by 2)
        pgeoh = 2.0 * (zo.at(K=0) - topt) * driver_constants.g  # m+2 s-2
        # convective-scale velocity w*
        # in the future, change 0.001 by ustar^3
        # TODO is 3d, swap to 2d when 2d temporaries are available
        zws = max(0.0, 0.001 - 1.5 * 0.41 * zkhvfl * pgeoh / pten)  # m+3 s-3

        if zws > constants.FLOAT_TINY:
            # convective-scale velocity w*
            zws = 1.2 * zws**0.3333
            # temperature excess
            ztexec = max(0.0, -1.5 * pahfs / (zrho * zws * 1004.64))  # K
            # moisture  excess
            zqexec = max(0.0, -1.5 * pqhfl / (zrho * zws))  # kg kg-1
        # zws for shallow convection closure (Grant 2001)
        # depth of the pbl
        pgeoh = pbl * driver_constants.g
        # convective-scale velocity W* (m/s)
        zws = max(0.0, 0.001 - 1.5 * 0.41 * zkhvfl * pgeoh / pten)
        zws = 1.2 * zws**0.3333


class DriverSetup:
    def __init__(self, stencil_factory: StencilFactory, GF_2020_config: GF2020Config):
        # Construct stencils
        self.setup = stencil_factory.from_dims_halo(
            func=setup,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "C1": GF_2020_config.C1,
                "ADV_TRIGGER": GF_2020_config.ADV_TRIGGER,
                "AUTOCONV": GF_2020_config.AUTOCONV,
                "USE_TRACER_TRANSP": GF_2020_config.USE_TRACER_TRANSP,
            },
        )

    def __call__(self, quantity_factory: QuantityFactory):
        self.driver_temporaries = DriverTemporaries.make(quantity_factory)
