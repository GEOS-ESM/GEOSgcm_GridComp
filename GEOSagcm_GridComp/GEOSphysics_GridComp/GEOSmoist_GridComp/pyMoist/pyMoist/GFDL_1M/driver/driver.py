"""GFDL_1M driver"""

from ndsl import QuantityFactory, StencilFactory, orchestrate
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ
from pyMoist.GFDL_1M.driver.config import config
from pyMoist.GFDL_1M.driver.ice_cloud.main import IceCloud
from pyMoist.GFDL_1M.driver.sat_tables import get_tables
from pyMoist.GFDL_1M.driver.support import (
    check_flags,
    ConfigConstants,
)
from pyMoist.GFDL_1M.driver.temporaries import Temporaries
from pyMoist.GFDL_1M.driver.outputs import Outputs
from pyMoist.GFDL_1M.driver.masks import Masks
from pyMoist.GFDL_1M.driver.terminal_fall.main import TerminalFall
from pyMoist.GFDL_1M.driver.warm_rain.main import WarmRain
from pyMoist.GFDL_1M.driver.fall_speed.main import FallSpeed
from pyMoist.GFDL_1M.driver.setup.main import Setup
from pyMoist.GFDL_1M.driver.finish.main import Finish


class driver:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        GFDL_1M_config: config,
    ):

        self.config_dependent_constants = ConfigConstants(GFDL_1M_config)

        # Check values for untested code paths
        check_flags(
            GFDL_1M_config,
            self.config_dependent_constants.DTS,
        )

        # -----------------------------------------------------------------------
        # initialize precipitation outputs
        # -----------------------------------------------------------------------

        self.outputs = Outputs(quantity_factory)

        # -----------------------------------------------------------------------
        # initialize temporaries
        # -----------------------------------------------------------------------

        self.temporaries = Temporaries(quantity_factory)

        # -----------------------------------------------------------------------
        # initialize masks
        # -----------------------------------------------------------------------

        self.masks = Masks(quantity_factory)

        # -----------------------------------------------------------------------
        # generate saturation specific humidity tables
        # -----------------------------------------------------------------------

        self.sat_tables = get_tables(stencil_factory.backend)

        # -----------------------------------------------------------------------
        # initialize stencils
        # -----------------------------------------------------------------------

        orchestrate(obj=self, config=stencil_factory.config.dace_config)

        self._setup = Setup(
            stencil_factory,
            quantity_factory,
            GFDL_1M_config,
            self.config_dependent_constants,
        )

        self._FallSpeed = FallSpeed(
            stencil_factory,
            quantity_factory,
            GFDL_1M_config,
            self.config_dependent_constants,
        )

        self._TerminalFall = TerminalFall(
            stencil_factory,
            quantity_factory,
            GFDL_1M_config,
            self.config_dependent_constants,
        )

        self._WarmRain = WarmRain(
            stencil_factory,
            quantity_factory,
            GFDL_1M_config,
            self.config_dependent_constants,
        )

        self._IceCloud = IceCloud(
            stencil_factory,
            quantity_factory,
            GFDL_1M_config,
            self.config_dependent_constants,
        )

        self._Finish = Finish(
            stencil_factory,
            quantity_factory,
            GFDL_1M_config,
            self.config_dependent_constants,
        )

    def __call__(
        self,
        GFDL_1M_config: config,
        qv: FloatField,
        ql: FloatField,
        qr: FloatField,
        qi: FloatField,
        qs: FloatField,
        qg: FloatField,
        qa: FloatField,
        qn: FloatField,  # NACTL + NACTI
        qv_dt: FloatField,
        ql_dt: FloatField,
        qr_dt: FloatField,
        qi_dt: FloatField,
        qs_dt: FloatField,
        qg_dt: FloatField,
        qa_dt: FloatField,
        t_dt: FloatField,
        t: FloatField,
        w: FloatField,
        uin: FloatField,
        vin: FloatField,
        udt: FloatField,
        vdt: FloatField,
        dz: FloatField,
        dp: FloatField,
        area: FloatFieldIJ,
        fr_land: FloatFieldIJ,
        cnv_frc: FloatFieldIJ,
        srf_type: FloatFieldIJ,
        eis: FloatFieldIJ,
        rhcrit3d: FloatField,
        anv_icefall: Float,
        ls_icefall: Float,
    ):
        # The driver modifies a number of variables (t, p, qX) but does not pass
        # the changes back to the rest of the model. To replicate this behavior,
        # temporary copies of these variables are used throughout the driver.
        self._setup(
            t,
            dp,
            rhcrit3d,
            qv,
            ql,
            qi,
            qr,
            qs,
            qg,
            qa,
            qn,
            self.temporaries.qv0,
            self.temporaries.ql0,
            self.temporaries.qr0,
            self.temporaries.qi0,
            self.temporaries.qs0,
            self.temporaries.qg0,
            self.temporaries.qa0,
            self.temporaries.qv1,
            self.temporaries.ql1,
            self.temporaries.qr1,
            self.temporaries.qi1,
            self.temporaries.qs1,
            self.temporaries.qg1,
            self.temporaries.qa1,
            dz,
            uin,
            vin,
            w,
            area,
            self.temporaries.t1,
            self.temporaries.dp1,
            self.temporaries.omq,
            self.temporaries.den,
            self.temporaries.p_dry,
            self.temporaries.m1,
            self.temporaries.u1,
            self.temporaries.v1,
            self.temporaries.w1,
            self.temporaries.onemsig,
            self.temporaries.ccn,
            self.temporaries.c_praut,
            self.temporaries.rh_limited,
            self.outputs.rain,
            self.outputs.snow,
            self.outputs.graupel,
            self.outputs.ice,
            self.outputs.m2_rain,
            self.outputs.m2_sol,
            self.outputs.revap,
            self.outputs.isubl,
        )

        for n in range(self.config_dependent_constants.NTIMES):
            self._FallSpeed(
                self.temporaries.ql1,
                self.temporaries.qi1,
                self.temporaries.qs1,
                self.temporaries.qg1,
                t,
                self.temporaries.t1,
                dz,
                self.temporaries.dz1,
                self.temporaries.den,
                self.temporaries.den1,
                self.temporaries.denfac,
                self.temporaries.p_dry,
                self.temporaries.vti,
                self.temporaries.vts,
                self.temporaries.vtg,
                cnv_frc,
                anv_icefall,
                ls_icefall,
            )

            self._TerminalFall(
                self.temporaries.t1,
                self.temporaries.qv1,
                self.temporaries.ql1,
                self.temporaries.qr1,
                self.temporaries.qg1,
                self.temporaries.qs1,
                self.temporaries.qi1,
                self.temporaries.dz1,
                self.temporaries.dp1,
                self.temporaries.den1,
                self.temporaries.vtg,
                self.temporaries.vts,
                self.temporaries.vti,
                self.outputs.rain,
                self.outputs.graupel,
                self.outputs.snow,
                self.outputs.ice,
                self.temporaries.rain1,
                self.temporaries.graupel1,
                self.temporaries.snow1,
                self.temporaries.ice1,
                self.temporaries.m1_sol,
                self.temporaries.w1,
                self.temporaries.ze,
                self.temporaries.zt,
                self.masks.is_frozen,
                self.masks.precip_fall,
            )

            self._WarmRain(
                self.temporaries.dz1,
                self.temporaries.t1,
                self.temporaries.qv1,
                self.temporaries.dp1,
                self.temporaries.ql1,
                self.temporaries.qr1,
                self.temporaries.qi1,
                self.temporaries.qs1,
                self.temporaries.qg1,
                self.temporaries.qa1,
                self.temporaries.ccn,
                self.temporaries.den,
                self.temporaries.denfac,
                self.temporaries.c_praut,
                self.temporaries.vtr,
                self.temporaries.evap1,
                self.temporaries.m1_rain,
                self.temporaries.w1,
                self.temporaries.rh_limited,
                eis,
                self.temporaries.onemsig,
                self.temporaries.rain1,
                self.temporaries.ze,
                self.temporaries.zt,
                self.temporaries.m1,
                self.temporaries.m1_sol,
                self.outputs.rain,
                self.outputs.revap,
                self.outputs.m2_rain,
                self.outputs.m2_sol,
                self.masks.precip_fall,
                self.sat_tables.table1,
                self.sat_tables.table2,
                self.sat_tables.table3,
                self.sat_tables.table4,
                self.sat_tables.des1,
                self.sat_tables.des2,
                self.sat_tables.des3,
                self.sat_tables.des4,
            )

            self._IceCloud(
                self.temporaries.t1,
                self.temporaries.p_dry,
                self.temporaries.dp1,
                self.temporaries.qv1,
                self.temporaries.ql1,
                self.temporaries.qr1,
                self.temporaries.qi1,
                self.temporaries.qs1,
                self.temporaries.qg1,
                self.temporaries.qa1,
                self.temporaries.den1,
                self.temporaries.denfac,
                self.temporaries.vts,
                self.temporaries.vtg,
                self.temporaries.vtr,
                self.temporaries.subl1,
                self.outputs.isubl,
                self.temporaries.rh_limited,
                self.temporaries.ccn,
                cnv_frc,
                srf_type,
                self.sat_tables.table1,
                self.sat_tables.table2,
                self.sat_tables.table3,
                self.sat_tables.table4,
                self.sat_tables.des1,
                self.sat_tables.des2,
                self.sat_tables.des3,
                self.sat_tables.des4,
            )

        self._Finish(
            self.temporaries.qv0,
            self.temporaries.ql0,
            self.temporaries.qr0,
            self.temporaries.qi0,
            self.temporaries.qs0,
            self.temporaries.qg0,
            self.temporaries.qa0,
            self.temporaries.qv1,
            self.temporaries.ql1,
            self.temporaries.qr1,
            self.temporaries.qi1,
            self.temporaries.qs1,
            self.temporaries.qg1,
            self.temporaries.qa1,
            self.temporaries.ccn,
            qv_dt,
            ql_dt,
            qr_dt,
            qi_dt,
            qs_dt,
            qg_dt,
            qa_dt,
            t,
            self.temporaries.t1,
            t_dt,
            w,
            self.temporaries.w1,
            uin,
            self.temporaries.u1,
            udt,
            vin,
            self.temporaries.v1,
            vdt,
            dz,
            dp,
            self.temporaries.dp1,
            self.temporaries.den,
            self.temporaries.p_dry,
            area,
            GFDL_1M_config.DT_MOIST,
            fr_land,
            cnv_frc,
            srf_type,
            eis,
            self.temporaries.rh_limited,
            self.temporaries.m1,
            anv_icefall,
            ls_icefall,
            self.outputs.revap,
            self.outputs.isubl,
            self.outputs.rain,
            self.outputs.snow,
            self.outputs.ice,
            self.outputs.graupel,
            self.outputs.m2_rain,
            self.outputs.m2_sol,
        )
