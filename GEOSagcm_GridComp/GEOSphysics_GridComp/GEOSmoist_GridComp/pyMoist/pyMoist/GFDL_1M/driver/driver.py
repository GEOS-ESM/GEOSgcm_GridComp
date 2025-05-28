"""GFDL_1M driver"""

from ndsl import QuantityFactory, StencilFactory, orchestrate
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ
from pyMoist.GFDL_1M.driver.check_flags import check_flags
from pyMoist.GFDL_1M.driver.config import MicrophysicsConfiguration
from pyMoist.GFDL_1M.driver.config_constants import ConfigConstants
from pyMoist.GFDL_1M.driver.fall_speed.main import FallSpeed
from pyMoist.GFDL_1M.driver.finish.main import Finish
from pyMoist.GFDL_1M.driver.ice_cloud.main import IceCloud
from pyMoist.GFDL_1M.driver.masks import Masks
from pyMoist.GFDL_1M.driver.outputs import Outputs
from pyMoist.GFDL_1M.driver.sat_tables import get_tables
from pyMoist.GFDL_1M.driver.setup.main import Setup
from pyMoist.GFDL_1M.driver.temporaries import Temporaries
from pyMoist.GFDL_1M.driver.terminal_fall.main import TerminalFall
from pyMoist.GFDL_1M.driver.warm_rain.main import WarmRain


class MicrophysicsDriver:
    """
    This class contains the GFDL single moment microphysics driver. The driver is broken
    into six components: Setup, FallSpeed, TerminalFall, WarmRain, IceCloud, and Finish.

    __init__:
        - checks validity of constants and trigger parameters for unimplemented options
        - initalizes internal fields
        - constructs stencils
        Arguments: StencilFactory, QuantityFactory, MicrophysicsConfiguration

    __call__:
        - evaluates stencils
        Arguments: various state fields (pressure, temperature, wind) and mixing ratios
    """

    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        GFDL_1M_config: MicrophysicsConfiguration,
    ):
        """
        Perform setup for the microphysics driver. Check flags for unimplemented options,
        initalize internal fields, and compile stencils.

        Arguments:
            stencil_factory: StencilFactory with model domain information
            quantity_factory: QuantityFactory with model domain information
            GFDL_1M_config: driver configuration
        """
        self.config_dependent_constants = ConfigConstants.make(GFDL_1M_config)

        # Check values for untested code paths
        check_flags(
            GFDL_1M_config,
            self.config_dependent_constants.DTS,
        )

        # -----------------------------------------------------------------------
        # initialize precipitation outputs
        # -----------------------------------------------------------------------

        self.outputs = Outputs.make(quantity_factory)

        # -----------------------------------------------------------------------
        # initialize temporaries
        # -----------------------------------------------------------------------

        self.temporaries = Temporaries.make(quantity_factory)

        # -----------------------------------------------------------------------
        # initialize masks
        # -----------------------------------------------------------------------

        self.masks = Masks.make(quantity_factory)

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
            GFDL_1M_config,
            self.config_dependent_constants,
        )

        self._fall_speed = FallSpeed(
            stencil_factory,
            quantity_factory,
            GFDL_1M_config,
            self.config_dependent_constants,
        )

        self._terminal_fall = TerminalFall(
            stencil_factory,
            quantity_factory,
            GFDL_1M_config,
            self.config_dependent_constants,
        )

        self._warm_rain = WarmRain(
            stencil_factory,
            quantity_factory,
            GFDL_1M_config,
            self.config_dependent_constants,
        )

        self._ice_cloud = IceCloud(
            stencil_factory,
            GFDL_1M_config,
            self.config_dependent_constants,
        )

        self._finish = Finish(
            stencil_factory,
            GFDL_1M_config,
            self.config_dependent_constants,
        )

    def __call__(
        self,
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
        u: FloatField,
        v: FloatField,
        u_dt: FloatField,
        v_dt: FloatField,
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
        """
        Evaluate the microphysics driver. The driver call is broken into six parts:
            - Setup: fill temporaries, compute required intermediary fields from inputs
            - FallSpeed: compute real fall speed of precipitates
            - TerminalFall: compute terminal fall speed of precipitates
            - WarmRain: warm rain cloud microphysics
            - IceCloud: ice cloud microphysical processes
            - Finish: compute output tendencies

        Arguments:
            GFDL_1M_config (in): driver configuration
            qv: (in): water vapor mixing ratio (kg/kg)
            ql (in): in cloud liquid mixing radio (kg/kg)
            qr (in): falling rain (kg/kg)
            qi (in): in cloud ice mixing radio (kg/kg)
            qs (in): in cloud snow mixing radio (kg/kg)
            qg (in): in cloud graupel mixing radio (kg/kg)
            qa (in): cloud fraction (convective + large scale)
            qn (in): liquid + ice concentration (m^-3)
            qv_dt (out): water vapor tendency
            ql_dt (out): in cloud liquid water tendency
            qr_dt (out): falling rain tendency
            qi_dt (out): in cloud frozen water tendency
            qs_dt (out): in cloud snow tendency
            qg_dt (out): in cloud graupel tendency
            qa_dt (out): cloud fraction (convective + large scale) tendency
            t_dt (out): atmospheric temperature tendency
            t (in): atmospheric temperature (K)
            w (in): vertical velocity (m/s)
            u (in): eastward winds (m/s)
            v (in): northward winds (m/s)
            u_dt (out): eastward wind tendency
            v_dt (out): northward wind tendency
            dz (in): layer thickness (m)
            dp (in): change in pressure between model levels (mb)
            area (in): grid cell area
            fr_land (in): land fraction
            cnv_frc (in): convection fraction
            srf_type (in): surface type
            eis (in): estimated inversion strength
            rhcrit3d (out): details unknown
            anv_icefall (in): internal parameter related to convective cloud icefall
            ls_icefall (in): internal parameter related to large scale cloud icefall
        """
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
            u,
            v,
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
            self._fall_speed(
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

            self._terminal_fall(
                self.temporaries.t1,
                self.temporaries.qv1,
                self.temporaries.ql1,
                self.temporaries.qr1,
                self.temporaries.qg1,
                self.temporaries.qs1,
                self.temporaries.qi1,
                self.temporaries.dz1,
                self.temporaries.dp1,
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

            self._warm_rain(
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

            self._ice_cloud(
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
                self.sat_tables.table2,
                self.sat_tables.table3,
                self.sat_tables.des2,
                self.sat_tables.des3,
            )

        self._finish(
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
            u,
            self.temporaries.u1,
            u_dt,
            v,
            self.temporaries.v1,
            v_dt,
            dp,
            self.temporaries.dp1,
            self.temporaries.m1,
            self.outputs.rain,
            self.outputs.snow,
            self.outputs.ice,
            self.outputs.graupel,
        )
