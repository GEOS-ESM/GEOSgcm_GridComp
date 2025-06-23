"""GFDL_1M driver"""

from ndsl import QuantityFactory, StencilFactory, orchestrate
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ
from pyMoist.GFDL_1M.driver.check_flags import check_flags
from pyMoist.GFDL_1M.config import GFDL1MConfig
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
        GFDL_1M_config: GFDL1MConfig,
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
        t: FloatField,
        w: FloatField,
        u: FloatField,
        v: FloatField,
        dz: FloatField,
        dp: FloatField,
        area: FloatFieldIJ,
        land_fraction: FloatFieldIJ,
        convection_fraction: FloatFieldIJ,
        surface_type: FloatFieldIJ,
        estimated_inversion_strength: FloatFieldIJ,
        rh_crit: FloatField,
        vapor: FloatField,
        liquid: FloatField,
        rain: FloatField,
        ice: FloatField,
        snow: FloatField,
        graupel: FloatField,
        cloud_fraction: FloatField,
        ice_concentration: FloatField,
        liquid_concentration: FloatField,
        dvapor_dt: FloatField,
        dliquid_dt: FloatField,
        drain_dt: FloatField,
        dice_dt: FloatField,
        dsnow_dt: FloatField,
        dgraupel_dt: FloatField,
        dcloud_fraction_dt: FloatField,
        dt_dt: FloatField,
        du_dt: FloatField,
        dv_dt: FloatField,
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
            t (in): atmospheric temperature (K)
            w (in): vertical velocity (m/s)
            u (in): eastward winds (m/s)
            v (in): northward winds (m/s)
            dz (in): layer thickness (m)
            dp (in): change in pressure between model levels (mb)
            area (in): grid cell area
            land_fraction (in): land fraction
            convection_fraction (in): convection fraction
            surface_type (in): surface type
            estimated_inversion_strength (in): estimated inversion strength
            rh_crit (in): critical relative humidity for pdf
            vapor: (in): water vapor mixing ratio (kg/kg)
            liquid (in): in cloud liquid mixing radio (kg/kg)
            rain (in): falling rain (kg/kg)
            ice (in): in cloud ice mixing radio (kg/kg)
            snow (in): in cloud snow mixing radio (kg/kg)
            graupel (in): in cloud graupel mixing radio (kg/kg)
            cloud_fraction (in): cloud fraction (convective + large scale)
            ice_concentration (in): ice concentration (m^-3)
            liquid_concentration (in): liquid concentration (m^-3)
            dvapor_dt (out): water vapor tendency
            dliquid_dt (out): in cloud liquid water tendency
            drain_dt (out): falling rain tendency
            dice_dt (out): in cloud frozen water tendency
            dsnow_dt (out): in cloud snow tendency
            dgraupel_dt (out): in cloud graupel tendency
            dcloud_fraction_dt (out): cloud fraction (convective + large scale) tendency
            dt_dt (out): atmospheric temperature tendency
            du_dt (out): eastward wind tendency
            dv_dt (out): northward wind tendency
            anv_icefall (in): internal parameter related to convective cloud icefall, details unknown
            ls_icefall (in): internal parameter related to large scale cloud icefall, details unknown
        """
        self._setup(
            t,
            dp,
            rh_crit,
            vapor,
            liquid,
            ice,
            rain,
            snow,
            graupel,
            cloud_fraction,
            ice_concentration,
            liquid_concentration,
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
                convection_fraction,
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
                estimated_inversion_strength,
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
                convection_fraction,
                surface_type,
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
            dvapor_dt,
            dliquid_dt,
            drain_dt,
            dice_dt,
            dsnow_dt,
            dgraupel_dt,
            dcloud_fraction_dt,
            t,
            self.temporaries.t1,
            dt_dt,
            w,
            self.temporaries.w1,
            u,
            self.temporaries.u1,
            du_dt,
            v,
            self.temporaries.v1,
            dv_dt,
            dp,
            self.temporaries.dp1,
            self.temporaries.m1,
            self.outputs.rain,
            self.outputs.snow,
            self.outputs.ice,
            self.outputs.graupel,
        )
