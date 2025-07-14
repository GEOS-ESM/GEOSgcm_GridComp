from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.GFDL_1M.driver.driver import MicrophysicsDriver
from pyMoist.GFDL_1M.finalize import Finalize
from pyMoist.GFDL_1M.masks import Masks
from pyMoist.GFDL_1M.outputs import Outputs
from pyMoist.GFDL_1M.PhaseChange.phase_change import PhaseChange
from pyMoist.GFDL_1M.setup import Setup
from pyMoist.GFDL_1M.state import CloudFractions, MixingRatios, VerticalMotion
from pyMoist.GFDL_1M.stencils import (
    _prepare_radiation,
    prepare_tendencies,
    _update_radiation,
    update_tendencies,
)
from pyMoist.GFDL_1M.temporaries import Temporaries
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable


class GFDL1M:
    """
    GFDL Single Moment microphysics

    The primary purpose of this code is to compute macro/microphysical tendencies to be applied to state
    variables (p, t, wind, etc.). This code requires all fields to be preloaded with Fortran memory or
    otherwise supplied between the __init__ and __call__ steps.

    Performs the following functions to achieve this goal:
    __init__
        - initalize saturaiton vapor pressure tables, initalize temporary/output fields, construct stencils
        Arguments: StencilFactory, QuantityFactory, GFDL1MConfig

    __call__
        - setup: compute additional required fields, create pristine copies of input variables
        - phase_change: create new condensates, perform phase change operations
        - driver: precipitate condensates
        - finalize: compute tendencies, prepare fields to be returned to the larger model
        Arguments: none (data needs to be pre-loaded)
    """

    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        GFDL_1M_config: GFDL1MConfig,
    ):
        self.stencil_factory = stencil_factory
        if self.stencil_factory.grid_indexing.n_halo != 0:
            raise ValueError(
                "halo needs to be zero for GFDL Single Moment microphysics"
            )
        self.quantity_factory = quantity_factory
        self.GFDL_1M_config = GFDL_1M_config

        # Initalize saturation tables
        self.saturation_tables = SaturationVaporPressureTable(
            self.stencil_factory.backend
        )

        # Initalize internal fields
        self.masks = Masks.make(quantity_factory=quantity_factory)
        self.temporaries = Temporaries.make(quantity_factory=quantity_factory)

        # Construct stencils

        self.prepare_tendencies = stencil_factory.from_dims_halo(
            func=prepare_tendencies,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self.setup = Setup(
            stencil_factory=stencil_factory,
            GFDL_1M_config=self.GFDL_1M_config,
            saturation_tables=self.saturation_tables,
            prepare_tendencies=self.prepare_tendencies,
        )

        self.phase_change = PhaseChange(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            GFDL_1M_config=self.GFDL_1M_config,
        )

        self.update_tendencies = stencil_factory.from_dims_halo(
            func=update_tendencies,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "DT_MOIST": self.GFDL_1M_config.DT_MOIST,
            },
        )

        self.prepare_radiation = stencil_factory.from_dims_halo(
            func=_prepare_radiation,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self.driver = MicrophysicsDriver(
            self.stencil_factory,
            self.quantity_factory,
            self.GFDL_1M_config,
        )

        self.update_radiation = stencil_factory.from_dims_halo(
            func=_update_radiation,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "DT_MOIST": GFDL_1M_config.DT_MOIST,
            },
        )

        self.finalize = Finalize(
            stencil_factory=stencil_factory,
            GFDL_1M_config=self.GFDL_1M_config,
            saturation_tables=self.saturation_tables,
            update_tendencies=self.update_tendencies,
        )

    def make_ouputs(self) -> Outputs:
        """Helper function to allocate empty outputs"""
        return Outputs.zeros(quantity_factory=self.quantity_factory)

    def __call__(
        self,
        geopotential_height_interface,
        p_interface,
        t,
        u,
        v,
        shallow_convective_rain,
        shallow_convective_snow,
        mixing_ratios: MixingRatios,
        cloud_fractions: CloudFractions,
        area,
        convection_fraction,
        surface_type,
        liquid_concentration,
        ice_concentration,
        vertical_motion: VerticalMotion,
        land_fraction,
        outputs: Outputs,
        evapc,
        sublc,
        rh_crit=None,
    ):
        """
        Args:
            (?) geopotential_height_interface: ?
            (?) p_interface: ?
            (?) t: ?
            (?) u: ?
            (?) v: ?
            (?) shallow_convective_rain: ?
            (?) shallow_convective_snow: ?
            (?) mixing_ratios: ?
            (?) cloud_fractions: ?
            (?) area: ?
            (?) convection_fraction: ?
            (?) surface_type: ?
            (?) liquid_concentration: ?
            (?) ice_concentration: ?
            (?) vertical_motion: ?
            (?) land_fraction: ?
            (out) outputs: ?
            (out) evapc: ?
            (out) sublc: ?
            (out | optional) rh_crit: ?
        """
        self.setup(
            geopotential_height_interface=geopotential_height_interface,
            p_interface=p_interface,
            t=t,
            u=u,
            v=v,
            shallow_convective_rain=shallow_convective_rain,
            shallow_convective_snow=shallow_convective_snow,
            mixing_ratios=mixing_ratios,
            cloud_fractions=cloud_fractions,
            masks=self.masks,
            outputs=outputs,
            temporaries=self.temporaries,
        )

        self.phase_change(
            estimated_inversion_strength=outputs.estimated_inversion_strength,
            p_mb=self.temporaries.p_mb,
            k_lcl=self.temporaries.k_lcl,
            p_interface_mb=self.temporaries.p_interface_mb,
            area=area,
            convection_fraction=convection_fraction,
            surface_type=surface_type,
            t=t,
            convective_liquid=mixing_ratios.convective_liquid,
            convective_ice=mixing_ratios.convective_ice,
            large_scale_liquid=mixing_ratios.large_scale_liquid,
            large_scale_ice=mixing_ratios.large_scale_ice,
            vapor=mixing_ratios.vapor,
            large_scale_cloud_fraction=cloud_fractions.large_scale,
            convective_cloud_fraction=cloud_fractions.convective,
            nactl=liquid_concentration,
            nacti=ice_concentration,
            qsat=self.temporaries.qsat,
            rhx=outputs.relative_humidity_after_pdf,
            rh_crit=rh_crit,
            evapc=evapc,
            sublc=sublc,
        )

        self.update_tendencies(
            u=u,
            v=v,
            t=t,
            vapor=mixing_ratios.vapor,
            rain=mixing_ratios.rain,
            snow=mixing_ratios.snow,
            graupel=mixing_ratios.graupel,
            convective_liquid=mixing_ratios.convective_liquid,
            convective_ice=mixing_ratios.convective_ice,
            large_scale_liquid=mixing_ratios.large_scale_liquid,
            large_scale_ice=mixing_ratios.large_scale_ice,
            convective_cloud_fraction=cloud_fractions.convective,
            large_scale_cloud_fraction=cloud_fractions.large_scale,
            du_dt=outputs.du_dt_macro,
            dv_dt=outputs.dv_dt_macro,
            dt_dt=outputs.dt_dt_macro,
            dvapor_dt=outputs.dvapor_dt_macro,
            dliquid_dt=outputs.dliquid_dt_macro,
            dice_dt=outputs.dice_dt_macro,
            dcloud_fraction_dt=outputs.dcloud_fraction_dt_macro,
            drain_dt=outputs.drain_dt_macro,
            dsnow_dt=outputs.dsnow_dt_macro,
            dgraupel_dt=outputs.dgraupel_dt_macro,
        )

        # prepare microphysics tendencies
        self.prepare_tendencies(
            u=u,
            v=v,
            t=t,
            vapor=mixing_ratios.vapor,
            rain=mixing_ratios.rain,
            snow=mixing_ratios.snow,
            graupel=mixing_ratios.graupel,
            convective_liquid=mixing_ratios.convective_liquid,
            convective_ice=mixing_ratios.convective_ice,
            large_scale_liquid=mixing_ratios.large_scale_liquid,
            large_scale_ice=mixing_ratios.large_scale_ice,
            convective_cloud_fraction=cloud_fractions.convective,
            large_scale_cloud_fraction=cloud_fractions.large_scale,
            du_dt=outputs.du_dt_micro,
            dv_dt=outputs.dv_dt_micro,
            dt_dt=outputs.dt_dt_micro,
            dvapor_dt=outputs.dvapor_dt_micro,
            dliquid_dt=outputs.dliquid_dt_micro,
            dice_dt=outputs.dice_dt_micro,
            dcloud_fraction_dt=outputs.dcloud_fraction_dt_micro,
            drain_dt=outputs.drain_dt_micro,
            dsnow_dt=outputs.dsnow_dt_micro,
            dgraupel_dt=outputs.dgraupel_dt_micro,
        )

        self.prepare_radiation(
            convective_cloud_fraction=cloud_fractions.convective,
            large_scale_cloud_fraction=cloud_fractions.large_scale,
            radiation_cloud_fraction=outputs.radiation_cloud_fraction,
            convective_liquid=mixing_ratios.convective_liquid,
            large_scale_liquid=mixing_ratios.large_scale_liquid,
            radiation_liquid=outputs.radiation_liquid,
            convective_ice=mixing_ratios.convective_ice,
            large_scale_ice=mixing_ratios.large_scale_ice,
            radiation_ice=outputs.radiation_ice,
            vapor=mixing_ratios.vapor,
            radiation_vapor=outputs.radiation_vapor,
            rain=mixing_ratios.rain,
            radiation_rain=outputs.radiation_rain,
            snow=mixing_ratios.snow,
            radiation_snow=outputs.radiation_snow,
            graupel=mixing_ratios.graupel,
            radiation_graupel=outputs.radiation_graupel,
        )

        self.driver(
            t=t,
            w=vertical_motion.velocity,
            u=u,
            v=v,
            dz=self.temporaries.layer_thickness_negative,
            dp=self.temporaries.dp,
            area=area,
            land_fraction=land_fraction,
            convection_fraction=convection_fraction,
            surface_type=surface_type,
            estimated_inversion_strength=outputs.estimated_inversion_strength,
            rh_crit=rh_crit,
            vapor=outputs.radiation_vapor,
            liquid=outputs.radiation_liquid,
            rain=outputs.radiation_rain,
            ice=outputs.radiation_ice,
            snow=outputs.radiation_snow,
            graupel=outputs.radiation_graupel,
            cloud_fraction=outputs.radiation_cloud_fraction,
            ice_concentration=ice_concentration,
            liquid_concentration=liquid_concentration,
            dvapor_dt=self.temporaries.dvapor_dt,
            dliquid_dt=self.temporaries.dliquid_dt,
            drain_dt=self.temporaries.drain_dt,
            dice_dt=self.temporaries.dice_dt,
            dsnow_dt=self.temporaries.dsnow_dt,
            dgraupel_dt=self.temporaries.dgraupel_dt,
            dcloud_fraction_dt=self.temporaries.dcloud_fraction_dt,
            dt_dt=self.temporaries.dt_dt,
            du_dt=self.temporaries.du_dt,
            dv_dt=self.temporaries.dv_dt,
        )

        self.update_radiation(
            t=t,
            u=u,
            v=v,
            radiation_cloud_fraction=outputs.radiation_cloud_fraction,
            radiation_ice=outputs.radiation_ice,
            radiation_liquid=outputs.radiation_liquid,
            radiation_vapor=outputs.radiation_vapor,
            radiation_rain=outputs.radiation_rain,
            radiation_snow=outputs.radiation_snow,
            radiation_graupel=outputs.radiation_graupel,
            dcloud_fraction_dt=self.temporaries.dcloud_fraction_dt,
            dt_dt=self.temporaries.dt_dt,
            du_dt=self.temporaries.du_dt,
            dv_dt=self.temporaries.dv_dt,
            dice_dt=self.temporaries.dice_dt,
            dliquid_dt=self.temporaries.dliquid_dt,
            dvapor_dt=self.temporaries.dvapor_dt,
            drain_dt=self.temporaries.drain_dt,
            dsnow_dt=self.temporaries.dsnow_dt,
            dgraupel_dt=self.temporaries.dgraupel_dt,
        )

        self.finalize(
            t=t,
            u=u,
            v=v,
            ice_concentration=ice_concentration,
            liquid_concentration=liquid_concentration,
            mixing_ratios=mixing_ratios,
            cloud_fractions=cloud_fractions,
            masks=self.masks,
            outputs=outputs,
            temporaries=self.temporaries,
            driver=self.driver,
        )
