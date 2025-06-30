from pyMoist.GFDL_1M.driver.driver import MicrophysicsDriver
from pyMoist.GFDL_1M.PhaseChange.phase_change import PhaseChange
from ndsl import QuantityFactory, StencilFactory
from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.GFDL_1M.getters_temporary import (
    mapl_get_resource_placeholder,
    esmf_placeholder,
    mapl_get_pointer_placeholder,
    mapl_get_pointer_placeholder,
    associated_checker,
)
from pyMoist.GFDL_1M.state import (
    LiquidWaterStaticEnergy,
    TotalWater,
    VericalMotion,
    MixingRatios,
    CloudFractions,
)
from pyMoist.GFDL_1M.stencils import (
    prepare_tendencies,
    prepare_radiation_quantities,
    update_tendencies,
    update_radiation_quantities,
)
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.GFDL_1M.outputs import Outputs
from pyMoist.GFDL_1M.temporaries import Temporaries
from pyMoist.GFDL_1M.masks import Masks
from pyMoist.GFDL_1M.finalize import Finalize
from pyMoist.GFDL_1M.setup import Setup
from pyMoist.interface.mapl.memory_factory import MAPLMemoryRepository


class GFDL1M:
    def __init__(
        self, stencil_factory: StencilFactory, quantity_factory: QuantityFactory, GFDL_1M_config: GFDL1MConfig
    ):
        self.stencil_factory = stencil_factory
        if self.stencil_factory.grid_indexing.n_halo != 0:
            raise ValueError("halo needs to be zero for GFDL Single Moment microphysics")
        self.quantity_factory = quantity_factory
        self.GFDL_1M_config = GFDL_1M_config

        # Initalize saturation tables
        self.saturation_tables = SaturationVaporPressureTable(self.stencil_factory.backend)

        # Initalize internal fields
        self.masks = Masks.make(quantity_factory=quantity_factory)
        self.outputs = Outputs.make(quantity_factory=quantity_factory)
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

        self.prepare_radiation_quantities = stencil_factory.from_dims_halo(
            func=prepare_radiation_quantities,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self.driver = MicrophysicsDriver(
            self.stencil_factory,
            self.quantity_factory,
            self.GFDL_1M_config,
        )

        self.update_radiation_quantities = stencil_factory.from_dims_halo(
            func=update_radiation_quantities,
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

    def __call__(self):
        self.setup(
            geopotential_height_interface=self.geopotential_height_interface,
            p_interface=self.p_interface,
            t=self.t,
            u=self.u,
            v=self.v,
            shallow_convective_rain=self.shallow_convective_rain,
            shallow_convective_snow=self.shallow_convective_snow,
            mixing_ratios=self.mixing_ratios,
            cloud_fractions=self.cloud_fractions,
            masks=self.masks,
            outputs=self.outputs,
            temporaries=self.temporaries,
        )

        self.phase_change(
            estimated_inversion_strength=self.outputs.estimated_inversion_strength,
            p_mb=self.temporaries.p_mb,
            k_lcl=self.temporaries.k_lcl,
            p_interface_mb=self.temporaries.p_interface_mb,
            area=self.area,
            convection_fraction=self.convection_fraction,
            surface_type=self.surface_type,
            t=self.t,
            convective_liquid=self.mixing_ratios.convective_liquid,
            convective_ice=self.mixing_ratios.convective_ice,
            large_scale_liquid=self.mixing_ratios.large_scale_liquid,
            large_scale_ice=self.mixing_ratios.large_scale_ice,
            vapor=self.mixing_ratios.vapor,
            large_scale_cloud_fraction=self.cloud_fractions.large_scale,
            convective_cloud_fraction=self.cloud_fractions.convective,
            nactl=self.liquid_concentration,
            nacti=self.ice_concentration,
            qsat=self.temporaries.qsat,
        )

        # pull rh_crit and rhx out of phase change component and send it back to the rest of the model
        rh_crit = self.phase_change.outputs.rh_crit
        self.outputs.relative_humidity_after_pdf = self.phase_change.outputs.rhx

        self.update_tendencies(
            u=self.u,
            v=self.v,
            t=self.t,
            vapor=self.mixing_ratios.vapor,
            rain=self.mixing_ratios.rain,
            snow=self.mixing_ratios.snow,
            graupel=self.mixing_ratios.graupel,
            convective_liquid=self.mixing_ratios.convective_liquid,
            convective_ice=self.mixing_ratios.convective_ice,
            large_scale_liquid=self.mixing_ratios.large_scale_liquid,
            large_scale_ice=self.mixing_ratios.large_scale_ice,
            convective_cloud_fraction=self.cloud_fractions.convective,
            large_scale_cloud_fraction=self.cloud_fractions.large_scale,
            du_dt=self.outputs.du_dt_macro,
            dv_dt=self.outputs.dv_dt_macro,
            dt_dt=self.outputs.dt_dt_macro,
            dvapor_dt=self.outputs.dvapor_dt_macro,
            dliquid_dt=self.outputs.dliquid_dt_macro,
            dice_dt=self.outputs.dice_dt_macro,
            dcloud_fraction_dt=self.outputs.dcloud_fraction_dt_macro,
            drain_dt=self.outputs.drain_dt_macro,
            dsnow_dt=self.outputs.dsnow_dt_macro,
            dgraupel_dt=self.outputs.dgraupel_dt_macro,
        )

        # prepare microphysics tendencies
        self.prepare_tendencies(
            u=self.u,
            v=self.v,
            t=self.t,
            vapor=self.mixing_ratios.vapor,
            rain=self.mixing_ratios.rain,
            snow=self.mixing_ratios.snow,
            graupel=self.mixing_ratios.graupel,
            convective_liquid=self.mixing_ratios.convective_liquid,
            convective_ice=self.mixing_ratios.convective_ice,
            large_scale_liquid=self.mixing_ratios.large_scale_liquid,
            large_scale_ice=self.mixing_ratios.large_scale_ice,
            convective_cloud_fraction=self.cloud_fractions.convective,
            large_scale_cloud_fraction=self.cloud_fractions.large_scale,
            du_dt=self.outputs.du_dt_micro,
            dv_dt=self.outputs.dv_dt_micro,
            dt_dt=self.outputs.dt_dt_micro,
            dvapor_dt=self.outputs.dvapor_dt_micro,
            dliquid_dt=self.outputs.dliquid_dt_micro,
            dice_dt=self.outputs.dice_dt_micro,
            dcloud_fraction_dt=self.outputs.dcloud_fraction_dt_micro,
            drain_dt=self.outputs.drain_dt_micro,
            dsnow_dt=self.outputs.dsnow_dt_micro,
            dgraupel_dt=self.outputs.dgraupel_dt_micro,
        )

        self.prepare_radiation_quantities(
            convective_cloud_fraction=self.cloud_fractions.convective,
            large_scale_cloud_fraction=self.cloud_fractions.large_scale,
            radiation_cloud_fraction=self.outputs.radiation_cloud_fraction,
            convective_liquid=self.mixing_ratios.convective_liquid,
            large_scale_liquid=self.mixing_ratios.large_scale_liquid,
            radiation_liquid=self.outputs.radiation_liquid,
            convective_ice=self.mixing_ratios.convective_ice,
            large_scale_ice=self.mixing_ratios.large_scale_ice,
            radiation_ice=self.outputs.radiation_ice,
            vapor=self.mixing_ratios.vapor,
            radiation_vapor=self.outputs.radiation_vapor,
            rain=self.mixing_ratios.rain,
            radiation_rain=self.outputs.radiation_rain,
            snow=self.mixing_ratios.snow,
            radiation_snow=self.outputs.radiation_snow,
            graupel=self.mixing_ratios.graupel,
            radiation_graupel=self.outputs.radiation_graupel,
        )

        self.driver(
            t=self.t,
            w=self.vertical_motion.velocity,
            u=self.u,
            v=self.v,
            dz=self.temporaries.layer_thickness_negative,
            dp=self.temporaries.dp,
            area=self.area,
            land_fraction=self.land_fraction,
            convection_fraction=self.convection_fraction,
            surface_type=self.surface_type,
            estimated_inversion_strength=self.outputs.estimated_inversion_strength,
            rh_crit=rh_crit,
            vapor=self.outputs.radiation_vapor,
            liquid=self.outputs.radiation_liquid,
            rain=self.outputs.radiation_rain,
            ice=self.outputs.radiation_ice,
            snow=self.outputs.radiation_snow,
            graupel=self.outputs.radiation_graupel,
            cloud_fraction=self.outputs.radiation_cloud_fraction,
            ice_concentration=self.ice_concentration,
            liquid_concentration=self.liquid_concentration,
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
            anv_icefall=self.GFDL_1M_config.ANV_ICEFALL,
            ls_icefall=self.GFDL_1M_config.LS_ICEFALL,
        )

        self.update_radiation_quantities(
            t=self.t,
            u=self.u,
            v=self.v,
            radiation_cloud_fraction=self.outputs.radiation_cloud_fraction,
            radiation_ice=self.outputs.radiation_ice,
            radiation_liquid=self.outputs.radiation_liquid,
            radiation_vapor=self.outputs.radiation_vapor,
            radiation_rain=self.outputs.radiation_rain,
            radiation_snow=self.outputs.radiation_snow,
            radiation_graupel=self.outputs.radiation_graupel,
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
            t=self.t,
            u=self.u,
            v=self.v,
            ice_concentration=self.ice_concentration,
            liquid_concentration=self.liquid_concentration,
            mixing_ratios=self.mixing_ratios,
            cloud_fractions=self.cloud_fractions,
            masks=self.masks,
            outputs=self.outputs,
            temporaries=self.temporaries,
            driver=self.driver,
        )
