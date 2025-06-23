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
from ndsl.dsl.typing import Float
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
            externals={
                "DT_MOIST": self.GFDL_1M_config.DT_MOIST,
            },
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
                "DT_MOIST": GFDL1MConfig.DT_MOIST,
            },
        )

        self.finalize = Finalize(
            stencil_factory=stencil_factory,
            GFDL_1M_config=self.GFDL_1M_config,
            saturation_tables=self.saturation_tables,
            FAC_RL=self.GFDL_1M_config.FAC_RL,
            MIN_RL=self.GFDL_1M_config.MIN_RL,
            MAX_RL=self.GFDL_1M_config.MAX_RL,
            FAC_RI=self.GFDL_1M_config.FAC_RI,
            MIN_RI=self.GFDL_1M_config.MIN_RI,
            MAX_RI=self.GFDL_1M_config.MAX_RI,
            update_tendencies=self.update_tendencies,
        )

    def get_fortran_data(self, *args, **kwds):
        ##### Link FORTRAN memory to Python memory #####
        ##### All links are live, so changes in Python will be reflected in FORTRAN #####
        ##### Not all "input" (model state) fields are modified #####

        # Get model state
        self.mixing_ratios = MixingRatios(
            vapor=mapl_get_pointer_placeholder("Q"),
            rain=mapl_get_pointer_placeholder("QRAIN"),
            snow=mapl_get_pointer_placeholder("QSNOW"),
            graupel=mapl_get_pointer_placeholder("QGRAUPEL"),
            convective_liquid=mapl_get_pointer_placeholder("QLCN"),
            convective_ice=mapl_get_pointer_placeholder("QICN"),
            large_scale_liquid=mapl_get_pointer_placeholder("QLLS"),
            large_scale_ice=mapl_get_pointer_placeholder("QILS"),
        )
        self.cloud_fractions = CloudFractions(
            convective=mapl_get_pointer_placeholder("CLCN"),
            large_scale=mapl_get_pointer_placeholder("CLLS"),
        )
        self.liquid_concentration = mapl_get_pointer_placeholder("NACTL")
        self.ice_concentration = mapl_get_pointer_placeholder("NACTI")
        self.area = mapl_get_pointer_placeholder("AREA")
        self.geopotential_height_interface = mapl_get_pointer_placeholder("ZLE")
        self.p_interface = mapl_get_pointer_placeholder("PLE")
        self.t = mapl_get_pointer_placeholder("T")
        self.u = mapl_get_pointer_placeholder("U")
        self.v = mapl_get_pointer_placeholder("V")
        self.land_fraction = mapl_get_pointer_placeholder("FRLAND")
        self.vertical_motion = VericalMotion(
            velocity=mapl_get_pointer_placeholder("W"),
            variance=mapl_get_pointer_placeholder("W2"),
            third_moment=mapl_get_pointer_placeholder("W3"),
        )
        self.liquid_water_static_energy = LiquidWaterStaticEnergy(
            flux=mapl_get_pointer_placeholder("WSL"),
            variance=mapl_get_pointer_placeholder("SL2"),
            third_moment=mapl_get_pointer_placeholder("SL3"),
        )
        self.total_water = TotalWater(
            flux=mapl_get_pointer_placeholder("WQT"),
            variance=mapl_get_pointer_placeholder("QT2"),
            third_moment=mapl_get_pointer_placeholder("QT3"),
        )
        self.convection_fraction = mapl_get_pointer_placeholder("CNV_FRC")
        self.surface_type = mapl_get_pointer_placeholder("SRF_TYPE")
        self.shallow_convective_rain = mapl_get_pointer_placeholder("SHLW_PRC3")
        self.shallow_convective_snow = mapl_get_pointer_placeholder("SHLW_SNO3")
        self.rh_crit = mapl_get_pointer_placeholder("RHCRIT")

        # Outputs: model fields originating from within GFDL
        self.outputs.liquid_radius = mapl_get_pointer_placeholder("RL")
        self.outputs.ice_radius = mapl_get_pointer_placeholder("RI")
        self.outputs.large_scale_nonanvil_precipitation_evaporation = mapl_get_resource_placeholder("EVAPC")
        self.outputs.large_scale_nonanvil_precipitation_sublimation = mapl_get_resource_placeholder("SUBLC")
        self.outputs.precipitated_rain = mapl_get_resource_placeholder("PRCP_RAIN")
        self.outputs.precipitated_snow = mapl_get_resource_placeholder("PRCP_SNOW")
        self.outputs.precipitated_ice = mapl_get_resource_placeholder("PRCP_ICE")
        self.outputs.precipitated_graupel = mapl_get_resource_placeholder("PRCP_GRAUPEL")

        # Outputs: model fields originating from within GFDL; radiation fields
        self.outputs.radiation_cloud_fraction = mapl_get_pointer_placeholder("FCLD")
        self.outputs.radiation_vapor = mapl_get_pointer_placeholder("QV")
        self.outputs.radiation_liquid = mapl_get_pointer_placeholder("QL")
        self.outputs.radiation_ice = mapl_get_pointer_placeholder("QI")
        self.outputs.radiation_rain = mapl_get_pointer_placeholder("QR")
        self.outputs.radiation_snow = mapl_get_pointer_placeholder("QS")
        self.outputs.radiation_graupel = mapl_get_pointer_placeholder("QG")
        self.outputs.lower_tropospheric_stability = mapl_get_pointer_placeholder("LTS")
        self.outputs.estimated_inversion_strength = mapl_get_pointer_placeholder("EIS")
        self.outputs.z_lcl = mapl_get_resource_placeholder("ZLCL")

        # Outputs: model fields originating from within GFDL; macrophysics/microphysics tendencies
        self.outputs.du_dt_macro = mapl_get_pointer_placeholder("DUDT_macro")
        self.outputs.dv_dt_macro = mapl_get_pointer_placeholder("DVDT_macro")
        self.outputs.dt_dt_macro = mapl_get_pointer_placeholder("DTDT_macro")
        self.outputs.dvapor_dt_macro = mapl_get_pointer_placeholder("DQVDT_macro")
        self.outputs.dliquid_dt_macro = mapl_get_pointer_placeholder("DQLDT_macro")
        self.outputs.dice_dt_macro = mapl_get_pointer_placeholder("DQIDT_macro")
        self.outputs.dcloud_fraction_dt_macro = mapl_get_pointer_placeholder("DQADT_macro")
        self.outputs.drain_dt_macro = mapl_get_pointer_placeholder("DQRDT_macro")
        self.outputs.dsnow_dt_macro = mapl_get_pointer_placeholder("DQSDT_macro")
        self.outputs.dgraupel_dt_macro = mapl_get_pointer_placeholder("DQGDT_macro")
        self.outputs.du_dt_micro = mapl_get_pointer_placeholder("DUDT_micro")
        self.outputs.dv_dt_micro = mapl_get_pointer_placeholder("DVDT_micro")
        self.outputs.dt_dt_micro = mapl_get_pointer_placeholder("DTDT_micro")
        self.outputs.dvapor_dt_micro = mapl_get_pointer_placeholder("DQVDT_micro")
        self.outputs.dliquid_dt_micro = mapl_get_pointer_placeholder("DQLDT_micro")
        self.outputs.dice_dt_micro = mapl_get_pointer_placeholder("DQIDT_micro")
        self.outputs.dcloud_fraction_dt_micro = mapl_get_pointer_placeholder("DQADT_micro")
        self.outputs.drain_dt_micro = mapl_get_pointer_placeholder("DQRDT_micro")
        self.outputs.dsnow_dt_micro = mapl_get_pointer_placeholder("DQSDT_micro")
        self.outputs.dgraupel_dt_micro = mapl_get_pointer_placeholder("DQGDT_micro")
        # Outputs: Exports to be filled
        self.outputs.large_scale_precip = mapl_get_pointer_placeholder("LS_PRCP")
        self.outputs.large_scale_snow = mapl_get_pointer_placeholder("LS_SNR")
        self.outputs.icefall = mapl_get_pointer_placeholder("ICE")
        self.outputs.freezing_rainfall = mapl_get_pointer_placeholder("FRZR")
        self.outputs.relative_humidity_after_pdf = mapl_get_pointer_placeholder("RHX")
        self.outputs.large_scale_nonanvil_precipitation_evaporation = mapl_get_pointer_placeholder("REV_LS")
        self.outputs.large_scale_nonanvil_precipitation_sublimation = mapl_get_pointer_placeholder("RSU_LS")
        self.outputs.large_scale_nonanvil_liquid_flux = mapl_get_pointer_placeholder("PFL_LS")
        self.outputs.large_scale_nonanvil_ice_flux = mapl_get_pointer_placeholder("PFI_LS")
        self.outputs.anvil_liquid_flux = mapl_get_pointer_placeholder("PFL_AN")
        self.outputs.anvil_ice_flux = mapl_get_pointer_placeholder("PFI_AN")
        self.outputs.large_scale_rainwater_source = mapl_get_pointer_placeholder("DQRL")
        self.outputs.simulated_reflectivity = mapl_get_pointer_placeholder("DBZ")
        self.outputs.maximum_reflectivity = mapl_get_pointer_placeholder("DBZ_MAX")
        self.outputs.one_km_agl_reflectivity = mapl_get_pointer_placeholder("DBZ_1KM")
        self.outputs.echo_top_reflectivity = mapl_get_pointer_placeholder("DBZ_TOP")
        self.outputs.minus_10c_reflectivity = mapl_get_pointer_placeholder("DBZ_M10C")
        # Unused fields, force to zero
        self.temporaries.all_zeros_3d = mapl_get_pointer_placeholder("CN_PRCP")
        self.temporaries.all_zeros_3d = mapl_get_pointer_placeholder("AN_PRCP")
        self.temporaries.all_zeros_3d = mapl_get_pointer_placeholder("SC_PRCP")
        self.temporaries.all_zeros_3d = mapl_get_pointer_placeholder("CN_SNR")
        self.temporaries.all_zeros_3d = mapl_get_pointer_placeholder("AN_SNR")
        self.temporaries.all_zeros_3d = mapl_get_pointer_placeholder("SC_SNR")

    def __call__(self, *args, **kwds):
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
            anv_icefall=self.ANV_ICEFALL,
            ls_icefall=self.LS_ICEFALL,
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
