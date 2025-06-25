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

    def get_fortran_data(
        self,
        mapl_internal: MAPLMemoryRepository,
        mapl_import: MAPLMemoryRepository,
        mapl_export: MAPLMemoryRepository,
        mapl_comp: MAPLMemoryRepository,
    ):
        ##### Link Fortran memory to Python memory #####
        ##### Fortran memory will only be modified if the __call__ function
        ##### is called from within a "with MAPLManagedMemory" statement #####
        ##### Not all linked fields are modified #####

        from numpy import float32

        # Get model state
        self.mixing_ratios = MixingRatios(
            vapor=(mapl_import.register("Q", float32, [X_DIM, Y_DIM, Z_DIM])),
            rain=(mapl_import.register("QRAIN", float32, [X_DIM, Y_DIM, Z_DIM])),
            snow=(mapl_import.register("QSNOW", float32, [X_DIM, Y_DIM, Z_DIM])),
            graupel=(mapl_import.register("QGRAUPEL", float32, [X_DIM, Y_DIM, Z_DIM])),
            convective_liquid=(mapl_import.register("QLCN", float32, [X_DIM, Y_DIM, Z_DIM])),
            convective_ice=(mapl_import.register("QICN", float32, [X_DIM, Y_DIM, Z_DIM])),
            large_scale_liquid=(mapl_import.register("QLLS", float32, [X_DIM, Y_DIM, Z_DIM])),
            large_scale_ice=(mapl_import.register("QILS", float32, [X_DIM, Y_DIM, Z_DIM])),
        )
        self.cloud_fractions = CloudFractions(
            convective=(mapl_import.register("CLCN", float32, [X_DIM, Y_DIM, Z_DIM])),
            large_scale=(mapl_import.register("CLLS", float32, [X_DIM, Y_DIM, Z_DIM])),
        )
        self.liquid_concentration = mapl_import.register("NACTL", float32, [X_DIM, Y_DIM, Z_DIM])
        self.ice_concentration = mapl_import.register("NACTI", float32, [X_DIM, Y_DIM, Z_DIM])
        self.area = mapl_import.register("AREA", float32, [X_DIM, Y_DIM])
        self.geopotential_height_interface = mapl_import.register(
            "ZLE", float32, [X_DIM, Y_DIM, Z_INTERFACE_DIM]
        )
        self.p_interface = mapl_import.register("PLE", float32, [X_DIM, Y_DIM])
        self.t = mapl_import.register("T", float32, [X_DIM, Y_DIM, Z_DIM])
        self.u = mapl_import.register("U", float32, [X_DIM, Y_DIM, Z_DIM])
        self.v = mapl_import.register("V", float32, [X_DIM, Y_DIM, Z_DIM])
        self.land_fraction = mapl_import.register("FRLAND", float32, [X_DIM, Y_DIM])
        self.vertical_motion = VericalMotion(
            velocity=(mapl_import.register("W", float32, [X_DIM, Y_DIM, Z_DIM])),
            variance=(mapl_import.register("W2", float32, [X_DIM, Y_DIM, Z_DIM])),
            third_moment=(mapl_import.register("W3", float32, [X_DIM, Y_DIM, Z_DIM])),
        )
        self.liquid_water_static_energy = LiquidWaterStaticEnergy(
            flux=(mapl_import.register("WSL", float32, [X_DIM, Y_DIM, Z_DIM])),
            variance=(mapl_import.register("SL2", float32, [X_DIM, Y_DIM, Z_DIM])),
            third_moment=(mapl_import.register("SL3", float32, [X_DIM, Y_DIM, Z_DIM])),
        )
        self.total_water = TotalWater(
            flux=(mapl_import.register("WQT", float32, [X_DIM, Y_DIM, Z_DIM])),
            variance=(mapl_import.register("QT2", float32, [X_DIM, Y_DIM, Z_DIM])),
            third_moment=(mapl_import.register("QT3", float32, [X_DIM, Y_DIM, Z_DIM])),
        )
        self.convection_fraction = mapl_export.register("CNV_FRC", float32, [X_DIM, Y_DIM])
        self.surface_type = mapl_export.register("SRF_TYPE", float32, [X_DIM, Y_DIM])
        self.shallow_convective_rain = mapl_get_pointer_placeholder("SHLW_PRC3")
        self.shallow_convective_snow = mapl_get_pointer_placeholder("SHLW_SNO3")
        self.rh_crit = mapl_get_pointer_placeholder("RHCRIT")

        # Outputs: model fields originating from within GFDL
        self.outputs.liquid_radius = mapl_export.register("RL", float32, [X_DIM, Y_DIM, Z_DIM])
        self.outputs.ice_radius = mapl_export.reiregister("RI", float32, [X_DIM, Y_DIM, Z_DIM])
        self.outputs.large_scale_nonanvil_precipitation_evaporation = mapl_export.register(
            "EVAPC", float32, [X_DIM, Y_DIM, Z_DIM]
        )
        self.outputs.large_scale_nonanvil_precipitation_sublimation = mapl_export.register(
            "SUBLC", float32, [X_DIM, Y_DIM, Z_DIM]
        )
        self.outputs.precipitated_rain = mapl_export.register("PRCP_RAIN", float32, [X_DIM, Y_DIM])
        self.outputs.precipitated_snow = mapl_export.register("PRCP_SNOW", float32, [X_DIM, Y_DIM])
        self.outputs.precipitated_ice = mapl_export.register("PRCP_ICE", float32, [X_DIM, Y_DIM])
        self.outputs.precipitated_graupel = mapl_export.register("PRCP_GRAUPEL", float32, [X_DIM, Y_DIM])

        # Outputs: model fields originating from within GFDL; radiation fields
        self.outputs.radiation_cloud_fraction = mapl_export.register("FCLD", float32, [X_DIM, Y_DIM, Z_DIM])
        self.outputs.radiation_vapor = mapl_export.register("QV", float32, [X_DIM, Y_DIM, Z_DIM])
        self.outputs.radiation_liquid = mapl_export.register("QL", float32, [X_DIM, Y_DIM, Z_DIM])
        self.outputs.radiation_ice = mapl_export.register("QI", float32, [X_DIM, Y_DIM, Z_DIM])
        self.outputs.radiation_rain = mapl_export.register("QR", float32, [X_DIM, Y_DIM, Z_DIM])
        self.outputs.radiation_snow = mapl_export.register("QS", float32, [X_DIM, Y_DIM, Z_DIM])
        self.outputs.radiation_graupel = mapl_export.register("QG", float32, [X_DIM, Y_DIM, Z_DIM])
        self.outputs.lower_tropospheric_stability = (
            mapl_export.register("LTS"),
            float32,
            [X_DIM, Y_DIM, Z_DIM],
        )
        self.outputs.estimated_inversion_strength = (
            mapl_export.register("EIS"),
            float32,
            [X_DIM, Y_DIM, Z_DIM],
        )
        self.outputs.z_lcl = (mapl_export.register("ZLCL"), float32, [X_DIM, Y_DIM])

        # Outputs: model fields originating from within GFDL; macrophysics/microphysics tendencies
        self.outputs.du_dt_macro = mapl_export.register("DUDT_macro", float32, [X_DIM, Y_DIM, Z_DIM])
        self.outputs.dv_dt_macro = mapl_export.register("DVDT_macro", float32, [X_DIM, Y_DIM, Z_DIM])
        self.outputs.dt_dt_macro = mapl_export.register("DTDT_macro", float32, [X_DIM, Y_DIM, Z_DIM])
        self.outputs.dvapor_dt_macro = mapl_export.register("DQVDT_macro", float32, [X_DIM, Y_DIM, Z_DIM])
        self.outputs.dliquid_dt_macro = mapl_export.register("DQLDT_macro", float32, [X_DIM, Y_DIM, Z_DIM])
        self.outputs.dice_dt_macro = mapl_export.register("DQIDT_macro", float32, [X_DIM, Y_DIM, Z_DIM])
        self.outputs.dcloud_fraction_dt_macro = mapl_export.register(
            "DQADT_macro", float32, [X_DIM, Y_DIM, Z_DIM]
        )
        self.outputs.drain_dt_macro = mapl_export.register("DQRDT_macro", float32, [X_DIM, Y_DIM, Z_DIM])
        self.outputs.dsnow_dt_macro = mapl_export.register("DQSDT_macro", float32, [X_DIM, Y_DIM, Z_DIM])
        self.outputs.dgraupel_dt_macro = mapl_export.register("DQGDT_macro", float32, [X_DIM, Y_DIM, Z_DIM])
        self.outputs.du_dt_micro = mapl_export.register("DUDT_micro", float32, [X_DIM, Y_DIM, Z_DIM])
        self.outputs.dv_dt_micro = mapl_export.register("DVDT_micro", float32, [X_DIM, Y_DIM, Z_DIM])
        self.outputs.dt_dt_micro = mapl_export.register("DTDT_micro", float32, [X_DIM, Y_DIM, Z_DIM])
        self.outputs.dvapor_dt_micro = mapl_export.register("DQVDT_micro", float32, [X_DIM, Y_DIM, Z_DIM])
        self.outputs.dliquid_dt_micro = mapl_export.register("DQLDT_micro", float32, [X_DIM, Y_DIM, Z_DIM])
        self.outputs.dice_dt_micro = mapl_export.register("DQIDT_micro", float32, [X_DIM, Y_DIM, Z_DIM])
        self.outputs.dcloud_fraction_dt_micro = mapl_export.register(
            "DQADT_micro", float32, [X_DIM, Y_DIM, Z_DIM]
        )
        self.outputs.drain_dt_micro = mapl_export.register("DQRDT_micro", float32, [X_DIM, Y_DIM, Z_DIM])
        self.outputs.dsnow_dt_micro = mapl_export.register("DQSDT_micro", float32, [X_DIM, Y_DIM, Z_DIM])
        self.outputs.dgraupel_dt_micro = mapl_export.register("DQGDT_micro", float32, [X_DIM, Y_DIM, Z_DIM])
        # Outputs: Exports to be filled
        self.outputs.large_scale_precip = mapl_export.register("LS_PRCP", float32, [X_DIM, Y_DIM, Z_DIM])
        self.outputs.large_scale_snow = mapl_export.register("LS_SNR", float32, [X_DIM, Y_DIM, Z_DIM])
        self.outputs.icefall = mapl_export.register("ICE", float32, [X_DIM, Y_DIM])
        self.outputs.freezing_rainfall = mapl_export.register("FRZR", float32, [X_DIM, Y_DIM])
        self.outputs.relative_humidity_after_pdf = mapl_export.register("RHX", float32, [X_DIM, Y_DIM, Z_DIM])
        self.outputs.large_scale_nonanvil_precipitation_evaporation = mapl_export.register(
            "REV_LS", float32, [X_DIM, Y_DIM, Z_DIM]
        )
        self.outputs.large_scale_nonanvil_precipitation_sublimation = mapl_export.register(
            "RSU_LS", float32, [X_DIM, Y_DIM, Z_DIM]
        )
        self.outputs.large_scale_nonanvil_liquid_flux = mapl_export.register(
            "PFL_LS", float32, [X_DIM, Y_DIM, Z_INTERFACE_DIM]
        )
        self.outputs.large_scale_nonanvil_ice_flux = mapl_export.register(
            "PFI_LS", float32, [X_DIM, Y_DIM, Z_INTERFACE_DIM]
        )
        self.outputs.anvil_liquid_flux = mapl_export.register(
            "PFL_AN", float32, [X_DIM, Y_DIM, Z_INTERFACE_DIM]
        )
        self.outputs.anvil_ice_flux = mapl_export.register("PFI_AN", float32, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        self.outputs.large_scale_rainwater_source = mapl_export.register(
            "DQRL", float32, [X_DIM, Y_DIM, Z_DIM]
        )
        self.outputs.simulated_reflectivity = mapl_export.register("DBZ", float32, [X_DIM, Y_DIM, Z_DIM])
        self.outputs.maximum_reflectivity = mapl_export.register("DBZ_MAX", float32, [X_DIM, Y_DIM])
        self.outputs.one_km_agl_reflectivity = mapl_export.register("DBZ_1KM", float32, [X_DIM, Y_DIM])
        self.outputs.echo_top_reflectivity = mapl_export.register("DBZ_TOP", float32, [X_DIM, Y_DIM])
        self.outputs.minus_10c_reflectivity = mapl_export.register("DBZ_M10C", float32, [X_DIM, Y_DIM])
        # Unused fields, force to zero
        self.temporaries.all_zeros_3d = self.quantity_factory.zeros(
            mapl_export.register("CN_PRCP", float32, [X_DIM, Y_DIM])
        )
        self.temporaries.all_zeros_3d = self.quantity_factory.zeros(
            mapl_export.register("AN_PRCP", float32, [X_DIM, Y_DIM])
        )
        self.temporaries.all_zeros_3d = self.quantity_factory.zeros(
            mapl_export.register("SC_PRCP", float32, [X_DIM, Y_DIM])
        )
        self.temporaries.all_zeros_3d = self.quantity_factory.zeros(
            mapl_export.register("CN_SNR", float32, [X_DIM, Y_DIM])
        )
        self.temporaries.all_zeros_3d = self.quantity_factory.zeros(
            mapl_export.register("AN_SNR", float32, [X_DIM, Y_DIM])
        )
        self.temporaries.all_zeros_3d = self.quantity_factory.zeros(
            mapl_export.register("SC_SNR", float32, [X_DIM, Y_DIM])
        )

    def __call__(self, mapl_import, mapl_export):
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
