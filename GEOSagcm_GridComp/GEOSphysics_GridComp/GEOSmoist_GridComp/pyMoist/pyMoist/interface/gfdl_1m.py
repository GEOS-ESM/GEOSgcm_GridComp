from typing import Any

from mpi4py import MPI

from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_INTERFACE_DIM
from pyMoist.GFDL_1M import GFDL1M, GFDL1MConfig, GFDL1MState
from pyMoist.interface import (
    InterfaceTransferType,
    MAPLManagedState,
    MAPLMemoryRepository,
    StencilBackendCompilerOverride,
    TimedCUDAProfiler,
)


class GFDL1MInterface:
    def __init__(
        self,
        quantity_factory: QuantityFactory,
        stencil_factory: StencilFactory,
        transfer_type: InterfaceTransferType,
    ) -> None:
        self._stencil_factory = stencil_factory
        self._quantity_factory = quantity_factory
        self._gfdl_1m: GFDL1M | None = None
        self._managed_state = MAPLManagedState(
            GFDL1MState.zeros(quantity_factory),
            transfer_type,
        )

    def init(
        self,
        config: GFDL1MConfig,
        mapl_internal: MAPLMemoryRepository,
        mapl_import: MAPLMemoryRepository,
        mapl_export: MAPLMemoryRepository,
    ):
        # Initialize the module
        with StencilBackendCompilerOverride(
            MPI.COMM_WORLD,
            self._stencil_factory.config.dace_config,
        ):
            self._gfdl_1m = GFDL1M(self._stencil_factory, self._quantity_factory, config)

        self._managed_state.register("mixing_ratio.vapor", "Q", mapl_internal)
        self._managed_state.register("mixing_ratio.rain", "QRAIN", mapl_internal)
        self._managed_state.register("mixing_ratio.snow", "QSNOW", mapl_internal)
        self._managed_state.register("mixing_ratio.graupel", "QGRAUPEL", mapl_internal)
        self._managed_state.register("mixing_ratio.large_scale_liquid", "QLLS", mapl_internal)
        self._managed_state.register("mixing_ratio.large_scale_ice", "QILS", mapl_internal)
        self._managed_state.register("mixing_ratio.convective_liquid", "QLCN", mapl_internal)
        self._managed_state.register("mixing_ratio.convective_ice", "QICN", mapl_internal)
        self._managed_state.register("cloud_fraction.large_scale", "CLLS", mapl_internal)
        self._managed_state.register("cloud_fraction.convective", "CLCN", mapl_internal)
        self._managed_state.register("concentration.liquid", "NACTL", mapl_internal)
        self._managed_state.register("concentration.ice", "NACTI", mapl_internal)

        self._managed_state.register_2D("area", "AREA", mapl_import)
        self._managed_state.register("p_interface", "PLE", mapl_import, dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM])
        self._managed_state.register("z_interface", "ZLE", mapl_import, dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM])
        self._managed_state.register("t", "T", mapl_import)
        self._managed_state.register("u", "U", mapl_import)
        self._managed_state.register("v", "V", mapl_import)
        self._managed_state.register_2D("land_fraction", "FRLAND", mapl_import)
        self._managed_state.register(
            "covariance_liquid_water_static_energy_and_total_water_specific_humidity", "SLQT", mapl_import
        )
        self._managed_state.register("omega", "OMEGA", mapl_import)
        self._managed_state.register("pdf_first_plume_fractional_area", "PDF_A", mapl_import)

        self._managed_state.register("vertical_motion.velocity", "W", mapl_import)
        self._managed_state.register("vertical_motion.variance", "W2", mapl_import)
        self._managed_state.register("vertical_motion.third_moment", "W3", mapl_import)

        self._managed_state.register("liquid_water_static_energy.flux", "WSL", mapl_import)
        self._managed_state.register("liquid_water_static_energy.variance", "SL2", mapl_import)
        self._managed_state.register("liquid_water_static_energy.third_moment", "SL3", mapl_import)

        self._managed_state.register("total_water.flux", "WQT", mapl_import)
        self._managed_state.register("total_water.variance", "QT2", mapl_import)
        self._managed_state.register("total_water.third_moment", "QT3", mapl_import)

        # MAPL_GetPointer fails on:
        # self_manage_state.register("surface_temperature", "TS", mapl_import)
        # self._manage_state.register(
        #     "scalar_diffusivity_interface",
        #     "KH",
        #     mapl_import,
        #     dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM],
        # )
        # self_manage_state.register("sensible_heat_flux", "SH", mapl_import)

        self._managed_state.register_2D("convection_fraction", "CNV_FRC", mapl_export, alloc=True)
        self._managed_state.register_2D("surface_type", "SRF_TYPE", mapl_export, alloc=True)
        self._managed_state.register("cloud_liquid_evaporation", "EVAPC", mapl_export, alloc=True)
        self._managed_state.register("cloud_ice_sublimation", "SUBLC", mapl_export, alloc=True)
        self._managed_state.register_2D("icefall", "ICE", mapl_export, alloc=True)
        self._managed_state.register_2D("freezing_rainfall", "FRZR", mapl_export, alloc=True)
        self._managed_state.register("relative_humidity_after_pdf", "RHX", mapl_export, alloc=True)
        self._managed_state.register("buoyancy_flux", "WTHV2", mapl_export, alloc=True)
        self._managed_state.register("liquid_water_flux", "WQL", mapl_export, alloc=True)
        self._managed_state.register("hydrostatic_pdf_iterations", "PDFITERS", mapl_export, alloc=True)
        self._managed_state.register_2D("lower_tropospheric_stability", "LTS", mapl_export, alloc=True)
        self._managed_state.register_2D("estimated_inversion_strength", "EIS", mapl_export, alloc=True)
        self._managed_state.register_2D("lcl_height", "ZLCL", mapl_export)
        self._managed_state.register("shallow_convection_rain", "SHLW_PRC3", mapl_export, alloc=True)
        self._managed_state.register("shallow_convection_snow", "SHLW_SNO3", mapl_export, alloc=True)
        self._managed_state.register("critical_relative_humidity_for_pdf", "RHCRIT", mapl_export, alloc=True)
        self._managed_state.register("large_scale_rainwater_source", "DQRL", mapl_export)

        self._managed_state.register("radiation_field.cloud_fraction", "FCLD", mapl_export, alloc=True)
        self._managed_state.register("radiation_field.vapor", "QV", mapl_export, alloc=True)
        self._managed_state.register("radiation_field.liquid", "QL", mapl_export, alloc=True)
        self._managed_state.register("radiation_field.ice", "QI", mapl_export, alloc=True)
        self._managed_state.register("radiation_field.rain", "QR", mapl_export, alloc=True)
        self._managed_state.register("radiation_field.snow", "QS", mapl_export, alloc=True)
        self._managed_state.register("radiation_field.graupel", "QG", mapl_export, alloc=True)

        self._managed_state.register("cloud_particle_effective_radius.liquid", "RL", mapl_export, alloc=True)
        self._managed_state.register("cloud_particle_effective_radius.ice", "RI", mapl_export, alloc=True)

        self._managed_state.register_2D("precipitation_at_surface.rain", "PRCP_RAIN", mapl_export, alloc=True)
        self._managed_state.register_2D("precipitation_at_surface.snow", "PRCP_SNOW", mapl_export, alloc=True)
        self._managed_state.register_2D("precipitation_at_surface.ice", "PRCP_ICE", mapl_export, alloc=True)
        self._managed_state.register_2D(
            "precipitation_at_surface.graupel", "PRCP_GRAUPEL", mapl_export, alloc=True
        )
        self._managed_state.register_2D(
            "precipitation_at_surface.shallow_convective_precipitation", "SC_PRCP", mapl_export, alloc=True
        )
        self._managed_state.register_2D(
            "precipitation_at_surface.deep_convective_precipitation", "CN_PRCP", mapl_export, alloc=True
        )
        self._managed_state.register_2D(
            "precipitation_at_surface.anvil_precipitation", "AN_PRCP", mapl_export, alloc=True
        )
        self._managed_state.register_2D(
            "precipitation_at_surface.shallow_convective_snow", "SC_SNR", mapl_export, alloc=True
        )
        self._managed_state.register_2D(
            "precipitation_at_surface.deep_convective_snow", "CN_SNR", mapl_export, alloc=True
        )
        self._managed_state.register_2D(
            "precipitation_at_surface.anvil_snow", "AN_SNR", mapl_export, alloc=True
        )

        self._managed_state.register_2D("non_anvil_large_scale.precip", "LS_PRCP", mapl_export, alloc=True)
        self._managed_state.register_2D("non_anvil_large_scale.snow", "LS_SNR", mapl_export, alloc=True)
        self._managed_state.register("non_anvil_large_scale.evaporation", "REV_LS", mapl_export, alloc=True)
        self._managed_state.register("non_anvil_large_scale.sublimation", "RSU_LS", mapl_export, alloc=True)
        self._managed_state.register(
            "non_anvil_large_scale.liquid_precip_flux",
            "PFL_LS",
            mapl_export,
            dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM],
            alloc=True,
        )
        self._managed_state.register(
            "non_anvil_large_scale.ice_precip_flux",
            "PFI_LS",
            mapl_export,
            dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM],
            alloc=True,
        )

        self._managed_state.register(
            "anvil.liquid_precip_flux",
            "PFL_AN",
            mapl_export,
            dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM],
            alloc=True,
        )
        self._managed_state.register(
            "anvil.ice_precip_flux",
            "PFI_AN",
            mapl_export,
            dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM],
            alloc=True,
        )

        self._managed_state.register(
            "tendencies.dcloud_fractiondt_macro", "DQADT_macro", mapl_export, alloc=True
        )
        self._managed_state.register("tendencies.dvapordt_macro", "DQVDT_macro", mapl_export, alloc=True)
        self._managed_state.register("tendencies.dicedt_macro", "DQIDT_macro", mapl_export, alloc=True)
        self._managed_state.register("tendencies.dliquiddt_macro", "DQLDT_macro", mapl_export, alloc=True)
        self._managed_state.register("tendencies.draindt_macro", "DQRDT_macro", mapl_export, alloc=True)
        self._managed_state.register("tendencies.dgraupeldt_macro", "DQGDT_macro", mapl_export, alloc=True)
        self._managed_state.register("tendencies.dsnowdt_macro", "DQSDT_macro", mapl_export, alloc=True)
        self._managed_state.register("tendencies.dudt_macro", "DUDT_macro", mapl_export, alloc=True)
        self._managed_state.register("tendencies.dvdt_macro", "DVDT_macro", mapl_export, alloc=True)
        self._managed_state.register("tendencies.dtdt_macro", "DTDT_macro", mapl_export, alloc=True)
        self._managed_state.register(
            "tendencies.dcloud_fractiondt_micro", "DQADT_micro", mapl_export, alloc=True
        )
        self._managed_state.register("tendencies.dvapordt_micro", "DQVDT_micro", mapl_export, alloc=True)
        self._managed_state.register("tendencies.dicedt_micro", "DQIDT_micro", mapl_export, alloc=True)
        self._managed_state.register("tendencies.dliquiddt_micro", "DQLDT_micro", mapl_export, alloc=True)
        self._managed_state.register("tendencies.draindt_micro", "DQRDT_micro", mapl_export, alloc=True)
        self._managed_state.register("tendencies.dgraupeldt_micro", "DQGDT_micro", mapl_export, alloc=True)
        self._managed_state.register("tendencies.dsnowdt_micro", "DQSDT_micro", mapl_export, alloc=True)
        self._managed_state.register("tendencies.dudt_micro", "DUDT_micro", mapl_export, alloc=True)
        self._managed_state.register("tendencies.dvdt_micro", "DVDT_micro", mapl_export, alloc=True)
        self._managed_state.register("tendencies.dtdt_micro", "DTDT_micro", mapl_export, alloc=True)
        self._managed_state.register("tendencies.dtdt_friction_pressure_weighted", "DTDTFRIC", mapl_export)

        self._managed_state.register("radar.simulated_reflectivity", "DBZ", mapl_export)
        self._managed_state.register_2D("radar.maximum_composite_reflectivity", "DBZ_MAX", mapl_export)
        self._managed_state.register_2D("radar.base_1km_agl_reflectivity", "DBZ_1KM", mapl_export)
        self._managed_state.register_2D("radar.echo_top_reflectivity", "DBZ_TOP", mapl_export)
        self._managed_state.register_2D("radar.minus_10c_reflectivity", "DBZ_M10C", mapl_export)

    def run(
        self,
        timings: dict[str, Any],
    ):
        if self._gfdl_1m is None:
            raise RuntimeError("GFDL1M Runtime called before initialization was done. Abort.")

        with TimedCUDAProfiler("GFDL 1M", timings):
            with TimedCUDAProfiler("GFDL 1M - State copy", timings):
                self._managed_state.fortran_to_ndsl()

            with TimedCUDAProfiler("GFDL 1M Numerics", timings):
                self._gfdl_1m(self._managed_state.ndsl_state)

            with TimedCUDAProfiler("GFDL 1M - State copy-back", timings):
                self._managed_state.ndsl_to_fortran()
