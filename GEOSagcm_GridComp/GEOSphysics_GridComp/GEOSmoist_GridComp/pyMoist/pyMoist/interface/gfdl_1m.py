from typing import Any
from mpi4py import MPI
from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from pyMoist.GFDL_1M import GFDL1M, GFDL1MState, GFDL1MConfig
from pyMoist.interface import (
    MAPLMemoryRepository,
    MAPLManagedMemory,
    StencilBackendCompilerOverride,
    TimedCUDAProfiler,
)
import numpy as np


class GFDL1MInterface:
    def __init__(self, quantity_factory: QuantityFactory, stencil_factory: StencilFactory) -> None:
        self._GFDL_1M_ready: bool = False
        self._stencil_factory = stencil_factory
        self._GFDL_1M_state = GFDL1MState.zeros(quantity_factory)
        self._quantity_factory = quantity_factory
        self._gfdl_1m: GFDL1M | None = None

    def init(
        self,
        config: GFDL1MConfig,
        mapl_internal: MAPLMemoryRepository,
        mapl_import: MAPLMemoryRepository,
        mapl_export: MAPLMemoryRepository,
    ):
        # Initalize the module
        with StencilBackendCompilerOverride(
            MPI.COMM_WORLD,
            self._stencil_factory.config.dace_config,
        ):
            self._gfdl_1m = GFDL1M(self._stencil_factory, self._quantity_factory, config)

        # Link Fortran memory to Python memory #####
        # Fortran memory will only be modified if GFDL1M.__call__
        # is called from within a "with MAPLManagedMemory" statement #####
        # Not all linked fields are modified #####

        # manage_state.register("area", "AREA", mapl_import, dims=[X_DIM, Y_DIM])
        # manage_state.register("total_water.flux", "WQT", mapl_import)

        mapl_internal.register("Q", np.float32, [X_DIM, Y_DIM, Z_DIM])
        mapl_internal.register("QRAIN", np.float32, [X_DIM, Y_DIM, Z_DIM])
        mapl_internal.register("QSNOW", np.float32, [X_DIM, Y_DIM, Z_DIM])
        mapl_internal.register("QGRAUPEL", np.float32, [X_DIM, Y_DIM, Z_DIM])
        mapl_internal.register("QLCN", np.float32, [X_DIM, Y_DIM, Z_DIM])
        mapl_internal.register("QICN", np.float32, [X_DIM, Y_DIM, Z_DIM])
        mapl_internal.register("QLLS", np.float32, [X_DIM, Y_DIM, Z_DIM])
        mapl_internal.register("QILS", np.float32, [X_DIM, Y_DIM, Z_DIM])
        mapl_internal.register("CLCN", np.float32, [X_DIM, Y_DIM, Z_DIM])
        mapl_internal.register("CLLS", np.float32, [X_DIM, Y_DIM, Z_DIM])
        mapl_internal.register("NACTL", np.float32, [X_DIM, Y_DIM, Z_DIM])
        mapl_internal.register("NACTI", np.float32, [X_DIM, Y_DIM, Z_DIM])

        mapl_import.register("AREA", np.float32, [X_DIM, Y_DIM])
        mapl_import.register("PLE", np.float32, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        mapl_import.register("ZLE", np.float32, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        mapl_import.register("T", np.float32, [X_DIM, Y_DIM, Z_DIM])
        mapl_import.register("U", np.float32, [X_DIM, Y_DIM, Z_DIM])
        mapl_import.register("V", np.float32, [X_DIM, Y_DIM, Z_DIM])
        mapl_import.register("FRLAND", np.float32, [X_DIM, Y_DIM])
        mapl_import.register("W", np.float32, [X_DIM, Y_DIM, Z_DIM])
        mapl_import.register("W2", np.float32, [X_DIM, Y_DIM, Z_DIM])
        mapl_import.register("W3", np.float32, [X_DIM, Y_DIM, Z_DIM])
        mapl_import.register("WSL", np.float32, [X_DIM, Y_DIM, Z_DIM])
        mapl_import.register("SL2", np.float32, [X_DIM, Y_DIM, Z_DIM])
        mapl_import.register("SL3", np.float32, [X_DIM, Y_DIM, Z_DIM])
        mapl_import.register("WQT", np.float32, [X_DIM, Y_DIM, Z_DIM])
        mapl_import.register("QT2", np.float32, [X_DIM, Y_DIM, Z_DIM])
        mapl_import.register("QT3", np.float32, [X_DIM, Y_DIM, Z_DIM])
        mapl_import.register("OMEGA", np.float32, [X_DIM, Y_DIM, Z_DIM])
        mapl_import.register("PDF_A", np.float32, [X_DIM, Y_DIM, Z_DIM])
        mapl_import.register("SLQT", np.float32, [X_DIM, Y_DIM, Z_DIM])
        # MAPL_GetPointer fails on:
        # mapl_import.register("TS", np.float32, [X_DIM, Y_DIM, Z_DIM])
        # mapl_import.register("KH", np.float32, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        # mapl_import.register("SH", np.float32, [X_DIM, Y_DIM, Z_DIM])

        mapl_export.register("CNV_FRC", np.float32, [X_DIM, Y_DIM], True)
        mapl_export.register("SRF_TYPE", np.float32, [X_DIM, Y_DIM], True)
        mapl_export.register("SHLW_PRC3", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("SHLW_SNO3", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("RHCRIT", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("RL", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("RI", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("EVAPC", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("SUBLC", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("PRCP_RAIN", np.float32, [X_DIM, Y_DIM], True)
        mapl_export.register("PRCP_SNOW", np.float32, [X_DIM, Y_DIM], True)
        mapl_export.register("PRCP_ICE", np.float32, [X_DIM, Y_DIM], True)
        mapl_export.register("PRCP_GRAUPEL", np.float32, [X_DIM, Y_DIM], True)
        mapl_export.register("FCLD", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("QV", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("QL", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("QI", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("QR", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("QS", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("QG", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("LTS", np.float32, [X_DIM, Y_DIM], True)
        mapl_export.register("EIS", np.float32, [X_DIM, Y_DIM], True)
        mapl_export.register("ZLCL", np.float32, [X_DIM, Y_DIM])
        mapl_export.register("DUDT_macro", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("DVDT_macro", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("DTDT_macro", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("DQVDT_macro", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("DQLDT_macro", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("DQIDT_macro", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("DQADT_macro", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("DQRDT_macro", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("DQSDT_macro", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("DQGDT_macro", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("DUDT_micro", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("DVDT_micro", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("DTDT_micro", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("DQVDT_micro", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("DQLDT_micro", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("DQIDT_micro", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("DQADT_micro", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("DQRDT_micro", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("DQSDT_micro", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("DQGDT_micro", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("LS_PRCP", np.float32, [X_DIM, Y_DIM], True)
        mapl_export.register("LS_SNR", np.float32, [X_DIM, Y_DIM], True)
        mapl_export.register("ICE", np.float32, [X_DIM, Y_DIM], True)
        mapl_export.register("FRZR", np.float32, [X_DIM, Y_DIM], True)
        mapl_export.register("RHX", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("REV_LS", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("RSU_LS", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("PFL_LS", np.float32, [X_DIM, Y_DIM, Z_INTERFACE_DIM], True)
        mapl_export.register("PFI_LS", np.float32, [X_DIM, Y_DIM, Z_INTERFACE_DIM], True)
        mapl_export.register("PFL_AN", np.float32, [X_DIM, Y_DIM, Z_INTERFACE_DIM], True)
        mapl_export.register("PFI_AN", np.float32, [X_DIM, Y_DIM, Z_INTERFACE_DIM], True)
        mapl_export.register("WTHV2", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("WQL", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("PDFITERS", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        mapl_export.register("DQRL", np.float32, [X_DIM, Y_DIM, Z_DIM])
        mapl_export.register("DTDTFRIC", np.float32, [X_DIM, Y_DIM, Z_DIM])
        mapl_export.register("DBZ", np.float32, [X_DIM, Y_DIM, Z_DIM])
        mapl_export.register("DBZ_MAX", np.float32, [X_DIM, Y_DIM])
        mapl_export.register("DBZ_1KM", np.float32, [X_DIM, Y_DIM])
        mapl_export.register("DBZ_TOP", np.float32, [X_DIM, Y_DIM])
        mapl_export.register("DBZ_M10C", np.float32, [X_DIM, Y_DIM])
        mapl_export.register("CN_PRCP", np.float32, [X_DIM, Y_DIM], True)
        mapl_export.register("AN_PRCP", np.float32, [X_DIM, Y_DIM], True)
        mapl_export.register("SC_PRCP", np.float32, [X_DIM, Y_DIM], True)
        mapl_export.register("CN_SNR", np.float32, [X_DIM, Y_DIM], True)
        mapl_export.register("AN_SNR", np.float32, [X_DIM, Y_DIM], True)
        mapl_export.register("SC_SNR", np.float32, [X_DIM, Y_DIM], True)

    def run(
        self,
        timings: dict[str, Any],
        mapl_internal_state: MAPLMemoryRepository,
        mapl_import_state: MAPLMemoryRepository,
        mapl_export_state: MAPLMemoryRepository,
    ):
        if self._gfdl_1m is None:
            raise RuntimeError("GFDL1M Runtime called before initialization was done. Abort.")
        with TimedCUDAProfiler("GFDL 1M", timings):
            with MAPLManagedMemory(mapl_internal_state) as mapl_internal:
                with MAPLManagedMemory(mapl_import_state) as mapl_import:
                    with MAPLManagedMemory(mapl_export_state) as mapl_export:
                        # Outputs: model fields originating from within GFDL
                        with TimedCUDAProfiler("GFDL 1M - State copy", timings):
                            cov_key = (
                                "covariance_liquid_water_static_energy_and_total_water_specific_humidity"
                            )
                            self._GFDL_1M_state.update_copy_memory(
                                {
                                    "area": mapl_import.AREA,
                                    "p_interface": mapl_import.PLE,
                                    "z_interface": mapl_import.ZLE,
                                    "t": mapl_import.T,
                                    "u": mapl_import.U,
                                    "v": mapl_import.V,
                                    "land_fraction": mapl_import.FRLAND,
                                    # See mapl_import.KH in registers
                                    "scalar_diffusivity_interface": None,
                                    "pdf_first_plume_fractional_area": mapl_import.PDF_A,
                                    cov_key: mapl_import.SLQT,
                                    # See mapl_import.TS in registers
                                    "surface_temperature": None,
                                    # See mapl_import.SH in registers
                                    "sensible_heat_flux": None,
                                    "omega": mapl_import.OMEGA,
                                    "convection_fraction": mapl_export.CNV_FRC,
                                    "surface_type": mapl_export.SRF_TYPE,
                                    "cloud_liquid_evaporation": mapl_export.EVAPC,
                                    "cloud_ice_sublimation": mapl_export.SUBLC,
                                    "icefall": mapl_export.ICE,
                                    "freezing_rainfall": mapl_export.FRZR,
                                    "relative_humidity_after_pdf": mapl_export.RHX,
                                    "buoyancy_flux": mapl_export.WTHV2,
                                    "liquid_water_flux": mapl_export.WQL,
                                    "hydrostatic_pdf_iterations": mapl_export.PDFITERS,
                                    "lower_tropospheric_stability": mapl_export.LTS,
                                    "estimated_inversion_strength": mapl_export.EIS,
                                    "lcl_height": mapl_export.ZLCL,
                                    "shallow_convection_rain": mapl_export.SHLW_PRC3,
                                    "shallow_convection_snow": mapl_export.SHLW_SNO3,
                                    "critical_relative_humidity_for_pdf": mapl_export.RHCRIT,
                                    "large_scale_rainwater_source": mapl_export.DQRL,
                                    "vertical_motion": {
                                        "velocity": mapl_import.W,
                                        "variance": mapl_import.W2,
                                        "third_moment": mapl_import.W3,
                                    },
                                    "mixing_ratio": {
                                        "vapor": mapl_internal.Q,
                                        "rain": mapl_internal.QRAIN,
                                        "snow": mapl_internal.QSNOW,
                                        "graupel": mapl_internal.QGRAUPEL,
                                        "large_scale_liquid": mapl_internal.QLLS,
                                        "large_scale_ice": mapl_internal.QILS,
                                        "convective_liquid": mapl_internal.QLCN,
                                        "convective_ice": mapl_internal.QICN,
                                    },
                                    "cloud_fraction": {
                                        "large_scale": mapl_internal.CLLS,
                                        "convective": mapl_internal.CLCN,
                                    },
                                    "concentration": {
                                        "liquid": mapl_internal.NACTL,
                                        "ice": mapl_internal.NACTI,
                                    },
                                    "liquid_water_static_energy": {
                                        "flux": mapl_import.WSL,
                                        "variance": mapl_import.SL2,
                                        "third_moment": mapl_import.SL3,
                                    },
                                    "total_water": {
                                        "flux": mapl_import.WQT,
                                        "variance": mapl_import.QT2,
                                        "third_moment": mapl_import.QT3,
                                    },
                                    "radiation_field": {
                                        "cloud_fraction": mapl_export.FCLD,
                                        "vapor": mapl_export.QV,
                                        "liquid": mapl_export.QL,
                                        "ice": mapl_export.QI,
                                        "rain": mapl_export.QR,
                                        "snow": mapl_export.QS,
                                        "graupel": mapl_export.QG,
                                    },
                                    "cloud_particle_effective_radius": {
                                        "liquid": mapl_export.RL,
                                        "ice": mapl_export.RI,
                                    },
                                    "precipitation_at_surface": {
                                        "rain": mapl_export.PRCP_RAIN,
                                        "snow": mapl_export.PRCP_SNOW,
                                        "ice": mapl_export.PRCP_ICE,
                                        "graupel": mapl_export.PRCP_GRAUPEL,
                                        "shallow_convective_precipitation": mapl_export.SC_PRCP,
                                        "deep_convective_precipitation": mapl_export.CN_PRCP,
                                        "anvil_precipitation": mapl_export.AN_PRCP,
                                        "shallow_convective_snow": mapl_export.SC_SNR,
                                        "deep_convective_snow": mapl_export.CN_SNR,
                                        "anvil_snow": mapl_export.AN_SNR,
                                    },
                                    "non_anvil_large_scale": {
                                        "precip": mapl_export.LS_PRCP,
                                        "snow": mapl_export.LS_SNR,
                                        "evaporation": mapl_export.REV_LS,
                                        "sublimation": mapl_export.RSU_LS,
                                        "liquid_precip_flux": mapl_export.PFL_LS,
                                        "ice_precip_flux": mapl_export.PFI_LS,
                                    },
                                    "anvil": {
                                        "liquid_precip_flux": mapl_export.PFL_AN,
                                        "ice_precip_flux": mapl_export.PFI_AN,
                                    },
                                    "tendencies": {
                                        "dcloud_fractiondt_macro": mapl_export.DQADT_macro,
                                        "dvapordt_macro": mapl_export.DQVDT_macro,
                                        "dicedt_macro": mapl_export.DQIDT_macro,
                                        "dliquiddt_macro": mapl_export.DQLDT_macro,
                                        "draindt_macro": mapl_export.DQRDT_macro,
                                        "dgraupeldt_macro": mapl_export.DQGDT_macro,
                                        "dsnowdt_macro": mapl_export.DQSDT_micro,
                                        "dudt_macro": mapl_export.DUDT_micro,
                                        "dvdt_macro": mapl_export.DVDT_micro,
                                        "dtdt_macro": mapl_export.DTDT_micro,
                                        "dcloud_fractiondt_micro": mapl_export.DQADT_micro,
                                        "dvapordt_micro": mapl_export.DQVDT_micro,
                                        "dicedt_micro": mapl_export.DQIDT_micro,
                                        "dliquiddt_micro": mapl_export.DQLDT_micro,
                                        "draindt_micro": mapl_export.DQRDT_micro,
                                        "dgraupeldt_micro": mapl_export.DQGDT_micro,
                                        "dsnowdt_micro": mapl_export.DQSDT_micro,
                                        "dudt_micro": mapl_export.DUDT_micro,
                                        "dvdt_micro": mapl_export.DVDT_micro,
                                        "dtdt_micro": mapl_export.DTDT_micro,
                                        "dtdt_friction_pressure_weighted": mapl_export.DTDTFRIC,
                                    },
                                    "radar": {
                                        "simulated_reflectivity": mapl_export.DBZ,
                                        "maximum_composite_reflectivity": mapl_export.DBZ_MAX,
                                        "base_1km_agl_reflectivity": mapl_export.DBZ_1KM,
                                        "echo_top_reflectivity": mapl_export.DBZ_TOP,
                                        "minus_10c_reflectivity": mapl_export.DBZ_M10C,
                                    },
                                },
                            )

                        # Call the module
                        with TimedCUDAProfiler("GFDL 1M Numerics", timings):
                            self._gfdl_1m(self._GFDL_1M_state)
