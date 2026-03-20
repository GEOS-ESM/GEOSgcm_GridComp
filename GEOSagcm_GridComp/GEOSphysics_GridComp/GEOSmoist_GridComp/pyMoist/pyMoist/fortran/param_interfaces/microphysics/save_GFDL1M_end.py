import os
from typing import Any

import f90nml
from MAPL_PythonBridge import UserCode, get_MAPLPy
from MAPL_PythonBridge.types import CVoidPointer
from mpi4py import MPI
from ndsl.constants import I_DIM, J_DIM, K_INTERFACE_DIM
from ndsl.dsl.typing import Float, Int

from pyMoist.fortran import get_NDSL_physics
from pyMoist.fortran.build_helper import InterfaceTransferType
from pyMoist.fortran.managed_state import MAPLManagedState
from pyMoist.fortran.memory_factory import MAPLMemoryRepository
from pyMoist.microphysics.GFDL_1M import GFDL1MConfig, GFDL1MState


def _default_or_get_from_namelist(default, name_in_namelist: str, namelist: dict[str, Any]) -> Any:
    return default if name_in_namelist not in namelist else namelist[name_in_namelist]


class GFDL1M_save_end(UserCode):
    def __init__(self) -> None:
        pass

    def init(self, mapl_state: CVoidPointer, import_state: CVoidPointer, export_state: CVoidPointer):
        maplpy = get_MAPLPy()
        ndsl_stack = get_NDSL_physics(mapl_state)
        self._managed_state = MAPLManagedState(
            GFDL1MState.empty(ndsl_stack.quantity_factory),
            InterfaceTransferType.CPU_MAP,
        )

        try:
            with open("input.nml", "r") as f:
                namelist = f90nml.read(f)
                if "gfdl_cloud_microphysics_nml" in namelist:
                    gfdl1m_namelist = namelist["gfdl_cloud_microphysics_nml"]
                else:
                    gfdl1m_namelist = {}
        except FileNotFoundError:
            raise FileNotFoundError(
                f"[DSL GFDL 1M] Could not find input.nml in curent directory ({os.getcwd()})"
            )

        do_evap = _default_or_get_from_namelist(False, "do_evap", gfdl1m_namelist)
        if do_evap:
            ccw_evap_eff_default = Float(0.0)
        else:
            ccw_evap_eff_default = Float(1e-2)
        do_subl = _default_or_get_from_namelist(False, "do_subl", gfdl1m_namelist)
        if do_subl:
            cci_evap_eff_default = Float(0.0)
        else:
            cci_evap_eff_default = Float(1e-2)
        use_aerosol_nn = maplpy.get_resource("USE_AEROSOL_NN:", mapl_state, default=True)

        config = GFDL1MConfig(
            USE_BERGERON=maplpy.get_resource("USE_BERGERON:", mapl_state, default=use_aerosol_nn),
            LPHYS_HYDROSTATIC=maplpy.get_resource("PHYS_HYDROSTATIC", mapl_state, default=True),
            LHYDROSTATIC=maplpy.get_resource("HYDROSTATIC", mapl_state, default=True),
            DT_MOIST=maplpy.get_resource("DSL__GFLD1M_DT", mapl_state, default=Float(0.0)),
            LMELTFRZ=maplpy.get_resource("MELTFRZ", mapl_state, default=True),
            TURNRHCRIT_PARAM=maplpy.get_resource("TURNRHCRIT:", mapl_state, default=Float(-9999.0)),
            PDFSHAPE=maplpy.get_resource("PDFSHAPE:", mapl_state, default=Int(1)),
            ANV_ICEFALL=maplpy.get_resource("ANV_ICEFALL:", mapl_state, default=Float(1.0)),
            LS_ICEFALL=maplpy.get_resource("LS_ICEFALL:", mapl_state, default=Float(1.0)),
            LIQ_RADII_PARAM=maplpy.get_resource("LIQ_RADII_PARAM:", mapl_state, default=Int(2)),
            ICE_RADII_PARAM=maplpy.get_resource("ICE_RADII_PARAM:", mapl_state, default=Int(1)),
            FAC_RI=maplpy.get_resource("FAC_RI:", mapl_state, default=Float(1.0)),
            MIN_RI=maplpy.get_resource("MIN_RI:", mapl_state, default=Float(5e-6)),
            MAX_RI=maplpy.get_resource("MAX_RI:", mapl_state, default=Float(100.0e-6)),
            FAC_RL=maplpy.get_resource("FAC_RL:", mapl_state, default=Float(1.0)),
            MIN_RL=maplpy.get_resource("MIN_RL:", mapl_state, default=Float(2.5e-6)),
            MAX_RL=maplpy.get_resource("MAX_RL:", mapl_state, default=Float(60e-6)),
            CCW_EVAP_EFF=maplpy.get_resource("CCW_EVAP_EFF:", mapl_state, default=ccw_evap_eff_default),
            CCI_EVAP_EFF=maplpy.get_resource("CCI_EVAP_EFF:", mapl_state, default=cci_evap_eff_default),
            CLD_MIN=_default_or_get_from_namelist(Float(0.5), "cld_min", gfdl1m_namelist),
            T_MIN=_default_or_get_from_namelist(Float(178.0), "t_min", gfdl1m_namelist),
            T_SUB=_default_or_get_from_namelist(Float(184.0), "t_sub", gfdl1m_namelist),
            MP_TIME=_default_or_get_from_namelist(Float(150.0), "mp_time", gfdl1m_namelist),
            RH_INC=_default_or_get_from_namelist(Float(0.25), "rh_inc", gfdl1m_namelist),
            RH_INS=_default_or_get_from_namelist(Float(0.25), "rh_inr", gfdl1m_namelist),
            RH_INR=_default_or_get_from_namelist(Float(0.25), "rh_ins", gfdl1m_namelist),
            TAU_R2G=_default_or_get_from_namelist(Float(900.0), "tau_r2g", gfdl1m_namelist),
            TAU_SMLT=_default_or_get_from_namelist(Float(900.0), "tau_smlt", gfdl1m_namelist),
            TAU_G2R=_default_or_get_from_namelist(Float(600.0), "tau_g2r", gfdl1m_namelist),
            TAU_IMLT=_default_or_get_from_namelist(Float(600.0), "tau_imlt", gfdl1m_namelist),
            TAU_I2S=_default_or_get_from_namelist(Float(1000.0), "tau_i2s", gfdl1m_namelist),
            TAU_L2R=_default_or_get_from_namelist(Float(900.0), "tau_l2r", gfdl1m_namelist),
            TAU_V2L=_default_or_get_from_namelist(Float(150.0), "tau_v2l", gfdl1m_namelist),
            TAU_L2V=_default_or_get_from_namelist(Float(300.0), "tau_l2v", gfdl1m_namelist),
            TAU_I2V=_default_or_get_from_namelist(Float(300.0), "tau_i2v", gfdl1m_namelist),
            TAU_S2V=_default_or_get_from_namelist(Float(600.0), "tau_s2v", gfdl1m_namelist),
            TAU_V2S=_default_or_get_from_namelist(Float(21600.0), "tau_v2s", gfdl1m_namelist),
            TAU_G2V=_default_or_get_from_namelist(Float(900.0), "tau_g2v", gfdl1m_namelist),
            TAU_V2G=_default_or_get_from_namelist(Float(21600.0), "tau_v2g", gfdl1m_namelist),
            TAU_REVP=_default_or_get_from_namelist(Float(600.0), "tau_revp", gfdl1m_namelist),
            TAU_FRZ=_default_or_get_from_namelist(Float(450.0), "tau_revp", gfdl1m_namelist),
            DW_LAND=_default_or_get_from_namelist(Float(0.20), "dw_land", gfdl1m_namelist),
            DW_OCEAN=_default_or_get_from_namelist(Float(0.10), "dw_ocean", gfdl1m_namelist),
            CCN_O=_default_or_get_from_namelist(Float(90.0), "ccn_o", gfdl1m_namelist),
            CCN_L=_default_or_get_from_namelist(Float(270.0), "ccn_l", gfdl1m_namelist),
            RTHRESHU=_default_or_get_from_namelist(Float(7.0e-6), "rthreshu", gfdl1m_namelist),
            RTHRESHS=_default_or_get_from_namelist(Float(10.0e-6), "rthreshs", gfdl1m_namelist),
            SAT_ADJ0=_default_or_get_from_namelist(Float(0.90), "sat_adj0", gfdl1m_namelist),
            QC_CRT=_default_or_get_from_namelist(Float(5.0e-8), "qc_crt", gfdl1m_namelist),
            QI_LIM=_default_or_get_from_namelist(Float(1.0), "qi_lim", gfdl1m_namelist),
            QL_MLT=_default_or_get_from_namelist(Float(2.0e-3), "ql_mlt", gfdl1m_namelist),
            QS_MLT=_default_or_get_from_namelist(Float(1.0e-6), "qs_mlt", gfdl1m_namelist),
            QL_GEN=_default_or_get_from_namelist(Float(1.0e-3), "ql_gen", gfdl1m_namelist),
            QI_GEN=_default_or_get_from_namelist(Float(9.82679e-5), "qi_gen", gfdl1m_namelist),
            QL0_MAX=_default_or_get_from_namelist(Float(2.0e-3), "ql0_max", gfdl1m_namelist),
            QI0_MAX=_default_or_get_from_namelist(Float(1.0e-4), "qi0_max", gfdl1m_namelist),
            QI0_CRT=_default_or_get_from_namelist(Float(1.0e-4), "qi0_crt", gfdl1m_namelist),
            QR0_CRT=_default_or_get_from_namelist(Float(1.0e-4), "qr0_crt", gfdl1m_namelist),
            QS0_CRT=_default_or_get_from_namelist(Float(1.0e-3), "qs0_crt", gfdl1m_namelist),
            C_PAUT=_default_or_get_from_namelist(Float(0.55), "c_paut", gfdl1m_namelist),
            C_PSACI=_default_or_get_from_namelist(Float(0.02), "c_psaci", gfdl1m_namelist),
            C_PIACR=_default_or_get_from_namelist(Float(5.0), "c_piacr", gfdl1m_namelist),
            C_CRACW=_default_or_get_from_namelist(Float(0.9), "c_cracw", gfdl1m_namelist),
            C_PGACS=_default_or_get_from_namelist(Float(2.0e-3), "c_pgacs", gfdl1m_namelist),
            C_PGACI=_default_or_get_from_namelist(Float(0.05), "c_pgaci", gfdl1m_namelist),
            ALIN=_default_or_get_from_namelist(Float(842.0), "alin", gfdl1m_namelist),
            CLIN=_default_or_get_from_namelist(Float(4.8), "clin", gfdl1m_namelist),
            CONST_VI=_default_or_get_from_namelist(False, "const_vi", gfdl1m_namelist),
            CONST_VS=_default_or_get_from_namelist(False, "const_vs", gfdl1m_namelist),
            CONST_VG=_default_or_get_from_namelist(False, "const_vg", gfdl1m_namelist),
            CONST_VR=_default_or_get_from_namelist(False, "const_vr", gfdl1m_namelist),
            VI_FAC=_default_or_get_from_namelist(Float(1.0), "vi_fac", gfdl1m_namelist),
            VS_FAC=_default_or_get_from_namelist(Float(1.0), "vs_fac", gfdl1m_namelist),
            VR_FAC=_default_or_get_from_namelist(Float(1.0), "vr_fac", gfdl1m_namelist),
            VG_FAC=_default_or_get_from_namelist(Float(1.0), "vg_fac", gfdl1m_namelist),
            VI_MAX=_default_or_get_from_namelist(Float(1.0), "vi_max", gfdl1m_namelist),
            VS_MAX=_default_or_get_from_namelist(Float(2.0), "vs_max", gfdl1m_namelist),
            VG_MAX=_default_or_get_from_namelist(Float(12.0), "vg_max", gfdl1m_namelist),
            VR_MAX=_default_or_get_from_namelist(Float(12.0), "vr_max", gfdl1m_namelist),
            FAST_SAT_ADJ=_default_or_get_from_namelist(False, "fast_sat_adj", gfdl1m_namelist),
            Z_SLOPE_LIQ=_default_or_get_from_namelist(True, "z_slope_liq", gfdl1m_namelist),
            Z_SLOPE_ICE=_default_or_get_from_namelist(False, "z_slope_ice", gfdl1m_namelist),
            USE_CCN=_default_or_get_from_namelist(False, "use_ccn", gfdl1m_namelist),
            USE_PPM=_default_or_get_from_namelist(False, "use_ppm", gfdl1m_namelist),
            MONO_PROF=_default_or_get_from_namelist(True, "mono_prof", gfdl1m_namelist),
            MP_PRINT=_default_or_get_from_namelist(False, "mp_print", gfdl1m_namelist),
            DE_ICE=_default_or_get_from_namelist(False, "de_ice", gfdl1m_namelist),
            SEDI_TRANSPORT=_default_or_get_from_namelist(False, "sedi_transport", gfdl1m_namelist),
            DO_SEDI_W=_default_or_get_from_namelist(False, "do_sedi_w", gfdl1m_namelist),
            DO_SEDI_HEAT=_default_or_get_from_namelist(False, "de_sedi_heat", gfdl1m_namelist),
            PROG_CCN=_default_or_get_from_namelist(False, "prog_ccn", gfdl1m_namelist),
            DO_BIGG=_default_or_get_from_namelist(False, "do_bigg", gfdl1m_namelist),
            DO_EVAP=do_evap,
            DO_SUBL=do_subl,
            DO_QA=_default_or_get_from_namelist(False, "do_qa", gfdl1m_namelist),
            PRECIPRAD=_default_or_get_from_namelist(True, "preciprad", gfdl1m_namelist),
            FIX_NEGATIVE=_default_or_get_from_namelist(False, "fix_negative", gfdl1m_namelist),
            ICLOUD_F=_default_or_get_from_namelist(Int(0), "icloud_f", gfdl1m_namelist),
            IRAIN_F=_default_or_get_from_namelist(Int(0), "irain_f", gfdl1m_namelist),
        )

    def run(
        self,
        mapl_state: CVoidPointer,
        import_state: CVoidPointer,
        export_state: CVoidPointer,
    ):
        pass

    def run_with_internal(
        self,
        mapl_state: CVoidPointer,
        import_state: CVoidPointer,
        export_state: CVoidPointer,
        internal_state: CVoidPointer,
    ):
        ndsl_stack = get_NDSL_physics(mapl_state)
        import_repository = MAPLMemoryRepository(import_state, ndsl_stack.quantity_factory)
        internal_repository = MAPLMemoryRepository(internal_state, ndsl_stack.quantity_factory)
        export_repository = MAPLMemoryRepository(export_state, ndsl_stack.quantity_factory)

        self._managed_state.register("mixing_ratio.vapor", "Q", internal_repository)
        self._managed_state.register("mixing_ratio.rain", "QRAIN", internal_repository)
        self._managed_state.register("mixing_ratio.snow", "QSNOW", internal_repository)
        self._managed_state.register("mixing_ratio.graupel", "QGRAUPEL", internal_repository)
        self._managed_state.register("mixing_ratio.large_scale_liquid", "QLLS", internal_repository)
        self._managed_state.register("mixing_ratio.large_scale_ice", "QILS", internal_repository)
        self._managed_state.register("mixing_ratio.convective_liquid", "QLCN", internal_repository)
        self._managed_state.register("mixing_ratio.convective_ice", "QICN", internal_repository)
        self._managed_state.register("cloud_fraction.large_scale", "CLLS", internal_repository)
        self._managed_state.register("cloud_fraction.convective", "CLCN", internal_repository)
        self._managed_state.register("concentration.liquid", "NACTL", internal_repository)
        self._managed_state.register("concentration.ice", "NACTI", internal_repository)

        self._managed_state.register_2D("area", "AREA", import_repository)
        self._managed_state.register(
            "p_interface", "PLE", import_repository, dims=[I_DIM, J_DIM, K_INTERFACE_DIM]
        )
        self._managed_state.register(
            "z_interface", "ZLE", import_repository, dims=[I_DIM, J_DIM, K_INTERFACE_DIM]
        )
        self._managed_state.register("t", "T", import_repository)
        self._managed_state.register("u", "U", import_repository)
        self._managed_state.register("v", "V", import_repository)
        self._managed_state.register_2D("land_fraction", "FRLAND", import_repository)
        self._managed_state.register(
            "covariance_liquid_water_static_energy_and_total_water_specific_humidity",
            "SLQT",
            import_repository,
        )
        self._managed_state.register("omega", "OMEGA", import_repository)
        self._managed_state.register("pdf_first_plume_fractional_area", "PDF_A", import_repository)

        self._managed_state.register("vertical_motion.velocity", "W", import_repository)
        self._managed_state.register("vertical_motion.variance", "W2", import_repository)
        self._managed_state.register("vertical_motion.third_moment", "W3", import_repository)

        self._managed_state.register("liquid_water_static_energy.flux", "WSL", import_repository)
        self._managed_state.register("liquid_water_static_energy.variance", "SL2", import_repository)
        self._managed_state.register("liquid_water_static_energy.third_moment", "SL3", import_repository)

        self._managed_state.register("total_water.flux", "WQT", import_repository)
        self._managed_state.register("total_water.variance", "QT2", import_repository)
        self._managed_state.register("total_water.third_moment", "QT3", import_repository)

        # MAPL_GetPointer fails on:
        # self_manage_state.register("surface_temperature", "TS", import_repository)
        # self._manage_state.register(
        #     "scalar_diffusivity_interface",
        #     "KH",
        #     import_repository,
        #     dims=[I_DIM, J_DIM, K_INTERFACE_DIM],
        # )
        # self_manage_state.register("sensible_heat_flux", "SH", import_repository)

        self._managed_state.register_2D("convection_fraction", "CNV_FRC", export_repository, alloc=True)
        self._managed_state.register_2D("surface_type", "SRF_TYPE", export_repository, alloc=True)
        self._managed_state.register("cloud_liquid_evaporation", "EVAPC", export_repository, alloc=True)
        self._managed_state.register("cloud_ice_sublimation", "SUBLC", export_repository, alloc=True)
        self._managed_state.register_2D("icefall", "ICE", export_repository, alloc=True)
        self._managed_state.register_2D("freezing_rainfall", "FRZR", export_repository, alloc=True)
        self._managed_state.register("relative_humidity_after_pdf", "RHX", export_repository, alloc=True)
        self._managed_state.register("buoyancy_flux", "WTHV2", export_repository, alloc=True)
        self._managed_state.register("liquid_water_flux", "WQL", export_repository, alloc=True)
        self._managed_state.register("hydrostatic_pdf_iterations", "PDFITERS", export_repository, alloc=True)
        self._managed_state.register_2D("lower_tropospheric_stability", "LTS", export_repository, alloc=True)
        self._managed_state.register_2D("estimated_inversion_strength", "EIS", export_repository, alloc=True)
        self._managed_state.register_2D("lcl_height", "ZLCL", export_repository)
        self._managed_state.register("shallow_convection_rain", "SHLW_PRC3", export_repository, alloc=True)
        self._managed_state.register("shallow_convection_snow", "SHLW_SNO3", export_repository, alloc=True)
        self._managed_state.register(
            "critical_relative_humidity_for_pdf", "RHCRIT", export_repository, alloc=True
        )
        self._managed_state.register("large_scale_rainwater_source", "DQRL", export_repository)

        self._managed_state.register("radiation_field.cloud_fraction", "FCLD", export_repository, alloc=True)
        self._managed_state.register("radiation_field.vapor", "QV", export_repository, alloc=True)
        self._managed_state.register("radiation_field.liquid", "QL", export_repository, alloc=True)
        self._managed_state.register("radiation_field.ice", "QI", export_repository, alloc=True)
        self._managed_state.register("radiation_field.rain", "QR", export_repository, alloc=True)
        self._managed_state.register("radiation_field.snow", "QS", export_repository, alloc=True)
        self._managed_state.register("radiation_field.graupel", "QG", export_repository, alloc=True)

        self._managed_state.register(
            "cloud_particle_effective_radius.liquid", "RL", export_repository, alloc=True
        )
        self._managed_state.register(
            "cloud_particle_effective_radius.ice", "RI", export_repository, alloc=True
        )

        self._managed_state.register_2D(
            "precipitation_at_surface.rain", "PRCP_RAIN", export_repository, alloc=True
        )
        self._managed_state.register_2D(
            "precipitation_at_surface.snow", "PRCP_SNOW", export_repository, alloc=True
        )
        self._managed_state.register_2D(
            "precipitation_at_surface.ice", "PRCP_ICE", export_repository, alloc=True
        )
        self._managed_state.register_2D(
            "precipitation_at_surface.graupel", "PRCP_GRAUPEL", export_repository, alloc=True
        )
        self._managed_state.register_2D(
            "precipitation_at_surface.shallow_convective_precipitation",
            "SC_PRCP",
            export_repository,
            alloc=True,
        )
        self._managed_state.register_2D(
            "precipitation_at_surface.deep_convective_precipitation", "CN_PRCP", export_repository, alloc=True
        )
        self._managed_state.register_2D(
            "precipitation_at_surface.anvil_precipitation", "AN_PRCP", export_repository, alloc=True
        )
        self._managed_state.register_2D(
            "precipitation_at_surface.shallow_convective_snow", "SC_SNR", export_repository, alloc=True
        )
        self._managed_state.register_2D(
            "precipitation_at_surface.deep_convective_snow", "CN_SNR", export_repository, alloc=True
        )
        self._managed_state.register_2D(
            "precipitation_at_surface.anvil_snow", "AN_SNR", export_repository, alloc=True
        )

        self._managed_state.register_2D(
            "non_anvil_large_scale.precip", "LS_PRCP", export_repository, alloc=True
        )
        self._managed_state.register_2D("non_anvil_large_scale.snow", "LS_SNR", export_repository, alloc=True)
        self._managed_state.register(
            "non_anvil_large_scale.evaporation", "REV_LS", export_repository, alloc=True
        )
        self._managed_state.register(
            "non_anvil_large_scale.sublimation", "RSU_LS", export_repository, alloc=True
        )
        self._managed_state.register(
            "non_anvil_large_scale.liquid_precip_flux",
            "PFL_LS",
            export_repository,
            dims=[I_DIM, J_DIM, K_INTERFACE_DIM],
            alloc=True,
        )
        self._managed_state.register(
            "non_anvil_large_scale.ice_precip_flux",
            "PFI_LS",
            export_repository,
            dims=[I_DIM, J_DIM, K_INTERFACE_DIM],
            alloc=True,
        )

        self._managed_state.register(
            "anvil.liquid_precip_flux",
            "PFL_AN",
            export_repository,
            dims=[I_DIM, J_DIM, K_INTERFACE_DIM],
            alloc=True,
        )
        self._managed_state.register(
            "anvil.ice_precip_flux",
            "PFI_AN",
            export_repository,
            dims=[I_DIM, J_DIM, K_INTERFACE_DIM],
            alloc=True,
        )

        self._managed_state.register(
            "tendencies.dcloud_fractiondt_macro", "DQADT_macro", export_repository, alloc=True
        )
        self._managed_state.register(
            "tendencies.dvapordt_macro", "DQVDT_macro", export_repository, alloc=True
        )
        self._managed_state.register("tendencies.dicedt_macro", "DQIDT_macro", export_repository, alloc=True)
        self._managed_state.register(
            "tendencies.dliquiddt_macro", "DQLDT_macro", export_repository, alloc=True
        )
        self._managed_state.register("tendencies.draindt_macro", "DQRDT_macro", export_repository, alloc=True)
        self._managed_state.register(
            "tendencies.dgraupeldt_macro", "DQGDT_macro", export_repository, alloc=True
        )
        self._managed_state.register("tendencies.dsnowdt_macro", "DQSDT_macro", export_repository, alloc=True)
        self._managed_state.register("tendencies.dudt_macro", "DUDT_macro", export_repository, alloc=True)
        self._managed_state.register("tendencies.dvdt_macro", "DVDT_macro", export_repository, alloc=True)
        self._managed_state.register("tendencies.dtdt_macro", "DTDT_macro", export_repository, alloc=True)
        self._managed_state.register(
            "tendencies.dcloud_fractiondt_micro", "DQADT_micro", export_repository, alloc=True
        )
        self._managed_state.register(
            "tendencies.dvapordt_micro", "DQVDT_micro", export_repository, alloc=True
        )
        self._managed_state.register("tendencies.dicedt_micro", "DQIDT_micro", export_repository, alloc=True)
        self._managed_state.register(
            "tendencies.dliquiddt_micro", "DQLDT_micro", export_repository, alloc=True
        )
        self._managed_state.register("tendencies.draindt_micro", "DQRDT_micro", export_repository, alloc=True)
        self._managed_state.register(
            "tendencies.dgraupeldt_micro", "DQGDT_micro", export_repository, alloc=True
        )
        self._managed_state.register("tendencies.dsnowdt_micro", "DQSDT_micro", export_repository, alloc=True)
        self._managed_state.register("tendencies.dudt_micro", "DUDT_micro", export_repository, alloc=True)
        self._managed_state.register("tendencies.dvdt_micro", "DVDT_micro", export_repository, alloc=True)
        self._managed_state.register("tendencies.dtdt_micro", "DTDT_micro", export_repository, alloc=True)
        self._managed_state.register(
            "tendencies.dtdt_friction_pressure_weighted", "DTDTFRIC", export_repository
        )

        self._managed_state.register("radar.simulated_reflectivity", "DBZ", export_repository)
        self._managed_state.register_2D("radar.maximum_composite_reflectivity", "DBZ_MAX", export_repository)
        self._managed_state.register_2D("radar.base_1km_agl_reflectivity", "DBZ_1KM", export_repository)
        self._managed_state.register_2D("radar.echo_top_reflectivity", "DBZ_TOP", export_repository)
        self._managed_state.register_2D("radar.minus_10c_reflectivity", "DBZ_M10C", export_repository)

        print(f"WRITING DATA TO GFDL1M_Out")
        self._managed_state.fortran_to_ndsl()
        self._managed_state.record("GFDL1M-Out_fortran")

    def finalize(
        self,
        mapl_state: CVoidPointer,
        import_state: CVoidPointer,
        export_state: CVoidPointer,
    ):
        print(f"SAVING DATA to {os.getcwd()}")
        self._managed_state.save_recorded()


CODE = GFDL1M_save_end()
