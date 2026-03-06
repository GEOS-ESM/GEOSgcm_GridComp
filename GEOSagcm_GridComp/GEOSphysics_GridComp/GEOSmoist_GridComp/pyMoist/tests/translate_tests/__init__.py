"""All translate tests must be imported here to be automatically discoverable by `pytest"""

from .GFDL_1M.driver.translate_GFDL_1M_Driver import TranslateGFDL_1M_Driver
from .GFDL_1M.driver.translate_GFDL_1M_DriverFinish import TranslateGFDL_1M_DriverFinish
from .GFDL_1M.driver.translate_GFDL_1M_DriverSetup import TranslateGFDL_1M_DriverSetup
from .GFDL_1M.driver.translate_GFDL_1M_FallSpeed import TranslateGFDL_1M_FallSpeed
from .GFDL_1M.driver.translate_GFDL_1M_IceCloud import TranslateGFDL_1M_IceCloud
from .GFDL_1M.driver.translate_GFDL_1M_TerminalFall import TranslateGFDL_1M_TerminalFall
from .GFDL_1M.driver.translate_GFDL_1M_WarmRain import TranslateGFDL_1M_WarmRain
from .GFDL_1M.driver.translate_GFDL_driver_tables import TranslateGFDL_driver_tables
from .GFDL_1M.PhaseChange.translate_GFDL_1M_Evaporate import TranslateGFDL_1M_Evaporate
from .GFDL_1M.PhaseChange.translate_GFDL_1M_HydrostaticPDF import TranslateGFDL_1M_HydrostaticPDF
from .GFDL_1M.PhaseChange.translate_GFDL_1M_MeltFreeze import TranslateGFDL_1M_MeltFreeze
from .GFDL_1M.PhaseChange.translate_GFDL_1M_PhaseChange import TranslateGFDL_1M_PhaseChange
from .GFDL_1M.PhaseChange.translate_GFDL_1M_RHCalculations import TranslateGFDL_1M_RHCalculations
from .GFDL_1M.PhaseChange.translate_GFDL_1M_Sublimate import TranslateGFDL_1M_Sublimate
from .GFDL_1M.translate_GFDL_1M_Finalize import TranslateGFDL_1M_Finalize
from .GFDL_1M.translate_GFDL_1M_RadiationCoupling import TranslateGFDL_1M_RadiationCoupling
from .GFDL_1M.translate_GFDL_1M_RedistributeClouds import TranslateGFDL_1M_RedistributeClouds
from .GFDL_1M.translate_GFDL_1M_Setup import TranslateGFDL_1M_Setup
from .UW.translate_compute_uwshcu import TranslateComputeUwshcuInv

from .UW.UW_translate_tests.translate_setup_inputs import TranslateSetupInputs
from .UW.UW_translate_tests.translate_prepare_inputs import TranslatePrepareInputs
from .UW.UW_translate_tests.translate_find_pbl import TranslateFindPbl
from .UW.UW_translate_tests.translate_find_klcl import TranslateFindKlcl
from .UW.UW_translate_tests.translate_compute_cin_cinlcl import TranslateComputeCinCinlcl
from .UW.UW_translate_tests.translate_compute_del_CIN import TranslateComputeDelCIN
from .UW.UW_translate_tests.translate_average_initial_and_final_CIN1 import TranslateAverageInitialFinalCIN1
from .UW.UW_translate_tests.translate_average_initial_and_final_CIN3 import TranslateAverageInitialFinalCIN3
from .UW.UW_translate_tests.translate_define_prel_cbmf import TranslateDefinePrelCbmf
from .UW.UW_translate_tests.translate_define_updraft_properties import TranslateDefineUpdraftProperties
from .UW.UW_translate_tests.translate_define_env_properties import TranslateDefineEnvProperties
from .UW.UW_translate_tests.translate_buoyancy_sorting import TranslateBuoyancySorting
from .UW.UW_translate_tests.translate_calc_ppen import TranslateCalcPpen
from .UW.UW_translate_tests.translate_recalc_condensate import TranslateRecalcCondensate
from .UW.UW_translate_tests.translate_calc_entrainment_mass_flux import TranslateCalcEntrainmentMassFlux
from .UW.UW_translate_tests.translate_calc_pbl_fluxes import TranslateCalcPblFluxes
from .UW.UW_translate_tests.translate_buoyancy_sorting_fluxes import TranslateBuoyancySortingFluxes
from .UW.UW_translate_tests.translate_penetrative_entrainment_fluxes import (
    TranslatePenetrativeEntrainmentFluxes,
)
from .UW.UW_translate_tests.translate_calc_momentum_tendency import TranslateMomentumTendency
from .UW.UW_translate_tests.translate_prevent_negative_condensate import TranslatePreventNegativeCondensate
from .UW.UW_translate_tests.translate_calc_thermodynamic_tendencies import TranslateThermodynamicTendencies
from .UW.UW_translate_tests.translate_tracer_tendencies import TranslateTracerTendencies
from .UW.UW_translate_tests.translate_compute_diagnostic_outputs import TranslateComputeDiagnosticOutputs
from .UW.UW_translate_tests.translate_calc_cumulus_condensate_at_interface import (
    TranslateCalcCumulusCondensate,
)
from .UW.UW_translate_tests.translate_adjust_implicit_CIN_inputs1 import TranslateAdjustImplicitCINInputs1
from .UW.UW_translate_tests.translate_adjust_implicit_CIN_inputs2 import TranslateAdjustImplicitCINInputs2
from .UW.UW_translate_tests.translate_recalc_environmental_variables import TranslateRecalcEnvVariables
from .UW.UW_translate_tests.translate_update_output_variables1 import TranslateUpdateOutputVars1
from .UW.UW_translate_tests.translate_compute_uwshcu_invert_after import TranslateComputeUwshcuInvertAfter

__all__ = [
    "TranslateGFDL_1M_Driver",
    "TranslateGFDL_1M_DriverFinish",
    "TranslateGFDL_1M_DriverSetup",
    "TranslateGFDL_1M_FallSpeed",
    "TranslateGFDL_1M_IceCloud",
    "TranslateGFDL_1M_TerminalFall",
    "TranslateGFDL_1M_WarmRain",
    "TranslateGFDL_driver_tables",
    "TranslateGFDL_1M_Evaporate",
    "TranslateGFDL_1M_HydrostaticPDF",
    "TranslateGFDL_1M_MeltFreeze",
    "TranslateGFDL_1M_PhaseChange",
    "TranslateGFDL_1M_RHCalculations",
    "TranslateGFDL_1M_Sublimate",
    "TranslateGFDL_1M_Finalize",
    "TranslateGFDL_1M_RadiationCoupling",
    "TranslateGFDL_1M_RedistributeClouds",
    "TranslateGFDL_1M_Setup",
    "TranslateComputeUwshcuInv",
]
