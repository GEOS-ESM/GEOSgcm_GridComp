"""All translate tests must be imported here to be automatically discoverable by `pytest"""

from .UW.translate_compute_uwshcu import TranslateComputeUwshcuInv
from .UW.UW_translate_tests.translate_setup_inputs import TranslateSetupInputs
from .UW.UW_translate_tests.translate_prepare_inputs import TranslatePrepareInputs
from .UW.UW_translate_tests.translate_find_pbl import TranslateFindPbl
from .UW.UW_translate_tests.translate_find_klcl import TranslateFindKlcl
from .UW.UW_translate_tests.translate_compute_cin_cinlcl import TranslateComputeCinCinlcl
from .UW.UW_translate_tests.translate_define_prel_cbmf import TranslateDefinePrelCbmf
from .UW.UW_translate_tests.translate_define_updraft_properties import TranslateDefineUpdraftProperties
from .UW.UW_translate_tests.translate_define_env_properties import TranslateDefineEnvProperties
from .UW.UW_translate_tests.translate_buoyancy_sorting import TranslateBuoyancySorting

__all__ = ["TranslateComputeUwshcuInv"]
