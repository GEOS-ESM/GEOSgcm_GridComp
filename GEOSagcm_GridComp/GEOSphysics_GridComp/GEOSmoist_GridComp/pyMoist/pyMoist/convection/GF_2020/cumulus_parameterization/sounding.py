from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.constants import WRTGRADS

class Sounding:
    def __init__(self, cumulus_parameterization_config: GF2020CumulusParameterizationConfig):
        if cumulus_parameterization_config.OUTPUT_SOUNDING == 1:
            raise NotImplementedError(
                "[NDSL] GF2020-->CumulusParameterization-->Sounding: Output Sounding capabilities have not"
                "been implemented. You should have been caught before getting here by the config checker."
                "Beware, something likely failing in the config checker as well - you may be unknowingly"
                "calling other untested/unimplemented sections."
            )

    def __call__(self, *args, **kwds):
        pass


class GATESounding:
    def __init__(self, cumulus_parameterization_config: GF2020CumulusParameterizationConfig):
        if WRTGRADS == 1:
            raise NotImplementedError(
                "[NDSL] GF2020-->CumulusParameterization-->GATESounding: GATE Sounding capabilities have not"
                "been implemented. You should have been caught before getting here by the config checker."
                "Beware, something likely failing in the config checker as well - you may be unknowingly"
                "calling other untested/unimplemented sections."
            )

    def __call__(self, *args, **kwds):
        pass
