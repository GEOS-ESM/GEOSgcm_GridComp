from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField, Int, Float 
from ndsl import Quantity, QuantityFactory, StencilFactory
import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import computation, interval, PARALLEL, log10, exp  # type: ignore
import pyMoist.constants as aerconstants

def _aer_activation_stencil():
    with computation(PARALLEL), interval(...):
        pass
def _get_act_frac_stencil():
    with computation(PARALLEL), interval(...):
        pass
def _act_frac_mat_stencil():
    with computation(PARALLEL), interval(...):
        pass
def _gser_stencil():
    with computation(PARALLEL), interval(...):
        pass
def _gcf_matrix_stencil():
    with computation(PARALLEL), interval(...):
        pass
class AerActivation:
    def __init__(
            self,
            stencil_factory: StencilFactory,
            quantity_factory: QuantityFactory,
            do_qa: bool,
    ):
        self._get_act_frac = stencil_factory.from_dims_halo(
                func=_get_act_frac_stencil,
                compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )
    def __call__(
        self,
        IM,
        JM,
        LM,
        q,
        t,
        plo,
        ple,
        zlo,
        zle,
        qlcn,
        qicn,
        qlls,
        qils,
        sh,
        evap,
        kpbl,
        tke,
        vvel,
        FRLAND,
        USE_AERO_BUFFER,
        AeroProps,
        aero_aci,
        NACTL,
        NACTI,
        NWFA,
        NN_LAND,
        NN_OCEAN

    ):
        self._get_act_frac_stencil(

        )
        