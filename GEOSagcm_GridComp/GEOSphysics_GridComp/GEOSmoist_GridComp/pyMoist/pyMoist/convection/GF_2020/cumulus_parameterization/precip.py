from ndsl import StencilFactory, QuantityFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import (
    GF2020CumulusParameterizationConfig,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.state import (
    GF2020CumulusParameterizationState,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.locals import (
    GF2020CumulusParameterizationLocals,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)
from ndsl.logging import ndsl_log
from ndsl.dsl.gt4py import (
    PARALLEL,
    computation,
    interval,
    FORWARD,
    function,
    BACKWARD,
    K,
)
from ndsl.dsl.typing import FloatField, FloatFieldIJ, IntField, Int, IntFieldIJ
import pyMoist.constants as constants
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import (
    FloatField_Plume,
    IntFieldIJ_Plume,
)
from pyMoist.shared_incloud_processes import ice_fraction


@function
def liquid_fraction(
    t,
    convection_fraction,
    surface_type,
    FRAC_MODIS,
):
    """
    Get the fraction of liquid condensates

    Args:
        t (in): temperature
        convection_fraction (in)
        surface_type (in)
        FRAC_MODIS (in): use fraction liq/ice content derived from MODIS/CALIPO sensors
    """
    if FRAC_MODIS == 1:
        liquid_fraction = 1.0 - ice_fraction(t, convection_fraction, surface_type)
    else:
        liquid_fraction = min(
            1.0,
            (
                max(0.0, (t - cumulus_parameterization_constants.T_ICE))
                / (cumulus_parameterization_constants.T_0 - cumulus_parameterization_constants.T_ICE)
            )
            ** 2,
        )

    return liquid_fraction


def partition_liquid_ice(
    t: FloatField,
    p: FloatField_Plume,
    geopotential_height: FloatField,
    topography_height_no_negative: FloatFieldIJ,
    surface_type: FloatFieldIJ,
    convection_fraction: FloatFieldIJ,
    error_code: IntFieldIJ_Plume,
    melting_layer: FloatField,
    part_liquid_ice: FloatField,
    plume: Int,
):
    """
    Partition total condensate into liquid and ice phases

    Args:

    """
    from __externals__ import MELT_GLAC, FRAC_MODIS, k_end

    with computation(PARALLEL), interval(...):
        # constants, set internally because they may differ from global constants
        # and need to only exist inside this stencil
        t1 = 276.16
        z_meltlayer1 = 4000.0
        z_meltlayer2 = 6000.0
        del_t = 3.0

        # prefill some fields
        part_liquid_ice = 1.0
        melting_layer = 0.0

    with computation(PARALLEL), interval(0, -1):
        if MELT_GLAC == True and plume == 2:
            if error_code[0, 0][plume] == 0:
                # get function of T for partition of total condensate into liq and ice phases
                part_liquid_ice = liquid_fraction(t, convection_fraction, surface_type, FRAC_MODIS)

    with computation(PARALLEL), interval(0, -1):
        if MELT_GLAC == True and plume == 2:
            if error_code[0, 0][plume] == 0:
                # define the melting layer (the layer will be between T_0+1 < TEMP < T_1
                if t <= (cumulus_parameterization_constants.T_0 - del_t):
                    melting_layer = 0.0

                elif t < (cumulus_parameterization_constants.T_0 + del_t) and t > (
                    cumulus_parameterization_constants.T_0 - del_t
                ):
                    melting_layer = (
                        (t - (cumulus_parameterization_constants.T_0 - del_t)) / (2.0 * del_t)
                    ) ** 2

                else:
                    melting_layer = 1.0

                melting_layer = melting_layer * (1.0 - melting_layer)

    with computation(FORWARD), interval(0, 1):
        if MELT_GLAC == True and plume == 2:
            # normalize vertical integral of melting_layer to 1
            norm: FloatFieldIJ = 0.0

    with computation(FORWARD), interval(...):
        if MELT_GLAC == True and plume == 2:
            if error_code[0, 0][plume] == 0:
                # normalize vertical integral of melting_layer to 1
                dp = 100.0 * (p[0, 0, 0][plume] - p[0, 0, 1][plume])
                norm = norm + melting_layer * dp / constants.MAPL_GRAV

    with computation(PARALLEL), interval(...):
        if MELT_GLAC == True and plume == 2:
            if error_code[0, 0][plume] == 0:
                # normalize vertical integral of melting_layer to 1
                melting_layer = (
                    melting_layer
                    / (norm + 1.0e-6)
                    * (100 * (p[0, 0, 0][plume] - p.at(K=k_end - 1, ddim=[plume])) / constants.MAPL_GRAV)
                )


def get_precip_fluxes(
    edto: FloatField,
    error_code: IntFieldIJ_Plume,
    plume: Int,
    ktop: IntField,
    pwdo: FloatField,
    pwo: FloatField,
    xmb: FloatField,
    prec_flx: FloatField,
    evap_flx: FloatField,
):
    with computation(BACKWARD), interval(...):
        if error_code[0, 0][plume] == 0:
            if K < ktop:
                prec_flx = prec_flx[0, 0, 1] + xmb * (pwo + edto * pwdo)
                prec_flx = max(0.0, prec_flx)

                evap_flx = evap_flx[0, 0, 1] - xmb * edto * pwdo
                evap_flx = max(0.0, evap_flx)


def output_evaporation_flux(
    error_code: IntFieldIJ_Plume,
    plume: Int,
    ktop: IntFieldIJ,
    po_cup: FloatField,
    evap_flx: FloatField,
    revsu_gf: FloatField,
):
    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0:
            if K <= ktop:
                dp: FloatFieldIJ = 100.0 * (po_cup - po_cup.at(K=K + 1))
                revsu_gf = revsu_gf + evap_flx * constants.MAPL_GRAV / dp


def output_deep_precipitation(
    cumulus: Int,
    error_code: IntFieldIJ_Plume,
    plume: Int,
    ktop: IntFieldIJ,
    prec_flx: FloatField,
    prfil_gf: FloatField,
):
    with computation(PARALLEL), interval(...):
        if cumulus == cumulus_parameterization_constants.deep:
            if error_code[0, 0][plume] == 0:
                if K <= ktop + 1:
                    prfil_gf = prec_flx


class PrecipFactor:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GF2020Config,
        cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
    ):
        pass
        # # make configuration visible at runtime
        # self.config = config
        # self.cumulus_parameterization_config = cumulus_parameterization_config

        # # construct stencils and functions
        # self._ = stencil_factory.from_dims_halo(
        #     func=,
        #     compute_dims=[X_DIM, Y_DIM, Z_DIM],
        #     externals={
        #         "MELT_GLAC": cumulus_parameterization_config.MELT_GLAC,
        #         "FRAC_MODIS": cumulus_parameterization_config.FRAC_MODIS,
        #     },
        # )

    def __call__(self, *args, **kwds):
        pass


class PrecipitationFlux:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GF2020Config,
        cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
    ):
        # make configuration visible at runtime
        self.config = config
        self.cumulus_parameterization_config = cumulus_parameterization_config

        # construct stencils and functions
        self._get_precip_fluxes = stencil_factory.from_dims_halo(
            func=get_precip_fluxes,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._get_precip_fluxes(
            # edto=,
            # error_code=,
            # plume=,
            # ktop=,
            # pwdo=,
            # pwo=,
            # xmb=,
            # prec_flx=,
            # evap_flx=,
        )


class RainEvapBelowCloudBase:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass


class CloudDissipation:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass


class OutputEvaporationFlux:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GF2020Config,
        cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
    ):
        # make configuration visible at runtime
        self.config = config
        self.cumulus_parameterization_config = cumulus_parameterization_config

        # construct stencils and functions
        self._output_evaporation_flux = stencil_factory.from_dims_halo(
            func=output_evaporation_flux,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._output_evaporation_flux(
            # error_code=,
            # plume=,
            # ktop=,
            # po_cup=,
            # evap_flx=,
            # revsu_gf=,
        )


class LightningFlassDensity:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass


class OutputDeepPrecipitation:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GF2020Config,
        cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
    ):
        # make configuration visible at runtime
        self.config = config
        self.cumulus_parameterization_config = cumulus_parameterization_config

        # construct stencils and functions
        self._output_deep_precipitation = stencil_factory.from_dims_halo(
            func=output_deep_precipitation,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._output_deep_precipitation(
            # cumulus=,
            # error_code=,
            # plume=,
            # ktop=,
            # prec_flx=,
            # prfil_gf=,
        )


class UpdateWorkfunctionsAndCondensates:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass
