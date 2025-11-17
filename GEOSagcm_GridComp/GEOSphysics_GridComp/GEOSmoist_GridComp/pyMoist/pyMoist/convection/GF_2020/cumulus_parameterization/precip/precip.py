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
from pyMoist.convection.GF_2020.cumulus_parameterization.precip.stencils import (
    partition_liquid_ice,
    get_precip_fluxes,
    output_evaporation_flux,
    output_deep_precipitation,
)


class PartitionLiquidIce:
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
        self._partition_liquid_ice = stencil_factory.from_dims_halo(
            func=partition_liquid_ice,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "MELT_ICE": cumulus_parameterization_config.MELT_ICE,
                "MODIS_FRACTION": cumulus_parameterization_config.MODIS_FRACTION,
            },
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._partition_liquid_ice(
            t=locals.t_new,
            p=state.output.p_cloud_levels_forced,
            geopotential_height=locals.geopotential_height_cloud_levels_forced,
            topography_height_no_negative=state.input_output.topography_height_no_negative,
            surface_type=state.input.surface_type,
            convection_fraction=state.input.convection_fraction,
            error_code=state.output.error_code,
            melting_layer=locals.melting_layer,
            part_liquid_ice=locals.partition_liquid_ice,
            plume=plume_dependent_constants.PLUME_INDEX,
        )


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
        #         "MELT_ICE": cumulus_parameterization_config.MELT_ICE,
        #         "MODIS_FRACTION": cumulus_parameterization_config.MODIS_FRACTION,
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
