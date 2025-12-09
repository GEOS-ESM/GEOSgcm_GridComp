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
from pyMoist.convection.GF_2020.cumulus_parameterization.moist_static_energy.stencils import (
    parcel_moist_static_energy,
    first_guess_moist_static_energy,
    moist_static_energy_inside_cloud,
)
from ndsl.logging import ndsl_log

# from pyMoist.convection.GF_2020.cumulus_parameterization.buoyancy.stencils import (
#     get_buoyancy,
# )


class ParcelMoistStaticEnergy:
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
        self._parcel_moist_static_energy = stencil_factory.from_dims_halo(
            func=parcel_moist_static_energy,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "BOUNDARY_CONDITION_METHOD": cumulus_parameterization_config.BOUNDARY_CONDITION_METHOD,
            },
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):

        if self.cumulus_parameterization_config.FRAC_MODIS != 1:
            ndsl_log.warning(
                " GF2020 cumulus parameterization called ParcelMoistStaticEnergy with "
                "untested BOUNDARY_CONDITION_METHOD option. Running untested code... proceed with caution"
            )

        self._parcel_moist_static_energy(
            error_code=state.output.error_code,
            t_excess=locals.t_excess,
            vapor_excess=locals.vapor_excess,
            add_buoyancy=locals.add_buoyancy,
            ocean_fraction=state.input.ocean_fraction,
            updraft_origin_level=state.output.updraft_origin_level,
            p=state.input_output.p_forced,
            environmenet_moist_static_energy=locals.environment_moist_static_energy_cloud_levels,
            environmenet_moist_static_energy_forced=locals.environment_moist_static_energy_cloud_levels_forced,
            t_perturbation=state.output.t_perturbation,
            moist_static_energy_origin_level=locals.moist_static_energy_origin_level,
            moist_static_energy_origin_level_forced=locals.moist_static_energy_origin_level_forced,
            AVERAGE_LAYER_DEPTH=plume_dependent_constants.AVERAGE_LAYER_DEPTH,
            plume=plume_dependent_constants.PLUME_INDEX,
        )


class UpdateMoistStaticEnergy:
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
        self._parcel_moist_static_energy = stencil_factory.from_dims_halo(
            func=parcel_moist_static_energy,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "BOUNDARY_CONDITION_METHOD": cumulus_parameterization_config.BOUNDARY_CONDITION_METHOD,
            },
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._parcel_moist_static_energy(
            error_code=state.output.error_code,
            t_excess=locals.t_excess,
            vapor_excess=locals.vapor_excess,
            add_buoyancy=locals.add_buoyancy,
            ocean_fraction=state.input.ocean_fraction,
            updraft_origin_level=state.output.updraft_origin_level,
            p=state.input_output.p_forced,
            environmenet_moist_static_energy=locals.environment_moist_static_energy_cloud_levels,
            environmenet_moist_static_energy_forced=locals.environment_moist_static_energy_cloud_levels_forced,
            t_perturbation=state.output.t_perturbation,
            moist_static_energy_origin_level=locals.moist_static_energy_origin_level,
            moist_static_energy_origin_level_forced=locals.moist_static_energy_origin_level_forced,
            AVERAGE_LAYER_DEPTH=plume_dependent_constants.AVERAGE_LAYER_DEPTH,
            plume=plume_dependent_constants.PLUME_INDEX,
        )


class FirstGuessMoistStaticEnergy:
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
        self._first_guess_mse = stencil_factory.from_dims_halo(
            func=first_guess_moist_static_energy,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._first_guess_mse(
            error_code=state.output.error_code,
            start_level=locals.start_level,
            cloud_top_level=state.output.cloud_top_level,
            mass_detrainment_updraft_forced=state.output.mass_detrainment_updraft_forced,
            mass_entrainment_updraft_forced=state.output.mass_entrainment_updraft_forced,
            normalized_massflux_updraft=locals.normalized_massflux_updraft,
            normalized_massflux_updraft_forced=state.output.normalized_massflux_updraft_forced,
            environment_moist_static_energy_forced=locals.environment_moist_static_energy_forced,
            environment_saturation_moist_static_energy_cloud_levels_forced=locals.environment_saturation_moist_static_energy_cloud_levels_forced,
            cloud_moist_static_energy_forced=locals.cloud_moist_static_energy_forced,
            vapor_excess=locals.vapor_excess,
            t_excess=locals.t_excess,
            add_buoyancy=locals.add_buoyancy,
            plume=plume_dependent_constants.PLUME_INDEX,
        )


class MoistStaticEnergyInsideCloud:
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
        self._moist_static_energy_inside_cloud = stencil_factory.from_dims_halo(
            func=moist_static_energy_inside_cloud,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._moist_static_energy_inside_cloud(
            error_code=state.output.error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
            # start_level=,
            # xhc=,
            # xhkb=,
            # ktop=,
            # up_massdetro=,
            # up_massentro=,
            # xzu=,
            # xhe=,
            # p_liq_ice=,
            # zqexec=,
            # ztexec=,
            # x_add_buoy=,
            # qrco=,
            # xhes_cup=,
        )

        # self._get_buoyancy(
        # hc: FloatField,
        # he_cup: FloatField,
        # hes_cup: FloatField,
        # error_code: IntField,
        # kbcon: IntField,
        # klcl: IntField,
        # ktop: IntField,
        # dby: FloatField,
        # )
