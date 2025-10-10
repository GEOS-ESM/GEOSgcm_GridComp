from ndsl import StencilFactory, QuantityFactory, ndsl_log
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.state import GF2020State
from pyMoist.convection.GF_2020.locals import GF2020Locals
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.convection.GF_2020.cumulus_parameterization.stencils import prefil_excess, set_new_t_vapor
from ndsl.constants import X_DIM, Y_DIM, Z_DIM


class CumulusParameterization:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        GF_2020_config: GF2020Config,
        cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
    ):
        self.GF_2020_config = GF_2020_config
        self.cu_param_config = cumulus_parameterization_config

        self._prefil_excess = stencil_factory.from_dims_halo(
            func=prefil_excess,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._set_new_t_vapor = stencil_factory.from_dims_halo(
            func=set_new_t_vapor,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"DT_MOIST": GF_2020_config.DT_MOIST},
        )

    def __call__(
        self,
        state: GF2020State,
        locals: GF2020Locals,
        saturation_tables: SaturationVaporPressureTable,
    ):
        if self.cu_param_config.PLUME_ORDER == 0:
            plume_types = ["shallow", "mid", "deep"]
            # ndsl_log.log(msg="running cumulus parameterization with order shallow, mid, deep")
        elif self.cu_param_config.PLUME_ORDER == 1:
            plume_types = ["shallow", "deep", "mid"]
            # ndsl_log.log(msg="running cumulus parameterization with order shallow, deep, mid")
        else:
            raise NotImplementedError("plume order not impelemented")

        for plume in plume_types:
            # temporary implementation, need a more elegant solution
            if plume == "shallow":
                downdraft_max_height_land = self.cu_param_config.DOWNDRAFT_MAX_HEIGHT_LAND_SHALLOW
                downdraft_max_height_ocean = self.cu_param_config.DOWNDRAFT_MAX_HEIGHT_OCEAN_SHALLOW
                updraft_max_height_land = self.cu_param_config.UPDRAFT_MAX_HEIGHT_LAND_SHALLOW
                updraft_max_height_ocean = self.cu_param_config.UPDRAFT_MAX_HEIGHT_OCEAN_SHALLOW
                minimum_evap_fraction_land = self.cu_param_config.MINIMUM_EVAP_FRACTION_LAND_SHALLOW
                minimum_evap_fraction_ocean = self.cu_param_config.MINIMUM_EVAP_FRACTION_OCEAN_SHALLOW
                maximum_evap_fraction_land = self.cu_param_config.MAXIMUM_EVAP_FRACTION_LAND_SHALLOW
                maximum_evap_fraction_ocean = self.cu_param_config.MAXIMUM_EVAP_FRACTION_OCEAN_SHALLOW
                cloud_base_mass_flux_factor = self.cu_param_config.CLOUD_BASE_MASS_FLUX_FACTOR_SHALLOW
                use_excess = self.cu_param_config.USE_EXCESS_SHALLOW
                enable_plume = self.cu_param_config.ENABLE_SHALLOW
            if plume == "mid":
                downdraft_max_height_land = self.cu_param_config.DOWNDRAFT_MAX_HEIGHT_LAND_MID
                downdraft_max_height_ocean = self.cu_param_config.DOWNDRAFT_MAX_HEIGHT_OCEAN_MID
                updraft_max_height_land = self.cu_param_config.UPDRAFT_MAX_HEIGHT_LAND_MID
                updraft_max_height_ocean = self.cu_param_config.UPDRAFT_MAX_HEIGHT_OCEAN_MID
                minimum_evap_fraction_land = self.cu_param_config.MINIMUM_EVAP_FRACTION_LAND_MID
                minimum_evap_fraction_ocean = self.cu_param_config.MINIMUM_EVAP_FRACTION_OCEAN_MID
                maximum_evap_fraction_land = self.cu_param_config.MAXIMUM_EVAP_FRACTION_LAND_MID
                maximum_evap_fraction_ocean = self.cu_param_config.MAXIMUM_EVAP_FRACTION_OCEAN_MID
                cloud_base_mass_flux_factor = self.cu_param_config.CLOUD_BASE_MASS_FLUX_FACTOR_MID
                use_excess_t_q = self.cu_param_config.USE_EXCESS_MID
                enable_plume = self.cu_param_config.ENABLE_MID
            if plume == "deep":
                downdraft_max_height_land = self.cu_param_config.DOWNDRAFT_MAX_HEIGHT_LAND_DEEP
                downdraft_max_height_ocean = self.cu_param_config.DOWNDRAFT_MAX_HEIGHT_OCEAN_DEEP
                updraft_max_height_land = self.cu_param_config.UPDRAFT_MAX_HEIGHT_LAND_DEEP
                updraft_max_height_ocean = self.cu_param_config.UPDRAFT_MAX_HEIGHT_OCEAN_DEEP
                minimum_evap_fraction_land = self.cu_param_config.MINIMUM_EVAP_FRACTION_LAND_DEEP
                minimum_evap_fraction_ocean = self.cu_param_config.MINIMUM_EVAP_FRACTION_OCEAN_DEEP
                maximum_evap_fraction_land = self.cu_param_config.MAXIMUM_EVAP_FRACTION_LAND_DEEP
                maximum_evap_fraction_ocean = self.cu_param_config.MAXIMUM_EVAP_FRACTION_OCEAN_DEEP
                cloud_base_mass_flux_factor = self.cu_param_config.CLOUD_BASE_MASS_FLUX_FACTOR_DEEP
                use_excess = self.cu_param_config.USE_EXCESS_DEEP
                enable_plume = self.cu_param_config.ENABLE_DEEP

            self._prefil_excess(
                locals.cumulus_parameterization_input.t_excess,
                locals.cumulus_parameterization_internal.t_excess,
                locals.cumulus_parameterization_input.vapor_excess,
                locals.cumulus_parameterization_internal.vapor_excess,
                locals.cumulus_parameterization_input.ocean_fraction,
                use_excess,
            )

            self._set_new_t_vapor(
                t=locals.cumulus_parameterization_input.t,
                vapor=locals.cumulus_parameterization_input.vapor_timestep_start,
                grid_scale_forcing_t=locals.cumulus_parameterization_input.grid_scale_forcing_t,
                grid_scale_forcing_vapor=locals.cumulus_parameterization_input.grid_scale_forcing_vapor,
                subgrid_scale_forcing_t=locals.cumulus_parameterization_input.subgrid_scale_forcing_t,
                subgrid_scale_forcing_vapor=locals.cumulus_parameterization_input.subgrid_scale_forcing_vapor,
                t_cu_param_internal=locals.cumulus_parameterization_internal.t,
                vapor_cu_param_internal=locals.cumulus_parameterization_internal.vapor,
                t_pbl_cu_param_internal=locals.cumulus_parameterization_internal.t_pbl,
                vapor_pbl_cu_param_internal=locals.cumulus_parameterization_internal.vapor_pbl,
                moist_static_energy=locals.cumulus_parameterization_internal.moist_static_energy,
            )
