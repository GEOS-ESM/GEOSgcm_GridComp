from f90nml import Namelist
from ndsl import StencilFactory
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.convection.GF_2020.state import GF2020State
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.cumulus_parameterization import (
    GF2020_CumulusParameterization,
)
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.convection.GF_2020.cumulus_parameterization.constants import (
    NUMBER_OF_PLUMES,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.state import (
    GF2020CumulusParameterizationState,
)
from pyMoist.convection_tracers import ConvectionTracers
import numpy as np


class TranslateGF2020_CumulusParameterization(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, namelist, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        self.in_vars["data_vars"] = {}

        self.out_vars = self.in_vars["data_vars"].copy()
        self.out_vars.update(
            {
                "t_excess": {},
                "vapor_excess": {},
                "grid_scale_forcing_t": {},
                "grid_scale_forcing_vapor": {},
                "subgrid_scale_forcing_t": {},
                "subgrid_scale_forcing_vapor": {},
                "seed_convection": {},
                "saturation_water_vapor": {},
                "ocean_fraction": {},
                "convection_fraction": {},
                "surface_type": {},
                "lateral_entrainment_rate": {},
                "last_error_code": {},
                "dtdt": {},
                "dvapordt": {},
                "dcloudicedt": {},
                "dudt": {},
                "dvdt": {},
                "dnliquiddt": {},
                "dnicedt": {},
                "dbuoyancydt": {},
                "dconvectiveicedt": {},
                "dlargescaleicedt": {},
                "dconvectiveliquiddt": {},
                "dlargescaleliquiddt": {},
                "dconvectivecloudfractiondt": {},
                "dlargescalecloudfractiondt": {},
                "error_code": {},
                "downdraft_origin_level": {},
                "lcl_level": {},
                "updraft_origin_level": {},
                "updraft_lfc_level": {},
                "cloud_top_level": {},
                "kstabi": {},
                "kstabm": {},
                "precip": {},
                "cloud_base_mass_flux_modified": {},
                "epsilon_forced": {},
                "total_normalized_integrated_condensate_forced": {},
                "scale_dependence_factor": {},
                "p_cloud_levels_forced": {},
                "entrainment_rate": {},
                "mass_entrainment_updraft_forced": {},
                "mass_entrainment_downdraft_forced": {},
                "mass_detrainment_updraft_forced": {},
                "mass_detrainment_downdraft_forced": {},
                "normalized_massflux_updraft_forced": {},
                "normalized_massflux_downdraft_forced": {},
                "condensate_to_fall_forced": {},
                "evaporate_in_downdraft_forced": {},
                "cloud_liquid_after_rain_forced": {},
                "t_updraft": {},
                "convective_cloud_fraction_output": {},
                "cloud_workfunction_0": {},
                "cloud_workfunction_1": {},
                "cloud_workfunction_2": {},
                "cloud_workfunction_3": {},
                "cloud_workfunction_1_pbl": {},
                "cloud_workfunction_1_cin": {},
                "cape_removal_time_scale": {},
                "pbl_time_scale": {},
                "lightning_density": {},
                "evaporation_sublimation_tendency": {},
                "convective_precip_flux": {},
                "t_perturbation": {},
                "grid_length": {},
                "pbl_level": {},
                "ccn": {},
                "air_density": {},
                "omega": {},
                "topography_height_no_negative": {},
                "sensible_heat_flux": {},
                "latent_heat_flux": {},
                "longitude_degrees": {},
                "latitude_degrees": {},
                "t_old": {},
                "vapor_old": {},
                "t_modified_by_advection": {},
                "vapor_modified_by_advection": {},
                "geopotential_height_forced": {},
                "p_forced": {},
                "p_surface": {},
                "t_surface": {},
                "u": {},
                "v": {},
                "w": {},
                "mass": {},
                "convective_scale_velocity": {},
                "buoyancy_excess": {},
                "large_scale_ice": {},
                "convective_ice": {},
                "large_scale_liquid": {},
                "convective_liquid": {},
                "large_scale_cloud_fraction": {},
                "convective_cloud_fraction": {},
                # "chemistry_tracers": {},
                # "chemistry_tracers_output": {},
            }
        )

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF2020-constants")
        self.cu_param_constants = data_loader.load("GF2020_CumulusParameterization-constants")
        self.convection_tracers_input = data_loader.load("GF2020_ConvectionTracers")

        # workaround because translate test cannot read in 4d fields
        self.manual_inputs = data_loader.load("GF2020_CumulusParameterization-In")

    def compute_func(self, **inputs):
        # initialize constants
        config = GF2020Config(SINGLE_COLUMN_MODE=False, **self.constants)
        cumulus_parameterization_config = GF2020CumulusParameterizationConfig(**self.cu_param_constants)

        # initialize convection tracers
        convection_tracers = ConvectionTracers.ones(
            self.quantity_factory,
            data_dimensions={
                "convection_tracers": config.NUMBER_OF_TRACERS,
                "size_three_dimension": 3,
                "size_four_dimension": 4,
            },
        )

        convection_tracers.tracers.field[:] = np.moveaxis(self.convection_tracers_input["tracers"], 0, 3)
        convection_tracers.vect_hcts.field[:] = self.convection_tracers_input["vect_hcts"]
        convection_tracers.kc_scal.field[:] = self.convection_tracers_input["kc_scal"]
        convection_tracers.fscav.field[:] = self.convection_tracers_input["fscav"]
        convection_tracers.convfaci2g.field[:] = self.convection_tracers_input["convfaci2g"]
        convection_tracers.retfactor.field[:] = self.convection_tracers_input["retfactor"]
        convection_tracers.liq_and_gas.field[:] = self.convection_tracers_input["liq_and_gas"]
        convection_tracers.online_cldliq.field[:] = self.convection_tracers_input["online_cldliq"]
        convection_tracers.online_vud.field[:] = self.convection_tracers_input["online_vud"]
        convection_tracers.ftemp_threshold.field[:] = self.convection_tracers_input["ftemp_threshold"]
        convection_tracers.use_gcc_washout.field[:] = self.convection_tracers_input["use_gcc_washout"]
        convection_tracers.use_gocart.field[:] = self.convection_tracers_input["use_gocart"]
        convection_tracers.is_wetdep.field[:] = self.convection_tracers_input["is_wetdep"]

        # initialize pyMoist saturation tables
        saturation_tables = SaturationVaporPressureTable(self.stencil_factory.backend)

        # initialize state
        state = GF2020CumulusParameterizationState.zeros(
            self.quantity_factory,
            data_dimensions={
                "plumes": NUMBER_OF_PLUMES,
                "convection_tracers": config.NUMBER_OF_TRACERS,
            },
        )

        # fill state with input data
        if cumulus_parameterization_config.PLUME_ORDER == 0:
            # fortran data is saved [deep, shallow, mid] in plume dimension
            # must be reordered to [shallow, mid, deep] for python
            # input fields
            state.input.t_excess.field[:] = self.manual_inputs["t_excess"]
            state.input.vapor_excess.field[:] = self.manual_inputs["vapor_excess"]
            state.input.grid_scale_forcing_t.field[:] = self.manual_inputs["grid_scale_forcing_t"]
            state.input.grid_scale_forcing_vapor.field[:] = self.manual_inputs["grid_scale_forcing_vapor"]
            state.input.subgrid_scale_forcing_t.field[:] = self.manual_inputs["subgrid_scale_forcing_t"]
            state.input.subgrid_scale_forcing_vapor.field[:] = self.manual_inputs[
                "subgrid_scale_forcing_vapor"
            ]
            state.input.seed_convection.field[:] = self.manual_inputs["seed_convection"]
            state.input.saturation_water_vapor.field[:] = self.manual_inputs["saturation_water_vapor"]
            state.input.ocean_fraction.field[:] = self.manual_inputs["ocean_fraction"]
            state.input.convection_fraction.field[:] = self.manual_inputs["convection_fraction"]
            state.input.surface_type.field[:] = self.manual_inputs["surface_type"]
            state.input.lateral_entrainment_rate.field[:] = self.manual_inputs["lateral_entrainment_rate"]
            state.input.last_error_code.field[:] = self.manual_inputs["last_error_code"]
            # output fields
            state.output.dtdt.field[:] = self.manual_inputs["dtdt"][:, :, :, [1, 2, 0]]
            state.output.dvapordt.field[:] = self.manual_inputs["dvapordt"][:, :, :, [1, 2, 0]]
            state.output.dcloudicedt.field[:] = self.manual_inputs["dcloudicedt"][:, :, :, [1, 2, 0]]
            state.output.dudt.field[:] = self.manual_inputs["dudt"][:, :, :, [1, 2, 0]]
            state.output.dvdt.field[:] = self.manual_inputs["dvdt"][:, :, :, [1, 2, 0]]
            state.output.dnliquiddt.field[:] = self.manual_inputs["dnliquiddt"][:, :, :, [1, 2, 0]]
            state.output.dnicedt.field[:] = self.manual_inputs["dnicedt"][:, :, :, [1, 2, 0]]
            state.output.dbuoyancydt.field[:] = self.manual_inputs["dbuoyancydt"][:, :, :, [1, 2, 0]]
            state.output.dconvectiveicedt.field[:] = self.manual_inputs["dconvectiveicedt"][
                :, :, :, [1, 2, 0]
            ]
            state.output.dlargescaleicedt.field[:] = self.manual_inputs["dlargescaleicedt"][
                :, :, :, [1, 2, 0]
            ]
            state.output.dconvectiveliquiddt.field[:] = self.manual_inputs["dconvectiveliquiddt"][
                :, :, :, [1, 2, 0]
            ]
            state.output.dlargescaleliquiddt.field[:] = self.manual_inputs["dlargescaleliquiddt"][
                :, :, :, [1, 2, 0]
            ]
            state.output.dconvectivecloudfractiondt.field[:] = self.manual_inputs[
                "dconvectivecloudfractiondt"
            ][:, :, :, [1, 2, 0]]
            state.output.dlargescalecloudfractiondt.field[:] = self.manual_inputs[
                "dlargescalecloudfractiondt"
            ][:, :, :, [1, 2, 0]]
            state.output.error_code.field[:] = self.manual_inputs["error_code"][:, :, [1, 2, 0]]
            state.output.downdraft_origin_level.field[:] = (
                self.manual_inputs["downdraft_origin_level"][:, :, [1, 2, 0]] - 1
            )
            state.output.lcl_level.field[:] = self.manual_inputs["lcl_level"][:, :, [1, 2, 0]] - 1
            state.output.updraft_origin_level.field[:] = (
                self.manual_inputs["updraft_origin_level"][:, :, [1, 2, 0]] - 1
            )
            state.output.updraft_lfc_level.field[:] = (
                self.manual_inputs["updraft_lfc_level"][:, :, [1, 2, 0]] - 1
            )
            state.output.cloud_top_level.field[:] = self.manual_inputs["cloud_top_level"][:, :, [1, 2, 0]] - 1
            state.output.kstabi.field[:] = self.manual_inputs["kstabi"][:, :, [1, 2, 0]] - 1
            state.output.kstabm.field[:] = self.manual_inputs["kstabm"][:, :, [1, 2, 0]] - 1
            state.output.precip.field[:] = self.manual_inputs["precip"][:, :, [1, 2, 0]]
            state.output.cloud_base_mass_flux_modified.field[:] = self.manual_inputs[
                "cloud_base_mass_flux_modified"
            ][:, :, [1, 2, 0]]
            state.output.epsilon_forced.field[:] = self.manual_inputs["epsilon_forced"][:, :, [1, 2, 0]]
            state.output.total_normalized_integrated_condensate_forced.field[:] = self.manual_inputs[
                "total_normalized_integrated_condensate_forced"
            ][:, :, [1, 2, 0]]
            state.output.scale_dependence_factor.field[:] = self.manual_inputs["scale_dependence_factor"][
                :, :, [1, 2, 0]
            ]
            state.output.p_cloud_levels_forced.field[:] = self.manual_inputs["p_cloud_levels_forced"][
                :, :, :, [1, 2, 0]
            ]
            state.output.entrainment_rate.field[:] = self.manual_inputs["entrainment_rate"]
            state.output.mass_entrainment_updraft_forced.field[:] = self.manual_inputs[
                "mass_entrainment_updraft_forced"
            ][:, :, :, [1, 2, 0]]
            state.output.mass_entrainment_downdraft_forced.field[:] = self.manual_inputs[
                "mass_entrainment_downdraft_forced"
            ][:, :, :, [1, 2, 0]]
            state.output.mass_detrainment_updraft_forced.field[:] = self.manual_inputs[
                "mass_detrainment_updraft_forced"
            ][:, :, :, [1, 2, 0]]
            state.output.mass_detrainment_downdraft_forced.field[:] = self.manual_inputs[
                "mass_detrainment_downdraft_forced"
            ][:, :, :, [1, 2, 0]]
            state.output.normalized_massflux_updraft_forced.field[:] = self.manual_inputs[
                "normalized_massflux_updraft_forced"
            ][:, :, :, [1, 2, 0]]
            state.output.normalized_massflux_downdraft_forced.field[:] = self.manual_inputs[
                "normalized_massflux_downdraft_forced"
            ][:, :, :, [1, 2, 0]]
            state.output.condensate_to_fall_forced.field[:] = self.manual_inputs["condensate_to_fall_forced"][
                :, :, :, [1, 2, 0]
            ]
            state.output.evaporate_in_downdraft_forced.field[:] = self.manual_inputs[
                "evaporate_in_downdraft_forced"
            ][:, :, :, [1, 2, 0]]
            state.output.cloud_liquid_after_rain_forced.field[:] = self.manual_inputs[
                "cloud_liquid_after_rain_forced"
            ][:, :, :, [1, 2, 0]]
            state.output.t_updraft.field[:] = self.manual_inputs["t_updraft"][:, :, :, [1, 2, 0]]
            state.output.convective_cloud_fraction.field[:] = self.manual_inputs[
                "convective_cloud_fraction_output"
            ][:, :, :, [1, 2, 0]]
            state.output.cloud_workfunction_0.field[:] = self.manual_inputs["cloud_workfunction_0"]
            state.output.cloud_workfunction_1.field[:] = self.manual_inputs["cloud_workfunction_1"]
            state.output.cloud_workfunction_2.field[:] = self.manual_inputs["cloud_workfunction_2"]
            state.output.cloud_workfunction_3.field[:] = self.manual_inputs["cloud_workfunction_3"]
            state.output.cloud_workfunction_1_pbl.field[:] = self.manual_inputs["cloud_workfunction_1_pbl"]
            state.output.cloud_workfunction_1_cin.field[:] = self.manual_inputs["cloud_workfunction_1_cin"]
            state.output.cape_removal_time_scale.field[:] = self.manual_inputs["cape_removal_time_scale"]
            state.output.pbl_time_scale.field[:] = self.manual_inputs["pbl_time_scale"]
            state.output.lightning_density.field[:] = self.manual_inputs["lightning_density"]
            state.output.evaporation_sublimation_tendency.field[:] = self.manual_inputs[
                "evaporation_sublimation_tendency"
            ]
            state.output.convective_precip_flux.field[:] = self.manual_inputs["convective_precip_flux"]
            state.output.t_perturbation.field[:] = self.manual_inputs["t_perturbation"]
            # input/output fields
            state.input_output.grid_length.field[:] = self.manual_inputs["grid_length"]
            state.input_output.pbl_level.field[:] = self.manual_inputs["pbl_level"] - 1
            state.input_output.ccn.field[:] = self.manual_inputs["ccn"]
            state.input_output.air_density.field[:] = self.manual_inputs["air_density"]
            state.input_output.omega.field[:] = self.manual_inputs["omega"]
            state.input_output.topography_height_no_negative.field[:] = self.manual_inputs[
                "topography_height_no_negative"
            ]
            state.input_output.sensible_heat_flux.field[:] = self.manual_inputs["sensible_heat_flux"]
            state.input_output.latent_heat_flux.field[:] = self.manual_inputs["latent_heat_flux"]
            state.input_output.longitude_degrees.field[:] = self.manual_inputs["longitude_degrees"]
            state.input_output.latitude_degrees.field[:] = self.manual_inputs["latitude_degrees"]
            state.input_output.t_old.field[:] = self.manual_inputs["t_old"]
            state.input_output.vapor_old.field[:] = self.manual_inputs["vapor_old"]
            state.input_output.t_modified_by_advection.field[:] = self.manual_inputs[
                "t_modified_by_advection"
            ]
            state.input_output.vapor_modified_by_advection.field[:] = self.manual_inputs[
                "vapor_modified_by_advection"
            ]
            state.input_output.geopotential_height_forced.field[:] = self.manual_inputs[
                "geopotential_height_forced"
            ]
            state.input_output.p_forced.field[:] = self.manual_inputs["p_forced"]
            state.input_output.p_surface.field[:] = self.manual_inputs["p_surface"]
            state.input_output.t_surface.field[:] = self.manual_inputs["t_surface"]
            state.input_output.u.field[:] = self.manual_inputs["u"]
            state.input_output.v.field[:] = self.manual_inputs["v"]
            state.input_output.w.field[:] = self.manual_inputs["w"]
            state.input_output.mass.field[:] = self.manual_inputs["mass"]
            state.input_output.convective_scale_velocity.field[:] = self.manual_inputs[
                "convective_scale_velocity"
            ]
            state.input_output.buoyancy_excess.field[:] = self.manual_inputs["buoyancy_excess"]
            state.input_output.large_scale_ice.field[:] = self.manual_inputs["large_scale_ice"]
            state.input_output.convective_ice.field[:] = self.manual_inputs["convective_ice"]
            state.input_output.large_scale_liquid.field[:] = self.manual_inputs["large_scale_liquid"]
            state.input_output.convective_liquid.field[:] = self.manual_inputs["convective_liquid"]
            state.input_output.large_scale_cloud_fraction.field[:] = self.manual_inputs[
                "large_scale_cloud_fraction"
            ]
            state.input_output.convective_cloud_fraction.field[:] = self.manual_inputs[
                "convective_cloud_fraction"
            ]
            state.input_output.chemistry_tracers.field[:] = self.manual_inputs["chemistry_tracers"]
            chemistry_tracers_input_5d = np.full(
                state.input_output.chemistry_tracers_output.field[:].shape, np.nan
            )
            for plume in range(NUMBER_OF_PLUMES):
                chemistry_tracers_input_5d[:, :, :, plume, :] = self.manual_inputs[
                    "chemistry_tracers_output"
                ][
                    :,
                    :,
                    :,
                    plume * config.NUMBER_OF_TRACERS : plume * config.NUMBER_OF_TRACERS
                    + config.NUMBER_OF_TRACERS,
                ]
            state.input_output.chemistry_tracers_output.field[:] = chemistry_tracers_input_5d[
                :, :, :, [1, 2, 0], :
            ]

            # initialize the test subject
            code = GF2020_CumulusParameterization(
                stencil_factory=self.stencil_factory,
                quantity_factory=self.quantity_factory,
                config=config,
                cumulus_parameterization_config=cumulus_parameterization_config,
            )

            # call the test subject
            code(state, saturation_tables, convection_tracers)

            # collapse plume dim for chemistry_tracers_output
            # NOTE ideally this has no numpy dependency
            chemistry_tracers_output_5d_reordered = state.input_output.chemistry_tracers_output.field[
                :, :, :, [2, 0, 1], :
            ]
            chemistry_tracers_output_4d = np.full(
                self.manual_inputs["chemistry_tracers_output"].shape, np.nan
            )
            for plume in range(NUMBER_OF_PLUMES):
                chemistry_tracers_output_4d[
                    :,
                    :,
                    :,
                    plume * config.NUMBER_OF_TRACERS : plume * config.NUMBER_OF_TRACERS
                    + config.NUMBER_OF_TRACERS,
                ] = chemistry_tracers_output_5d_reordered[:, :, :, plume, :]

            outputs = {
                # input fields
                "t_excess": state.input.t_excess.field[:],
                "vapor_excess": state.input.vapor_excess.field[:],
                "grid_scale_forcing_t": state.input.grid_scale_forcing_t.field[:],
                "grid_scale_forcing_vapor": state.input.grid_scale_forcing_vapor.field[:],
                "subgrid_scale_forcing_t": state.input.subgrid_scale_forcing_t.field[:],
                "subgrid_scale_forcing_vapor": state.input.subgrid_scale_forcing_vapor.field[:],
                "seed_convection": state.input.seed_convection.field[:],
                "saturation_water_vapor": state.input.saturation_water_vapor.field[:],
                "ocean_fraction": state.input.ocean_fraction.field[:],
                "convection_fraction": state.input.convection_fraction.field[:],
                "surface_type": state.input.surface_type.field[:],
                "lateral_entrainment_rate": state.input.lateral_entrainment_rate.field[:],
                "last_error_code": state.input.last_error_code.field[:],
                # output fields
                "dtdt": state.output.dtdt.field[:, :, :, [2, 0, 1]],
                "dvapordt": state.output.dvapordt.field[:, :, :, [2, 0, 1]],
                "dcloudicedt": state.output.dcloudicedt.field[:, :, :, [2, 0, 1]],
                "dudt": state.output.dudt.field[:, :, :, [2, 0, 1]],
                "dvdt": state.output.dvdt.field[:, :, :, [2, 0, 1]],
                "dnliquiddt": state.output.dnliquiddt.field[:, :, :, [2, 0, 1]],
                "dnicedt": state.output.dnicedt.field[:, :, :, [2, 0, 1]],
                "dbuoyancydt": state.output.dbuoyancydt.field[:, :, :, [2, 0, 1]],
                "dconvectiveicedt": state.output.dconvectiveicedt.field[:, :, :, [2, 0, 1]],
                "dlargescaleicedt": state.output.dlargescaleicedt.field[:, :, :, [2, 0, 1]],
                "dconvectiveliquiddt": state.output.dconvectiveliquiddt.field[:, :, :, [2, 0, 1]],
                "dlargescaleliquiddt": state.output.dlargescaleliquiddt.field[:, :, :, [2, 0, 1]],
                "dconvectivecloudfractiondt": state.output.dconvectivecloudfractiondt.field[
                    :, :, :, [2, 0, 1]
                ],
                "dlargescalecloudfractiondt": state.output.dlargescalecloudfractiondt.field[
                    :, :, :, [2, 0, 1]
                ],
                "error_code": state.output.error_code.field[:, :, [2, 0, 1]],
                "downdraft_origin_level": state.output.downdraft_origin_level.field[:, :, [2, 0, 1]] + 1,
                "lcl_level": state.output.lcl_level.field[:, :, [2, 0, 1]] + 1,
                "updraft_origin_level": state.output.updraft_origin_level.field[:, :, [2, 0, 1]] + 1,
                "updraft_lfc_level": state.output.updraft_lfc_level.field[:, :, [2, 0, 1]] + 1,
                "cloud_top_level": state.output.cloud_top_level.field[:, :, [2, 0, 1]] + 1,
                "kstabi": state.output.kstabi.field[:, :, [2, 0, 1]] + 1,
                "kstabm": state.output.kstabm.field[:, :, [2, 0, 1]] + 1,
                "precip": state.output.precip.field[:, :, [2, 0, 1]],
                "cloud_base_mass_flux_modified": state.output.cloud_base_mass_flux_modified.field[
                    :, :, [2, 0, 1]
                ],
                "epsilon_forced": state.output.epsilon_forced.field[:, :, [2, 0, 1]],
                "total_normalized_integrated_condensate_forced": state.output.total_normalized_integrated_condensate_forced.field[
                    :, :, [2, 0, 1]
                ],
                "scale_dependence_factor": state.output.scale_dependence_factor.field[:, :, [2, 0, 1]],
                "p_cloud_levels_forced": state.output.p_cloud_levels_forced.field[:, :, :, [2, 0, 1]],
                "entrainment_rate": state.output.entrainment_rate.field[:, :, :, [2, 0, 1]],
                "mass_entrainment_updraft_forced": state.output.mass_entrainment_updraft_forced.field[
                    :, :, :, [2, 0, 1]
                ],
                "mass_entrainment_downdraft_forced": state.output.mass_entrainment_downdraft_forced.field[
                    :, :, :, [2, 0, 1]
                ],
                "mass_detrainment_updraft_forced": state.output.mass_detrainment_updraft_forced.field[
                    :, :, :, [2, 0, 1]
                ],
                "mass_detrainment_downdraft_forced": state.output.mass_detrainment_downdraft_forced.field[
                    :, :, :, [2, 0, 1]
                ],
                "normalized_massflux_updraft_forced": state.output.normalized_massflux_updraft_forced.field[
                    :, :, :, [2, 0, 1]
                ],
                "normalized_massflux_downdraft_forced": state.output.normalized_massflux_downdraft_forced.field[
                    :, :, :, [2, 0, 1]
                ],
                "condensate_to_fall_forced": state.output.condensate_to_fall_forced.field[:, :, :, [2, 0, 1]],
                "evaporate_in_downdraft_forced": state.output.evaporate_in_downdraft_forced.field[
                    :, :, :, [2, 0, 1]
                ],
                "cloud_liquid_after_rain_forced": state.output.cloud_liquid_after_rain_forced.field[
                    :, :, :, [2, 0, 1]
                ],
                "t_updraft": state.output.t_updraft.field[:, :, :, [2, 0, 1]],
                "convective_cloud_fraction_output": state.output.convective_cloud_fraction.field[
                    :, :, :, [2, 0, 1]
                ],
                "cloud_workfunction_0": state.output.cloud_workfunction_0.field[:],
                "cloud_workfunction_1": state.output.cloud_workfunction_1.field[:],
                "cloud_workfunction_2": state.output.cloud_workfunction_2.field[:],
                "cloud_workfunction_3": state.output.cloud_workfunction_3.field[:],
                "cloud_workfunction_1_pbl": state.output.cloud_workfunction_1_pbl.field[:],
                "cloud_workfunction_1_cin": state.output.cloud_workfunction_1_cin.field[:],
                "cape_removal_time_scale": state.output.cape_removal_time_scale.field[:],
                "pbl_time_scale": state.output.pbl_time_scale.field[:],
                "lightning_density": state.output.lightning_density.field[:],
                "evaporation_sublimation_tendency": state.output.evaporation_sublimation_tendency.field[:],
                "convective_precip_flux": state.output.convective_precip_flux.field[:],
                "t_perturbation": state.output.t_perturbation.field[:],
                # input/output fields
                "grid_length": state.input_output.grid_length.field[:],
                "pbl_level": state.input_output.pbl_level.field[:] + 1,
                "ccn": state.input_output.ccn.field[:],
                "air_density": state.input_output.air_density.field[:],
                "omega": state.input_output.omega.field[:],
                "topography_height_no_negative": state.input_output.topography_height_no_negative.field[:],
                "sensible_heat_flux": state.input_output.sensible_heat_flux.field[:],
                "latent_heat_flux": state.input_output.latent_heat_flux.field[:],
                "longitude_degrees": state.input_output.longitude_degrees.field[:],
                "latitude_degrees": state.input_output.latitude_degrees.field[:],
                "t_old": state.input_output.t_old.field[:],
                "vapor_old": state.input_output.vapor_old.field[:],
                "t_modified_by_advection": state.input_output.t_modified_by_advection.field[:],
                "vapor_modified_by_advection": state.input_output.vapor_modified_by_advection.field[:],
                "geopotential_height_forced": state.input_output.geopotential_height_forced.field[:],
                "p_forced": state.input_output.p_forced.field[:],
                "p_surface": state.input_output.p_surface.field[:],
                "t_surface": state.input_output.t_surface.field[:],
                "u": state.input_output.u.field[:],
                "v": state.input_output.v.field[:],
                "w": state.input_output.w.field[:],
                "mass": state.input_output.mass.field[:],
                "convective_scale_velocity": state.input_output.convective_scale_velocity.field[:],
                "buoyancy_excess": state.input_output.buoyancy_excess.field[:],
                "large_scale_ice": state.input_output.large_scale_ice.field[:],
                "convective_ice": state.input_output.convective_ice.field[:],
                "large_scale_liquid": state.input_output.large_scale_liquid.field[:],
                "convective_liquid": state.input_output.convective_liquid.field[:],
                "large_scale_cloud_fraction": state.input_output.large_scale_cloud_fraction.field[:],
                "convective_cloud_fraction": state.input_output.convective_cloud_fraction.field[:],
                "chemistry_tracers": state.input_output.chemistry_tracers.field[:],
                "chemistry_tracers_output": chemistry_tracers_output_4d,
            }
        else:
            raise NotImplementedError(
                "Plume order unsupported. Please implement a way to transfer the fortran shape into python shape."
            )

        return outputs
