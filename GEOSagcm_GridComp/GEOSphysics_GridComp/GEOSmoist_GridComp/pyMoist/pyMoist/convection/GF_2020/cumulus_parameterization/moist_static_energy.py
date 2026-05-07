from ndsl import NDSLRuntime, Quantity, StencilFactory
from ndsl.constants import I_DIM, J_DIM, K_DIM
from ndsl.dsl.gt4py import FORWARD, PARALLEL, K, computation, interval
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ, Int, IntFieldIJ

import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.buoyancy import get_buoyancy
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import FloatField_Plume, IntFieldIJ_Plume
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import GF2020PlumeDependentConstants
from pyMoist.convection.GF_2020.cumulus_parameterization.shared_functions import get_cloud_boundary_conditions


def parcel_moist_static_energy(
    error_code: IntFieldIJ_Plume,
    t_excess: FloatFieldIJ,
    vapor_excess: FloatFieldIJ,
    add_buoyancy: FloatFieldIJ,
    ocean_fraction: FloatFieldIJ,
    updraft_origin_level: IntFieldIJ_Plume,
    p: FloatField,
    environment_moist_static_energy: FloatField,
    environment_moist_static_energy_forced: FloatField,
    t_perturbation: FloatField,
    moist_static_energy_origin_level: FloatFieldIJ,
    moist_static_energy_origin_level_forced: FloatFieldIJ,
    AVERAGE_LAYER_DEPTH: Float,
    plume: Int,
):
    """Determine the moist static energy of the forced and unforced parcel.

    Args:
        error_code (IntFieldIJ_Plume)
        t_excess (FloatFieldIJ)
        vapor_excess (FloatFieldIJ)
        add_buoyancy (FloatFieldIJ)
        ocean_fraction (FloatFieldIJ)
        updraft_origin_level (IntFieldIJ_Plume)
        p (FloatField)
        environment_moist_static_energy (FloatField)
        environment_moist_static_energy_forced (FloatField)
        t_perturbation (FloatField)
        moist_static_energy_origin_level (FloatFieldIJ)
        moist_static_energy_origin_level_forced (FloatFieldIJ)
        AVERAGE_LAYER_DEPTH (Float)
        plume (Int)
    """
    from __externals__ import BOUNDARY_CONDITION_METHOD, k_end

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            modification = (cumulus_parameterization_constants.XLV * vapor_excess + cumulus_parameterization_constants.CP * t_excess) + add_buoyancy

            moist_static_energy_origin_level = get_cloud_boundary_conditions(
                field=environment_moist_static_energy,
                scalar_perturbation=modification,
                p=p,
                updraft_origin_level=updraft_origin_level[0, 0][plume],
                ocean_fraction=ocean_fraction,
                BOUNDARY_CONDITION_METHOD=BOUNDARY_CONDITION_METHOD,
                AVERAGE_LAYER_DEPTH=AVERAGE_LAYER_DEPTH,
                k_end=k_end,
                compute_perturbation=False,
                perturbation_field=t_perturbation,
            )
            moist_static_energy_origin_level_forced = get_cloud_boundary_conditions(
                field=environment_moist_static_energy_forced,
                scalar_perturbation=modification,
                p=p,
                updraft_origin_level=updraft_origin_level[0, 0][plume],
                ocean_fraction=ocean_fraction,
                BOUNDARY_CONDITION_METHOD=BOUNDARY_CONDITION_METHOD,
                AVERAGE_LAYER_DEPTH=AVERAGE_LAYER_DEPTH,
                k_end=k_end,
                compute_perturbation=False,
                perturbation_field=t_perturbation,
            )


def first_guess_moist_static_energy(
    error_code: IntFieldIJ_Plume,
    start_level: IntFieldIJ,
    cloud_top_level: IntFieldIJ_Plume,
    mass_detrainment_updraft_forced: FloatField_Plume,
    mass_entrainment_updraft_forced: FloatField_Plume,
    normalized_massflux_updraft: FloatField,
    normalized_massflux_updraft_forced: FloatField_Plume,
    environment_moist_static_energy_forced: FloatField,
    environment_saturation_moist_static_energy_cloud_levels_forced: FloatField,
    cloud_moist_static_energy_forced: FloatField,
    vapor_excess: FloatFieldIJ,
    t_excess: FloatFieldIJ,
    add_buoyancy: FloatFieldIJ,
    plume: Int,
):
    """Compute moist static energy within the cloud (between LCL/start_level
    and equilibrium level/cloud_top_level)

    Args:
        error_code (IntFieldIJ_Plume)
        start_level (IntFieldIJ)
        cloud_top_level (IntFieldIJ_Plume)
        mass_detrainment_updraft_forced (FloatField_Plume)
        mass_entrainment_updraft_forced (FloatField_Plume)
        normalized_massflux_updraft (FloatField)
        normalized_massflux_updraft_forced (FloatField_Plume)
        environment_moist_static_energy_forced (FloatField)
        environment_saturation_moist_static_energy_cloud_levels_forced (FloatField)
        cloud_moist_static_energy_forced (FloatField)
        vapor_excess (FloatFieldIJ)
        t_excess (FloatFieldIJ)
        add_buoyancy (FloatFieldIJ)
        plume (Int)
    """
    with computation(FORWARD), interval(1, None):
        if error_code[0, 0][plume] == 0:
            if K >= start_level + 1 and K <= cloud_top_level[0, 0][plume] + 1:  # mass cons option
                denom: FloatFieldIJ = (
                    normalized_massflux_updraft[0, 0, -1] - 0.5 * mass_detrainment_updraft_forced[0, 0, -1][plume] + mass_entrainment_updraft_forced[0, 0, -1][plume]
                )
                if denom > 0.0:
                    cloud_moist_static_energy_forced = (
                        cloud_moist_static_energy_forced[0, 0, -1] * normalized_massflux_updraft_forced[0, 0, -1][plume]
                        - 0.5 * mass_detrainment_updraft_forced[0, 0, -1][plume] * cloud_moist_static_energy_forced[0, 0, -1]
                        + mass_entrainment_updraft_forced[0, 0, -1][plume] * environment_moist_static_energy_forced[0, 0, -1]
                    ) / denom
                    if K == start_level + 1:
                        perturbation: FloatFieldIJ = (
                            cumulus_parameterization_constants.XLV * vapor_excess + cumulus_parameterization_constants.CP * t_excess
                        ) + add_buoyancy
                        cloud_moist_static_energy_forced = cloud_moist_static_energy_forced + perturbation * mass_entrainment_updraft_forced[0, 0, -1][plume] / denom
                else:
                    cloud_moist_static_energy_forced = cloud_moist_static_energy_forced[0, 0, -1]

    with computation(PARALLEL), interval(0, -1):
        if error_code[0, 0][plume] == 0:
            if K >= cloud_top_level[0, 0][plume] + 2:
                cloud_moist_static_energy_forced = environment_saturation_moist_static_energy_cloud_levels_forced


def moist_static_energy_inside_cloud(
    error_code: IntFieldIJ_Plume,
    cloud_top_level: IntFieldIJ_Plume,
    start_level: IntFieldIJ,
    moist_static_energy_origin_level: FloatFieldIJ,
    cloud_moist_static_energy: FloatField,
    environment_moist_static_energy_modified: FloatField,
    environment_saturation_moist_static_energy_cloud_levels_modified: FloatField,
    mass_detrainment_updraft_forced: FloatField_Plume,
    mass_entrainment_updraft_forced: FloatField_Plume,
    normalized_massflux_updraft_modified: FloatField,
    partition_liquid_ice: FloatField,
    vapor_excess: FloatFieldIJ,
    t_excess: FloatFieldIJ,
    add_buoyancy: FloatFieldIJ,
    cloud_liquid_after_rain_forced: FloatField_Plume,
    plume: Int,
):
    """Compute moist static energy within the cloud.

    Args:
        error_code (IntFieldIJ_Plume)
        cloud_top_level (IntFieldIJ_Plume)
        start_level (IntFieldIJ)
        moist_static_energy_origin_level (FloatFieldIJ)
        cloud_moist_static_energy (FloatField)
        environment_moist_static_energy_modified (FloatField)
        environment_saturation_moist_static_energy_cloud_levels_modified (FloatField)
        mass_detrainment_updraft_forced (FloatField_Plume)
        mass_entrainment_updraft_forced (FloatField_Plume)
        normalized_massflux_updraft_modified (FloatField)
        partition_liquid_ice (FloatField)
        vapor_excess (FloatFieldIJ)
        t_excess (FloatFieldIJ)
        add_buoyancy (FloatFieldIJ)
        cloud_liquid_after_rain_forced (FloatField_Plume)
        plume (Int)
    """
    with computation(PARALLEL), interval(...):
        cloud_moist_static_energy = 0.0
        if error_code[0, 0][plume] == 0:
            if K <= start_level:
                cloud_moist_static_energy = moist_static_energy_origin_level

    with computation(FORWARD), interval(...):
        if error_code[0, 0][plume] == 0:
            if K >= start_level + 1 and K <= cloud_top_level[0, 0][plume] + 1:
                denom: FloatFieldIJ = (
                    normalized_massflux_updraft_modified.at(K=K - 1)
                    - 0.5 * mass_detrainment_updraft_forced.at(K=K - 1, ddim=[plume])
                    + mass_entrainment_updraft_forced.at(K=K - 1, ddim=[plume])
                )
                if denom == 0.0:
                    cloud_moist_static_energy = cloud_moist_static_energy.at(K=K - 1)
                else:
                    cloud_moist_static_energy = (
                        cloud_moist_static_energy.at(K=K - 1) * normalized_massflux_updraft_modified.at(K=K - 1)
                        - 0.5 * mass_detrainment_updraft_forced.at(K=K - 1, ddim=[plume]) * cloud_moist_static_energy.at(K=K - 1)
                        + mass_entrainment_updraft_forced.at(K=K - 1, ddim=[plume]) * environment_moist_static_energy_modified.at(K=K - 1)
                    ) / denom
                    if K == start_level + 1:
                        x_add: FloatFieldIJ = (cumulus_parameterization_constants.XLV * vapor_excess + cumulus_parameterization_constants.CP * t_excess) + add_buoyancy
                        cloud_moist_static_energy = cloud_moist_static_energy + x_add * mass_entrainment_updraft_forced.at(K=K - 1, ddim=[plume]) / denom

                # include glaciation effects on cloud_moist_static_energy
                cloud_moist_static_energy = (
                    cloud_moist_static_energy + cumulus_parameterization_constants.XLF * (1.0 - partition_liquid_ice) * cloud_liquid_after_rain_forced[0, 0, 0][plume]
                )

    with computation(PARALLEL), interval(0, -1):
        if error_code[0, 0][plume] == 0:
            if K >= cloud_top_level[0, 0][plume] + 2:
                cloud_moist_static_energy = environment_saturation_moist_static_energy_cloud_levels_modified
                normalized_massflux_updraft_modified = 0.0


class StaticControl(NDSLRuntime):
    """Update cloud moist static energy and compute buoyancy remaining after convection."""

    def __init__(
        self,
        stencil_factory: StencilFactory,
        config: GF2020Config,
        cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
    ):
        # init NDSLRuntime
        super().__init__(stencil_factory)

        # make configuration visible at runtime
        self.config = config
        self.cumulus_parameterization_config = cumulus_parameterization_config

        # construct stencils and functions
        self._moist_static_energy_inside_cloud = stencil_factory.from_dims_halo(
            func=moist_static_energy_inside_cloud,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )

        self._get_buoyancy = stencil_factory.from_dims_halo(
            func=get_buoyancy,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )

    def __call__(
        self,
        error_code: Quantity,
        start_level: Quantity,
        lcl_level: Quantity,
        updraft_lfc_level: Quantity,
        cloud_top_level: Quantity,
        cloud_moist_static_energy_modified: Quantity,
        moist_static_energy_origin_level_modified: Quantity,
        environment_moist_static_energy_modified: Quantity,
        environment_moist_static_energy_cloud_levels_modified: Quantity,
        environment_saturation_moist_static_energy_cloud_levels_modified: Quantity,
        mass_detrainment_updraft_forced: Quantity,
        mass_entrainment_updraft_forced: Quantity,
        normalized_massflux_updraft_modified: Quantity,
        partition_liquid_ice: Quantity,
        vapor_excess: Quantity,
        t_excess: Quantity,
        add_buoyancy: Quantity,
        cloud_liquid_after_rain_forced: Quantity,
        d_buoyancy_modified: Quantity,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._moist_static_energy_inside_cloud(
            error_code=error_code,
            cloud_top_level=cloud_top_level,
            start_level=start_level,
            moist_static_energy_origin_level=moist_static_energy_origin_level_modified,
            cloud_moist_static_energy=cloud_moist_static_energy_modified,
            environment_moist_static_energy_modified=environment_moist_static_energy_modified,
            environment_saturation_moist_static_energy_cloud_levels_modified=environment_saturation_moist_static_energy_cloud_levels_modified,
            mass_detrainment_updraft_forced=mass_detrainment_updraft_forced,
            mass_entrainment_updraft_forced=mass_entrainment_updraft_forced,
            normalized_massflux_updraft_modified=normalized_massflux_updraft_modified,
            partition_liquid_ice=partition_liquid_ice,
            vapor_excess=vapor_excess,
            t_excess=t_excess,
            add_buoyancy=add_buoyancy,
            cloud_liquid_after_rain_forced=cloud_liquid_after_rain_forced,
            plume=plume_dependent_constants.PLUME_INDEX,
        )

        self._get_buoyancy(
            lcl_level=lcl_level,
            updraft_lfc_level=updraft_lfc_level,
            cloud_top_level=cloud_top_level,
            cloud_moist_static_energy=cloud_moist_static_energy_modified,
            environment_moist_static_energy=environment_moist_static_energy_cloud_levels_modified,
            environment_saturation_moist_static_energy=environment_saturation_moist_static_energy_cloud_levels_modified,
            d_buoyancy=d_buoyancy_modified,
            error_code=error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
        )
