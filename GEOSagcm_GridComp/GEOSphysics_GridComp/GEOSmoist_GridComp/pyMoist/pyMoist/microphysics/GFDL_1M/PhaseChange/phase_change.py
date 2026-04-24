"""This module is the wrapper for the GFDL_1M microphysics scheme (in progress).
I/O and error handling is performed here.
Calculations can be found in deeper functions."""

from ndsl import Local, NDSLRuntime, Quantity, QuantityFactory, StencilFactory
from ndsl.constants import I_DIM, J_DIM, K_DIM
from ndsl.dsl.typing import Float

from pyMoist.constants import FLOAT_TINY
from pyMoist.microphysics.GFDL_1M.config import GFDL1MConfig
from pyMoist.microphysics.GFDL_1M.PhaseChange.evaporate import evaporate
from pyMoist.microphysics.GFDL_1M.PhaseChange.hydrostatic_pdf import hydrostatic_pdf
from pyMoist.microphysics.GFDL_1M.PhaseChange.melt_freeze import melt_freeze
from pyMoist.microphysics.GFDL_1M.PhaseChange.rh_calculations import fill_rh_crit_export, rh_calculations
from pyMoist.microphysics.GFDL_1M.PhaseChange.sublimate import sublimate
from pyMoist.saturation_tables import SaturationVaporPressureTable
from pyMoist.shared.incloud_processes import fix_up_clouds


class PhaseChange(NDSLRuntime):
    """This class is the wrapper for the GFDL_1M microphysics scheme. I/O and error handling
    are performed at this level, all calculations are performed within deeper functions.
    """

    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GFDL1MConfig,
        saturation_tables: SaturationVaporPressureTable,
    ):
        # init NDSLRuntime
        super().__init__(stencil_factory)

        if not config.USE_BERGERON:
            raise NotImplementedError(
                "Untested option for use_bergeron. Code may be missing or incomplete. "
                "Disable this error manually to continue."
            )

        if config.PDFSHAPE >= 5:
            raise NotImplementedError(f"PDF_SHAPE={config.PDFSHAPE} hasn't been ported from the Fortran")

        if config.PDFSHAPE > 1 and config.PDFSHAPE < 5:
            raise NotImplementedError(
                f"PDF_SHAPE={config.PDFSHAPE} is ported but untested. "
                "Disable this error manually to continue."
            )

        # make config and tables visible at runtime
        self.config = config
        self.saturation_tables = saturation_tables

        # innitalize locals
        self._alpha: Local = self.make_local(quantity_factory, [I_DIM, J_DIM, K_DIM], Float)

        # construct stencils
        self._rh_calculations = stencil_factory.from_dims_halo(
            func=rh_calculations,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "DW_LAND": config.DW_LAND,
                "DW_OCEAN": config.DW_OCEAN,
                "TURNRHCRIT_PARAM": config.TURNRHCRIT_PARAM,
            },
        )

        self._fill_rh_crit_export = stencil_factory.from_dims_halo(
            func=fill_rh_crit_export,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )

        self._hydrostatic_pdf = stencil_factory.from_dims_halo(
            func=hydrostatic_pdf,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "DT_MOIST": config.DT_MOIST,
                "PDF_SHAPE": config.PDFSHAPE,
                "USE_BERGERON": config.USE_BERGERON,
                "FLOAT_TINY": FLOAT_TINY,
            },
        )

        self._melt_freeze = stencil_factory.from_dims_halo(
            func=melt_freeze,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "DT_MOIST": config.DT_MOIST,
            },
        )
        self._evaporate = stencil_factory.from_dims_halo(
            func=evaporate,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "DT_MOIST": config.DT_MOIST,
                "CCW_EVAP_EFF": config.CCW_EVAP_EFF,
            },
        )
        self._sublimate = stencil_factory.from_dims_halo(
            func=sublimate,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "DT_MOIST": config.DT_MOIST,
                "CCI_EVAP_EFF": config.CCW_EVAP_EFF,
            },
        )
        self._fix_up_clouds = stencil_factory.from_dims_halo(
            func=fix_up_clouds,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )

        # Dev NOTE: this is an orchestration workaround. Direct call to
        #           `self.tables.X` fails closure capture for
        #           argument reconstruction at call time
        self._ese = self.saturation_tables.ese
        self._esw = self.saturation_tables.esw
        self._esx = self.saturation_tables.esx
        self._estfrz = self.saturation_tables.frz
        self._estlqu = self.saturation_tables.lqu

    def __call__(
        self,
        t: Quantity,
        mixing_ratio_vapor: Quantity,
        mixing_ratio_large_scale_liquid: Quantity,
        mixing_ratio_convective_liquid: Quantity,
        mixing_ratio_large_scale_ice: Quantity,
        mixing_ratio_convective_ice: Quantity,
        cloud_fraction_large_scale: Quantity,
        cloud_fraction_convective: Quantity,
        concentration_ice: Quantity,
        concentration_liquid: Quantity,
        relative_humidity_after_pdf: Quantity,
        estimated_inversion_strength: Quantity,
        area: Quantity,
        critical_relative_humidity_for_pdf: Quantity,
        pdf_iters: Quantity,
        cloud_liquid_evaporation: Quantity,
        cloud_ice_sublimation: Quantity,
        convection_fraction: Quantity,
        surface_type: Quantity,
        local_lcl_level: Quantity,
        local_p_mb: Quantity,
        local_p_interface_mb: Quantity,
        local_saturation_specific_humidity: Quantity,
    ):
        """Allow for phase change of excess vapor/liquid/ice before calling the microphysics driver

        Args:
            t (Quantity): temperature (Kelvin)
            mixing_ratio_vapor (Quantity): water vapor mixing ratio (kg/kg)
            mixing_ratio_large_scale_liquid (Quantity): large scale (non-convective) cloud liquid (kg/kg)
            mixing_ratio_convective_liquid (Quantity): convective liquid (kg/kg)
            mixing_ratio_large_scale_ice (Quantity): large scale (non-convective) cloud ice (kg/kg)
            mixing_ratio_convective_ice (Quantity): convective ice (kg/kg)
            cloud_fraction_large_scale (Quantity): (unitless)
            cloud_fraction_convective (Quantity): (unitless)
            concentration_ice (Quantity): cloud ice particle concentration (m^-3)
            concentration_liquid (Quantity): cloud liquid particle concentration (m^-3)
            relative_humidity_after_pdf (Quantity): (unitless)
            estimated_inversion_strength (Quantity): K
            area (Quantity): grid cell area (m^2)
            critical_relative_humidity_for_pdf (Quantity): (unitless)
            pdf_iters (Quantity): number of iterations in the hydrostatic pdf before exit
            cloud_liquid_evaporation (Quantity): (kg kg-1 s-1)
            cloud_ice_sublimation (Quantity): (kg kg-1 s-1)
            convection_fraction (Quantity)
            surface_type (Quantity)
            local_lcl_level (Quantity)
            local_p_mb (Quantity): grid center pressure (mb)
            local_p_interface_mb (Quantity): grid edge pressure (mb)
            local_saturation_specific_humidity (Quantity)

        """
        self._rh_calculations(
            estimated_inversion_strength=estimated_inversion_strength,
            p_mb=local_p_mb,
            p_interface_mb=local_p_interface_mb,
            area=area,
            lcl_level=local_lcl_level,
            alpha=self._alpha,
        )

        if critical_relative_humidity_for_pdf is not None:
            self._fill_rh_crit_export(self._alpha, critical_relative_humidity_for_pdf)

        self._hydrostatic_pdf(
            alpha=self._alpha,
            convection_fraction=convection_fraction,
            surface_type=surface_type,
            p_mb=local_p_mb,
            mixing_ratio_vapor=mixing_ratio_vapor,
            mixing_ratio_large_scale_liquid=mixing_ratio_large_scale_liquid,
            mixing_ratio_convective_liquid=mixing_ratio_convective_liquid,
            mixing_ratio_large_scale_ice=mixing_ratio_large_scale_ice,
            mixing_ratio_convective_ice=mixing_ratio_convective_ice,
            t=t,
            large_scale_cloud_fraction=cloud_fraction_large_scale,
            convective_cloud_fraction=cloud_fraction_convective,
            ice_concentration=concentration_ice,
            relative_humidity=relative_humidity_after_pdf,
            pdf_iters=pdf_iters,
            ese=self._ese,
            esw=self._esw,
            esx=self._esx,
            estfrz=self._estfrz,
            estlqu=self._estlqu,
        )

        if self.config.LMELTFRZ:
            self._melt_freeze(
                convection_fraction=convection_fraction,
                surface_type=surface_type,
                t=t,
                mixing_ratio_liquid=mixing_ratio_convective_liquid,
                mixing_ratio_ice=mixing_ratio_convective_ice,
            )
            self._melt_freeze(
                convection_fraction=convection_fraction,
                surface_type=surface_type,
                t=t,
                mixing_ratio_liquid=mixing_ratio_large_scale_liquid,
                mixing_ratio_ice=mixing_ratio_large_scale_ice,
            )

        if self.config.CCW_EVAP_EFF > 0.0 and not self.config.DO_EVAP:
            self._evaporate(
                p_mb=local_p_mb,
                t=t,
                mixing_ratio_vapor=mixing_ratio_vapor,
                mixing_ratio_convective_liquid=mixing_ratio_convective_liquid,
                mixing_ratio_convective_ice=mixing_ratio_convective_ice,
                convective_cloud_fraction=cloud_fraction_convective,
                liquid_concentration=concentration_liquid,
                ice_concentration=concentration_ice,
                saturation_specific_humidity=local_saturation_specific_humidity,
                evaporation=cloud_liquid_evaporation,
            )

        if self.config.CCI_EVAP_EFF > 0.0 and not self.config.DO_SUBL:
            self._sublimate(
                p_mb=local_p_mb,
                t=t,
                mixing_ratio_vapor=mixing_ratio_vapor,
                mixing_ratio_convective_liquid=mixing_ratio_convective_liquid,
                mixing_ratio_convective_ice=mixing_ratio_convective_ice,
                convective_cloud_fraction=cloud_fraction_convective,
                liquid_concentration=concentration_liquid,
                ice_concentration=concentration_ice,
                saturation_specific_humidity=local_saturation_specific_humidity,
                sublimation=cloud_ice_sublimation,
            )

        self._fix_up_clouds(
            mixing_ratio_vapor=mixing_ratio_vapor,
            t=t,
            mixing_ratio_large_scale_liquid=mixing_ratio_large_scale_liquid,
            mixing_ratio_large_scale_ice=mixing_ratio_large_scale_ice,
            large_scale_cloud_fraction=cloud_fraction_large_scale,
            mixing_ratio_convective_liquid=mixing_ratio_convective_liquid,
            mixing_ratio_convective_ice=mixing_ratio_convective_ice,
            convective_cloud_fraction=cloud_fraction_convective,
        )
