"""This module is the wrapper for the GFDL_1M microphysics scheme (in progress).
I/O and error handling is performed here.
Calculations can be found in deeper functions."""

from ndsl import Local, NDSLRuntime, QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float
from pyMoist.constants import FLOAT_TINY
from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.GFDL_1M.PhaseChange.evaporate import evaporate
from pyMoist.GFDL_1M.PhaseChange.hydrostatic_pdf import hydrostatic_pdf
from pyMoist.GFDL_1M.PhaseChange.melt_freeze import melt_freeze
from pyMoist.GFDL_1M.PhaseChange.rh_calculations import fill_rh_crit_export, rh_calculations
from pyMoist.GFDL_1M.PhaseChange.sublimate import sublimate
from pyMoist.saturation_tables import SaturationVaporPressureTable
from pyMoist.shared_incloud_processes import fix_up_clouds


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
        self._alpha: Local = self.make_local(quantity_factory, [X_DIM, Y_DIM, Z_DIM], Float)

        # construct stencils
        self._rh_calculations = stencil_factory.from_dims_halo(
            func=rh_calculations,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "DW_LAND": config.DW_LAND,
                "DW_OCEAN": config.DW_OCEAN,
                "TURNRHCRIT_PARAM": config.TURNRHCRIT_PARAM,
            },
        )

        self._fill_rh_crit_export = stencil_factory.from_dims_halo(
            func=fill_rh_crit_export,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._hydrostatic_pdf = stencil_factory.from_dims_halo(
            func=hydrostatic_pdf,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "DT_MOIST": config.DT_MOIST,
                "PDF_SHAPE": config.PDFSHAPE,
                "USE_BERGERON": config.USE_BERGERON,
                "FLOAT_TINY": FLOAT_TINY,
            },
        )

        self._meltfrz = stencil_factory.from_dims_halo(
            func=melt_freeze,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "DT_MOIST": config.DT_MOIST,
            },
        )
        self._evap = stencil_factory.from_dims_halo(
            func=evaporate,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "DT_MOIST": config.DT_MOIST,
                "CCW_EVAP_EFF": config.CCW_EVAP_EFF,
            },
        )
        self._subl = stencil_factory.from_dims_halo(
            func=sublimate,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "DT_MOIST": config.DT_MOIST,
                "CCI_EVAP_EFF": config.CCW_EVAP_EFF,
            },
        )
        self._fix_up_clouds = stencil_factory.from_dims_halo(
            func=fix_up_clouds,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
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
        t,
        mixing_ratio_vapor,
        mixing_ratio_large_scale_liquid,
        mixing_ratio_convective_liquid,
        mixing_ratio_large_scale_ice,
        mixing_ratio_convective_ice,
        cloud_fraction_large_scale,
        cloud_fraction_convective,
        concentration_ice,
        concentration_liquid,
        relative_humidity_after_pdf,
        estimated_inversion_strength,
        area,
        critical_relative_humidity_for_pdf,
        pdf_iters,
        cloud_liquid_evaporation,
        cloud_ice_sublimation,
        convection_fraction,
        surface_type,
        local_lcl_level,
        local_p_mb,
        local_p_interface_mb,
        local_saturation_specific_humidity,
    ):
        """

        Args:

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
            vapor=mixing_ratio_vapor,
            large_scale_liquid=mixing_ratio_large_scale_liquid,
            convective_liquid=mixing_ratio_convective_liquid,
            large_scale_ice=mixing_ratio_large_scale_ice,
            convective_ice=mixing_ratio_convective_ice,
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
            self._meltfrz(
                convection_fraction=convection_fraction,
                surface_type=surface_type,
                t=t,
                liquid=mixing_ratio_convective_liquid,
                ice=mixing_ratio_convective_ice,
            )
            self._meltfrz(
                convection_fraction=convection_fraction,
                surface_type=surface_type,
                t=t,
                liquid=mixing_ratio_large_scale_liquid,
                ice=mixing_ratio_large_scale_ice,
            )

        if self.config.CCW_EVAP_EFF > 0.0 and not self.config.DO_EVAP:
            self._evap(
                p_mb=local_p_mb,
                t=t,
                vapor=mixing_ratio_vapor,
                convective_liquid=mixing_ratio_convective_liquid,
                convective_ice=mixing_ratio_convective_ice,
                convective_cloud_fraction=cloud_fraction_convective,
                liquid_concentration=concentration_liquid,
                ice_concentration=concentration_ice,
                saturation_specific_humidity=local_saturation_specific_humidity,
                evaporation=cloud_liquid_evaporation,
            )

        if self.config.CCI_EVAP_EFF > 0.0 and not self.config.DO_SUBL:
            self._subl(
                p_mb=local_p_mb,
                t=t,
                vapor=mixing_ratio_vapor,
                convective_liquid=mixing_ratio_convective_liquid,
                convective_ice=mixing_ratio_convective_ice,
                convective_cloud_fraction=cloud_fraction_convective,
                liquid_concentration=concentration_liquid,
                ice_concentration=concentration_ice,
                saturation_specific_humidity=local_saturation_specific_humidity,
                sublimation=cloud_ice_sublimation,
            )

        self._fix_up_clouds(
            vapor=mixing_ratio_vapor,
            t=t,
            large_scale_liquid=mixing_ratio_large_scale_liquid,
            large_scale_ice=mixing_ratio_large_scale_ice,
            large_scale_cloud_fraction=cloud_fraction_large_scale,
            convective_liquid=mixing_ratio_convective_liquid,
            convective_ice=mixing_ratio_convective_ice,
            convective_cloud_fraction=cloud_fraction_convective,
        )
