"""This module is the wrapper for the GFDL_1M microphysics scheme (in progress).
I/O and errorhandling is performed here.
Calculations can be found in deeper functions."""

from ndsl import QuantityFactory, StencilFactory, orchestrate
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField, FloatFieldIJ
from pyMoist.constants import FLOAT_TINY
from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.GFDL_1M.PhaseChange.evaporate import evaporate
from pyMoist.GFDL_1M.PhaseChange.hydrostatic_pdf import hydrostatic_pdf
from pyMoist.GFDL_1M.PhaseChange.masks import Masks
from pyMoist.GFDL_1M.PhaseChange.melt_freeze import melt_freeze
from pyMoist.GFDL_1M.PhaseChange.outputs import Outputs
from pyMoist.GFDL_1M.PhaseChange.rh_calculations import rh_calculations
from pyMoist.GFDL_1M.PhaseChange.sublimate import sublimate
from pyMoist.GFDL_1M.PhaseChange.temporaries import Temporaries
from pyMoist.saturation_tables.formulation import SaturationFormulation
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.shared_incloud_processes import fix_up_clouds


class PhaseChange:
    """This class is the wrapper for the GFDL_1M microphysics scheme. I/O and error handling
    are perfromed at this level, all calculations are performed within deeper functions.
    """

    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        GFDL_1M_config: GFDL1MConfig,
        formulation: SaturationFormulation = SaturationFormulation.Staars,
    ):
        if GFDL_1M_config.USE_BERGERON is not True:
            raise NotImplementedError(
                "Untested option for use_bergeron. Code may be missing or incomplete. \
                    Disable this error manually to continue."
            )

        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory
        self.GFDL_1M_config = GFDL_1M_config

        # -----------------------------------------------------------------------
        # initialize precipitation outputs
        # -----------------------------------------------------------------------

        self.outputs = Outputs.make(quantity_factory)

        # -----------------------------------------------------------------------
        # initialize temporaries
        # -----------------------------------------------------------------------

        self.temporaries = Temporaries.make(quantity_factory)

        # -----------------------------------------------------------------------
        # initialize masks
        # -----------------------------------------------------------------------

        self.masks = Masks.make(quantity_factory)

        # -----------------------------------------------------------------------
        # Initalize QSat tables
        # -----------------------------------------------------------------------
        self.tables = SaturationVaporPressureTable(
            self.stencil_factory.backend,
            formulation=formulation,
        )

        # -----------------------------------------------------------------------
        # Initalizse stencils
        # -----------------------------------------------------------------------
        orchestrate(obj=self, config=stencil_factory.config.dace_config)

        self._rh_calculations = self.stencil_factory.from_dims_halo(
            func=rh_calculations,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "DW_LAND": GFDL_1M_config.DW_LAND,
                "DW_OCEAN": GFDL_1M_config.DW_OCEAN,
                "TURNRHCRIT_PARAM": GFDL_1M_config.TURNRHCRIT_PARAM,
            },
        )

        self._hydrostatic_pdf = self.stencil_factory.from_dims_halo(
            func=hydrostatic_pdf,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "DT_MOIST": GFDL_1M_config.DT_MOIST,
                "PDF_SHAPE": GFDL_1M_config.PDF_SHAPE,
                "USE_BERGERON": GFDL_1M_config.USE_BERGERON,
                "FLOAT_TINY": FLOAT_TINY,
            },
        )

        self._meltfrz = self.stencil_factory.from_dims_halo(
            func=melt_freeze,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "DT_MOIST": GFDL_1M_config.DT_MOIST,
            },
        )
        self._evap = self.stencil_factory.from_dims_halo(
            func=evaporate,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "DT_MOIST": GFDL_1M_config.DT_MOIST,
                "CCW_EVAP_EFF": GFDL_1M_config.CCW_EVAP_EFF,
            },
        )
        self._subl = self.stencil_factory.from_dims_halo(
            func=sublimate,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "DT_MOIST": GFDL_1M_config.DT_MOIST,
                "CCI_EVAP_EFF": GFDL_1M_config.CCW_EVAP_EFF,
            },
        )
        self._fix_up_clouds = self.stencil_factory.from_dims_halo(
            func=fix_up_clouds,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        estimated_inversion_strength: FloatFieldIJ,
        p_mb: FloatField,
        k_lcl: FloatFieldIJ,
        p_interface_mb: FloatField,
        area: FloatFieldIJ,
        convection_fraction: FloatFieldIJ,
        surface_type: FloatFieldIJ,
        t: FloatField,
        convective_liquid: FloatField,
        convective_ice: FloatField,
        large_scale_liquid: FloatField,
        large_scale_ice: FloatField,
        vapor: FloatField,
        large_scale_cloud_fraction: FloatField,
        convective_cloud_fraction: FloatField,
        nactl: FloatField,
        nacti: FloatField,
        qsat: FloatField,
    ):
        self._rh_calculations(
            estimated_inversion_strength=estimated_inversion_strength,
            minrhcrit=self.temporaries.minrhcrit,
            p_mb=p_mb,
            p_interface_mb=p_interface_mb,
            area=area,
            alpha=self.temporaries.alpha,
            k_lcl=k_lcl,
            rh_crit_3d=self.outputs.rh_crit,
        )

        self._hydrostatic_pdf(
            alpha=self.temporaries.alpha,
            convection_fraction=convection_fraction,
            surface_type=surface_type,
            p_mb=p_mb,
            vapor=vapor,
            large_scale_liquid=large_scale_liquid,
            convective_liquid=convective_liquid,
            large_scale_ice=large_scale_ice,
            convective_ice=convective_ice,
            t=t,
            large_scale_cloud_fraction=large_scale_cloud_fraction,
            convective_cloud_fraction=convective_cloud_fraction,
            nacti=nacti,
            rhx=self.outputs.rhx,
            ese=self.tables.ese,
            esw=self.tables.esw,
            esx=self.tables.esx,
            estfrz=self.tables.frz,
            estlqu=self.tables.lqu,
        )

        if self.GFDL_1M_config.MELTFRZ:
            self._meltfrz(convection_fraction, surface_type, t, convective_liquid, convective_ice)
            self._meltfrz(convection_fraction, surface_type, t, large_scale_liquid, large_scale_ice)

        if self.GFDL_1M_config.CCW_EVAP_EFF > 0.0:
            self._evap(
                p_mb,
                t,
                vapor,
                convective_liquid,
                convective_ice,
                convective_cloud_fraction,
                nactl,
                nacti,
                qsat,
                self.outputs.evapc,
            )

        if self.GFDL_1M_config.CCI_EVAP_EFF > 0.0:
            self._subl(
                p_mb,
                t,
                vapor,
                convective_liquid,
                convective_ice,
                convective_cloud_fraction,
                nactl,
                nacti,
                qsat,
                self.outputs.sublc,
            )

        self._fix_up_clouds(
            vapor,
            t,
            large_scale_liquid,
            large_scale_ice,
            large_scale_cloud_fraction,
            convective_liquid,
            convective_ice,
            convective_cloud_fraction,
        )
