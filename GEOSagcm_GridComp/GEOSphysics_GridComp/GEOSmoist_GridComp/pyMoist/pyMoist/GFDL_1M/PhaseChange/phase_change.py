"""This module is the wrapper for the GFDL_1M microphysics scheme (in progress).
I/O and errorhandling is performed here.
Calculations can be found in deeper functions."""

from ndsl import QuantityFactory, StencilFactory, orchestrate
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ
from pyMoist.saturation_tables.formulation import SaturationFormulation
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.shared_incloud_processes import fix_up_clouds
from pyMoist.GFDL_1M.PhaseChange.config import PhaseChangeConfiguration
from pyMoist.GFDL_1M.PhaseChange.rh_calculations import rh_calculations
from pyMoist.GFDL_1M.PhaseChange.hydrostatic_pdf import hydrostatic_pdf
from pyMoist.GFDL_1M.PhaseChange.melt_freeze import melt_freeze
from pyMoist.GFDL_1M.PhaseChange.evaporate import evaporate
from pyMoist.GFDL_1M.PhaseChange.sublimate import sublimate
from pyMoist.GFDL_1M.PhaseChange.outputs import Outputs
from pyMoist.GFDL_1M.PhaseChange.masks import Masks
from pyMoist.GFDL_1M.PhaseChange.temporaries import Temporaries
from pyMoist.constants import FLOAT_TINY


class PhaseChange:
    """This class is the wrapper for the GFDL_1M microphysics scheme. I/O and error handling
    are perfromed at this level, all calculations are performed within deeper functions.
    """

    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        phase_change_config: PhaseChangeConfiguration,
        formulation: SaturationFormulation = SaturationFormulation.Staars,
    ):
        if phase_change_config.USE_BERGERON is not True:
            raise NotImplementedError(
                "Untested option for use_bergeron. Code may be missing or incomplete. \
                    Disable this error manually to continue."
            )

        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory
        self.phase_change_config = phase_change_config

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
                "DW_LAND": phase_change_config.DW_LAND,
                "DW_OCEAN": phase_change_config.DW_OCEAN,
                "TURNRHCRIT_PARAM": phase_change_config.TURNRHCRIT_PARAM,
            },
        )

        self._hydrostatic_pdf = self.stencil_factory.from_dims_halo(
            func=hydrostatic_pdf,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "DT_MOIST": phase_change_config.DT_MOIST,
                "PDF_SHAPE": phase_change_config.PDF_SHAPE,
                "USE_BERGERON": phase_change_config.USE_BERGERON,
                "FLOAT_TINY": FLOAT_TINY,
            },
        )

        self._meltfrz = self.stencil_factory.from_dims_halo(
            func=melt_freeze,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "DT_MOIST": phase_change_config.DT_MOIST,
            },
        )
        self._evap = self.stencil_factory.from_dims_halo(
            func=evaporate,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "DT_MOIST": phase_change_config.DT_MOIST,
                "CCW_EVAP_EFF": phase_change_config.CCW_EVAP_EFF,
            },
        )
        self._subl = self.stencil_factory.from_dims_halo(
            func=sublimate,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "DT_MOIST": phase_change_config.DT_MOIST,
                "CCI_EVAP_EFF": phase_change_config.CCW_EVAP_EFF,
            },
        )
        self._fix_up_clouds = self.stencil_factory.from_dims_halo(
            func=fix_up_clouds,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        phase_change_config: PhaseChangeConfiguration,
        estimated_inversion_strength: FloatFieldIJ,
        p_mb: FloatField,
        klcl: FloatFieldIJ,
        p_interface_mb: FloatField,
        area: FloatFieldIJ,
        cnv_frc: FloatFieldIJ,
        srf_type: FloatFieldIJ,
        t: FloatField,
        qlcn: FloatField,
        qicn: FloatField,
        qlls: FloatField,
        qils: FloatField,
        q: FloatField,
        clls: FloatField,
        clcn: FloatField,
        nactl: FloatField,
        nacti: FloatField,
        qst: FloatField,
    ):
        self._rh_calculations(
            eis=estimated_inversion_strength,
            minrhcrit=self.temporaries.minrhcrit,
            p_mb=p_mb,
            p_interface_mb=p_interface_mb,
            area=area,
            alpha=self.temporaries.alpha,
            klcl=klcl,
        )

        self._hydrostatic_pdf(
            alpha=self.temporaries.alpha,
            cnv_frc=cnv_frc,
            srf_type=srf_type,
            p_mb=p_mb,
            q=q,
            qlls=qlls,
            qlcn=qlcn,
            qils=qils,
            qicn=qicn,
            t=t,
            clls=clls,
            clcn=clcn,
            nacti=nacti,
            rhx=self.outputs.rhx,
            ese=self.tables.ese,
            esw=self.tables.esw,
            esx=self.tables.esx,
            estfrz=self.tables.frz,
            estlqu=self.tables.lqu,
        )

        if self.phase_change_config.DO_MELT_FREEZE:
            self._meltfrz(cnv_frc, srf_type, t, qlcn, qicn)
            self._meltfrz(cnv_frc, srf_type, t, qlls, qils)

        if self.phase_change_config.CCW_EVAP_EFF > 0.0:
            self._evap(
                p_mb,
                t,
                q,
                qlcn,
                qicn,
                clcn,
                nactl,
                nacti,
                qst,
                self.outputs.evapc,
            )

        if self.phase_change_config.CCI_EVAP_EFF > 0.0:
            self._subl(
                p_mb,
                t,
                q,
                qlcn,
                qicn,
                clcn,
                nactl,
                nacti,
                qst,
                self.outputs.sublc,
            )

        self._fix_up_clouds(q, t, qlls, qils, clls, qlcn, qicn, clcn)
