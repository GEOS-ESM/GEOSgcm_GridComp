from ndsl import QuantityFactory, StencilFactory, orchestrate
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ
from pyMoist.saturation.formulation import SaturationFormulation
from pyMoist.saturation.qsat import QSat

from .GFDL_1M_util import (
    evap,
    fix_up_clouds,
    get_last,
    hybrid_index_2dout,
    hystpdf,
    initial_calc,
    meltfrz,
    subl,
)


class evap_subl_pdf:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        formulation: SaturationFormulation = SaturationFormulation.Staars,
        USE_BERGERON: bool = True,
    ):
        if USE_BERGERON is not True:
            raise NotImplementedError(
                "Untested option for USE_BERGERON. Code may be missing or incomplete. \
                    Disable this error manually to continue."
            )

        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory

        # Initalize QSat tables for later calculations
        self.qsat = QSat(
            self.stencil_factory,
            self.quantity_factory,
            formulation=formulation,
        )

        orchestrate(obj=self, config=stencil_factory.config.dace_config)
        self._get_last = self.stencil_factory.from_dims_halo(
            func=get_last,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )
        self._hybrid_index_2dout = self.stencil_factory.from_dims_halo(
            func=hybrid_index_2dout,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )
        self._initial_calc = self.stencil_factory.from_dims_halo(
            func=initial_calc,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )
        self._hystpdf = self.stencil_factory.from_dims_halo(
            func=hystpdf,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "USE_BERGERON": USE_BERGERON,
            },
        )
        self._meltfrz = self.stencil_factory.from_dims_halo(
            func=meltfrz,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )
        self._evap = self.stencil_factory.from_dims_halo(
            func=evap,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )
        self._subl = self.stencil_factory.from_dims_halo(
            func=subl,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )
        self._fix_up_clouds = self.stencil_factory.from_dims_halo(
            func=fix_up_clouds,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )
        self._tmp = self.quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        self._minrhcrit = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._PLEmb_top = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._PLmb_at_klcl = self.quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        self._halo = self.stencil_factory.grid_indexing.n_halo
        self._k_mask = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._alpha = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.rhx = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.evapc = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.sublc = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        for i in range(0, self._k_mask.view[:].shape[0]):
            for j in range(0, self._k_mask.view[:].shape[1]):
                for k in range(0, self._k_mask.view[:].shape[2]):
                    self._k_mask.view[i, j, k] = k + 1

    def __call__(
        self,
        EIS: FloatFieldIJ,
        dw_land: Float,
        dw_ocean: Float,
        PDFSHAPE: Float,
        TURNRHCRIT_PARAM: Float,
        PLmb: FloatField,
        KLCL: FloatFieldIJ,
        PLEmb: FloatField,
        AREA: FloatFieldIJ,
        DT_MOIST: Float,
        CNV_FRC: FloatFieldIJ,
        SRF_TYPE: FloatFieldIJ,
        T: FloatField,
        QLCN: FloatField,
        QICN: FloatField,
        QLLS: FloatField,
        QILS: FloatField,
        CCW_EVAP_EFF: Float,
        CCI_EVAP_EFF: Float,
        Q: FloatField,
        CLLS: FloatField,
        CLCN: FloatField,
        NACTL: FloatField,
        NACTI: FloatField,
        QST: FloatField,
        LMELTFRZ: bool = True,
    ):
        if LMELTFRZ is not True:
            raise NotImplementedError(
                "Untested option for LMELTFRZ. Code may be missing or incomplete. \
                    Disable this error manually to continue."
            )
        if PDFSHAPE != 1:
            raise NotImplementedError(
                "Untested option for PDFSHAPE. Code may be missing or incomplete. \
                    Disable this error manually to continue."
            )
        if CCW_EVAP_EFF <= 0:
            raise NotImplementedError(
                "Untested option for CCW_EVAP_EFF. Code may be missing or incomplete. \
                    Disable this error manually to continue."
            )
        if CCI_EVAP_EFF <= 0:
            raise NotImplementedError(
                "Untested option for CCI_EVAP_EFF. Code may be missing or incomplete. \
                    Disable this error manually to continue."
            )

        self._get_last(PLEmb, self._tmp, self._PLEmb_top)

        # Temporary implementation of hybrid_index_2dout.py, perhaps not working as
        # indended (backend issue), will need to be addressed at later date
        self._hybrid_index_2dout(PLmb, self._k_mask, KLCL, self._PLmb_at_klcl)

        self._initial_calc(
            EIS,
            dw_land,
            dw_ocean,
            TURNRHCRIT_PARAM,
            self._minrhcrit,
            self._PLmb_at_klcl,
            PLmb,
            self._PLEmb_top,
            AREA,
            self._alpha,
        )

        self._hystpdf(
            DT_MOIST,
            self._alpha,
            PDFSHAPE,
            CNV_FRC,
            SRF_TYPE,
            PLmb,
            Q,
            QLLS,
            QLCN,
            QILS,
            QICN,
            T,
            CLLS,
            CLCN,
            NACTL,
            NACTI,
            self.rhx,
            self.qsat.ese,
            self.qsat.esw,
            self.qsat.esx,
            self.qsat.esw.view[0][12316],
            self.qsat.esw.view[0][8316],
        )

        if LMELTFRZ is True:
            self._meltfrz(DT_MOIST, CNV_FRC, SRF_TYPE, T, QLCN, QICN)
            self._meltfrz(DT_MOIST, CNV_FRC, SRF_TYPE, T, QLLS, QILS)

        if CCW_EVAP_EFF > 0.0:
            self._evap(
                DT_MOIST,
                CCW_EVAP_EFF,
                PLmb,
                T,
                Q,
                QLCN,
                QICN,
                CLCN,
                NACTL,
                NACTI,
                QST,
                self.evapc,
            )

        if CCI_EVAP_EFF > 0.0:
            self._subl(
                DT_MOIST,
                CCW_EVAP_EFF,
                PLmb,
                T,
                Q,
                QLCN,
                QICN,
                CLCN,
                NACTL,
                NACTI,
                QST,
                self.sublc,
            )

        self._fix_up_clouds(Q, T, QLLS, QILS, CLLS, QLCN, QICN, CLCN)
