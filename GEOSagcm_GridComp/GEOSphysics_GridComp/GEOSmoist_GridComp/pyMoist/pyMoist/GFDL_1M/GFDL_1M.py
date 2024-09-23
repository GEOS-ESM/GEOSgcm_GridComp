import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import (
    computation,
    interval,
    PARALLEL,
    FORWARD,
    atan,
    sin,
    tan,
    sqrt,
    tanh,
    exp,
    log10,
)
from ndsl.boilerplate import get_factories_single_tile_numpy
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField, FloatFieldIJ, Float, IntField, IntFieldIJ
from ndsl import StencilFactory, QuantityFactory, orchestrate
import numpy as np
import pyMoist.pyMoist_constants as constants
from pyMoist.saturation.formulation import SaturationFormulation
from pyMoist.saturation.qsat import QSat
from .GFDL_1M_util import (
    get_last,
    hybrid_index_2dout,
    initial_calc,
    hystpdf,
    meltfrz,
    evap,
    subl,
)


class evap_subl_pdf:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        formulation: SaturationFormulation = SaturationFormulation.Staars,
        use_table_lookup: bool = True,
    ):
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
        self._tmp = self.quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        self._minrhcrit = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._PLEmb_top = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._PLmb_at_klcl = self.quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        self._halo = self.stencil_factory.grid_indexing.n_halo
        self._k_mask = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._alpha = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._evapc = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        for i in range(0, self._k_mask.view[:].shape[0]):
            for j in range(0, self._k_mask.view[:].shape[1]):
                for k in range(0, self._k_mask.view[:].shape[2]):
                    self._k_mask.view[i, j, k] = k + 1

        # TESTING
        self.TESTVAR_1 = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.TESTVAR_2 = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.TESTVAR_3 = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.TESTVAR_4 = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.TESTVAR_5 = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.TESTVAR_6 = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        self.TESTVAR_1.view[:] = -999.0
        self.TESTVAR_2.view[:] = -999.0
        self.TESTVAR_3.view[:] = -999.0
        self.TESTVAR_4.view[:] = -999.0
        self.TESTVAR_5.view[:] = -999.0
        self.TESTVAR_6.view[:] = -999.0
        # END TESTING

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
        QCm: FloatField,
    ):
        # Theoretically, the following stencil and for loop should provide the same (correct) result. However, they both provide different incorrect results
        # The for loop is currently closer to being correct, with only the k level incorrect, so I am using that.
        self._get_last(PLEmb, self._tmp, self._PLEmb_top)

        # Temporary implementation of hybrid_index_2dout.py, perhaps not working as indended (backend issue), will need to be addressed at later date
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
            self.qsat.ese,
            self.qsat.esw,
            self.qsat.esx,
            self.qsat.esw.view[0][12316],
            self.qsat.esw.view[0][8316],
            self.TESTVAR_1,
            self.TESTVAR_2,
            self.TESTVAR_3,
            self.TESTVAR_4,
            self.TESTVAR_5,
            self.TESTVAR_6,
        )

        # self._meltfrz(DT_MOIST, CNV_FRC, SRF_TYPE, T, QLCN, QICN)
        # self._meltfrz(DT_MOIST, CNV_FRC, SRF_TYPE, T, QLLS, QILS)

        # RHCRIT = Float(1.0)
        # self._evap(DT_MOIST, CCW_EVAP_EFF, RHCRIT, PLmb, T, Q, QLCN, QICN, CLCN, NACTL, NACTI, QST, self._evapc, QCm)

        # self._subl(DT_MOIST, CCW_EVAP_EFF, RHCRIT, PLmb, T, Q, QLCN, QICN, CLCN, NACTL, NACTI, QST, self._evapc)
