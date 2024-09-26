from ndsl import Namelist, Quantity, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.GFDL_1M.GFDL_1M import evap_subl_pdf


class TranslateGFDL_1M(TranslateFortranData2Py):
    def __init__(self, grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory
        self._grid = grid
        self.max_error = 1e-9

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "EIS": {},
            "PLmb": {},
            "KLCL": {},
            "PLEmb": {},
            "AREA": {},
            "CNV_FRC": {},
            "SRF_TYPE": {},
            "T": {},
            "QLCN": {},
            "QICN": {},
            "QLLS": {},
            "QILS": {},
            "Q": {},
            "CLLS": {},
            "CLCN": {},
            "NACTL": {},
            "NACTI": {},
            "QST": {},
        }

        # Float Inputs
        self.in_vars["parameters"] = [
            "dw_land",
            "dw_ocean",
            "TURNRHCRIT_PARAM",
            "DT_MOIST_Float",
            "CCW_EVAP_EFF",
            "CCI_EVAP_EFF",
            "PDFSHAPE",
        ]

        # FloatField Outputs
        self.out_vars = {
            "T": self.grid.compute_dict(),
            "Q": self.grid.compute_dict(),
            "QLLS": self.grid.compute_dict(),
            "QILS": self.grid.compute_dict(),
            "CLLS": self.grid.compute_dict(),
            "QLCN": self.grid.compute_dict(),
            "QICN": self.grid.compute_dict(),
            "CLCN": self.grid.compute_dict(),
            "RHX": self.grid.compute_dict(),
            "EVAPC": self.grid.compute_dict(),
            "SUBLC": self.grid.compute_dict(),
        }

    def make_ij_field(self, data) -> Quantity:
        qty = self.quantity_factory.empty(
            [X_DIM, Y_DIM],
            "n/a",
        )
        qty.view[:, :] = qty.np.asarray(data[:, :])
        return qty

    def make_ijk_field(self, data) -> Quantity:
        qty = self.quantity_factory.empty(
            [X_DIM, Y_DIM, Z_DIM],
            "n/a",
        )
        qty.view[:, :, :] = qty.np.asarray(data[:, :, :])
        return qty

    def compute(self, inputs):
        # FloatField Variables
        EIS = self.make_ij_field(inputs["EIS"])
        PLmb = self.make_ijk_field(inputs["PLmb"])
        KLCL = self.make_ij_field(inputs["KLCL"])
        PLEmb = self.make_ijk_field(inputs["PLEmb"])
        AREA = self.make_ij_field(inputs["AREA"])
        CNV_FRC = self.make_ij_field(inputs["CNV_FRC"])
        SRF_TYPE = self.make_ij_field(inputs["SRF_TYPE"])
        T = self.make_ijk_field(inputs["T"])
        QLCN = self.make_ijk_field(inputs["QLCN"])
        QICN = self.make_ijk_field(inputs["QICN"])
        QLLS = self.make_ijk_field(inputs["QLLS"])
        QILS = self.make_ijk_field(inputs["QILS"])
        Q = self.make_ijk_field(inputs["Q"])
        CLLS = self.make_ijk_field(inputs["CLLS"])
        CLCN = self.make_ijk_field(inputs["CLCN"])
        NACTL = self.make_ijk_field(inputs["NACTL"])
        NACTI = self.make_ijk_field(inputs["NACTI"])
        QST = self.make_ijk_field(inputs["QST"])

        # Float Variables
        dw_land = Float(inputs["dw_land"])
        dw_ocean = Float(inputs["dw_ocean"])
        TURNRHCRIT_PARAM = Float(inputs["TURNRHCRIT_PARAM"])
        DT_MOIST = Float(inputs["DT_MOIST_Float"])
        CCW_EVAP_EFF = Float(inputs["CCW_EVAP_EFF"])
        CCI_EVAP_EFF = Float(inputs["CCI_EVAP_EFF"])
        PDFSHAPE = Float(inputs["PDFSHAPE"])

        code = evap_subl_pdf(self.stencil_factory, self.quantity_factory)

        code(
            EIS,
            dw_land,
            dw_ocean,
            PDFSHAPE,
            TURNRHCRIT_PARAM,
            PLmb,
            KLCL,
            PLEmb,
            AREA,
            DT_MOIST,
            CNV_FRC,
            SRF_TYPE,
            T,
            QLCN,
            QICN,
            QLLS,
            QILS,
            CCW_EVAP_EFF,
            CCI_EVAP_EFF,
            Q,
            CLLS,
            CLCN,
            NACTL,
            NACTI,
            QST,
        )

        return {
            "T": T.view[:],
            "Q": Q.view[:],
            "QLLS": QLLS.view[:],
            "QILS": QILS.view[:],
            "CLLS": CLLS.view[:],
            "QLCN": QLCN.view[:],
            "QICN": QICN.view[:],
            "CLCN": CLCN.view[:],
            "RHX": code.rhx.view[:],
            "EVAPC": code.evapc.view[:],
            "SUBLC": code.sublc.view[:],
        }
