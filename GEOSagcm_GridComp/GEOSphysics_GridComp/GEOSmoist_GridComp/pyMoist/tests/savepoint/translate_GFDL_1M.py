from ndsl import Namelist, Quantity, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.GFDL_1M.GFDL_1M import evap_subl_pdf
from pyMoist.saturation.formulation import SaturationFormulation
import os
import xarray as xr
import numpy as np


class TranslateGFDL_1M(TranslateFortranData2Py):
    def __init__(self, grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory
        self._grid = grid
        self.max_error = 1e-9

        # # TESTING
        # # print(os.path.dirname(os.path.dirname(__file__)))
        # with xr.open_dataset(
        #     os.path.join(
        #         os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
        #         "pyMoist",
        #         "saturation",
        #         "netCDFs",
        #         "QSat_Tables.nc",
        #     )
        # ) as ds:
        #     ese_array = ds.data_vars["ese"].values[0, 0, :]
        #     esw_array = ds.data_vars["esw"].values[0, 0, :]
        #     esx_array = ds.data_vars["esx"].values[0, 0, :]

        # with xr.open_dataset("/home/charleskrop/netcdfs/QSat-Out.nc") as ds:
        #     estfrz = ds.data_vars["ESTFRZ_TEST"].values[0, 0, 0]
        #     estlqu = ds.data_vars["ESTLQU_TEST"].values[0, 0, 0]

        # esefrz_index = -999
        # eswfrz_index = -999
        # esxfrz_index = -999
        # eselqu_index = -999
        # eswlqu_index = -999
        # esxlqu_index = -999
        # for i in range(len(ese_array)):
        #     if ese_array[i] == estfrz:
        #         esefrz_index = i
        #     if esw_array[i] == estfrz:
        #         eswfrz_index = i
        #     if esx_array[i] == estfrz:
        #         esxfrz_index = i
        #     if ese_array[i] == estlqu:
        #         eselqu_index = i
        #     if esw_array[i] == estlqu:
        #         eswlqu_index = i
        #     if esx_array[i] == estlqu:
        #         esxlqu_index = i

        # print("esefrz_index: ", esefrz_index)
        # print("eswfrz_index: ", eswfrz_index)
        # print("esxfrz_index: ", esxfrz_index)
        # print("eselqu_index: ", eselqu_index)
        # print("eswlqu_index: ", eswlqu_index)
        # print("esxlqu_index: ", esxlqu_index)

        # with xr.open_dataset("/home/charleskrop/netcdfs/GFDL_1M-Out.nc") as ds:
        #     TESTVAR_1 = ds.data_vars["TESTVAR_1"].values[0, 0, :, :, :]
        #     TESTVAR_2 = ds.data_vars["TESTVAR_2"].values[0, 0, :, :, :]
        #     TESTVAR_3 = ds.data_vars["TESTVAR_3"].values[0, 0, :, :, :]
        #     TESTVAR_4 = ds.data_vars["TESTVAR_4"].values[0, 0, :, :, :]
        #     TESTVAR_5 = ds.data_vars["TESTVAR_5"].values[0, 0, :, :, :]
        #     if np.all(TESTVAR_2 == TESTVAR_5):
        #         print("YUH")
        #     else:
        #         print("NAH")
        # # END TESTING

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
            "RADIUS": {},
            "QCm": {},
        }

        # Float Inputs
        self.in_vars["parameters"] = [
            "dw_land",
            "dw_ocean",
            "TURNRHCRIT_PARAM",
            "DT_MOIST",
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
            # "TESTVAR_1": self.grid.compute_dict(),
            # "TESTVAR_2": self.grid.compute_dict(),
            # "TESTVAR_3": self.grid.compute_dict(),
            # "TESTVAR_4": self.grid.compute_dict(),
            # "TESTVAR_5": self.grid.compute_dict(),
            # "TESTVAR_6": self.grid.compute_dict(),
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
        QCm = self.make_ijk_field(inputs["QCm"])

        # Float Variables
        dw_land = Float(inputs["dw_land"])
        dw_ocean = Float(inputs["dw_ocean"])
        TURNRHCRIT_PARAM = Float(inputs["TURNRHCRIT_PARAM"])
        DT_MOIST = Float(inputs["DT_MOIST"])
        CCW_EVAP_EFF = Float(inputs["CCW_EVAP_EFF"])
        CCI_EVAP_EFF = Float(inputs["CCI_EVAP_EFF"])
        PDFSHAPE = Float(inputs["PDFSHAPE"])

        code = evap_subl_pdf(
            self.stencil_factory,
            self.quantity_factory,
        )

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
            QCm,
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
            # "TESTVAR_1": code.TESTVAR_1.view[:],
            # "TESTVAR_2": code.TESTVAR_2.view[:],
            # "TESTVAR_3": code.TESTVAR_3.view[:],
            # "TESTVAR_4": code.TESTVAR_4.view[:],
            # "TESTVAR_5": code.TESTVAR_5.view[:],
            # "TESTVAR_6": code.TESTVAR_6.view[:],
        }
