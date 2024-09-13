from ndsl import Namelist, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.evap_subl_pdf import evap_subl_pdf


class Translateevap_subl_pdf(TranslateFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.compute_func = evap_subl_pdf(
            self.stencil_factory,
            self.grid.quantity_factory,
        )
        print(grid.__dict__)
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
        ]

        # FloatField Outputs
        self.out_vars = {
            "T": self.grid.compute_dict(),
            "QLCN": self.grid.compute_dict(),
            "QICN": self.grid.compute_dict(),
            "QLLS": self.grid.compute_dict(),
            "QILS": self.grid.compute_dict(),
            "Q": self.grid.compute_dict(),
            "CLCN": self.grid.compute_dict(),
            "NACTL": self.grid.compute_dict(),
            "NACTI": self.grid.compute_dict(),
            "QST": self.grid.compute_dict(),
            "RADIUS": self.grid.compute_dict(),
            "QCm": self.grid.compute_dict(),
        }

    # Calculated Outputs
    def compute_from_storage(self, inputs):
        self.compute_func(**inputs)
        return {
            "T": inputs["T"],
            "QLCN": inputs["QLCN"],
            "QICN": inputs["QICN"],
            "QLLS": inputs["QLLS"],
            "QILS": inputs["QILS"],
            "Q": inputs["Q"],
            "CLCN": inputs["CLCN"],
            "NACTL": inputs["NACTL"],
            "NACTI": inputs["NACTI"],
            "QST": inputs["QST"],
            "RADIUS": inputs["RADIUS"],
            "QCm": inputs["QCm"],
        }
