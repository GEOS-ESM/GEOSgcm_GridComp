from ndsl import Namelist, StencilFactory
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.constants import X_DIM, Y_DIM, Z_DIM

from ndsl import Namelist, StencilFactory
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.evap_subl_pdf import evap_subl_pdf
from ndsl.constants import X_DIM, Y_DIM, Z_DIM

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
        #FloatField Inputs
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
            "testminrhcrit": {},
            "testPLEmb": {},
            "testturnrhcrit": {},
            "testrhcritpart1": {},
            "testrhcritpart2": {},
            "testrhcrit": {},
            "testalpha": {},
            "testtancalc": {},
        }

        #Float Inputs
        self.in_vars["parameters"] = ["dw_land", "dw_ocean", "TURNRHCRIT_PARAM", "DT_MOIST"]

        #FloatField Outputs
        self.out_vars = {
            "testminrhcrit": self.grid.compute_dict(),
            "testPLEmb": self.grid.compute_dict(),
            "testturnrhcrit": self.grid.compute_dict(),
            "testrhcritpart1": self.grid.compute_dict(),
            "testrhcritpart2": self.grid.compute_dict(),
            "testrhcrit": self.grid.compute_dict(),
            "testalpha": self.grid.compute_dict(),
            "testtancalc": self.grid.compute_dict(),
            "T": self.grid.compute_dict(),
            "QLCN": self.grid.compute_dict(),
            "QICN": self.grid.compute_dict(),
        }

    #Calculated Outputs
    def compute_from_storage(self, inputs):
        self.compute_func(**inputs)
        return {"testminrhcrit": inputs["testminrhcrit"], "testPLEmb": inputs["testPLEmb"], "testturnrhcrit": inputs["testturnrhcrit"],
                "testrhcritpart1": inputs["testrhcritpart1"], "testrhcritpart2": inputs["testrhcritpart2"],
                "testrhcrit": inputs["testrhcrit"], "testalpha": inputs["testalpha"], "testtancalc": inputs["testtancalc"],
                "T": inputs["T"], "QLCN": inputs["QLCN"], "QICN": inputs["QICN"]}
