from ndsl import Namelist, StencilFactory
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.radiation_coupling import RadiationCoupling


class TranslateRadCouple(TranslateFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.compute_func = RadiationCoupling(  # type: ignore
            self.stencil_factory,
            self.grid.quantity_factory,
            do_qa=False,  # Change to do_qa=namelist.do_qa,
            # if QSAT module procedures QSAT0 and QSAT3 are implemented
        )
        self._grid = grid

        self.max_error = 1e-10

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "Q": self.grid.compute_dict(),
            "T": self.grid.compute_dict(),
            "QLLS": self.grid.compute_dict(),
            "QILS": self.grid.compute_dict(),
            "CLLS": self.grid.compute_dict(),
            "QLCN": self.grid.compute_dict(),
            "QICN": self.grid.compute_dict(),
            "CLCN": self.grid.compute_dict(),
            "PLmb": self.grid.compute_dict(),
            "QRAIN": self.grid.compute_dict(),
            "QSNOW": self.grid.compute_dict(),
            "QGRAUPEL": self.grid.compute_dict(),
            "NACTL": self.grid.compute_dict(),
            "NACTI": self.grid.compute_dict(),
            "RAD_QV": self.grid.compute_dict(),
            "RAD_QL": self.grid.compute_dict(),
            "RAD_QI": self.grid.compute_dict(),
            "RAD_QR": self.grid.compute_dict(),
            "RAD_QS": self.grid.compute_dict(),
            "RAD_QG": self.grid.compute_dict(),
            "RAD_CF": self.grid.compute_dict(),
            "CLDREFFI": self.grid.compute_dict(),
            "CLDREFFL": self.grid.compute_dict(),
        }

        # Float Inputs
        self.in_vars["parameters"] = [
            "FAC_RL",
            "MIN_RL",
            "MAX_RL",
            "FAC_RI",
            "MIN_RI",
            "MAX_RI",
        ]

        # FloatField Outputs
        self.out_vars = {
            "RAD_QV": self.grid.compute_dict(),
            "RAD_QL": self.grid.compute_dict(),
            "RAD_QI": self.grid.compute_dict(),
            "RAD_QR": self.grid.compute_dict(),
            "RAD_QS": self.grid.compute_dict(),
            "RAD_QG": self.grid.compute_dict(),
            "RAD_CF": self.grid.compute_dict(),
            "CLDREFFI": self.grid.compute_dict(),
            "CLDREFFL": self.grid.compute_dict(),
        }

    def compute(self, inputs):
        self.make_storage_data_input_vars(inputs)
        self.compute_func(**inputs)
        return self.slice_output(inputs)
