from ndsl import Namelist, StencilFactory
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.radiation_coupling import RadiationCoupling
from ndsl.constants import X_DIM, Y_DIM, Z_DIM


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
            do_qa=namelist.do_qa,
        )
        self._grid = grid
        self.max_error = 1e-9

        # ADAPT BELOW TO INPUTS
        #
        # fillq_info = self.grid.compute_dict()
        # fillq_info["serialname"] = "fq"
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
        }
        #Inputs that are Floats
        self.in_vars["parameters"] = ['FAC_RL', 'MIN_RL', 'MAX_RL', 'FAC_RI', 'MIN_RI', 'MAX_RI']

        #Outputs
        self.out_vars = {
            "Q": self.grid.compute_dict(),
            "T": self.grid.compute_dict(),
            "QLLS": self.grid.compute_dict(),
            "QILS": self.grid.compute_dict(),
            "CLLS": self.grid.compute_dict(),
            "QLCN": self.grid.compute_dict(),
            "QICN": self.grid.compute_dict(),
            "CLCN": self.grid.compute_dict(),
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

    def compute_from_storage(self, inputs):
        #self._grid.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="unknown")
        outputs = {
            "RAD_QV": self._grid.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="unknown"),
            "RAD_QL": self._grid.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="unknown"),
            "RAD_QI": self._grid.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="unknown"),
            "RAD_QR": self._grid.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="unknown"),
            "RAD_QS": self._grid.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="unknown"),
            "RAD_QG": self._grid.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="unknown"),
            "RAD_CF": self._grid.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="unknown"),
            "CLDREFFI": self._grid.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="unknown"),
            "CLDREFFL": self._grid.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="unknown")
        }
        self.compute_func(**inputs, **outputs)
        return {**outputs, "T": inputs["T"], "Q": inputs["Q"], "QLLS": inputs["QLLS"], "QILS": inputs["QILS"], "CLLS": inputs["CLLS"], "QLCN": inputs["QLCN"], "QICN": inputs["QICN"], "CLCN": inputs["CLCN"]}
