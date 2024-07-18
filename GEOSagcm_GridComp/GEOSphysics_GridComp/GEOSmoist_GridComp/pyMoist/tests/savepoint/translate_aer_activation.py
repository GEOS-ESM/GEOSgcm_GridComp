from ndsl import Namelist, StencilFactory
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.aer_activation import AerActivation
from ndsl.constants import X_DIM, Y_DIM, Z_DIM

class TranslateAerActivation(TranslateFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: Namelist, 
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory
        self._grid = grid
        self.max_error = 1e-9

        #FloatField Inputs
        self.in_vars["data_vars"] = {
            "AERO_F_DUST": {"names_4d": ["n_modes"]},
            "AERO_F_ORGANIC": {"names_4d": ["n_modes"]},
            "TMP3D": {},
            "NACTL": {},
            "PLE": {},
            "AERO_F_SOOT": {"names_4d": ["n_modes"]},
            "T": {},
            "AERO_HYGROSCOPICITY": {"names_4d": ["n_modes"]},
            "ZLE0": {},
            "QLLS": {},
            "AERO_NUM": {"names_4d": ["n_modes"]},
            "ZL0": {},
            "PLmb": {},
            "QILS": {},
            "AERO_DGN": {"names_4d": ["n_modes"]},
            "TKE": {},
            "NACTI": {},
            "NWFA": {},
            "QLCN": {},
            "AERO_SIGMA": {"names_4d": ["n_modes"]},
            "AERO_DENSITY": {"names_4d": ["n_modes"]},
            "QICN": {},
        }


        #Float/Int Inputs
        self.in_vars["parameters"] = ["n_modes", "SH", "EVAP", "KPBL", "FRLAND", "CCN_LND", "CCN_OCN"]

        #FloatField Outputs
        self.out_vars["data_vars"] = {
            "NACTL": {},
            "NACTI": {},
            "NWFA": {},
        }

     #Calculated Outputs
    def compute_from_storage(self, inputs):
        aer_activation_call = AerActivation(
            self.stencil_factory,
            self.quantity_factory,
            inputs["n_modes"],
            USE_AERSOL_NN = True
        )
        outputs = {
            "NACTL": self._grid.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="unknown"),
            "NACTI": self._grid.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="unknown"),
            "NWFA": self._grid.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="unknown")
        }
        aer_activation_call(**inputs, **outputs)
        return {**outputs}