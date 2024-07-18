from ndsl import Namelist, StencilFactory
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.aer_activation import AerActivation
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float


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

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "AERO_F_DUST": {},
            "AERO_F_ORGANIC": {},
            "TMP3D": {},
            "NACTL": {},
            "PLE": {},
            "AERO_F_SOOT": {},
            "T": {},
            "AERO_HYGROSCOPICITY": {},
            "ZLE0": {},
            "QLLS": {},
            "AERO_NUM": {},
            "ZL0": {},
            "PLmb": {},
            "QILS": {},
            "AERO_DGN": {},
            "TKE": {},
            "NACTI": {},
            "NWFA": {},
            "QLCN": {},
            "AERO_SIGMA": {},
            "AERO_DENSITY": {},
            "QICN": {},
        }

        # Float/Int Inputs
        self.in_vars["parameters"] = [
            "n_modes",
            "SH",
            "EVAP",
            "KPBL",
            "FRLAND",
            "CCN_LND",
            "CCN_OCN",
        ]

        # FloatField Outputs
        self.out_vars = {
            "NACTL": self.grid.compute_dict(),
            "NACTI": self.grid.compute_dict(),
            "NWFA": self.grid.compute_dict(),
        }

    def compute(self, inputs):
        aer_activation = AerActivation(
            self.stencil_factory,
            self.quantity_factory,
            int(inputs["n_modes"]),
            USE_AERSOL_NN=True,
        )

        # Outputs
        nactl = inputs["NACTL"].astype(Float)
        nacti = inputs["NACTI"].astype(Float)
        nwfa = inputs["NWFA"].astype(Float)

        # Inputs
        aero_f_dust = inputs["AERO_F_DUST"].astype(Float)

        aer_activation.ddim_debug(aero_f_dust)

        return {
            "NACTL": nactl,
            "NACTI": nacti,
            "NWFA": nwfa,
        }
