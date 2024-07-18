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
        aero_f_organic = inputs["AERO_F_ORGANIC"].astype(Float)
        tmp3d = inputs["TMP3D"].astype(Float) 
        nactl = inputs["NACTL"].astype(Float)
        ple = inputs["PLE"].astype(Float)
        aero_f_soot = inputs["AERO_F_SOOT"].astype(Float)
        t = inputs["T"].astype(Float)
        aero_hygroscopicity = inputs["AERO_HYGROSCOPICITY"].astype(Float)
        zle0 = inputs["ZLE0"].astype(Float)
        qlls = inputs["QLLS"].astype(Float)
        aero_num = inputs["AERO_F_NUM"].astype(Float)
        zl0 = inputs["ZL0"].astype(Float)
        plmb = inputs["PLmb"].astype(Float)
        qils = inputs["QILS"].astype(Float)
        aero_dgn = inputs["AERO_F_DGN"].astype(Float)
        tke = inputs["TKE"].astype(Float)
        nacti = inputs["NACTI"].astype(Float)
        nwfa = inputs["NWFA"].astype(Float)
        qlcn = inputs["QLCN"].astype(Float)
        aero_sigma = inputs["AERO_F_SIGMA"].astype(Float)
        aero_density = inputs["AERO_F_DENSITY"].astype(Float)
        qicn = inputs["QICN"].astype(Float)

        aer_activation.ddim_debug(aero_f_dust)
        aer_activation.ddim_debug(aero_f_organic)
        aer_activation.ddim_debug(tmp3d)
        aer_activation.ddim_debug(nactl)
        aer_activation.ddim_debug(ple)
        aer_activation.ddim_debug(aero_f_soot)
        aer_activation.ddim_debug(t)
        aer_activation.ddim_debug(aero_hygroscopicity)
        aer_activation.ddim_debug(zle0)
        aer_activation.ddim_debug(qlls)
        aer_activation.ddim_debug(aero_num)
        aer_activation.ddim_debug(zl0)
        aer_activation.ddim_debug(plmb)
        aer_activation.ddim_debug(qils)
        aer_activation.ddim_debug(aero_dgn)
        aer_activation.ddim_debug(tke)
        aer_activation.ddim_debug(nacti)
        aer_activation.ddim_debug(nwfa)
        aer_activation.ddim_debug(qlcn)
        aer_activation.ddim_debug(aero_sigma)
        aer_activation.ddim_debug(aero_density)
        aer_activation.ddim_debug(qicn)

        return {
            "NACTL": nactl,
            "NACTI": nacti,
            "NWFA": nwfa,
        }

