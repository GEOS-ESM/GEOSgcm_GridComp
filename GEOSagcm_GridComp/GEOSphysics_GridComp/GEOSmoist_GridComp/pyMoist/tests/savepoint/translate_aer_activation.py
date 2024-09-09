from ndsl import Namelist, Quantity, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float
from ndsl.stencils.testing.translate import TranslateFortranData2Py

from pyMoist.aer_activation import AerActivation


class TranslateAerActivation(TranslateFortranData2Py):
    def __init__(self, grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory
        self._grid = grid

        self.nmodes_quantity_factory = AerActivation.make_nmodes_quantity_factory(
            self.quantity_factory
        )

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
            "FRLAND": {},
        }

        # Float/Int Inputs
        self.in_vars["parameters"] = [
            "n_modes",
            "SH",
            "EVAP",
            "KPBL",
            "CCN_LND",
            "CCN_OCN",
        ]

        # FloatField Outputs
        self.out_vars = {
            "NACTL": self.grid.compute_dict(),
            "NACTI": self.grid.compute_dict(),
            "NWFA": self.grid.compute_dict(),
        }

    def make_ijk_field(self, data) -> Quantity:
        qty = self.quantity_factory.empty(
            [X_DIM, Y_DIM, Z_DIM],
            "n/a",
        )
        qty.view[:, :, :] = qty.np.asarray(data[:, :, :])
        return qty

    def make_ij_field(self, data) -> Quantity:
        qty = self.quantity_factory.empty(
            [X_DIM, Y_DIM],
            "n/a",
        )
        qty.view[:, :] = qty.np.asarray(data[:, :])
        return qty

    def make_nmodes_ijk_field(self, data) -> Quantity:
        qty = self.nmodes_quantity_factory.empty(
            [X_DIM, Y_DIM, Z_DIM, "n_modes"],
            "n/a",
        )
        qty.view[:, :, :, :] = qty.np.asarray(data[:, :, :, :])
        return qty

    def compute(self, inputs):
        aer_activation = AerActivation(
            self.stencil_factory,
            self.quantity_factory,
            int(inputs["n_modes"]),
            USE_AERSOL_NN=True,
        )

        # Outputs
        nactl = self.make_ijk_field(inputs["NACTL"])
        nacti = self.make_ijk_field(inputs["NACTI"])
        nwfa = self.make_ijk_field(inputs["NWFA"])

        # Inputs
        aero_hygroscopicity = self.make_nmodes_ijk_field(inputs["AERO_HYGROSCOPICITY"])
        aero_sigma = self.make_nmodes_ijk_field(inputs["AERO_SIGMA"])
        aero_dgn = self.make_nmodes_ijk_field(inputs["AERO_DGN"])
        aero_num = self.make_nmodes_ijk_field(inputs["AERO_NUM"])

        frland = self.make_ij_field(inputs["FRLAND"])

        tmp3d = self.make_ijk_field(inputs["TMP3D"])
        t = self.make_ijk_field(inputs["T"])
        qlls = self.make_ijk_field(inputs["QLLS"])
        plmb = self.make_ijk_field(inputs["PLmb"])
        qils = self.make_ijk_field(inputs["QILS"])
        # TKE is a 73 level Field, but level 0 is never indexed
        tke = self.make_ijk_field(inputs["TKE"][:, :, 1:])
        nacti = self.make_ijk_field(inputs["NACTI"])
        qlcn = self.make_ijk_field(inputs["QLCN"])
        qicn = self.make_ijk_field(inputs["QICN"])

        n_modes = inputs["n_modes"]
        ccn_lnd = Float(inputs["CCN_LND"])
        ccn_ocn = Float(inputs["CCN_OCN"])

        # Unused - but present in the original Fortran
        # aero_f_dust = inputs["AERO_F_DUST"].astype(Float)
        # aero_f_organic = inputs["AERO_F_ORGANIC"].astype(Float)
        # aero_f_soot = inputs["AERO_F_SOOT"].astype(Float)
        # aero_density = inputs["AERO_DENSITY"].astype(Float)
        # ple = inputs["PLE"].astype(Float)
        # zle0 = inputs["ZLE0"].astype(Float)
        # zl0 = inputs["ZL0"].astype(Float)
        # sh = inputs["SH"].astype(Float)
        # evap = inputs["EVAP"].astype(Float)
        # kpbl = inputs["KPBL"].astype(Float)

        plmb.view[:, :, :] = plmb.view[:, :, :] * 100.0

        aer_activation(
            aero_dgn=aero_dgn,
            aero_num=aero_num,
            nacti=nacti,
            t=t,
            plo=plmb,
            qicn=qicn,
            qils=qils,
            qlcn=qlcn,
            qlls=qlls,
            nn_land=Float(ccn_lnd * 1.0e6),
            frland=frland,
            nn_ocean=Float(ccn_ocn * 1.0e6),
            aero_hygroscopicity=aero_hygroscopicity,
            nwfa=nwfa,
            nactl=nactl,
            vvel=tmp3d,
            tke=tke,
            aero_sigma=aero_sigma,
        )

        return {
            "NACTL": nactl.view[:, :, :],
            "NACTI": nacti.view[:, :, :],
            "NWFA": nwfa.view[:, :, :],
        }
