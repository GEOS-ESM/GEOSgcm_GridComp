from f90nml import Namelist
from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Float, Int
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array
from pyMoist.convection.GF_2020.GF_functions.get_melting_profile import GetMeltingProfile


class TranslateGetMeltingProfile(TranslateFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "MELT_GLAC": {},
            "cumulus": {},
            "edto": {},
            "ierr": {},
            "melting_layer": {},
            "p_liq_ice": {},
            "po_cup": {},
            "pwdo": {},
            "pwo": {},
            "qrco": {},
            "tn_cup": {},
        }

        # Float/Int Inputs
        self.in_vars["parameters"] = []

        # FloatField Outputs
        self.out_vars = {
            "melting": self.grid.compute_dict(),
        }

    def compute(self, inputs):

        get_melting_profile = GetMeltingProfile(
            self.stencil_factory,
            self.grid.quantity_factory,
        )

        # Field inputs
        MELT_GLAC = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int
        )
        safe_assign_array(MELT_GLAC.view[:, :, :], inputs["MELT_GLAC"])

        cumulus = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int
        )
        safe_assign_array(cumulus.view[:, :, :], inputs["cumulus"])

        edto = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(edto.view[:, :, :], inputs["edto"])

        ierr = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int
        )
        safe_assign_array(ierr.view[:, :, :], inputs["ierr"])

        melting_layer = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(melting_layer.view[:, :, :], inputs["melting_layer"])

        p_liq_ice = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(p_liq_ice.view[:, :, :], inputs["p_liq_ice"])

        po_cup = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(po_cup.view[:, :, :], inputs["po_cup"])

        pwdo = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(pwdo.view[:, :, :], inputs["pwdo"])

        pwo = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(pwo.view[:, :, :], inputs["pwo"])

        qrco = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(qrco.view[:, :, :], inputs["qrco"])

        tn_cup = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(tn_cup.view[:, :, :], inputs["tn_cup"])

        melting = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")

        get_melting_profile(
            # In
            MELT_GLAC=MELT_GLAC,
            cumulus=cumulus,
            edto=edto,
            ierr=ierr,
            melting_layer=melting_layer,
            p_liq_ice=p_liq_ice,
            po_cup=po_cup,
            pwdo=pwdo,
            pwo=pwo,
            qrco=qrco,
            tn_cup=tn_cup,
            # Out
            melting=melting,
        )

        return {
            "melting": melting.view[:],
        }
