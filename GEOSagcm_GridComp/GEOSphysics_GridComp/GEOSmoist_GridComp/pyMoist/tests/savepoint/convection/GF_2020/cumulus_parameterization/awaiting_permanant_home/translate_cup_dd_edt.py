from f90nml import Namelist
from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float, Int
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array
from pyMoist.convection.GF_2020.GF_functions.cup_dd_edt import CupDDEdt


class TranslateCupDDEdt(TranslateFortranData2Py):
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
            "ccn": {},
            "cumulus": {},
            "edtmax": {},
            "edtmin": {},
            "ierr": {},
            "kbcon": {},
            "ktop": {},
            "maxens2": {},
            "p": {},
            "psum2": {},
            "psumh": {},
            "pwav": {},
            "pwev": {},
            "us": {},
            "vs": {},
            "z": {},
            "aeroevap": {},
        }

        # Float/Int Inputs
        self.in_vars["parameters"] = []

        # FloatField Outputs
        self.out_vars = {
            "ierr": self.grid.compute_dict(),
            "edt": self.grid.compute_dict(),
            "edtc": self.grid.compute_dict(),
        }

    def compute(self, inputs):

        cup_dd_edt = CupDDEdt(
            self.stencil_factory,
            self.grid.quantity_factory,
        )

        # # Field inputs
        ccn = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(ccn.view[:, :, :], inputs["ccn"])

        cumulus = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=Int,
        )
        safe_assign_array(cumulus.view[:, :, :], inputs["cumulus"])

        edtmax = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(edtmax.view[:, :, :], inputs["edtmax"])

        edtmin = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(edtmin.view[:, :, :], inputs["edtmin"])

        ierr = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=Int,
        )
        safe_assign_array(ierr.view[:, :, :], inputs["ierr"])

        kbcon = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int
        )
        safe_assign_array(kbcon.view[:, :, :], inputs["kbcon"])

        ktop = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int
        )
        safe_assign_array(ktop.view[:, :, :], inputs["ktop"])

        maxens2 = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int
        )
        safe_assign_array(maxens2.view[:, :, :], inputs["maxens2"])

        p = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(p.view[:, :, :], inputs["p"])

        psum2 = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(psum2.view[:, :, :], inputs["psum2"])

        psumh = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(psumh.view[:, :, :], inputs["psumh"])

        pwav = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(pwav.view[:, :, :], inputs["pwav"])

        pwev = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(pwev.view[:, :, :], inputs["pwev"])

        us = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(us.view[:, :, :], inputs["us"])

        vs = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(vs.view[:, :, :], inputs["vs"])

        z = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(z.view[:, :, :], inputs["z"])

        aeroevap = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=Int,
        )
        safe_assign_array(aeroevap.view[:, :, :], inputs["aeroevap"])

        edt = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )

        edtc = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )

        cup_dd_edt(
            # In
            ccn=ccn,
            cumulus=cumulus,
            edtmax=edtmax,
            edtmin=edtmin,
            kbcon=kbcon,
            ktop=ktop,
            maxens2=maxens2,
            p=p,
            psum2=psum2,
            psumh=psumh,
            pwav=pwav,
            pwev=pwev,
            us=us,
            vs=vs,
            z=z,
            aeroevap=aeroevap,
            # Out
            edt=edt,
            edtc=edtc,
            ierr=ierr,
        )

        return {
            "edt": edt.view[:],
            "ierr": ierr.view[:],
            "edtc": edtc.view[:],
        }
