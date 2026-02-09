from f90nml import Namelist
from ndsl import Quantity, QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Float, FloatField, Int
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array
from pyMoist.UW.compute_uwshcu import setup_inputs, ComputeUwshcuInv
from pyMoist.UW.config import UWConfiguration


class TranslateSetupInputs(TranslateFortranData2Py):
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
            "PLE": {},
            "QICN": {},
            "QILS": {},
            "QLCN": {},
            "QLLS": {},
            "ZLE": {},
        }

        # FloatField Outputs
        self.out_vars = {
            "DP": self.grid.compute_dict(),
            "MASS": self.grid.compute_dict(),
            "PK": self.grid.compute_dict(),
            "PKE": self.grid.compute_dict(),
            "PL": self.grid.compute_dict(),
            "QITOT": self.grid.compute_dict(),
            "QLTOT": self.grid.compute_dict(),
            "RKFRE": self.grid.compute_dict(),
            "ZL0": self.grid.compute_dict(),
            "ZLE0": self.grid.compute_dict(),
        }

    def compute(self, inputs):

        _setup_inputs = self.stencil_factory.from_dims_halo(
            func=setup_inputs,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        # Inputs
        PLE = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(PLE.view[:, :, :], inputs["PLE"])
        QLLS = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(QLLS.view[:, :, :], inputs["QLLS"])
        QLCN = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(QLCN.view[:, :, :], inputs["QLCN"])
        QILS = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(QILS.view[:, :, :], inputs["QILS"])
        QICN = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(QICN.view[:, :, :], inputs["QICN"])
        QLLS = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(QLLS.view[:, :, :], inputs["QLLS"])
        ZLE = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(ZLE.view[:, :, :], inputs["ZLE"])

        # Outputs
        PKE = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        PL = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        PK = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        ZLE0 = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        ZL0 = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        DP = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        MASS = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        RKFRE = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM], units="n/a")
        QLTOT = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        QITOT = QuantityFactory.zeros(self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")

        _setup_inputs(
            PLE=PLE,
            QLLS=QLLS,
            QLCN=QLCN,
            QILS=QILS,
            QICN=QICN,
            ZLE=ZLE,
            PKE=PKE,
            PL=PL,
            PK=PK,
            ZLE0=ZLE0,
            ZL0=ZL0,
            DP=DP,
            MASS=MASS,
            RKFRE=RKFRE,
            QLTOT=QLTOT,
            QITOT=QITOT,
        )

        return {
            "DP": DP.view[:],
            "MASS": MASS.view[:],
            "PK": PK.view[:],
            "PKE": PKE.view[:],
            "PL": PL.view[:],
            "QITOT": QITOT.view[:],
            "QLTOT": QLTOT.view[:],
            "RKFRE": RKFRE.view[:],
            "ZL0": ZL0.view[:],
            "ZLE0": ZLE0.view[:],
        }
