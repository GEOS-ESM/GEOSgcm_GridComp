from ndsl import Namelist, QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Float, Int
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array
from pyMoist.convection.GF_2020.GF_functions.get_partition_liq_ice import (
    GetPartitionLiqIce,
)


class TranslateGetPartitionLiqIce(TranslateFortranData2Py):
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
            "ierr": {},
            "cnvfrc": {},
            "srftype": {},
            "tn": {},
            "z1": {},
            "zo_cup": {},
            "po_cup": {},
            "FRAC_MODIS": {},
        }

        # Float/Int Inputs
        self.in_vars["parameters"] = []

        # FloatField Outputs
        self.out_vars = {
            "melting_layer": self.grid.compute_dict(),
            "p_liq_ice": self.grid.compute_dict(),
        }

    def compute(self, inputs):

        get_partition_liq_ice = GetPartitionLiqIce(
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

        FRAC_MODIS = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int
        )
        safe_assign_array(FRAC_MODIS.view[:, :, :], inputs["FRAC_MODIS"])

        cnvfrc = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(cnvfrc.view[:, :, :], inputs["cnvfrc"])

        ierr = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int
        )
        safe_assign_array(ierr.view[:, :, :], inputs["ierr"])

        po_cup = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a"
        )
        safe_assign_array(po_cup.view[:, :, :], inputs["po_cup"])

        z1 = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a"
        )
        safe_assign_array(z1.view[:, :, :], inputs["z1"])

        zo_cup = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a"
        )
        safe_assign_array(zo_cup.view[:, :, :], inputs["zo_cup"])

        srftype = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a"
        )
        safe_assign_array(srftype.view[:, :, :], inputs["srftype"])

        tn = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a"
        )
        safe_assign_array(tn.view[:, :, :], inputs["tn"])

        melting_layer = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a"
        )
        p_liq_ice = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a"
        )

        get_partition_liq_ice(
            # In
            MELT_GLAC=MELT_GLAC,
            cumulus=cumulus,
            ierr=ierr,
            cnvfrc=cnvfrc,
            srftype=srftype,
            tn=tn,
            z1=z1,
            zo_cup=zo_cup,
            po_cup=po_cup,
            FRAC_MODIS=FRAC_MODIS,
            # Out
            melting_layer=melting_layer,
            p_liq_ice=p_liq_ice,
        )

        return {"melting_layer": melting_layer.view[:], "p_liq_ice": p_liq_ice.view[:]}
