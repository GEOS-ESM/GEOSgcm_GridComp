from f90nml import Namelist
from gt4py.cartesian.gtscript import int32

import pyMoist.constants as constants
from ndsl import StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Float, Int
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array
from pyMoist.saturation_tables import get_saturation_vapor_pressure_table
from pyMoist.UW.compute_uwshcu import calc_entrainment_mass_flux
from pyMoist.UW.config import UWConfiguration


class TranslateCalcEntrainmentMassFlux(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "condensation": {},
            "dp0": {},
            "exnifc0": {},
            "kbup": {},
            "pifc0": {},
            "pmid0": {},
            "ppen": {},
            "qt0": {},
            "qtu": {},
            "rei": {},
            "ssqt0": {},
            "ssthl0": {},
            "sstr0": {},
            "ssu0": {},
            "ssv0": {},
            "thl0": {},
            "thlu": {},
            "thv0bot": {},
            "thv0top": {},
            "tru": {},
            "u0": {},
            "v0": {},
            "uu": {},
            "vu": {},
            "umf": {},
            "kpen": {},
            "tr0_CalcEntrain": {},
        }

        # FloatField Outputs
        self.out_vars = {
            "qtu_emf": self.grid.compute_dict(),
            "thlu_emf": self.grid.compute_dict(),
            "tru_emf": self.grid.compute_dict(),
            "uu_emf": self.grid.compute_dict(),
            "vu_emf": self.grid.compute_dict(),
            "emf": self.grid.compute_dict(),
        }

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("ComputeUwshcuInv-constants")

    def compute(self, inputs):
        config = UWConfiguration(**self.constants)

        self.quantity_factory.add_data_dimensions(
            {
                "ntracers": constants.NCNST,
            }
        )

        self._calc_entrainment_mass_flux = self.stencil_factory.from_dims_halo(
            func=calc_entrainment_mass_flux,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "ncnst": config.NCNST,
                "dotransport": config.dotransport,
                "rpen": config.rpen,
                "dt": config.dt,
                "use_cumpenent": config.use_cumpenent,
            },
        )

        # Inputs
        condensation = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a", dtype=bool)

        for i in range(0, 24):
            for j in range(0, 24):
                if inputs["condensation"][i, j] == 1:
                    condensation.view[i, j] = False
                else:
                    condensation.view[i, j] = True

        exnifc0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(exnifc0.view[:], inputs["exnifc0"])
        kpen = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int)
        safe_assign_array(kpen.view[:], inputs["kpen"] - 1)
        pifc0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(pifc0.view[:], inputs["pifc0"])
        ppen = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        safe_assign_array(ppen.view[:], inputs["ppen"])
        qt0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(qt0.view[:], inputs["qt0"])
        qtu = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(qtu.view[:], inputs["qtu"])
        ssqt0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(ssqt0.view[:], inputs["ssqt0"])
        ssthl0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(ssthl0.view[:], inputs["ssthl0"])
        thl0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(thl0.view[:], inputs["thl0"])
        thlu = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(thlu.view[:], inputs["thlu"])
        thv0bot = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(thv0bot.view[:], inputs["thv0bot"])
        thv0top = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(thv0top.view[:], inputs["thv0top"])
        uu = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(uu.view[:], inputs["uu"])
        vu = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(vu.view[:], inputs["vu"])
        tru = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM, "ntracers"], units="n/a")
        safe_assign_array(tru.view[:], inputs["tru"])
        kbup = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int)
        safe_assign_array(kbup.view[:], inputs["kbup"] - 1)
        umf = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(umf.view[:], inputs["umf"])
        rei = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(rei.view[:], inputs["rei"])
        dp0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(dp0.view[:], inputs["dp0"])
        pmid0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(pmid0.view[:], inputs["pmid0"])
        u0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(u0.view[:], inputs["u0"])
        v0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(v0.view[:], inputs["v0"])
        ssu0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(ssu0.view[:], inputs["ssu0"])
        ssv0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(ssv0.view[:], inputs["ssv0"])
        tr0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM, "ntracers"], units="n/a")
        safe_assign_array(tr0.view[:], inputs["tr0_CalcEntrain"])
        sstr0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM, "ntracers"], units="n/a")
        safe_assign_array(sstr0.view[:], inputs["sstr0"])

        # Outputs
        qtu_emf = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        thlu_emf = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        tru_emf = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM, "ntracers"], units="n/a")
        uu_emf = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        vu_emf = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        emf = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")

        saturation_vapor_pressure_table = get_saturation_vapor_pressure_table(self.stencil_factory.backend)
        self.ese = saturation_vapor_pressure_table.ese
        self.esx = saturation_vapor_pressure_table.esx

        umf_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        qtflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        slflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        uflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        vflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        cush_inout = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")

        # The iteration you want to test
        iter_test = int32(0)

        # # Call stencils
        self._calc_entrainment_mass_flux(
            condensation=condensation,
            thlu=thlu,
            qtu=qtu,
            uu=uu,
            vu=vu,
            tru=tru,
            tru_emf=tru_emf,
            kpen=kpen,
            kbup=kbup,
            pifc0=pifc0,
            thv0bot=thv0bot,
            thv0top=thv0top,
            exnifc0=exnifc0,
            umf_zint=umf,
            ppen=ppen,
            rei=rei,
            dp0=dp0,
            thl0=thl0,
            ssthl0=ssthl0,
            pmid0=pmid0,
            qt0=qt0,
            ssqt0=ssqt0,
            u0=u0,
            ssu0=ssu0,
            v0=v0,
            ssv0=ssv0,
            tr0=tr0,
            sstr0=sstr0,
            thlu_emf=thlu_emf,
            qtu_emf=qtu_emf,
            uu_emf=uu_emf,
            vu_emf=vu_emf,
            emf=emf,
            iteration=iter_test,
        )

        return {
            "qtu_emf": qtu_emf.view[:],
            "thlu_emf": thlu_emf.view[:],
            "tru_emf": tru_emf.view[:],
            "uu_emf": uu_emf.view[:],
            "vu_emf": vu_emf.view[:],
            "emf": emf.view[:],
        }
