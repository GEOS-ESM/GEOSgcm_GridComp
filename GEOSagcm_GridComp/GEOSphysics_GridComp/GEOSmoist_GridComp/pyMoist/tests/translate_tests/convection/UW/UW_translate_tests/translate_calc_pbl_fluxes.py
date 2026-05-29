from f90nml import Namelist
from ndsl import StencilFactory
from ndsl.constants import I_DIM, J_DIM, K_DIM, K_INTERFACE_DIM
from ndsl.dsl.gt4py import int32
from ndsl.dsl.typing import Int
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array

from pyMoist.convection.UW.compute_uwshcu import calc_pbl_fluxes
from pyMoist.convection.UW.config import UWConfiguration
from pyMoist.saturation_tables import get_saturation_vapor_pressure_table


class TranslateCalcPblFluxes(TranslateFortranData2Py):
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
            "cbmf": {},
            "condensation": {},
            "kinv": {},
            "pifc0": {},
            "pmid0": {},
            "qt0": {},
            "qtsrc": {},
            "ssqt0": {},
            "ssthl0": {},
            "sstr0": {},
            "ssu0": {},
            "ssv0": {},
            "thl0": {},
            "thlsrc": {},
            "tr0_PblFlux": {},
            "trsrc": {},
            "u0": {},
            "v0": {},
            "vsrc": {},
            "usrc": {},
            "exnifc0": {},
        }

        # FloatField Outputs
        self.out_vars = {
            "qtflx": self.grid.compute_dict(),
            "slflx": self.grid.compute_dict(),
            "trflx": self.grid.compute_dict(),
            "uflx": self.grid.compute_dict(),
            "vflx": self.grid.compute_dict(),
            "xflx": self.grid.compute_dict(),
        }

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("ComputeUwshcuInv-constants")
        self.constants["JASON"] = True

    def compute(self, inputs):
        config = UWConfiguration(**self.constants)

        self.quantity_factory.add_data_dimensions(
            {
                "ntracers": config.NCNST,
            }
        )

        self._calc_pbl_fluxes = self.stencil_factory.from_dims_halo(
            func=calc_pbl_fluxes,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "ncnst": config.NCNST,
                "dt": config.dt,
                "dotransport": config.dotransport,
            },
        )

        # Inputs
        condensation = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a", dtype=bool)

        for i in range(0, 24):
            for j in range(0, 24):
                if inputs["condensation"][i, j] == 1:
                    condensation.view[i, j] = False
                else:
                    condensation.view[i, j] = True

        cbmf = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(cbmf.view[:], inputs["cbmf"])
        kinv = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a", dtype=Int)
        safe_assign_array(kinv.view[:], inputs["kinv"])
        pifc0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(pifc0.view[:], inputs["pifc0"])
        pmid0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(pmid0.view[:], inputs["pmid0"])
        qt0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(qt0.view[:], inputs["qt0"])
        qtsrc = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(qtsrc.view[:], inputs["qtsrc"])
        ssqt0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(ssqt0.view[:], inputs["ssqt0"])
        ssthl0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(ssthl0.view[:], inputs["ssthl0"])
        sstr0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM, "ntracers"], units="n/a")
        safe_assign_array(sstr0.view[:], inputs["sstr0"])
        ssu0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(ssu0.view[:], inputs["ssu0"])
        ssv0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(ssv0.view[:], inputs["ssv0"])
        thl0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(thl0.view[:], inputs["thl0"])
        thlsrc = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(thlsrc.view[:], inputs["thlsrc"])
        tr0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM, "ntracers"], units="n/a")
        safe_assign_array(tr0.view[:], inputs["tr0_PblFlux"])
        trsrc = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, "ntracers"], units="n/a")
        safe_assign_array(trsrc.view[:], inputs["trsrc"])
        u0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(u0.view[:], inputs["u0"])
        v0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(v0.view[:], inputs["v0"])
        vsrc = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(vsrc.view[:], inputs["vsrc"])
        usrc = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(usrc.view[:], inputs["usrc"])
        exnifc0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(exnifc0.view[:], inputs["exnifc0"])

        # Outputs
        qtflx = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        slflx = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        trflx = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM, "ntracers"], units="n/a")
        uflx = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        vflx = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        xflx = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        xflx_ndim = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM, "ntracers"], units="n/a")

        saturation_vapor_pressure_table = get_saturation_vapor_pressure_table(self.stencil_factory.backend)
        self.ese = saturation_vapor_pressure_table.ese
        self.esx = saturation_vapor_pressure_table.esx

        umf_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        qtflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        slflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        uflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        vflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        cush_inout = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")

        # The iteration you want to test
        iter_test = int32(0)

        # # Call stencils
        self._calc_pbl_fluxes(
            condensation=condensation,
            qtsrc=qtsrc,
            qt0=qt0,
            ssqt0=ssqt0,
            pifc0=pifc0,
            pmid0=pmid0,
            kinv=kinv,
            cbmf=cbmf,
            xflx=xflx,
            qtflx=qtflx,
            uflx=uflx,
            vflx=vflx,
            slflx=slflx,
            thlsrc=thlsrc,
            thl0=thl0,
            ssthl0=ssthl0,
            exnifc0=exnifc0,
            usrc=usrc,
            u0=u0,
            ssu0=ssu0,
            vsrc=vsrc,
            v0=v0,
            ssv0=ssv0,
            trsrc=trsrc,
            tr0=tr0,
            sstr0=sstr0,
            trflx=trflx,
            xflx_ndim=xflx_ndim,
            iteration=iter_test,
        )

        return {
            "qtflx": qtflx.view[:],
            "slflx": slflx.view[:],
            "trflx": trflx.view[:],
            "uflx": uflx.view[:],
            "vflx": vflx.view[:],
            "xflx": xflx.view[:],
        }
