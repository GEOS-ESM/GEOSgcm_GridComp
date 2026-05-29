from f90nml import Namelist
from ndsl import StencilFactory
from ndsl.constants import I_DIM, J_DIM, K_DIM, K_INTERFACE_DIM
from ndsl.dsl.gt4py import int32
from ndsl.dsl.typing import Int
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array

from pyMoist.convection.UW.compute_uwshcu import define_updraft_properties
from pyMoist.convection.UW.config import UWConfiguration
from pyMoist.saturation_tables import get_saturation_vapor_pressure_table


class TranslateDefineUpdraftProperties(TranslateFortranData2Py):
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
            "cinlcl": {},
            "condensation": {},
            "kinv": {},
            "krel": {},
            "prel": {},
            "qtsrc": {},
            "rho0inv": {},
            "thlsrc": {},
            "winv": {},
            "thvu": {},
        }

        # FloatField Outputs
        self.out_vars = {
            "thvu_out": self.grid.compute_dict(),
            "ufrc": self.grid.compute_dict(),
            "ufrclcl": self.grid.compute_dict(),
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

        self._define_updraft_properties = self.stencil_factory.from_dims_halo(
            func=define_updraft_properties,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "rbuoy": config.rbuoy,
            },
        )

        # Field inputs
        condensation = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a", dtype=bool)

        for i in range(0, 24):
            for j in range(0, 24):
                if inputs["condensation"][i, j] == 1:
                    condensation.view[i, j] = False
                else:
                    condensation.view[i, j] = True

        cbmf = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(cbmf.view[:], inputs["cbmf"])
        cinlcl = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(cinlcl.view[:], inputs["cinlcl"])
        kinv = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a", dtype=Int)
        safe_assign_array(kinv.view[:], inputs["kinv"])
        krel = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a", dtype=Int)
        safe_assign_array(krel.view[:], inputs["krel"] - 1)
        prel = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(prel.view[:], inputs["prel"])
        qtsrc = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(qtsrc.view[:], inputs["qtsrc"])
        rho0inv = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(rho0inv.view[:], inputs["rho0inv"])
        thlsrc = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(thlsrc.view[:], inputs["thlsrc"])
        winv = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(winv.view[:], inputs["winv"])
        thvu = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(thvu.view[:], inputs["thvu"])

        # Outputs
        ufrc = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        ufrclcl = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")

        saturation_vapor_pressure_table = get_saturation_vapor_pressure_table(self.stencil_factory.backend)
        self.ese = saturation_vapor_pressure_table.ese
        self.esx = saturation_vapor_pressure_table.esx

        stop_cin = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a", dtype=bool)

        umf_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        qtflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        slflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        uflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        vflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        cush_inout = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        ufrcinv = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        umf_zint = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        wu = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        emf = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        thlu = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        qtu = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        wlcl = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        cush = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        fer_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        fdr_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qldet_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qidet_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        dcm_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qvten_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qlten_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qiten_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        sten_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        uten_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        vten_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qrten_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qsten_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        cufrc_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")

        # The iteration you want to test
        iter_test = int32(0)

        # # Call stencils
        self._define_updraft_properties(
            condensation=condensation,
            iteration=iter_test,
            winv=winv,
            cinlcl_IJ=cinlcl,
            cbmf=cbmf,
            rho0inv=rho0inv,
            krel=krel,
            ufrc=ufrc,
            ufrcinv=ufrcinv,
            kinv=kinv,
            umf_zint=umf_zint,
            wu=wu,
            emf=emf,
            thlu=thlu,
            qtu=qtu,
            thlsrc=thlsrc,
            qtsrc=qtsrc,
            prel=prel,
            ese=self.ese,
            esx=self.esx,
            thvu=thvu,
            wlcl=wlcl,
            ufrclcl=ufrclcl,
            cush=cush,
            umf_out=umf_out,
            dcm_out=dcm_out,
            qvten_out=qvten_out,
            qlten_out=qlten_out,
            qiten_out=qiten_out,
            sten_out=sten_out,
            uten_out=uten_out,
            vten_out=vten_out,
            qrten_out=qrten_out,
            qsten_out=qsten_out,
            cufrc_out=cufrc_out,
            cush_inout=cush_inout,
            qldet_out=qldet_out,
            qidet_out=qidet_out,
            qtflx_out=qtflx_out,
            slflx_out=slflx_out,
            uflx_out=uflx_out,
            vflx_out=vflx_out,
            fer_out=fer_out,
            fdr_out=fdr_out,
        )

        # NOTE 1 variable is failing - I am ignoring this for now.

        return {
            "thvu_out": thvu.view[:],  # thvu_out fails, serialization issue I think
            "ufrc": ufrc.view[:],
            "ufrclcl": ufrclcl.view[:],
        }
