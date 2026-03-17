from f90nml import Namelist
from gt4py.cartesian.gtscript import int32
from ndsl import StencilFactory
from ndsl.constants import I_DIM, J_DIM, K_DIM, K_INTERFACE_DIM
from ndsl.dsl.typing import Int
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array

import pyMoist.constants as constants
from pyMoist.saturation_tables import get_saturation_vapor_pressure_table
from pyMoist.UW.compute_uwshcu import penetrative_entrainment_fluxes
from pyMoist.UW.config import UWConfiguration


class TranslatePenetrativeEntrainmentFluxes(TranslateFortranData2Py):
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
            "emf": {},
            "kbup": {},
            "kpen": {},
            "pifc0": {},
            "pmid0": {},
            "qt0": {},
            "qtflx": {},
            "qtu_emf": {},
            "slflx": {},
            "ssqt0": {},
            "ssthl0": {},
            "sstr0": {},
            "ssu0": {},
            "ssv0": {},
            "thl0": {},
            "thlu_emf": {},
            "tr0_PenEntrainFlux": {},
            "trflx": {},
            "tru_emf": {},
            "u0": {},
            "v0": {},
            "vflx": {},
            "vu_emf": {},
            "uu_emf": {},
            "uflx": {},
            "umf": {},
            "kinv": {},
            "cbmf": {},
            "ql0": {},
            "qi0": {},
            "exnifc0": {},
            "krel": {},
        }

        # FloatField Outputs
        self.out_vars = {
            "qtflx": self.grid.compute_dict(),
            "slflx": self.grid.compute_dict(),
            "trflx": self.grid.compute_dict(),
            "uflx": self.grid.compute_dict(),
            "vflx": self.grid.compute_dict(),
            "qlten_sink": self.grid.compute_dict(),
            "qiten_sink": self.grid.compute_dict(),
            "uemf": self.grid.compute_dict(),
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

        self._penetrative_entrainment_fluxes = self.stencil_factory.from_dims_halo(
            func=penetrative_entrainment_fluxes,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "ncnst": config.NCNST,
                "k0": config.k0,
                "dt": config.dt,
                "dotransport": config.dotransport,
                "use_momenflx": config.use_momenflx,
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

        kbup = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a", dtype=Int)
        safe_assign_array(kbup.view[:], inputs["kbup"] - 1)
        kpen = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a", dtype=Int)
        safe_assign_array(kpen.view[:], inputs["kpen"] - 1)
        pifc0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(pifc0.view[:], inputs["pifc0"])
        pmid0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(pmid0.view[:], inputs["pmid0"])
        qt0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(qt0.view[:], inputs["qt0"])
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
        tr0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM, "ntracers"], units="n/a")
        safe_assign_array(tr0.view[:], inputs["tr0_PenEntrainFlux"])
        u0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(u0.view[:], inputs["u0"])
        v0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(v0.view[:], inputs["v0"])
        exnifc0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(exnifc0.view[:], inputs["exnifc0"])
        umf = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(umf.view[:], inputs["umf"])
        thlu_emf = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(thlu_emf.view[:], inputs["thlu_emf"])
        qtu_emf = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(qtu_emf.view[:], inputs["qtu_emf"])
        uu_emf = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(uu_emf.view[:], inputs["uu_emf"])
        vu_emf = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(vu_emf.view[:], inputs["vu_emf"])
        tru_emf = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM, "ntracers"], units="n/a")
        safe_assign_array(tru_emf.view[:], inputs["tru_emf"])
        emf = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(emf.view[:], inputs["emf"])
        cbmf = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(cbmf.view[:], inputs["cbmf"])
        kinv = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a", dtype=Int)
        safe_assign_array(kinv.view[:], inputs["kinv"] - 1)
        krel = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a", dtype=Int)
        safe_assign_array(krel.view[:], inputs["krel"] - 1)
        qi0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(qi0.view[:], inputs["qi0"])
        ql0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(ql0.view[:], inputs["ql0"])

        # Outputs
        qtflx = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(qtflx.view[:], inputs["qtflx"])
        slflx = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(slflx.view[:], inputs["slflx"])
        trflx = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM, "ntracers"], units="n/a")
        safe_assign_array(trflx.view[:], inputs["trflx"])
        uflx = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(uflx.view[:], inputs["uflx"])
        vflx = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(vflx.view[:], inputs["vflx"])

        saturation_vapor_pressure_table = get_saturation_vapor_pressure_table(self.stencil_factory.backend)
        self.ese = saturation_vapor_pressure_table.ese
        self.esx = saturation_vapor_pressure_table.esx

        umf_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        qtflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        slflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        uflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        vflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        cush_inout = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        cush = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        fer_out= self.quantity_factory.zeros(dims=[I_DIM, J_DIM,K_DIM], units="n/a")
        fdr_out=self.quantity_factory.zeros(dims=[I_DIM, J_DIM,K_DIM], units="n/a")
        qldet_out=self.quantity_factory.zeros(dims=[I_DIM, J_DIM,K_DIM], units="n/a")
        qidet_out=self.quantity_factory.zeros(dims=[I_DIM, J_DIM,K_DIM], units="n/a")
        dcm_out=self.quantity_factory.zeros(dims=[I_DIM, J_DIM,K_DIM], units="n/a")
        qvten_out=self.quantity_factory.zeros(dims=[I_DIM, J_DIM,K_DIM], units="n/a")
        qlten_out=self.quantity_factory.zeros(dims=[I_DIM, J_DIM,K_DIM], units="n/a")
        qiten_out=self.quantity_factory.zeros(dims=[I_DIM, J_DIM,K_DIM], units="n/a")
        sten_out=self.quantity_factory.zeros(dims=[I_DIM, J_DIM,K_DIM], units="n/a")
        uten_out=self.quantity_factory.zeros(dims=[I_DIM, J_DIM,K_DIM], units="n/a")
        vten_out=self.quantity_factory.zeros(dims=[I_DIM, J_DIM,K_DIM], units="n/a")
        qrten_out=self.quantity_factory.zeros(dims=[I_DIM, J_DIM,K_DIM], units="n/a")
        qsten_out=self.quantity_factory.zeros(dims=[I_DIM, J_DIM,K_DIM], units="n/a")
        cufrc_out=self.quantity_factory.zeros(dims=[I_DIM, J_DIM,K_DIM], units="n/a")
        

        uemf = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        qlten_sink = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qiten_sink = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        # The iteration you want to test
        iter_test = int32(0)

        # # Call stencils
        self._penetrative_entrainment_fluxes(
            condensation=condensation,
            kbup=kbup,
            kpen=kpen,
            exnifc0=exnifc0,
            emf=emf,
            thlu_emf=thlu_emf,
            thl0=thl0,
            ssthl0=ssthl0,
            pifc0=pifc0,
            pmid0=pmid0,
            qtu_emf=qtu_emf,
            qt0=qt0,
            ssqt0=ssqt0,
            uu_emf=uu_emf,
            vu_emf=vu_emf,
            u0=u0,
            v0=v0,
            ssu0=ssu0,
            ssv0=ssv0,
            trflx=trflx,
            tru_emf=tru_emf,
            tr0=tr0,
            sstr0=sstr0,
            kinv=kinv,
            cbmf=cbmf,
            uflx=uflx,
            vflx=vflx,
            slflx=slflx,
            qtflx=qtflx,
            uemf=uemf,
            krel=krel,
            umf_zint=umf,
            ql0=ql0,
            qi0=qi0,
            ese=self.ese,
            esx=self.esx,
            qlten_sink=qlten_sink,
            qiten_sink=qiten_sink,
            iteration=iter_test,
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

        return {
            "qtflx": qtflx.view[:],
            "slflx": slflx.view[:],
            "trflx": trflx.view[:],
            "uflx": uflx.view[:],
            "vflx": vflx.view[:],
            "qlten_sink": qlten_sink.view[:],
            "qiten_sink": qiten_sink.view[:],
            "uemf": uemf.view[:],
        }
