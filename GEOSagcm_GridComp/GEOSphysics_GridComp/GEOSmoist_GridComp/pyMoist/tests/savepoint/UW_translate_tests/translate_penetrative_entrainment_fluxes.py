from f90nml import Namelist

import pyMoist.constants as constants
from ndsl import StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.gt4py import int32
from ndsl.dsl.typing import Float, Int
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array
from pyMoist.saturation_tables import get_saturation_vapor_pressure_table
from pyMoist.UW.compute_uwshcu import penetrative_entrainment_fluxes
from pyMoist.UW.config import UWConfiguration


class TranslatePenetrativeEntrainmentFluxes(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
        # UW_config: UWConfiguration,
    ):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory
        # self.UW_config = UW_config

        self._penetrative_entrainment_fluxes = self.stencil_factory.from_dims_halo(
            func=penetrative_entrainment_fluxes,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ncnst": 23},
        )

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

        # Float/Int Inputs
        self.in_vars["parameters"] = [
            "dotransport",
            "ncnst",
            "k0",
            "tr0",
            "windsrcavg",
            "qtsrchgt",
            "qtsrc_fac",
            "thlsrc_fac",
            "frc_rasn",
            "rbuoy",
            "epsvarw",
            "use_CINcin",
            "mumin1",
            "rmaxfrac",
            "PGFc",
            "niter_xc",
            "criqc",
            "rle",
            "cridist_opt",
            "mixscale",
            "rkm",
            "dt",
            "detrhgt",
            "rdrag",
            "use_self_detrain",
            "detrhgt",
            "use_cumpenent",
            "rpen",
            "use_momenflx",
            "rdrop",
            "iter_cin",
        ]

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

    def compute(self, inputs):
        self.UW_config = UWConfiguration(Int(inputs["ncnst"]), Int(inputs["k0"]), Int(inputs["windsrcavg"]))

        self.quantity_factory.add_data_dimensions(
            {
                "ntracers": constants.NCNST,
            }
        )

        # Float/Int Inputs
        dotransport = Int(inputs["dotransport"])
        k0 = Int(inputs["k0"])
        windsrcavg = Int(inputs["windsrcavg"])
        qtsrchgt = Float(inputs["qtsrchgt"])
        qtsrc_fac = Float(inputs["qtsrc_fac"])
        thlsrc_fac = Float(inputs["thlsrc_fac"])
        frc_rasn = Float(inputs["frc_rasn"])
        rbuoy = Float(inputs["rbuoy"])
        epsvarw = Float(inputs["epsvarw"])
        use_CINcin = Int(inputs["use_CINcin"])
        mumin1 = Float(inputs["mumin1"])
        rmaxfrac = Float(inputs["rmaxfrac"])
        PGFc = Float(inputs["PGFc"])
        dt = Float(inputs["dt"])
        niter_xc = Int(inputs["niter_xc"])
        criqc = Float(inputs["criqc"])
        rle = Float(inputs["rle"])
        cridist_opt = Int(inputs["cridist_opt"])
        mixscale = Float(inputs["mixscale"])
        rdrag = Float(inputs["rdrag"])
        rkm = Float(inputs["rkm"])
        use_self_detrain = Int(inputs["use_self_detrain"])
        detrhgt = Float(inputs["detrhgt"])
        use_cumpenent = Int(inputs["use_cumpenent"])
        rpen = Float(inputs["rpen"])
        use_momenflx = Int(inputs["use_momenflx"])
        rdrop = Float(inputs["rdrop"])
        iter_cin = Int(inputs["iter_cin"])

        # Inputs
        condensation = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a", dtype=bool)

        for i in range(0, 24):
            for j in range(0, 24):
                if inputs["condensation"][i, j] == 1:
                    condensation.view[i, j] = False
                else:
                    condensation.view[i, j] = True

        kbup = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int)
        safe_assign_array(kbup.view[:], inputs["kbup"] - 1)
        kpen = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int)
        safe_assign_array(kpen.view[:], inputs["kpen"] - 1)
        pifc0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(pifc0.view[:], inputs["pifc0"])
        pmid0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(pmid0.view[:], inputs["pmid0"])
        qt0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(qt0.view[:], inputs["qt0"])
        ssqt0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(ssqt0.view[:], inputs["ssqt0"])
        ssthl0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(ssthl0.view[:], inputs["ssthl0"])
        sstr0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM, "ntracers"], units="n/a")
        safe_assign_array(sstr0.view[:], inputs["sstr0"])
        ssu0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(ssu0.view[:], inputs["ssu0"])
        ssv0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(ssv0.view[:], inputs["ssv0"])
        thl0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(thl0.view[:], inputs["thl0"])
        tr0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM, "ntracers"], units="n/a")
        safe_assign_array(tr0.view[:], inputs["tr0_PenEntrainFlux"])
        u0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(u0.view[:], inputs["u0"])
        v0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(v0.view[:], inputs["v0"])
        exnifc0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(exnifc0.view[:], inputs["exnifc0"])
        umf = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(umf.view[:], inputs["umf"])
        thlu_emf = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(thlu_emf.view[:], inputs["thlu_emf"])
        qtu_emf = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(qtu_emf.view[:], inputs["qtu_emf"])
        uu_emf = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(uu_emf.view[:], inputs["uu_emf"])
        vu_emf = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(vu_emf.view[:], inputs["vu_emf"])
        tru_emf = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM, "ntracers"], units="n/a")
        safe_assign_array(tru_emf.view[:], inputs["tru_emf"])
        emf = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(emf.view[:], inputs["emf"])
        cbmf = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(cbmf.view[:], inputs["cbmf"])
        kinv = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int)
        safe_assign_array(kinv.view[:], inputs["kinv"] - 1)
        krel = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int)
        safe_assign_array(krel.view[:], inputs["krel"] - 1)
        qi0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(qi0.view[:], inputs["qi0"])
        ql0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(ql0.view[:], inputs["ql0"])

        # Outputs
        qtflx = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(qtflx.view[:], inputs["qtflx"])
        slflx = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(slflx.view[:], inputs["slflx"])
        trflx = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM, "ntracers"], units="n/a")
        safe_assign_array(trflx.view[:], inputs["trflx"])
        uflx = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(uflx.view[:], inputs["uflx"])
        vflx = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(vflx.view[:], inputs["vflx"])

        saturation_vapor_pressure_table = get_saturation_vapor_pressure_table(self.stencil_factory.backend)
        self.ese = saturation_vapor_pressure_table.ese
        self.esx = saturation_vapor_pressure_table.esx

        umf_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        qtflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        slflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        uflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        vflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        cush_inout = umf_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")

        uemf = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        qlten_sink = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        qiten_sink = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
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
            dotransport=dotransport,
            trflx=trflx,
            tru_emf=tru_emf,
            tr0=tr0,
            sstr0=sstr0,
            use_momenflx=use_momenflx,
            kinv=kinv,
            cbmf=cbmf,
            uflx=uflx,
            vflx=vflx,
            slflx=slflx,
            qtflx=qtflx,
            uemf=uemf,
            krel=krel,
            umf_zint=umf,
            k0=k0,
            ql0=ql0,
            qi0=qi0,
            dt=dt,
            ese=self.ese,
            esx=self.esx,
            umf_out=umf_out,
            qtflx_out=qtflx_out,
            slflx_out=slflx_out,
            uflx_out=uflx_out,
            vflx_out=vflx_out,
            qlten_sink=qlten_sink,
            qiten_sink=qiten_sink,
            cush_inout=cush_inout,
            iteration=iter_test,
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
