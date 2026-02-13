from f90nml import Namelist
from gt4py.cartesian.gtscript import int32

import pyMoist.constants as constants
from ndsl import StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Float, Int
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array
from pyMoist.saturation_tables import get_saturation_vapor_pressure_table
from pyMoist.UW.compute_uwshcu import buoyancy_sorting_fluxes
from pyMoist.UW.config import UWConfiguration


class TranslateBuoyancySortingFluxes(TranslateFortranData2Py):
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

        self._buoyancy_sorting_fluxes = self.stencil_factory.from_dims_halo(
            func=buoyancy_sorting_fluxes,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ncnst": 23},
        )

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "condensation": {},
            "exnifc0": {},
            "kbup": {},
            "krel": {},
            "pifc0": {},
            "pmid0": {},
            "qt0": {},
            "qtflx": {},
            "qtu": {},
            "slflx": {},
            "ssqt0": {},
            "ssthl0": {},
            "sstr0": {},
            "ssu0": {},
            "ssv0": {},
            "thl0": {},
            "thlu": {},
            "tr0_BuoySortFlux": {},
            "trflx": {},
            "tru": {},
            "u0": {},
            "v0": {},
            "vflx": {},
            "vu": {},
            "uflx": {},
            "umf": {},
            "uu": {},
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
        krel = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int)
        safe_assign_array(krel.view[:], inputs["krel"] - 1)
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
        safe_assign_array(tr0.view[:], inputs["tr0_BuoySortFlux"])
        u0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(u0.view[:], inputs["u0"])
        v0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(v0.view[:], inputs["v0"])
        exnifc0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(exnifc0.view[:], inputs["exnifc0"])
        umf = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(umf.view[:], inputs["umf"])
        thlu = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(thlu.view[:], inputs["thlu"])
        qtu = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(qtu.view[:], inputs["qtu"])
        uu = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(uu.view[:], inputs["uu"])
        vu = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(vu.view[:], inputs["vu"])
        tru = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM, "ntracers"], units="n/a")
        safe_assign_array(tru.view[:], inputs["tru"])

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

        # The iteration you want to test
        iter_test = int32(0)

        # # Call stencils
        self._buoyancy_sorting_fluxes(
            condensation=condensation,
            kbup=kbup,
            krel=krel,
            exnifc0=exnifc0,
            umf_zint=umf,
            thlu=thlu,
            thl0=thl0,
            ssthl0=ssthl0,
            pifc0=pifc0,
            pmid0=pmid0,
            qtu=qtu,
            qt0=qt0,
            ssqt0=ssqt0,
            uu=uu,
            u0=u0,
            v0=v0,
            vu=vu,
            ssu0=ssu0,
            ssv0=ssv0,
            dotransport=dotransport,
            trflx=trflx,
            tru=tru,
            tr0=tr0,
            sstr0=sstr0,
            qtflx=qtflx,
            uflx=uflx,
            vflx=vflx,
            slflx=slflx,
            iteration=iter_test,
        )

        return {
            "qtflx": qtflx.view[:],
            "slflx": slflx.view[:],
            "trflx": trflx.view[:],
            "uflx": uflx.view[:],
            "vflx": vflx.view[:],
        }
