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
from pyMoist.UW.compute_uwshcu import recalc_condensate
from pyMoist.UW.config import UWConfiguration


# Dev NOTE: The data for this translate test comes from combining two files in
#           a single nc file using the following NCO tool:
#           `ncks -A ComputeUwshcu-In.nc RecalcCondensate-In.nc`


class TranslateRecalcCondensate(TranslateFortranData2Py):
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

        self._recalc_condensate = self.stencil_factory.from_dims_halo(
            func=recalc_condensate,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ncnst": 23},
        )

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "condensation": {},
            "exnifc0": {},
            "fer": {},
            "kbup": {},
            "kpen": {},
            "pifc0": {},
            "ppen": {},
            "qt0": {},
            "qtu": {},
            "ssqt0": {},
            "ssthl0": {},
            "thl0": {},
            "thlu": {},
            "thv0bot": {},
            "thv0top": {},
            "zifc0": {},
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
            "cldhgt": self.grid.compute_dict(),
            "diten": self.grid.compute_dict(),
            "dwten": self.grid.compute_dict(),
            "emf": self.grid.compute_dict(),
            "fdr": self.grid.compute_dict(),
            "qtu_top": self.grid.compute_dict(),
            "thlu_top": self.grid.compute_dict(),
            "ufrc": self.grid.compute_dict(),
            "umf": self.grid.compute_dict(),
            "xco": self.grid.compute_dict(),
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

        exnifc0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(exnifc0.view[:], inputs["exnifc0"])
        fer = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(fer.view[:], inputs["fer"])
        kbup_IJ = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a", dtype=Int)
        safe_assign_array(kbup_IJ.view[:], inputs["kbup"] - 1)
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
        zifc0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(zifc0.view[:], inputs["zifc0"])
        krel = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int)
        safe_assign_array(krel.view[:], inputs["krel"] - 1)

        # Outputs

        cldhgt = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        diten = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        dwten = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        emf = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        fdr = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        qtu_top = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        thlu_top = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        ufrc = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        umf = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        xco = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        kbup = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int)
        dwten_temp = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        diten_temp = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        umf_temp = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        cush = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")

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
        self._recalc_condensate(
            condensation=condensation,
            fer=fer,
            kpen=kpen,
            ppen=ppen,
            thlu=thlu,
            thl0=thl0,
            ssthl0=ssthl0,
            qtu=qtu,
            qt0=qt0,
            ssqt0=ssqt0,
            pifc0=pifc0,
            ese=self.ese,
            esx=self.esx,
            umf_out=umf_out,
            qtflx_out=qtflx_out,
            slflx_out=slflx_out,
            uflx_out=uflx_out,
            vflx_out=vflx_out,
            criqc=criqc,
            thv0bot=thv0bot,
            thv0top=thv0top,
            exnifc0=exnifc0,
            zifc0=zifc0,
            kbup_IJ=kbup_IJ,
            kbup=kbup,
            krel=krel,
            k0=k0,
            umf_zint=umf,
            emf=emf,
            ufrc=ufrc,
            dwten=dwten,
            diten=diten,
            dwten_temp=dwten_temp,
            diten_temp=diten_temp,
            thlu_top=thlu_top,
            qtu_top=qtu_top,
            cldhgt=cldhgt,
            umf_temp=umf_temp,
            fdr=fdr,
            xco=xco,
            cush=cush,
            cush_inout=cush_inout,
            iteration=iter_test,
        )

        return {
            "cldhgt": cldhgt.view[:],
            "diten": diten.view[:],
            "dwten": dwten.view[:],
            "emf": emf.view[:],
            "fdr": fdr.view[:],
            "qtu_top": qtu_top.view[:],
            "thlu_top": thlu_top.view[:],
            "ufrc": ufrc.view[:],
            "umf": umf.view[:],
            "xco": xco.view[:],
        }
