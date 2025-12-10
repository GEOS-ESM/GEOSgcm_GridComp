from f90nml import Namelist

from ndsl import Quantity, QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Float, FloatField, Int, Bool
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array
from pyMoist.UW.compute_uwshcu import (
    define_env_properties,
)
from gt4py.cartesian.gtscript import int32
from pyMoist.UW.config import UWConfiguration
import pyMoist.constants as constants
from pyMoist.saturation_tables import (
    get_saturation_vapor_pressure_table,
)


# Dev NOTE: The data for this translate test comes from combining two files in
#           a single nc file using the following NCO tool:
#           `ncks -A ComputeUwshcu-In.nc DefineEnvProperties-In.nc`


class TranslateDefineEnvProperties(TranslateFortranData2Py):
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

        self._define_env_properties = self.stencil_factory.from_dims_halo(
            func=define_env_properties,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ncnst": 23},
        )

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "condensation": {},
            "kinv": {},
            "krel": {},
            "pifc0": {},
            "pmid0": {},
            "prel": {},
            "qt0": {},
            "ssqt0": {},
            "ssthl0": {},
            "sstr0": {},
            "ssu0": {},
            "ssv0": {},
            "thl0": {},
            "trsrc": {},
            "u0": {},
            "v0": {},
            "usrc": {},
            "vsrc": {},
            "thv0rel": {},
            "tr0_DefineEnvProperties": {},
            "uu": {},
            "vu": {},
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
            "dpe": self.grid.compute_dict(),
            "exne": self.grid.compute_dict(),
            "pe": self.grid.compute_dict(),
            "qsat_pe": self.grid.compute_dict(),
            "qte": self.grid.compute_dict(),
            "thle": self.grid.compute_dict(),
            "thvebot": self.grid.compute_dict(),
            "tre": self.grid.compute_dict(),
            "tru": self.grid.compute_dict(),
            "ue": self.grid.compute_dict(),
            "uplus": self.grid.compute_dict(),
            "uu": self.grid.compute_dict(),
            "ve": self.grid.compute_dict(),
            "vplus": self.grid.compute_dict(),
            "vu": self.grid.compute_dict(),
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

        # Field inputs
        condensation = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a", dtype=bool)

        for i in range(0, 24):
            for j in range(0, 24):
                if inputs["condensation"][i, j] == 1:
                    condensation.view[i, j] = False
                else:
                    condensation.view[i, j] = True

        kinv = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int)
        safe_assign_array(kinv.view[:], inputs["kinv"])
        krel = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int)
        safe_assign_array(krel.view[:], inputs["krel"] - 1)
        pifc0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(pifc0.view[:], inputs["pifc0"])
        pmid0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(pmid0.view[:], inputs["pmid0"])
        prel = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(prel.view[:], inputs["prel"])
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
        trsrc = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, "ntracers"], units="n/a")
        safe_assign_array(trsrc.view[:], inputs["trsrc"])
        u0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(u0.view[:], inputs["u0"])
        v0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(v0.view[:], inputs["v0"])
        usrc = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(usrc.view[:], inputs["usrc"])
        vsrc = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(vsrc.view[:], inputs["vsrc"])
        thv0rel = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(thv0rel.view[:], inputs["thv0rel"])
        tr0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM, "ntracers"], units="n/a")
        safe_assign_array(tr0.view[:], inputs["tr0_DefineEnvProperties"])
        uu = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(uu.view[:], inputs["uu"])
        vu = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(vu.view[:], inputs["vu"])
        # Outputs
        dpe = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        exne = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        pe = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        qsat_pe = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        qte = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        thle = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        thvebot = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        tre = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, "ntracers"], units="n/a")
        tru = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM, "ntracers"], units="n/a")
        ue = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        uplus = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")

        ve = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        vplus = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        uplus_3D = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        vplus_3D = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")

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
        self._define_env_properties(
            condensation=condensation,
            iteration=iter_test,
            krel=krel,
            kinv=kinv,
            PGFc=PGFc,
            ssu0=ssu0,
            ssv0=ssv0,
            prel=prel,
            pifc0=pifc0,
            uu=uu,
            vu=vu,
            usrc=usrc,
            vsrc=vsrc,
            dotransport=dotransport,
            tru=tru,
            trsrc=trsrc,
            thv0rel=thv0rel,
            thl0=thl0,
            ssthl0=ssthl0,
            pmid0=pmid0,
            qt0=qt0,
            ssqt0=ssqt0,
            u0=u0,
            v0=v0,
            tre=tre,
            tr0=tr0,
            sstr0=sstr0,
            uplus=uplus,
            vplus=vplus,
            uplus_3D=uplus_3D,
            vplus_3D=vplus_3D,
            qsat_pe=qsat_pe,
            pe=pe,
            thle=thle,
            qte=qte,
            dpe=dpe,
            exne=exne,
            thvebot=thvebot,
            ue=ue,
            ve=ve,
        )

        return {
            "dpe": dpe.view[:],
            "exne": exne.view[:],
            "pe": pe.view[:],
            "qsat_pe": qsat_pe.view[:],
            "qte": qte.view[:],
            "thle": thle.view[:],
            "thvebot": thvebot.view[:],
            "tre": tre.view[:],
            "tru": tru.view[:],
            "ue": ue.view[:],
            "uplus": uplus.view[:],
            "uu": uu.view[:],
            "ve": ve.view[:],
            "vplus": vplus.view[:],
            "vu": vu.view[:],
        }
