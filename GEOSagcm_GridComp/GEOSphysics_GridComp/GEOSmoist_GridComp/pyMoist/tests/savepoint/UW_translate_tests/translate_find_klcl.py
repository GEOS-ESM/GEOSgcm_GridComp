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
from pyMoist.UW.compute_uwshcu import find_cumulus_characteristics, find_klcl
from pyMoist.UW.config import UWConfiguration


# Dev NOTE: The data for this translate test comes from combining two files in
#           a single nc file using the following NCO tool:
#           `ncks -A ComputeUwshcu-In.nc FindKlcl-In.nc`


class TranslateFindKlcl(TranslateFortranData2Py):
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

        self._find_cumulus_characteristics = self.stencil_factory.from_dims_halo(
            func=find_cumulus_characteristics,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ncnst": 23},
        )
        self._find_klcl = self.stencil_factory.from_dims_halo(
            func=find_klcl,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ncnst": 23},
        )

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "condensation": {},
            "pifc0": {},
            "pmid0": {},
            "qt0": {},
            "ssu0": {},
            "ssv0": {},
            "thvlmin": {},
            "u0": {},
            "v0": {},
            "t0": {},
            "qv0": {},
            "qtavg": {},
            "uavg": {},
            "vavg": {},
            "kinv": {},
            "thl0": {},
            "ssthl0": {},
            "ssqt0": {},
            "tr0_FindKlcl": {},
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
            "qtsrc": self.grid.compute_dict(),
            "thlsrc": self.grid.compute_dict(),
            "thvlsrc": self.grid.compute_dict(),
            "usrc": self.grid.compute_dict(),
            "vsrc": self.grid.compute_dict(),
            "trsrc": self.grid.compute_dict(),
            "klcl": self.grid.compute_dict(),
            "plcl": self.grid.compute_dict(),
            "qt0lcl": self.grid.compute_dict(),
            "thl0lcl": self.grid.compute_dict(),
            "thv0lcl": self.grid.compute_dict(),
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

        pifc0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(pifc0.view[:], inputs["pifc0"])
        pmid0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(pmid0.view[:], inputs["pmid0"])
        qt0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(qt0.view[:], inputs["qt0"])
        ssu0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(ssu0.view[:], inputs["ssu0"])
        ssv0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(ssv0.view[:], inputs["ssv0"])
        thvlmin = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(thvlmin.view[:], inputs["thvlmin"])
        u0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(u0.view[:], inputs["u0"])
        v0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(v0.view[:], inputs["v0"])
        t0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(t0.view[:], inputs["t0"])
        qv0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(qv0.view[:], inputs["qv0"])
        qtavg = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(qtavg.view[:], inputs["qtavg"])
        uavg = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(uavg.view[:], inputs["uavg"])
        vavg = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(vavg.view[:], inputs["vavg"])
        kinv = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int)
        safe_assign_array(kinv.view[:], inputs["kinv"])
        thl0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(thl0.view[:], inputs["thl0"])
        ssthl0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(ssthl0.view[:], inputs["ssthl0"])
        ssqt0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(ssqt0.view[:], inputs["ssqt0"])
        tr0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM, "ntracers"], units="n/a")
        safe_assign_array(tr0.view[:], inputs["tr0_FindKlcl"])

        # Outputs
        klcl = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int)
        plcl = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        qt0lcl = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        qtsrc = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        thl0lcl = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        thlsrc = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        thv0lcl = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        thvlsrc = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        trsrc = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, "ntracers"], units="n/a")
        usrc = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        vsrc = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        tpert_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        qpert_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        shfx = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        evap = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        # The iteration you want to test
        iter_test = int32(0)

        # # Call stencils
        self._find_cumulus_characteristics(
            condensation=condensation,
            windsrcavg=windsrcavg,
            pifc0=pifc0,
            t0=t0,
            qv0=qv0,
            shfx=shfx,
            evap=evap,
            thlsrc_fac=thlsrc_fac,
            qtsrc_fac=qtsrc_fac,
            qt0=qt0,
            qtavg=qtavg,
            thvlmin=thvlmin,
            uavg=uavg,
            vavg=vavg,
            kinv=kinv,
            u0=u0,
            v0=v0,
            ssu0=ssu0,
            ssv0=ssv0,
            pmid0=pmid0,
            dotransport=dotransport,
            tr0=tr0,
            trsrc=trsrc,
            qtsrc=qtsrc,
            thvlsrc=thvlsrc,
            thlsrc=thlsrc,
            usrc=usrc,
            vsrc=vsrc,
            tpert_out=tpert_out,
            qpert_out=qpert_out,
            iteration=iter_test,
        )

        saturation_vapor_pressure_table = get_saturation_vapor_pressure_table(self.stencil_factory.backend)
        self.ese = saturation_vapor_pressure_table.ese
        self.esx = saturation_vapor_pressure_table.esx

        umf_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        qtflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        slflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        uflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        vflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        cush_inout = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")

        self._find_klcl(
            condensation=condensation,
            iteration=iter_test,
            pifc0=pifc0,
            umf_out=umf_out,
            qtflx_out=qtflx_out,
            slflx_out=slflx_out,
            uflx_out=uflx_out,
            vflx_out=vflx_out,
            qtsrc=qtsrc,
            thlsrc=thlsrc,
            ese=self.ese,
            esx=self.esx,
            k0=k0,
            thl0=thl0,
            ssthl0=ssthl0,
            pmid0=pmid0,
            qt0=qt0,
            ssqt0=ssqt0,
            plcl=plcl,
            klcl=klcl,
            thl0lcl=thl0lcl,
            qt0lcl=qt0lcl,
            thv0lcl=thv0lcl,
            cush_inout=cush_inout,
        )

        return {
            "qtsrc": qtsrc.view[:],
            "thlsrc": thlsrc.view[:],
            "thvlsrc": thvlsrc.view[:],
            "usrc": usrc.view[:],
            "vsrc": vsrc.view[:],
            "trsrc": trsrc.view[:],
            "klcl": klcl.view[:],  # klcl should fail by 1
            "plcl": plcl.view[:],
            "qt0lcl": qt0lcl.view[:],
            "thl0lcl": thl0lcl.view[:],
            "thv0lcl": thv0lcl.view[:],
        }
