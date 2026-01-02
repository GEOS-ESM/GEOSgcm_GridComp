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
from pyMoist.UW.compute_uwshcu import compute_cin_cinlcl
from pyMoist.UW.config import UWConfiguration


# Dev NOTE: The data for this translate test comes from combining two files in
#           a single nc file using the following NCO tool:
#           `ncks -A ComputeUwshcu-In.nc ComputeCinCinlcl-In.nc`


class TranslateComputeCinCinlcl(TranslateFortranData2Py):
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

        self._compute_cin_cinlcl = self.stencil_factory.from_dims_halo(
            func=compute_cin_cinlcl,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ncnst": 23},
        )

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "condensation": {},
            "kinv": {},
            "klcl": {},
            "pifc0": {},
            "plcl": {},
            "qtsrc": {},
            "thlsrc": {},
            "thv0lcl": {},
            "thv0top": {},
            "thv0bot": {},
            "thvlmin": {},
            "thvlsrc": {},
            "tkeavg": {},
            "trsrc": {},
            "usrc": {},
            "vsrc": {},
            "rkfre": {},
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
            "cin_i": self.grid.compute_dict(),
            "cinlcl_i": self.grid.compute_dict(),
            "ke": self.grid.compute_dict(),
            "kinv_o": self.grid.compute_dict(),
            "klfc_o": self.grid.compute_dict(),
            "plcl_o": self.grid.compute_dict(),
            "plfc_o": self.grid.compute_dict(),
            "qtsrc_o": self.grid.compute_dict(),
            "thlsrc_o": self.grid.compute_dict(),
            "thvlmin_o": self.grid.compute_dict(),
            "thv0lcl_o": self.grid.compute_dict(),
            "tkeavg_o": self.grid.compute_dict(),
            "trsrc_o": self.grid.compute_dict(),
            "usrc_o": self.grid.compute_dict(),
            "vsrc_o": self.grid.compute_dict(),
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
        klcl = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int)
        safe_assign_array(klcl.view[:], inputs["klcl"] - 1)
        pifc0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(pifc0.view[:], inputs["pifc0"])
        plcl = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(plcl.view[:], inputs["plcl"])
        qtsrc = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(qtsrc.view[:], inputs["qtsrc"])
        thlsrc = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(thlsrc.view[:], inputs["thlsrc"])
        thv0bot = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(thv0bot.view[:], inputs["thv0bot"])
        thv0lcl = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(thv0lcl.view[:], inputs["thv0lcl"])
        thv0top = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(thv0top.view[:], inputs["thv0top"])
        thvlmin = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(thvlmin.view[:], inputs["thvlmin"])
        thvlsrc = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(thvlsrc.view[:], inputs["thvlsrc"])
        tkeavg = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(tkeavg.view[:], inputs["tkeavg"])
        trsrc = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, "ntracers"], units="n/a")
        safe_assign_array(trsrc.view[:], inputs["trsrc"])
        usrc = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(usrc.view[:], inputs["usrc"])
        vsrc = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(vsrc.view[:], inputs["vsrc"])
        rkfre = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        safe_assign_array(rkfre.view[:], inputs["rkfre"])

        # Outputs
        cin_i = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        cinlcl_i = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        ke = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        kinv_o = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int)
        klfc_o = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int)
        plcl_o = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        plfc_o = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        qtsrc_o = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        thlsrc_o = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        thvlmin_o = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        thv0lcl_o = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        tkeavg_o = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        trsrc_o = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, "ntracers"], units="n/a")
        usrc_o = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        vsrc_o = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        cin_IJ = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        cinlcl_IJ = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        plfc_IJ = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        klfc_IJ = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a", dtype=Int)
        plfc = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        klfc = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int)
        cin = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        cinlcl = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        thvubot = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        thvutop = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        klcl_o = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int)
        thvlsrc_o = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")

        umf_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        qtflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        slflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        uflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        vflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        cush_inout = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")

        saturation_vapor_pressure_table = get_saturation_vapor_pressure_table(self.stencil_factory.backend)
        self.ese = saturation_vapor_pressure_table.ese
        self.esx = saturation_vapor_pressure_table.esx

        stop_cin = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a", dtype=bool)

        # The iteration you want to test
        iter_test = int32(0)

        # # Call stencils
        self._compute_cin_cinlcl(
            condensation=condensation,
            stop_cin=stop_cin,
            klcl=klcl,
            kinv=kinv,
            thvlsrc=thvlsrc,
            pifc0=pifc0,
            thv0bot=thv0bot,
            thv0top=thv0top,
            plcl=plcl,
            thv0lcl=thv0lcl,
            thlsrc=thlsrc,
            qtsrc=qtsrc,
            ese=self.ese,
            esx=self.esx,
            umf_out=umf_out,
            qtflx_out=qtflx_out,
            slflx_out=slflx_out,
            uflx_out=uflx_out,
            vflx_out=vflx_out,
            cin_IJ=cin_IJ,
            cinlcl_IJ=cinlcl_IJ,
            plfc_IJ=plfc_IJ,
            klfc_IJ=klfc_IJ,
            plfc=plfc,
            klfc=klfc,
            cin=cin,
            thvubot=thvubot,
            thvutop=thvutop,
            k0=k0,
            iteration=iter_test,
            rbuoy=rbuoy,
            rkfre=rkfre,
            tkeavg=tkeavg,
            epsvarw=epsvarw,
            thvlmin=thvlmin,
            usrc=usrc,
            vsrc=vsrc,
            dotransport=dotransport,
            trsrc=trsrc,
            trsrc_o=trsrc_o,
            cin_i=cin_i,
            cinlcl_i=cinlcl_i,
            ke=ke,
            kinv_o=kinv_o,
            klcl_o=klcl_o,
            klfc_o=klfc_o,
            plcl_o=plcl_o,
            plfc_o=plfc_o,
            tkeavg_o=tkeavg_o,
            thvlmin_o=thvlmin_o,
            qtsrc_o=qtsrc_o,
            thvlsrc_o=thvlsrc_o,
            thlsrc_o=thlsrc_o,
            usrc_o=usrc_o,
            vsrc_o=vsrc_o,
            thv0lcl_o=thv0lcl_o,
            cinlcl=cinlcl,
            cush_inout=cush_inout,
        )

        return {
            "cin_i": cin_i.view[:],
            "cinlcl_i": cinlcl_i.view[:],
            "ke": ke.view[:],
            "kinv_o": kinv_o.view[:],  # kinv_o will fail by 1
            "klfc_o": klfc_o.view[:],  # klfc_o will fail by 1
            "plcl_o": plcl_o.view[:],
            "plfc_o": plfc_o.view[:],
            "qtsrc_o": qtsrc_o.view[:],
            "thlsrc_o": thlsrc_o.view[:],
            "thvlmin_o": thvlmin_o.view[:],
            "thv0lcl_o": thv0lcl_o.view[:],
            "tkeavg_o": tkeavg_o.view[:],
            "trsrc_o": trsrc_o.view[:],
            "usrc_o": usrc_o.view[:],
            "vsrc_o": vsrc_o.view[:],
        }
