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
from pyMoist.convection.UW.compute_uwshcu import compute_cin_cinlcl
from pyMoist.convection.UW.config import UWConfiguration
from pyMoist.saturation_tables import get_saturation_vapor_pressure_table


class TranslateComputeCinCinlcl(TranslateFortranData2Py):
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

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("ComputeUwshcuInv-constants")

    def compute(self, inputs):
        config = UWConfiguration(**self.constants)

        self.quantity_factory.add_data_dimensions(
            {
                "ntracers": constants.NCNST,
            }
        )

        self._compute_cin_cinlcl = self.stencil_factory.from_dims_halo(
            func=compute_cin_cinlcl,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "ncnst": config.NCNST,
                "dotransport": config.dotransport,
                "k0": config.k0,
                "rbuoy": config.rbuoy,
                "epsvarw": config.epsvarw,
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

        kinv = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a", dtype=Int)
        safe_assign_array(kinv.view[:], inputs["kinv"])
        klcl = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a", dtype=Int)
        safe_assign_array(klcl.view[:], inputs["klcl"] - 1)
        pifc0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(pifc0.view[:], inputs["pifc0"])
        plcl = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(plcl.view[:], inputs["plcl"])
        qtsrc = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(qtsrc.view[:], inputs["qtsrc"])
        thlsrc = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(thlsrc.view[:], inputs["thlsrc"])
        thv0bot = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(thv0bot.view[:], inputs["thv0bot"])
        thv0lcl = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(thv0lcl.view[:], inputs["thv0lcl"])
        thv0top = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(thv0top.view[:], inputs["thv0top"])
        thvlmin = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(thvlmin.view[:], inputs["thvlmin"])
        thvlsrc = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(thvlsrc.view[:], inputs["thvlsrc"])
        tkeavg = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(tkeavg.view[:], inputs["tkeavg"])
        trsrc = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, "ntracers"], units="n/a")
        safe_assign_array(trsrc.view[:], inputs["trsrc"])
        usrc = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(usrc.view[:], inputs["usrc"])
        vsrc = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(vsrc.view[:], inputs["vsrc"])
        rkfre = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(rkfre.view[:], inputs["rkfre"])

        # Outputs
        cin_i = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        cinlcl_i = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        ke = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        kinv_o = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a", dtype=Int)
        klfc_o = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a", dtype=Int)
        plcl_o = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        plfc_o = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qtsrc_o = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        thlsrc_o = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        thvlmin_o = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        thv0lcl_o = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        tkeavg_o = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        trsrc_o = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, "ntracers"], units="n/a")
        usrc_o = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        vsrc_o = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        cin_IJ = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        cinlcl_IJ = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        plfc_IJ = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        klfc_IJ = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a", dtype=Int)
        plfc = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        klfc = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a", dtype=Int)
        cin = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        cinlcl = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        thvubot = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        thvutop = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        klcl_o = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a", dtype=Int)
        thvlsrc_o = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")

        umf_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        qtflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        slflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        uflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        vflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        cush_inout = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
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

        saturation_vapor_pressure_table = get_saturation_vapor_pressure_table(self.stencil_factory.backend)
        self.ese = saturation_vapor_pressure_table.ese
        self.esx = saturation_vapor_pressure_table.esx

        stop_cin = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a", dtype=bool)

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
            cin_IJ=cin_IJ,
            cinlcl_IJ=cinlcl_IJ,
            plfc_IJ=plfc_IJ,
            klfc_IJ=klfc_IJ,
            plfc=plfc,
            klfc=klfc,
            cin=cin,
            RKFRE=rkfre,
            thvubot=thvubot,
            thvutop=thvutop,
            iteration=iter_test,
            tkeavg=tkeavg,
            thvlmin=thvlmin,
            usrc=usrc,
            vsrc=vsrc,
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
