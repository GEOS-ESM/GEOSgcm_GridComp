from f90nml import Namelist
from ndsl import StencilFactory
from ndsl.constants import I_DIM, J_DIM, K_DIM, K_INTERFACE_DIM
from ndsl.dsl.gt4py import int32
from ndsl.dsl.typing import Int
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array

import pyMoist.constants as constants
from pyMoist.convection.UW.compute_uwshcu import calc_cumulus_base_mass_flux, define_prel_krel
from pyMoist.convection.UW.config import UWConfiguration
from pyMoist.saturation_tables import get_saturation_vapor_pressure_table


class TranslateDefinePrelCbmf(TranslateFortranData2Py):
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
            "cin": {},
            "cinlcl": {},
            "dp0": {},
            "exnifc0": {},
            "kinv": {},
            "klcl": {},
            "pifc0": {},
            "plcl": {},
            "thv0bot": {},
            "thv0lcl": {},
            "thv0top": {},
            "tkeavg": {},
            "rkfre": {},
        }

        # FloatField Outputs
        self.out_vars = {
            "cbmf": self.grid.compute_dict(),
            "krel": self.grid.compute_dict(),
            "prel": self.grid.compute_dict(),
            "thv0rel": self.grid.compute_dict(),
            "ufrcinv": self.grid.compute_dict(),
            "wcrit": self.grid.compute_dict(),
            "winv": self.grid.compute_dict(),
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

        self._define_prel_krel = self.stencil_factory.from_dims_halo(
            func=define_prel_krel,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )

        self._calc_cumulus_base_mass_flux = self.stencil_factory.from_dims_halo(
            func=calc_cumulus_base_mass_flux,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "use_CINcin": config.use_CINcin,
                "rbuoy": config.rbuoy,
                "epsvarw": config.epsvarw,
                "dt": config.dt,
                "mumin1": config.mumin1,
                "rmaxfrac": config.rmaxfrac,
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
        cin = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(cin.view[:], inputs["cin"])
        cinlcl = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(cinlcl.view[:], inputs["cinlcl"])
        dp0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(dp0.view[:], inputs["dp0"])
        exnifc0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(exnifc0.view[:], inputs["exnifc0"])
        pifc0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(pifc0.view[:], inputs["pifc0"])
        plcl = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(plcl.view[:], inputs["plcl"])
        thv0bot = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(thv0bot.view[:], inputs["thv0bot"])
        thv0lcl = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(thv0lcl.view[:], inputs["thv0lcl"])
        thv0top = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(thv0top.view[:], inputs["thv0top"])
        tkeavg = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(tkeavg.view[:], inputs["tkeavg"])
        rkfre = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(rkfre.view[:], inputs["rkfre"])

        # Outputs
        cbmf = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        krel = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a", dtype=Int)
        prel = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        thv0rel = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        ufrcinv = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        wcrit = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        winv = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        rho0inv = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
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
        self._define_prel_krel(
            condensation=condensation,
            iteration=iter_test,
            klcl=klcl,
            kinv=kinv,
            pifc0=pifc0,
            thv0bot=thv0bot,
            plcl=plcl,
            thv0lcl=thv0lcl,
            krel=krel,
            prel=prel,
            thv0rel=thv0rel,
        )

        self._calc_cumulus_base_mass_flux(
            condensation=condensation,
            iteration=iter_test,
            cin_IJ=cin,
            cinlcl_IJ=cinlcl,
            RKFRE=rkfre,
            tkeavg=tkeavg,
            kinv=kinv,
            pifc0=pifc0,
            thv0top=thv0top,
            exnifc0=exnifc0,
            dp0=dp0,
            winv=winv,
            cbmf=cbmf,
            rho0inv=rho0inv,
            ufrcinv=ufrcinv,
            wcrit=wcrit,
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
            "cbmf": cbmf.view[:],
            "krel": krel.view[:],  # krel will fail by 1
            "prel": prel.view[:],
            "thv0rel": thv0rel.view[:],
            "ufrcinv": ufrcinv.view[:],
            "wcrit": wcrit.view[:],
            "winv": winv.view[:],
        }
